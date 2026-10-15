// Copyright (c) 2026, Oliver Kohlbacher and the DIALibGen authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/LibraryExport.h>
#include <odia/AtomicFile.h>
#include <odia/DIANNLibraryFile.h>
#include <odia/ModifiedSequence.h>
#include <odia/Trace.h>

#include "Parallel.h"

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

namespace ODIA
{
  namespace
  {
    struct FormatName { LibraryFormat format; const char* name; const char* extension; };
    constexpr FormatName FORMATS[] = {
      {LibraryFormat::DIALibGenParquet, "dialibgen_parquet", ".parquet"},
      {LibraryFormat::DIANNTSV, "diann_tsv", ".tsv"},
      {LibraryFormat::DIANN2Parquet, "diann2_parquet", ".parquet"},
      {LibraryFormat::OpenSwathPQP, "openswath_pqp", ".pqp"},
      {LibraryFormat::OpenSwathTSV, "openswath_tsv", ".tsv"},
    };

    /// The extension as written. Case-sensitive, as -out always was: the
    /// streaming writer and refine's sidecars key off the exact suffix, and a
    /// `.PARQUET` must not be resolved to a layout that a later check does not
    /// recognise.
    std::string extensionOf(const std::string& filename)
    {
      return std::filesystem::path(filename).extension().string();
    }

    void ok(const arrow::Status& status)
    {
      if (!status.ok()) { throw std::runtime_error("Parquet: " + status.ToString()); }
    }

    /// A member of a protein group, "sp|P02768|ALBU_HUMAN", as accession and
    /// entry name; a member without the db|accession|name form is both.
    std::pair<std::string_view, std::string_view> accessionAndName(std::string_view member)
    {
      const auto bar1 = member.find('|');
      const auto bar2 = bar1 == std::string_view::npos ? bar1 : member.find('|', bar1 + 1);
      if (bar2 == std::string_view::npos) { return {member, member}; }
      return {member.substr(bar1 + 1, bar2 - bar1 - 1), member.substr(bar2 + 1)};
    }

    std::vector<std::string_view> members(std::string_view group)
    {
      std::vector<std::string_view> out;
      std::size_t pos = 0;
      while (pos < group.size())
      {
        std::size_t end = group.find(';', pos);
        if (end == std::string_view::npos) { end = group.size(); }
        if (end > pos) { out.push_back(group.substr(pos, end - pos)); }
        pos = end + 1;
      }
      return out;
    }

    std::string strippedSequence(std::string_view modified)
    {
      std::string out;
      int depth = 0;
      for (const char c : modified)
      {
        if (c == '(' || c == '[') { ++depth; continue; }
        if (c == ')' || c == ']') { --depth; continue; }
        if (depth == 0 && std::isupper(static_cast<unsigned char>(c))) { out.push_back(c); }
      }
      return out;
    }

    // ------------------------------------------------------- DIA-NN 2 Parquet

    std::shared_ptr<arrow::Schema> diann2Schema(const std::vector<std::pair<std::string, std::string>>& metadata)
    {
      const auto s = arrow::utf8();
      const auto i = arrow::int64();
      const auto f = arrow::float32();
      std::vector<std::shared_ptr<arrow::Field>> fields;
      for (const auto& [name, type] : std::vector<std::pair<const char*, std::shared_ptr<arrow::DataType>>>{
             {"Precursor.Id", s}, {"Modified.Sequence", s}, {"Stripped.Sequence", s},
             {"Precursor.Charge", i}, {"Proteotypic", i}, {"Decoy", i}, {"N.Term", i}, {"C.Term", i},
             {"RT", f}, {"IM", f}, {"Q.Value", f}, {"Peptidoform.Q.Value", f}, {"PTM.Site.Confidence", f},
             {"PG.Q.Value", f}, {"Precursor.Mz", f}, {"Product.Mz", f}, {"Relative.Intensity", f},
             {"Fragment.Type", s}, {"Fragment.Charge", i}, {"Fragment.Series.Number", i},
             {"Fragment.Loss.Type", s}, {"Exclude.From.Quant", i},
             {"Protein.Ids", s}, {"Protein.Group", s}, {"Protein.Names", s}, {"Genes", s}, {"Flags", i}})
      { fields.push_back(arrow::field(name, type, /*nullable=*/false)); }
      std::shared_ptr<arrow::KeyValueMetadata> kv;
      if (!metadata.empty())
      {
        kv = std::make_shared<arrow::KeyValueMetadata>();
        for (const auto& [key, value] : metadata) { kv->Append(key, value); }
      }
      return arrow::schema(fields, kv);
    }

    /// One row group: the transitions of precursors [first, last).
    std::shared_ptr<arrow::Table> diann2RowGroup(const Library& library, const std::shared_ptr<arrow::Schema>& schema,
                                                 std::size_t first, std::size_t last, std::size_t& unknown_mods)
    {
      const auto& p = library.precursors();
      const auto& t = library.transitions();
      arrow::StringBuilder id, seq, stripped, ftype, loss, ids, group, names, genes;
      arrow::Int64Builder charge, proteotypic, decoy, n_term, c_term, fcharge, ordinal, exclude, flags;
      arrow::FloatBuilder rt, im, q, pq, ptm, pgq, pmz, qmz, intensity;
      std::size_t rows = 0;
      for (std::size_t i = first; i < last; ++i) { rows += p.transition_count[i]; }
      for (auto* b : {&charge, &proteotypic, &decoy, &n_term, &c_term, &fcharge, &ordinal, &exclude, &flags})
      { ok(b->Reserve(static_cast<std::int64_t>(rows))); }
      for (auto* b : {&rt, &im, &q, &pq, &ptm, &pgq, &pmz, &qmz, &intensity})
      { ok(b->Reserve(static_cast<std::int64_t>(rows))); }

      std::string precursor_id, accessions, entry_names;
      for (std::size_t i = first; i < last; ++i)
      {
        const auto stored = library.strings().get(p.modified_sequence[i]);
        const std::string modified = diannModifiedSequence(stored, &unknown_mods);
        const std::string plain = strippedSequence(modified);
        const int z = p.charge[i];
        // As in the TSV: DIA-NN gives a decoy its target's sequence and charge,
        // so the id alone would merge the pair on reading.
        precursor_id = modified + std::to_string(z) + (p.decoy[i] ? "_decoy" : "");
        accessions.clear(); entry_names.clear();
        const auto list = members(library.strings().get(p.protein_group[i]));
        for (std::size_t k = 0; k < list.size(); ++k)
        {
          const auto [accession, entry] = accessionAndName(list[k]);
          if (k) { accessions += ';'; entry_names += ';'; }
          accessions += accession; entry_names += entry;
        }
        const bool has_rt = !std::isnan(p.irt[i]);
        const bool has_im = i < p.im.size() && !std::isnan(p.im[i]);
        for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
        {
          const std::uint32_t j = p.transition_begin[i] + k;
          ok(id.Append(precursor_id)); ok(seq.Append(modified)); ok(stripped.Append(plain));
          charge.UnsafeAppend(z); proteotypic.UnsafeAppend(list.size() == 1 ? 1 : 0);
          decoy.UnsafeAppend(p.decoy[i] ? 1 : 0); n_term.UnsafeAppend(0); c_term.UnsafeAppend(0);
          rt.UnsafeAppend(has_rt ? p.irt[i] : 0.0f); im.UnsafeAppend(has_im ? p.im[i] : 0.0f);
          q.UnsafeAppend(0.0f); pq.UnsafeAppend(0.0f); ptm.UnsafeAppend(1.0f); pgq.UnsafeAppend(0.0f);
          pmz.UnsafeAppend(static_cast<float>(fromFixed(p.mz[i])));
          qmz.UnsafeAppend(static_cast<float>(fromFixed(t.product_mz[j])));
          intensity.UnsafeAppend(std::isnan(t.library_intensity[j]) ? 0.0f : t.library_intensity[j]);
          ok(ftype.Append(std::string(toString(t.type[j]))));
          fcharge.UnsafeAppend(t.charge[j]); ordinal.UnsafeAppend(t.ordinal[j]);
          ok(loss.Append(std::string(toString(t.loss[j]))));
          exclude.UnsafeAppend(0);
          ok(ids.Append(accessions)); ok(group.Append(accessions)); ok(names.Append(entry_names));
          ok(genes.Append("")); flags.UnsafeAppend(0);
        }
      }
      std::vector<std::shared_ptr<arrow::Array>> arrays;
      const auto finish = [&](arrow::ArrayBuilder& b) {
        std::shared_ptr<arrow::Array> a;
        ok(b.Finish(&a));
        arrays.push_back(std::move(a));
      };
      for (arrow::ArrayBuilder* b : std::initializer_list<arrow::ArrayBuilder*>{
             &id, &seq, &stripped, &charge, &proteotypic, &decoy, &n_term, &c_term,
             &rt, &im, &q, &pq, &ptm, &pgq, &pmz, &qmz, &intensity,
             &ftype, &fcharge, &ordinal, &loss, &exclude, &ids, &group, &names, &genes, &flags})
      { finish(*b); }
      return arrow::Table::Make(schema, arrays, static_cast<std::int64_t>(rows));
    }

    // ------------------------------------------------------------- OpenSWATH

    std::string uniModString(const OpenMS::AASequence& sequence) { return sequence.toUniModString(); }

    double lossMass(LossType loss)
    {
      switch (loss)
      {
        case LossType::Water: return 18.0105646837;
        case LossType::Ammonia: return 17.0265491015;
        case LossType::Phospho: return 97.9768955;
        case LossType::Metaphosphate: return 79.9663304;
        case LossType::CO: return 27.9949146;
        case LossType::None: case LossType::Other: break;
      }
      return 0.0;
    }

    OpenMS::Residue::ResidueType ionType(FragmentType type, LossType loss)
    {
      using R = OpenMS::Residue;
      if (loss == LossType::Water && type == FragmentType::B) { return R::BIonMinusH20; }
      if (loss == LossType::Water && type == FragmentType::Y) { return R::YIonMinusH20; }
      if (loss == LossType::Ammonia && type == FragmentType::B) { return R::BIonMinusNH3; }
      if (loss == LossType::Ammonia && type == FragmentType::Y) { return R::YIonMinusNH3; }
      switch (type)
      {
        case FragmentType::A: return R::AIon;
        case FragmentType::B: return R::BIon;
        case FragmentType::C: return R::CIon;
        case FragmentType::X: return R::XIon;
        case FragmentType::Y: return R::YIon;
        case FragmentType::Z: return R::ZIon;
        case FragmentType::Precursor: return R::Precursor;
        case FragmentType::Unknown: break;
      }
      return R::NonIdentified;
    }

    void addModification(std::vector<OpenMS::TargetedExperiment::Peptide::Modification>& mods, int location,
                         const OpenMS::ResidueModification& rmod)
    {
      OpenMS::TargetedExperiment::Peptide::Modification mod;
      mod.location = location;
      mod.mono_mass_delta = rmod.getDiffMonoMass();
      mod.avg_mass_delta = rmod.getDiffAverageMass();
      mod.unimod_id = rmod.getUniModRecordId();
      mods.push_back(mod);
    }
  } // namespace

  std::vector<std::string> LibraryExport::formatNames()
  {
    std::vector<std::string> out{"auto"};
    for (const auto& f : FORMATS) { out.emplace_back(f.name); }
    return out;
  }

  const char* LibraryExport::name(LibraryFormat format)
  {
    for (const auto& f : FORMATS) { if (f.format == format) { return f.name; } }
    return "?";
  }

  LibraryFormat LibraryExport::resolve(const std::string& out, const std::string& requested)
  {
    const std::string extension = extensionOf(out);
    if (requested.empty() || requested == "auto")
    {
      if (extension == ".parquet") { return LibraryFormat::DIALibGenParquet; }
      if (extension == ".tsv") { return LibraryFormat::DIANNTSV; }
      if (extension == ".pqp") { return LibraryFormat::OpenSwathPQP; }
      throw std::invalid_argument("-out must end in .parquet, .tsv or .pqp (got '" + out + "')");
    }
    for (const auto& f : FORMATS)
    {
      if (requested != f.name) { continue; }
      if (extension != f.extension)
      {
        throw std::invalid_argument("-out_format " + requested + " writes a " + f.extension +
                                    " file; -out is '" + out + "'");
      }
      return f.format;
    }
    std::string known;
    for (const auto& n : formatNames()) { known += (known.empty() ? "" : ", ") + n; }
    throw std::invalid_argument("unknown -out_format '" + requested + "' (known: " + known + ")");
  }

  bool LibraryExport::carriesRecipe(LibraryFormat format)
  {
    return format == LibraryFormat::DIALibGenParquet || format == LibraryFormat::DIANN2Parquet;
  }

  void LibraryExport::storeDIANN2Parquet(const std::string& filename, const Library& library,
                                         const std::vector<std::pair<std::string, std::string>>& metadata,
                                         unsigned threads)
  {
    trace::Scope write_step("write_diann2_parquet", 1);
    write_step.fields().rows = static_cast<std::int64_t>(library.transitionCount());
    const auto schema = diann2Schema(metadata);
    AtomicFile output(filename);
    auto opened = arrow::io::FileOutputStream::Open(output.temporaryPath().string());
    if (!opened.ok()) { throw std::runtime_error("cannot write library: " + filename); }
    const auto outfile = *opened;
    // ZSTD and per-row-group dictionaries for the strings each precursor
    // repeats on every fragment row, as DIALibGen's own layout; DIA-NN 2.0
    // reads the file through Arrow either way.
    const auto properties = parquet::WriterProperties::Builder()
                              .compression(parquet::Compression::ZSTD)
                              ->enable_dictionary()
                              ->build();
    auto writer_result = parquet::arrow::FileWriter::Open(*schema, arrow::default_memory_pool(), outfile, properties,
                                                          parquet::default_arrow_writer_properties());
    if (!writer_result.ok()) { throw std::runtime_error("cannot write Parquet: " + writer_result.status().ToString()); }
    auto writer = std::move(*writer_result);
    // Without store_schema (DIA-NN's own file has no ARROW:schema), the
    // schema's metadata reaches the file only this way.
    if (schema->metadata()) { ok(writer->AddKeyValueMetadata(schema->metadata())); }

    // Row groups of 65,536 precursors (~800k rows), built a wave at a time in
    // parallel and written in order: the file does not depend on threads.
    constexpr std::size_t GROUP = std::size_t(1) << 16;
    const std::size_t n = library.precursorCount();
    const std::size_t groups = (n + GROUP - 1) / GROUP;
    const unsigned workers = std::max(1u, std::min(parallel::workers(threads), 8u));
    std::atomic<std::size_t> unknown_mods{0};
    for (std::size_t wave = 0; wave < groups; wave += workers)
    {
      const std::size_t count = std::min<std::size_t>(workers, groups - wave);
      std::vector<std::shared_ptr<arrow::Table>> tables(count);
      parallel::forEach(count, workers, [&](std::size_t k)
      {
        const std::size_t first = (wave + k) * GROUP;
        std::size_t unknown = 0;
        tables[k] = diann2RowGroup(library, schema, first, std::min(n, first + GROUP), unknown);
        unknown_mods += unknown;
      });
      for (const auto& table : tables)
      {
        if (table->num_rows() == 0) { continue; }
        ok(writer->WriteTable(*table, table->num_rows()));
      }
    }
    if (n == 0)
    {
      // An empty library is still a readable file with the schema.
      std::size_t unknown = 0;
      const auto empty = diann2RowGroup(library, schema, 0, 0, unknown);
      ok(writer->WriteTable(*empty, 1));
    }
    ok(writer->Close());
    ok(outfile->Close());
    output.commit();
    if (unknown_mods.load())
    {
      std::cerr << "WARNING: " << unknown_mods.load() << " modifications have no UniMod accession and were written by "
                   "name; DIA-NN may skip those precursors\n";
    }
  }

  void LibraryExport::toTargetedExperiment(const Library& library, OpenMS::TargetedExperiment& experiment)
  {
    using OpenMS::TargetedExperiment;
    using OpenMS::ReactionMonitoringTransition;
    const auto& p = library.precursors();
    const auto& t = library.transitions();

    std::vector<TargetedExperiment::Peptide> peptides;
    std::vector<TargetedExperiment::Protein> proteins;
    std::vector<ReactionMonitoringTransition> transitions;
    peptides.reserve(library.precursorCount());
    transitions.reserve(library.transitionCount());
    std::unordered_map<std::string, std::size_t> protein_index;
    // Parsed once per distinct sequence: precursors of one peptide share it.
    std::unordered_map<std::uint32_t, OpenMS::AASequence> parsed;

    for (std::size_t i = 0; i < library.precursorCount(); ++i)
    {
      auto found = parsed.find(p.modified_sequence[i]);
      if (found == parsed.end())
      {
        const std::string stored(library.strings().get(p.modified_sequence[i]));
        OpenMS::AASequence sequence;
        try { sequence = OpenMS::AASequence::fromString(stored); }
        catch (const std::exception& e)
        { throw std::runtime_error("OpenSWATH export: cannot parse '" + stored + "': " + e.what()); }
        found = parsed.emplace(p.modified_sequence[i], std::move(sequence)).first;
      }
      const OpenMS::AASequence& sequence = found->second;
      const std::string unimod = uniModString(sequence);
      const int z = p.charge[i];
      const bool decoy = p.decoy[i] != 0;
      const std::string prefix = decoy ? "DECOY_" : "";
      // Unique per precursor; the index keeps two copies of one precursor apart.
      const std::string peptide_id = prefix + std::to_string(i) + "_" + unimod + "_" + std::to_string(z);

      TargetedExperiment::Peptide peptide;
      peptide.id = peptide_id;
      peptide.sequence = sequence.toUnmodifiedString();
      peptide.setMetaValue("full_peptide_name", unimod);
      peptide.setPeptideGroupLabel(prefix + unimod);
      peptide.setChargeState(z);
      if (!std::isnan(p.irt[i]))
      {
        TargetedExperiment::RetentionTime rt;
        rt.setRT(p.irt[i]);
        rt.retention_time_type = OpenMS::TargetedExperimentHelper::RetentionTime::RTType::IRT;
        peptide.rts.push_back(rt);
      }
      if (i < p.im.size() && !std::isnan(p.im[i])) { peptide.setDriftTime(p.im[i]); }
      if (sequence.hasNTerminalModification()) { addModification(peptide.mods, -1, *sequence.getNTerminalModification()); }
      if (sequence.hasCTerminalModification())
      { addModification(peptide.mods, static_cast<int>(sequence.size()), *sequence.getCTerminalModification()); }
      for (OpenMS::Size k = 0; k < sequence.size(); ++k)
      {
        if (sequence[k].isModified()) { addModification(peptide.mods, static_cast<int>(k), *sequence[k].getModification()); }
      }
      for (const auto member : members(library.strings().get(p.protein_group[i])))
      {
        const std::string protein_id(member);
        peptide.protein_refs.push_back(protein_id);
        if (protein_index.emplace(protein_id, proteins.size()).second)
        {
          TargetedExperiment::Protein protein;
          protein.id = protein_id;
          OpenMS::CVTerm accession;
          accession.setCVIdentifierRef("MS");
          accession.setAccession("MS:1000885");
          accession.setName("protein accession");
          accession.setValue(OpenMS::DataValue(std::string(accessionAndName(member).first)));
          protein.addCVTerm(accession);
          proteins.push_back(std::move(protein));
        }
      }

      const double precursor_mz = fromFixed(p.mz[i]);
      for (std::uint32_t k = 0; k < p.transition_count[i]; ++k)
      {
        const std::uint32_t j = p.transition_begin[i] + k;
        const int fragment_charge = t.charge[j];
        std::string annotation = std::string(toString(t.type[j])) + std::to_string(t.ordinal[j]);
        if (const double loss = lossMass(t.loss[j]); loss > 0.0)
        { annotation += "-" + std::to_string(static_cast<long>(std::lround(loss))); }
        annotation += "^" + std::to_string(fragment_charge);

        ReactionMonitoringTransition transition;
        transition.setNativeID(prefix + std::to_string(j) + "_" + annotation + "_" + std::to_string(i));
        transition.setPeptideRef(peptide_id);
        transition.setPrecursorMZ(precursor_mz);
        transition.setProductMZ(fromFixed(t.product_mz[j]));
        transition.setLibraryIntensity(std::isnan(t.library_intensity[j]) ? 0.0 : t.library_intensity[j]);
        auto product = transition.getProduct();
        product.setChargeState(fragment_charge);
        TargetedExperiment::Interpretation interpretation;
        interpretation.rank = 1;
        interpretation.ordinal = t.ordinal[j];
        interpretation.iontype = ionType(t.type[j], t.loss[j]);
        if (const double loss = lossMass(t.loss[j]); loss > 0.0)
        {
          OpenMS::CVTerm neutral_loss;
          neutral_loss.setCVIdentifierRef("MS");
          neutral_loss.setAccession("MS:1001524");
          neutral_loss.setName("fragment neutral loss");
          neutral_loss.setValue(-loss);
          interpretation.addCVTerm(neutral_loss);
        }
        product.addInterpretation(interpretation);
        transition.setProduct(product);
        transition.setDecoyTransitionType(decoy ? ReactionMonitoringTransition::DECOY
                                                : ReactionMonitoringTransition::TARGET);
        transition.setMetaValue("annotation", annotation);
        transition.setDetectingTransition(true);
        transition.setIdentifyingTransition(false);
        transition.setQuantifyingTransition(true);
        transitions.push_back(std::move(transition));
      }
      peptides.push_back(std::move(peptide));
    }
    experiment.setProteins(std::move(proteins));
    experiment.setPeptides(std::move(peptides));
    experiment.setTransitions(std::move(transitions));
  }

  void LibraryExport::storeOpenSwath(const std::string& filename, const Library& library, bool pqp)
  {
    if (const std::size_t decoys = library.decoyCount(); decoys > 0)
    {
      // OpenMS keys PEPTIDE (and through it PROTEIN) by sequence, and a
      // DIALibGen decoy carries its target's sequence: the writer then marks
      // that peptide and its proteins as decoys. Precursors and transitions
      // keep their own flags.
      std::cerr << "WARNING: " << decoys << " decoys share their targets' sequences, so OpenSWATH would mark those "
                   "targets as decoys; for OpenSWATH, generate with -generation:decoys none and run "
                   "OpenSwathDecoyGenerator\n";
    }
    trace::Scope build_step("build_targeted_experiment", 1);
    build_step.fields().rows = static_cast<std::int64_t>(library.transitionCount());
    OpenMS::TargetedExperiment experiment;
    toTargetedExperiment(library, experiment);
    build_step.end();
    trace::Scope write_step(pqp ? "write_pqp" : "write_openswath_tsv", 1);
    AtomicFile output(filename);
    const std::string staged = output.temporaryPath().string();
    if (pqp) { OpenMS::TransitionPQPFile().convertTargetedExperimentToPQP(staged.c_str(), experiment); }
    else { OpenMS::TransitionTSVFile().convertTargetedExperimentToTSV(staged.c_str(), experiment); }
    if (!std::filesystem::exists(output.temporaryPath()))
    { throw std::runtime_error("OpenMS wrote no file for " + filename); }
    output.commit();
  }

  void LibraryExport::store(const std::string& filename, const Library& library, LibraryFormat format,
                            const std::vector<std::pair<std::string, std::string>>& metadata, unsigned threads)
  {
    switch (format)
    {
      case LibraryFormat::DIALibGenParquet:
        DIANNLibraryFile::storeParquetCompact(filename, library, DIANNLibraryFile::Fingerprint{}, "", metadata, threads);
        return;
      case LibraryFormat::DIANNTSV: DIANNLibraryFile::storeTSV(filename, library, threads); return;
      case LibraryFormat::DIANN2Parquet: storeDIANN2Parquet(filename, library, metadata, threads); return;
      case LibraryFormat::OpenSwathPQP: storeOpenSwath(filename, library, true); return;
      case LibraryFormat::OpenSwathTSV: storeOpenSwath(filename, library, false); return;
    }
  }

} // namespace ODIA
