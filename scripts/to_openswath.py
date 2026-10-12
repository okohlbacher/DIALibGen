#!/usr/bin/env python3
"""DIALibGen TSV -> OpenSWATH transition TSV (readable by OpenMS TargetedFileConverter).   usage: to_openswath.py library.tsv openswath.tsv"""
import sys, pandas as pd
t = pd.read_csv(sys.argv[1], sep="\t")
pd.DataFrame({
    "PrecursorMz": t["Precursor.Mz"], "ProductMz": t["Product.Mz"], "LibraryIntensity": t["Relative.Intensity"],
    "NormalizedRetentionTime": t["RT"], "PrecursorIonMobility": t["IM"],
    "PeptideSequence": t["Modified.Sequence"].str.replace(r"\(.*?\)", "", regex=True).str.strip("."),   # '.(Acetyl)PEPTIDE' -> 'PEPTIDE'
    "ModifiedPeptideSequence": t["Modified.Sequence"], "PrecursorCharge": t["Precursor.Charge"],
    "ProteinId": t["Protein.Group"], "TransitionGroupId": t["Precursor.Id"],
    "FragmentType": t["Fragment.Type"], "FragmentCharge": t["Fragment.Charge"],
    "FragmentSeriesNumber": t["Fragment.Series.Number"], "Decoy": t["Decoy"],
}).to_csv(sys.argv[2], sep="\t", index=False)
