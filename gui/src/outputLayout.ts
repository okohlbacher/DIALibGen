// The output layouts the CLI writes (-out_format), each with the extension
// -out must carry. `auto` (no -out_format) is what the extension alone picks:
// .parquet DIALibGen's layout, .tsv DIA-NN's TSV, .pqp OpenSWATH's PQP. Only
// the two layouts that share an extension with an automatic one are sent
// explicitly, so the default command line is unchanged.

export type Layout = 'dialibgen_parquet' | 'diann_tsv' | 'diann2_parquet' | 'openswath_pqp' | 'openswath_tsv'
export type Extension = 'parquet' | 'tsv' | 'pqp'

export interface LayoutInfo { id: Layout; label: string; ext: Extension; help: string }

export const LAYOUTS: readonly LayoutInfo[] = [
  { id: 'dialibgen_parquet', label: 'DIALibGen Parquet', ext: 'parquet',
    help: 'Carries the config verbatim in its schema metadata, so the library states how it was built. DIA-NN does not read it.' },
  { id: 'diann_tsv', label: 'DIA-NN TSV', ext: 'tsv',
    help: 'DIA-NN dialect, read by DIA-NN 1.7 and 2.x. A TSV cannot carry the recipe, so this library is not reproducible from itself.' },
  { id: 'diann2_parquet', label: 'DIA-NN 2 Parquet', ext: 'parquet',
    help: 'DIA-NN 2.0’s own Parquet library layout, which DIA-NN 2.0 loads. Carries the config in its metadata.' },
  { id: 'openswath_pqp', label: 'OpenSWATH PQP', ext: 'pqp',
    help: 'OpenSWATH’s PQP, written by OpenMS. Run OpenSwathDecoyGenerator on it when the library has no decoys.' },
  { id: 'openswath_tsv', label: 'OpenSWATH TSV', ext: 'tsv',
    help: 'OpenSWATH’s assay TSV, as TargetedFileConverter reads it. Run OpenSwathDecoyGenerator when the library has no decoys.' }
]

export const DEFAULT_LAYOUT: Layout = 'dialibgen_parquet'

export function layoutInfo(id: Layout): LayoutInfo {
  return LAYOUTS.find((l) => l.id === id) ?? LAYOUTS[0]
}

export function isLayout(value: unknown): value is Layout {
  return LAYOUTS.some((l) => l.id === value)
}

/// The extension of @p path, spelled exactly as the CLI requires, or null.
export function extensionOf(path: string): Extension | null {
  const m = /\.(parquet|tsv|pqp)$/.exec(path)
  return m ? m[1] as Extension : null
}

/// What the CLI picks from the extension alone.
export function autoLayout(ext: Extension): Layout {
  return ext === 'parquet' ? 'dialibgen_parquet' : ext === 'tsv' ? 'diann_tsv' : 'openswath_pqp'
}

/// The layout a run of @p path writes: the selected one where its extension
/// matches, else what the extension picks.
export function effectiveLayout(selected: Layout, path: string): Layout | null {
  const ext = extensionOf(path)
  if (!ext) return null
  return layoutInfo(selected).ext === ext ? selected : autoLayout(ext)
}

/// The -out_format value to send, or undefined when the extension says it.
export function outFormatFor(layout: Layout): Layout | undefined {
  const ext = layoutInfo(layout).ext
  return autoLayout(ext) === layout ? undefined : layout
}

export function withExtension(path: string, ext: Extension): string {
  return path.replace(/\.(parquet|tsv|pqp)$/i, '') + '.' + ext
}
