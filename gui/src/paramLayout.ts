// How the tool's own effective config becomes a form.
//
// The TOOL is the source of truth for which keys exist and what they default
// to: the backend runs `-write_config` and hands the JSON over, and every spec
// below is an OVERLAY on that -- prose, grouping, and the choice lists the JSON
// cannot express. Types are INFERRED from the tool's own default value, so a
// spec can never disagree with the CLI about whether something is a number.
//
// A key the tool grows and this file has not heard of still appears, inferred,
// under Advanced. That is deliberate: a new CLI parameter that the GUI silently
// hides is worse than one with a missing description, and it is what keeps the
// two from drifting without anybody noticing.

export type ParamKind =
  | 'bool'
  | 'int'
  | 'double'
  | 'string'
  | 'int-range'
  | 'double-range'
  | 'int-list'
  | 'string-list'
  | 'unsupported'

export interface ParamSpec {
  name: string
  kind: ParamKind
  /// Human label. Defaults to the key with underscores turned into spaces.
  label: string
  description: string
  group: 'core' | 'advanced'
  choices?: string[]
  min?: number
  max?: number
  /// Not shown in the form: managed elsewhere in the UI, or not a setting.
  hidden?: boolean
}

type Overlay = Omit<Partial<ParamSpec>, 'name' | 'kind'> & { description: string }

// Order matters: this is the order the core section is rendered in, and it runs
// roughly in pipeline order -- digest, then precursors, then prediction.
const OVERLAY: Record<string, Overlay> = {
  schema_version: {
    description: 'Config format version. Written for you; a config from a newer version is refused.',
    hidden: true
  },
  enzyme: {
    group: 'core',
    description:
      'Trypsin/P cuts after K and R regardless of a following proline, which is what DIA-NN does. ' +
      'Plain Trypsin cost 169,044 peptides on the human proteome.',
    choices: [
      'Trypsin/P',
      'Trypsin',
      'Lys-C',
      'Lys-C/P',
      'Lys-N',
      'Arg-C',
      'Chymotrypsin',
      'Glu-C',
      'Asp-N',
      'Asp-N/B',
      'PepsinA',
      'no cleavage',
      'unspecific cleavage'
    ]
  },
  missed_cleavages: { group: 'core', description: 'Missed cleavage sites allowed per peptide.', min: 0, max: 5 },
  peptide_length: { group: 'core', description: 'Peptide length window, in residues.', min: 1, max: 100 },
  precursor_charges: {
    group: 'core',
    description:
      'Charge states to enumerate. 2 and 3 alone cover 92.76% of observed precursors; a precursor the ' +
      'library cannot express is a ceiling on identifications, not a tuning knob.'
  },
  precursor_mz: { group: 'core', description: 'Precursor m/z window. Match it to the instrument method.' },
  fixed_modifications: {
    group: 'core',
    label: 'fixed modifications',
    description:
      'One OpenMS modification per line, e.g. "Carbamidomethyl (C)". EMPTY for a non-alkylated ' +
      'preparation. This is the single setting most likely to be wrong, and being wrong is expensive.'
  },
  variable_modifications: {
    group: 'core',
    description: 'One OpenMS modification per line, e.g. "Oxidation (M)". Each one multiplies the library.'
  },
  max_variable_modifications: {
    group: 'core',
    description: 'Most variable modifications on one peptide.',
    min: 0,
    max: 5
  },
  decoys: {
    group: 'core',
    description:
      '"none" is deliberate: a library is an interchange artefact and the consumer decides its own null. ' +
      'DIA-NN searches shipped decoys IN ADDITION to the ones it generates.',
    choices: ['none', 'mutate', 'pseudo_reverse', 'reverse', 'shuffle']
  },
  instrument: {
    group: 'core',
    description: 'Instrument the MS2 model conditions on.',
    choices: ['QE', 'Lumos', 'timsTOF', 'SciexTOF', 'Fusion', 'Eclipse', 'Velos', 'Elite', 'OrbitrapTribrid', 'ThermoTribrid']
  },
  nce: { group: 'core', description: 'Normalised collision energy the MS2 model conditions on.', min: 0, max: 100 },
  irt_rescale: {
    group: 'core',
    description:
      'OFF means the RT column is the model’s raw 0..1 output, NOT iRT, and is not interchangeable ' +
      'with another tool’s iRT library. Turn it on only to export.'
  },
  derive_ion_mobility: {
    group: 'core',
    description: 'Emit 1/K0 alongside CCS. Off costs a diaPASEF consumer the entire mobility dimension.'
  },

  fragment_mz: { group: 'advanced', description: 'Fragment m/z window.' },
  max_fragment_charge: { group: 'advanced', description: 'Highest fragment charge to enumerate.', min: 1, max: 3 },
  fragments: { group: 'advanced', description: 'Fragments kept per precursor, as [minimum, maximum].', min: 1, max: 50 },
  min_relative_intensity: {
    group: 'advanced',
    description: 'Drop fragments below this fraction of the base peak.',
    min: 0,
    max: 1
  },
  reserved_doubly_charged: {
    group: 'advanced',
    description: 'Slots reserved for doubly charged fragments, which the fragment selection otherwise under-picks.'
  },
  n_terminal_methionine_excision: {
    group: 'advanced',
    description: 'Also digest each protein with its initiator methionine removed.'
  },
  free_cysteine_rt_correction: {
    group: 'advanced',
    description:
      'Apply the measured RT offset for peptides with unalkylated cysteines. Self-gating: an alkylated ' +
      'library is identical whether this is on or off.'
  },
  recompute_decoy_mz: {
    group: 'advanced',
    description: 'Recompute decoy precursor and fragment m/z from the decoy sequence instead of inheriting the target’s.'
  },

  // The three model paths are set by the Models picker, which points the whole
  // trio at one directory. Leaving them in the form as well would give two
  // controls for one thing and let them disagree.
  rt_model: { description: 'Set by the model directory picker.', hidden: true },
  ms2_model: { description: 'Set by the model directory picker.', hidden: true },
  ccs_model: { description: 'Set by the model directory picker.', hidden: true }
}

const RANGE_KEYS = new Set(['peptide_length', 'precursor_mz', 'fragment_mz', 'fragments'])

function isInt(n: number): boolean {
  return Number.isInteger(n)
}

/// Infer the widget kind from the tool's own default value. A two-element
/// numeric array is a [min, max] range only for the keys the tool documents as
/// one -- everything else stays a list, so a future two-charge default does not
/// silently become a range.
export function inferKind(name: string, value: unknown): ParamKind {
  if (typeof value === 'boolean') return 'bool'
  if (typeof value === 'number') return isInt(value) ? 'int' : 'double'
  if (typeof value === 'string') return 'string'
  if (Array.isArray(value)) {
    const nums = value.filter((v) => typeof v === 'number') as number[]
    if (RANGE_KEYS.has(name) && value.length === 2 && nums.length === 2) {
      return nums.every(isInt) ? 'int-range' : 'double-range'
    }
    if (value.length > 0 && nums.length === value.length) return 'int-list'
    if (value.every((v) => typeof v === 'string')) return 'string-list'
    if (value.length === 0) return 'string-list' // an empty list is a list of strings until told otherwise
  }
  // null, objects, mixed arrays: shown read-only rather than dropped, so a key
  // the form cannot edit is still visible instead of vanishing.
  return 'unsupported'
}

/// Build the full spec list from the tool's effective config.
export function buildSpecs(defaults: Record<string, unknown>): ParamSpec[] {
  const order = Object.keys(OVERLAY)
  const keys = Object.keys(defaults).sort((a, b) => {
    const ia = order.indexOf(a)
    const ib = order.indexOf(b)
    if (ia === -1 && ib === -1) return a.localeCompare(b)
    if (ia === -1) return 1
    if (ib === -1) return -1
    return ia - ib
  })
  return keys.map((name) => {
    const o = OVERLAY[name]
    const kind = inferKind(name, defaults[name])
    return {
      name,
      kind,
      label: o?.label ?? name.replace(/_/g, ' '),
      // An undescribed key is a key this file has not caught up with. Say so
      // rather than render an empty hint: the blank looks like a finished field.
      description: o?.description ?? 'Not yet described here — see `DIALibraryGenerator -write_config`.',
      group: o?.group ?? 'advanced',
      choices: o?.choices,
      min: o?.min,
      max: o?.max,
      hidden: o?.hidden ?? false
    }
  })
}

/// Why @p name currently does nothing, or null when it is live. The field stays
/// visible -- hiding it would make the option undiscoverable and make the form
/// jump around as switches flip -- but it is disabled and says why.
export function inertBecause(name: string, values: Record<string, unknown>): string | null {
  if (name === 'max_variable_modifications') {
    const v = values.variable_modifications
    if (Array.isArray(v) && v.length === 0) {
      return 'No effect: no variable modifications are configured.'
    }
  }
  if (name === 'recompute_decoy_mz' && values.decoys === 'none') {
    return 'No effect: no decoys are being generated.'
  }
  return null
}
