// A window.dialibgen that answers from memory, so App can be rendered without
// Tauri. The defaults below are a REAL `-write_config` dump, trimmed: the whole
// point of the form is that it is built from whatever the tool emits, so a
// mock that invents its own shape would test nothing.

import type { DialibgenApi, ModelStatus, Progress, RunResult } from '../types'

export const SAMPLE_CONFIG: Record<string, unknown> = {
  schema_version: 1,
  enzyme: 'Trypsin/P',
  missed_cleavages: 1,
  peptide_length: [7, 30],
  precursor_charges: [1, 2, 3, 4],
  precursor_mz: [350.0, 1200.0],
  fragment_mz: [200.0, 1800.0],
  max_fragment_charge: 2,
  fragments: [4, 12],
  fixed_modifications: ['Carbamidomethyl (C)'],
  variable_modifications: [],
  max_variable_modifications: 2,
  n_terminal_methionine_excision: true,
  free_cysteine_rt_correction: true,
  derive_ion_mobility: true,
  min_relative_intensity: 0.0001,
  reserved_doubly_charged: 3,
  decoys: 'none',
  rt_model: '',
  ms2_model: '',
  ccs_model: '',
  instrument: 'QE',
  nce: 30.0,
  irt_rescale: false,
  recompute_decoy_mz: false
}

export interface MockOptions {
  models?: ModelStatus
  config?: Record<string, unknown>
  configError?: string
  /// Fails the probe, as an uninstalled or unrunnable binary does.
  probeFails?: boolean
}

export interface MockBridge {
  api: DialibgenApi
  emitLog: (line: string) => void
  emitProgress: (p: Progress) => void
  emitDone: (r: RunResult) => void
  runs: unknown[]
  saved: Record<string, unknown>[]
}

export function installMockBridge(opts: MockOptions = {}): MockBridge {
  const logs: ((s: string) => void)[] = []
  const progs: ((p: Progress) => void)[] = []
  const dones: ((r: RunResult) => void)[] = []
  const runs: unknown[] = []
  const saved: Record<string, unknown>[] = []
  let presets: Record<string, Record<string, unknown>> = {}

  const api: DialibgenApi = {
    probe: async () =>
      opts.probeFails
        ? { bin: 'DIALibraryGenerator', source: 'path', ok: false, detail: 'cannot execute: not found' }
        : { bin: '/opt/bin/DIALibraryGenerator', source: 'bundled', ok: true, version: '0.2.0',
            detail: 'DIALibraryGenerator 0.2.0 (bundled)' },
    models: async () => opts.models ?? { dir: '/models', missing: [] },
    defaultConfig: async () => {
      if (opts.configError) throw new Error(opts.configError)
      return opts.config ?? { ...SAMPLE_CONFIG }
    },

    pickFasta: async () => '/data/proteins.fasta',
    pickOutput: async () => '/data/library.parquet',
    pickDirectory: async () => '/models',
    pickConfig: async () => '/data/config.json',
    saveConfig: async () => '/data/config.json',
    readConfig: async () => ({ enzyme: 'Lys-C' }),
    writeConfig: async () => true,
    revealPath: async () => true,

    run: async (params) => {
      runs.push(params)
      return { started: true, runId: runs.length }
    },
    cancel: async () => ({ cancelled: true }),

    loadSettings: async () => ({ schemaVersion: 1, lastUsed: null, presets }),
    saveLast: async (v) => {
      saved.push(v)
      return true
    },
    savePreset: async (name, v) => {
      presets = { ...presets, [name]: v }
      return true
    },
    deletePreset: async (name) => {
      const next = { ...presets }
      delete next[name]
      presets = next
      return true
    },

    onLog: (cb) => {
      logs.push(cb)
      return () => logs.splice(logs.indexOf(cb), 1)
    },
    onProgress: (cb) => {
      progs.push(cb)
      return () => progs.splice(progs.indexOf(cb), 1)
    },
    onDone: (cb) => {
      dones.push(cb)
      return () => dones.splice(dones.indexOf(cb), 1)
    }
  }

  window.dialibgen = api
  return {
    api,
    emitLog: (s) => logs.forEach((f) => f(s)),
    emitProgress: (p) => progs.forEach((f) => f(p)),
    emitDone: (r) => dones.forEach((f) => f(r)),
    runs,
    saved
  }
}
