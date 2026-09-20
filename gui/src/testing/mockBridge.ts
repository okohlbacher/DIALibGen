// In-memory bridge using the native 0.11.0 `-mode generate -write_config` dump.
// Only the three machine-specific model paths are normalized to empty strings
// in generation-defaults.json. Native contract tests check this fixture for drift.

import type { DialibgenApi, ModelStatus, RunResult, TuningOption } from '../types'
import generationDefaults from './generation-defaults.json'
import refinementDefaults from './refinement-defaults.json'
import tuningOptions from './tuning-options.json'

export const SAMPLE_CONFIG: Record<string, unknown> = generationDefaults

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
  emitDone: (r: RunResult) => void
  runs: unknown[]
  saved: Record<string, unknown>[]
}

export function installMockBridge(opts: MockOptions = {}): MockBridge {
  const logs: ((s: string) => void)[] = []
  const dones: ((r: RunResult) => void)[] = []
  const runs: unknown[] = []
  const saved: Record<string, unknown>[] = []
  let presets: Record<string, Record<string, unknown>> = {}

  const api: DialibgenApi = {
    probe: async () =>
      opts.probeFails
        ? { bin: 'DIALibGen', source: 'path', ok: false, detail: 'cannot execute: not found' }
        : { bin: '/opt/bin/DIALibGen', source: 'bundled', ok: true, version: '0.11.0',
            detail: 'DIALibGen 0.11.0 (bundled)' },
    models: async () => opts.models ?? { dir: '/models', missing: [] },
    defaultConfig: async (mode = 'generate') => {
      if (opts.configError) throw new Error(opts.configError)
      return opts.config ?? (mode === 'generate' ? { ...SAMPLE_CONFIG } : { ...refinementDefaults, ...(mode === 'tune' ? { filter: false, write_rt: false } : {}) })
    },

    tuningOptions: async () => tuningOptions as TuningOption[],
    pickLibrary: async () => '/data/library.parquet',
    pickReport: async () => '/data/report.parquet',
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
      return { started: true }
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
    onDone: (cb) => {
      dones.push(cb)
      return () => dones.splice(dones.indexOf(cb), 1)
    }
  }

  window.dialibgen = api
  return {
    api,
    emitLog: (s) => logs.forEach((f) => f(s)),
    emitDone: (r) => dones.forEach((f) => f(r)),
    runs,
    saved
  }
}
