// The single bridge between the React frontend and the Rust backend. App.tsx
// speaks only to `window.dialibgen`, so the tests can install a mock without
// touching Tauri and a future host (a web build, say) needs no change above
// this file.

import { invoke } from '@tauri-apps/api/core'
import { listen, type UnlistenFn } from '@tauri-apps/api/event'
import { open, save } from '@tauri-apps/plugin-dialog'
import { revealItemInDir } from '@tauri-apps/plugin-opener'
import type {
  BinaryInfo,
  DialibgenApi,
  ModelStatus,
  RunResult,
  RunStarted,
  Settings
} from './types'

const FASTA = [{ name: 'Protein FASTA', extensions: ['fasta', 'fa', 'faa'] }]
const LIBRARY = [
  { name: 'Spectral library', extensions: ['parquet', 'tsv'] }
]
const JSON_FILTER = [{ name: 'Library config (JSON)', extensions: ['json'] }]

// listen() is async but the renderer wants a synchronous unsubscribe (it
// returns one from a useEffect). Hand back a function that unlistens once the
// async subscription resolves, and no-ops safely if it fires before then.
function subscribe<T>(event: string, cb: (payload: T) => void): () => void {
  let un: UnlistenFn | null = null
  let cancelled = false
  listen<T>(event, (e) => cb(e.payload)).then((u) => {
    if (cancelled) u()
    else un = u
  })
  return () => {
    cancelled = true
    if (un) un()
  }
}

const api: DialibgenApi = {
  probe: () => invoke<BinaryInfo>('probe'),
  models: (dir?: string) => invoke<ModelStatus>('models', { dir: dir ?? null }),
  defaultConfig: () => invoke<Record<string, unknown>>('default_config'),

  // Native dialogs live in the frontend via the dialog plugin; the backend only
  // ever receives already-chosen paths.
  pickFasta: async () => {
    const r = await open({ multiple: false, filters: FASTA })
    return typeof r === 'string' ? r : null
  },
  pickOutput: async (defaultPath?: string) => {
    const r = await save({ defaultPath, filters: LIBRARY })
    return typeof r === 'string' ? r : null
  },
  pickDirectory: async () => {
    const r = await open({ directory: true, multiple: false })
    return typeof r === 'string' ? r : null
  },
  pickConfig: async () => {
    const r = await open({ multiple: false, filters: JSON_FILTER })
    return typeof r === 'string' ? r : null
  },
  saveConfig: async (defaultPath?: string) => {
    const r = await save({ defaultPath, filters: JSON_FILTER })
    return typeof r === 'string' ? r : null
  },
  readConfig: (path: string) => invoke<Record<string, unknown>>('read_config', { path }),
  writeConfig: (path: string, value: Record<string, unknown>) =>
    invoke<boolean>('write_config', { path, value }),
  // Reveals in the OS file manager rather than opening: a 3 GB Parquet handed
  // to whatever claims the extension is not a kindness.
  revealPath: async (path: string) => {
    await revealItemInDir(path)
    return true
  },

  run: (params) => invoke<RunStarted>('run', { params }),
  cancel: () => invoke<{ cancelled: boolean }>('cancel'),

  loadSettings: () => invoke<Settings>('load_settings'),
  saveLast: (values) => invoke<boolean>('save_last', { values }),
  savePreset: (name, values) => invoke<boolean>('save_preset', { name, values }),
  deletePreset: (name) => invoke<boolean>('delete_preset', { name }),

  onLog: (cb) => subscribe<string>('dialibgen:log', cb),
  onDone: (cb) => subscribe<RunResult>('dialibgen:done', cb)
}

window.dialibgen = api
