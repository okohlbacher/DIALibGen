// Shared shapes between the Rust backend and the React frontend. These mirror
// the serde structs the Tauri commands return (see src-tauri/src/*.rs).

export interface BinaryInfo {
  bin: string
  source: 'env' | 'bundled' | 'path'
  ok: boolean
  version?: string
  detail: string
}

/// Which of the three PeptDeep ONNX exports are present. `missing` empty means
/// the tool can run; anything else is the reason it cannot.
export interface ModelStatus {
  dir: string | null
  missing: string[]
}

export interface RunStarted {
  started: boolean
  reason?: string
}

export interface RunResult {
  ok: boolean
  code: number | null
  /// Size of the library that was written, when one was. Absent on failure,
  /// and the difference between "exited 0" and "produced something".
  bytes?: number | null
}

export interface Settings {
  schemaVersion: number
  lastUsed: Record<string, unknown> | null
  presets: Record<string, Record<string, unknown>>
}

// Declared here (not in api.ts) so tests can install a mock bridge without
// importing api.ts, whose module body executes Tauri plumbing.
declare global {
  interface Window {
    dialibgen: DialibgenApi
  }
}

export interface DialibgenApi {
  probe: () => Promise<BinaryInfo>
  models: (dir?: string) => Promise<ModelStatus>
  /// The tool's own `-write_config` output: every default materialised. The
  /// form is built from this, so the GUI cannot offer a key the CLI lacks nor
  /// default one differently.
  defaultConfig: () => Promise<Record<string, unknown>>

  pickFasta: () => Promise<string | null>
  pickOutput: (defaultPath?: string) => Promise<string | null>
  pickDirectory: () => Promise<string | null>
  pickConfig: () => Promise<string | null>
  saveConfig: (defaultPath?: string) => Promise<string | null>
  readConfig: (path: string) => Promise<Record<string, unknown>>
  writeConfig: (path: string, value: Record<string, unknown>) => Promise<boolean>
  revealPath: (path: string) => Promise<boolean>

  run: (params: {
    in: string
    out: string
    config: Record<string, unknown>
    modelDir?: string | null
    threads?: number | null
  }) => Promise<RunStarted>
  cancel: () => Promise<{ cancelled: boolean }>

  loadSettings: () => Promise<Settings>
  saveLast: (values: Record<string, unknown>) => Promise<boolean>
  savePreset: (name: string, values: Record<string, unknown>) => Promise<boolean>
  deletePreset: (name: string) => Promise<boolean>

  onLog: (cb: (line: string) => void) => () => void
  onDone: (cb: (r: RunResult) => void) => () => void
}
