import { useCallback, useEffect, useMemo, useRef, useState, type JSX } from 'react'
import ParamField from './ParamField'
import { buildSpecs, inertBecause, type ParamSpec } from './paramLayout'
import type { BinaryInfo, ModelStatus, RunResult } from './types'

type Values = Record<string, unknown>

interface Outcome {
  ok: boolean
  text: string
  path?: string
}

const MAX_LOG_LINES = 2000

function humanBytes(n: number): string {
  const u = ['B', 'kB', 'MB', 'GB', 'TB']
  let i = 0
  let v = n
  while (v >= 1000 && i < u.length - 1) {
    v /= 1000
    i += 1
  }
  return `${v < 10 && i > 0 ? v.toFixed(1) : Math.round(v)} ${u[i]}`
}

/// Swap a library's extension while keeping the rest of the path, so the
/// format toggle does not make the user re-pick the file.
export function withExtension(path: string, ext: 'parquet' | 'tsv'): string {
  return path.replace(/\.(parquet|tsv)$/i, '') + '.' + ext
}

export default function App(): JSX.Element {
  const [info, setInfo] = useState<BinaryInfo | null>(null)
  const [models, setModels] = useState<ModelStatus | null>(null)
  const [modelDir, setModelDir] = useState<string>('')
  const [defaults, setDefaults] = useState<Values | null>(null)
  const [configError, setConfigError] = useState<string | null>(null)
  const [values, setValues] = useState<Values>({})
  const [fasta, setFasta] = useState('')
  const [out, setOut] = useState('')
  const [threads, setThreads] = useState(1)
  const [showAdvanced, setShowAdvanced] = useState(false)
  const [log, setLog] = useState<string[]>([])
  const [running, setRunning] = useState(false)
  const [outcome, setOutcome] = useState<Outcome | null>(null)
  const [presets, setPresets] = useState<Record<string, Values>>({})
  const logRef = useRef<HTMLPreElement | null>(null)
  // The done handler is installed once and so cannot close over `out`.
  const outRef = useRef(out)

  const specs: ParamSpec[] = useMemo(() => (defaults ? buildSpecs(defaults) : []), [defaults])
  const visible = specs.filter((s) => !s.hidden)
  const core = visible.filter((s) => s.group === 'core')
  const advanced = visible.filter((s) => s.group === 'advanced')

  // ---- startup: probe the binary, then ask IT for the parameter set --------
  useEffect(() => {
    let alive = true
    void (async () => {
      const i = await window.dialibgen.probe()
      if (!alive) return
      setInfo(i)
      try {
        const d = await window.dialibgen.defaultConfig()
        if (!alive) return
        setDefaults(d)
        setValues({ ...d, nce: null })
        setConfigError(null)
      } catch (e) {
        if (alive) setConfigError(String(e))
      }
      const s = await window.dialibgen.loadSettings()
      if (!alive) return
      setPresets(s.presets ?? {})
      const last = s.lastUsed as (Values & { __modelDir?: string }) | null
      if (last) {
        if (typeof last.__modelDir === 'string') setModelDir(last.__modelDir)
        setValues((v) => ({ ...v, ...stripPrivate(last) }))
      }
    })()
    return () => {
      alive = false
    }
  }, [])

  // ---- model availability, re-checked whenever the directory changes ------
  useEffect(() => {
    let alive = true
    void window.dialibgen.models(modelDir || undefined).then((m) => {
      if (alive) setModels(m)
    })
    return () => {
      alive = false
    }
  }, [modelDir])

  // ---- run event streams --------------------------------------------------
  useEffect(() => {
    const offLog = window.dialibgen.onLog((line) => {
      setLog((l) => (l.length >= MAX_LOG_LINES ? [...l.slice(1), line] : [...l, line]))
    })
    const offDone = window.dialibgen.onDone((r: RunResult) => {
      setRunning(false)
      if (r.ok && r.bytes) {
        setOutcome({ ok: true, text: `Wrote ${humanBytes(r.bytes)}`, path: outRef.current })
      } else if (r.ok) {
        // Exit 0 and nothing on disk is not success: it means the path the tool
        // was given is not the path anybody can read back.
        setOutcome({ ok: false, text: 'The tool exited cleanly but no library is on disk.' })
      } else {
        setOutcome({ ok: false, text: `Failed (exit ${r.code ?? 'killed'}). See the log.` })
      }
    })
    return () => {
      offLog()
      offDone()
    }
  }, [])

  useEffect(() => {
    if (logRef.current) logRef.current.scrollTop = logRef.current.scrollHeight
  }, [log])

  const set = useCallback((name: string, v: unknown) => {
    setValues((old) => ({ ...old, [name]: v }))
  }, [])

  const modelsReady = !!models && models.missing.length === 0
  const outExt = /\.parquet$/.test(out) ? 'parquet' : /\.tsv$/.test(out) ? 'tsv' : null
  const threadsValid = Number.isSafeInteger(threads) && threads >= 0
  const canRun = threadsValid && !running && !!fasta && !!out && !!outExt && modelsReady && !!defaults && !!info?.ok

  async function start(): Promise<void> {
    setRunning(true)
    outRef.current = out
    setLog([])
    setOutcome(null)
    try {
      if (!await window.dialibgen.saveLast({ ...values, __modelDir: modelDir })) {
        throw new Error('could not save settings; check that the app config directory is writable')
      }
      const r = await window.dialibgen.run({
        in: fasta,
        out,
        config: toolConfig(values),
        modelDir: modelDir || null,
        threads
      })
      if (!r.started) {
        setRunning(false)
        setOutcome({ ok: false, text: r.reason ?? 'could not start' })
      }
    } catch (e) {
      setRunning(false)
      setOutcome({ ok: false, text: `Could not start: ${String(e)}` })
    }
  }

  return (
    <div className="app">
      <header>
        <h1>DIALibGen</h1>
        <span className={`badge ${info?.ok ? 'ok' : 'bad'}`}>
          {info ? info.detail : 'probing…'}
        </span>
        <span className={`badge ${modelsReady ? 'ok' : 'bad'}`}>
          {models === null
            ? 'checking models…'
            : modelsReady
              ? 'models ready'
              : `${models.missing.length} model${models.missing.length === 1 ? '' : 's'} missing`}
        </span>
      </header>

      <main>
        <section className="controls">
          <div className="field">
            <label htmlFor="fasta">protein FASTA</label>
            <div className="row tight">
              <input id="fasta" type="text" value={fasta} placeholder="(not set)"
                     onChange={(e) => setFasta(e.target.value)} />
              <button type="button" className="secondary slim"
                      onClick={async () => {
                        const p = await window.dialibgen.pickFasta()
                        if (p) setFasta(p)
                      }}>
                Browse…
              </button>
            </div>
          </div>

          <div className="field">
            <label htmlFor="out">output library</label>
            <div className="row tight">
              <input id="out" type="text" value={out} placeholder="(not set)"
                     onChange={(e) => setOut(e.target.value)} />
              <button type="button" className="secondary slim"
                      onClick={async () => {
                        const p = await window.dialibgen.pickOutput(out || undefined)
                        if (p) setOut(p)
                      }}>
                Browse…
              </button>
            </div>
            <div className="row tight formats">
              <button type="button" className={`chip${outExt === 'parquet' ? ' on' : ''}`}
                      disabled={!out} onClick={() => setOut(withExtension(out, 'parquet'))}>
                .parquet
              </button>
              <button type="button" className={`chip${outExt === 'tsv' ? ' on' : ''}`}
                      disabled={!out} onClick={() => setOut(withExtension(out, 'tsv'))}>
                .tsv
              </button>
            </div>
            <p className="help">
              {outExt === 'tsv'
                ? 'DIA-NN dialect. A TSV cannot carry the recipe, so this library is not reproducible from itself.'
                : outExt === 'parquet'
                  ? 'Carries the config verbatim in its schema metadata, so the library states how it was built.'
                  : 'The name must end in .parquet or .tsv.'}
            </p>
          </div>

          <div className="field">
            <label htmlFor="models">model directory</label>
            <div className="row tight">
              <input id="models" type="text" value={modelDir} placeholder="(search the usual places)"
                     onChange={(e) => setModelDir(e.target.value)} />
              <button type="button" className="secondary slim"
                      onClick={async () => {
                        const p = await window.dialibgen.pickDirectory()
                        if (p) setModelDir(p)
                      }}>
                Browse…
              </button>
            </div>
            {models && models.missing.length > 0 ? (
              <p className="help bad">
                Not found: {models.missing.join(', ')}. Release builds include these models.
                For a custom installation, select a directory containing all three .onnx files.
              </p>
            ) : (
              <p className="help">{models?.dir ?? ''}</p>
            )}
          </div>

          <div className="field">
            <label htmlFor="threads">threads</label>
            <input id="threads" type="number" min={0} step={1} value={threads} aria-invalid={!threadsValid}
                   onChange={(e) => setThreads(Number(e.target.value) || 0)} />
            <p className={`help${threadsValid ? '' : ' bad'}`}>{threadsValid
              ? '1 is the default. 0 uses all available cores, with more memory.'
              : 'Enter a whole number of threads (0 or more).'}</p>
          </div>

          <div className="row">
            <button type="button" disabled={!canRun} onClick={() => void start()}>
              {running ? 'Running…' : 'Generate library'}
            </button>
            <button type="button" className="secondary" disabled={!running}
                    onClick={() => void window.dialibgen.cancel().catch((e: unknown) =>
                      setOutcome({ ok: false, text: `Could not cancel: ${String(e)}` }))}>
              Cancel
            </button>
          </div>

          <Presets
            presets={presets}
            values={values}
            onApply={(v) => setValues((old) => ({ ...old, ...v }))}
            onSaved={setPresets}
            onError={(text) => setOutcome({ ok: false, text })}
          />

          {outcome && (
            <div className={`outcome ${outcome.ok ? 'ok' : 'bad'}`}>
              <span>{outcome.text}</span>
              {outcome.ok && outcome.path && (
                <button type="button" className="secondary slim"
                        onClick={() => void window.dialibgen.revealPath(outcome.path as string)}>
                  Show in folder
                </button>
              )}
            </div>
          )}
        </section>

        <section className="right">
          <div className="params">
            <div className="section-head">
              <h2>Library configuration</h2>
              <div className="row tight actions">
                <button type="button" className="secondary slim"
                        disabled={!defaults} onClick={() => defaults && setValues({ ...defaults, nce: null })}>
                  Reset
                </button>
                <button type="button" className="secondary slim" onClick={() => void loadConfig(setValues).catch((e: unknown) =>
                  setOutcome({ ok: false, text: `Could not load config: ${String(e)}` }))}>
                  Load…
                </button>
                <button type="button" className="secondary slim"
                        onClick={() => void saveConfig(toolConfig(values)).catch((e: unknown) =>
                          setOutcome({ ok: false, text: `Could not save config: ${String(e)}` }))}>
                  Save…
                </button>
              </div>
            </div>

            {configError && (
              <p className="help bad">
                Could not read the tool’s own defaults: {configError}. The form is unavailable until the
                binary runs — check the path in the badge above.
              </p>
            )}

            {core.map((s) => (
              <ParamField key={s.name} spec={s} value={values[s.name]}
                          inert={inertBecause(s.name, values)}
                          onChange={(v) => set(s.name, v)} />
            ))}

            {advanced.length > 0 && (
              <>
                <button type="button" className="disclose"
                        aria-expanded={showAdvanced}
                        onClick={() => setShowAdvanced((x) => !x)}>
                  {showAdvanced ? '▾' : '▸'} Advanced ({advanced.length})
                </button>
                {showAdvanced &&
                  advanced.map((s) => (
                    <ParamField key={s.name} spec={s} value={values[s.name]}
                                inert={inertBecause(s.name, values)}
                                onChange={(v) => set(s.name, v)} />
                  ))}
              </>
            )}
          </div>

          <div className="logpane">
            <div className="section-head">
              <h2>Log</h2>
            </div>
            <pre ref={logRef} className="log">{log.join('\n')}</pre>
          </div>
        </section>
      </main>
    </div>
  )
}

/// Keys the GUI keeps in the same object for convenience but the tool must
/// never see: it refuses an unknown config key outright.
export function stripPrivate(v: Values): Values {
  const out: Values = {}
  for (const [k, val] of Object.entries(v)) {
    if (!k.startsWith('__')) out[k] = val
  }
  return out
}

// A blank NCE lets the CLI resolve the selected instrument's default. Presets
// retain null to preserve this choice; JSON sent to the CLI must omit the key.
function toolConfig(v: Values): Values {
  const config = stripPrivate(v)
  if (config.nce == null) delete config.nce
  return config
}

async function loadConfig(setValues: (f: (v: Values) => Values) => void): Promise<void> {
  const p = await window.dialibgen.pickConfig()
  if (!p) return
  const parsed = (await window.dialibgen.readConfig(p)) as Values
  // A partial config preserves the form's other values; a read failure leaves
  // the form unchanged and is reported by the caller.
  setValues((v) => ({ ...v, ...parsed, nce: parsed.nce ?? null }))
}

async function saveConfig(values: Values): Promise<void> {
  const p = await window.dialibgen.saveConfig('library-config.json')
  if (!p) return
  await window.dialibgen.writeConfig(p, values)
}

function Presets({
  presets,
  values,
  onApply,
  onSaved,
  onError
}: {
  presets: Record<string, Values>
  values: Values
  onApply: (v: Values) => void
  onSaved: (p: Record<string, Values>) => void
  onError: (message: string) => void
}): JSX.Element {
  const [name, setName] = useState('')
  const [selected, setSelected] = useState('')
  const names = Object.keys(presets).sort()
  return (
    <div className="field presets">
      <label htmlFor="preset">presets</label>
      <div className="row tight">
        <select id="preset" value={selected} onChange={(e) => {
          const key = e.target.value
          setSelected(key)
          const p = presets[key]
          if (p) onApply(p)
        }}>
          <option value="">{names.length ? 'apply a preset…' : '(none saved)'}</option>
          {names.map((n) => (
            <option key={n} value={n}>{n}</option>
          ))}
        </select>
        <button type="button" className="secondary slim" disabled={!selected}
                onClick={() => void (async () => {
                  if (!await window.dialibgen.deletePreset(selected)) throw new Error('settings could not be written')
                  const next = { ...presets }
                  delete next[selected]
                  onSaved(next)
                  setSelected('')
                })().catch((e: unknown) => onError(`Could not delete preset: ${String(e)}`))}>
          Delete
        </button>
      </div>
      <div className="row tight">
        <input type="text" value={name} placeholder="save current as…"
               aria-label="preset name" onChange={(e) => setName(e.target.value)} />
        <button type="button" className="secondary slim" disabled={!name.trim()}
                onClick={() => void (async () => {
                  const clean = name.trim()
                  const payload = stripPrivate(values)
                  if (!await window.dialibgen.savePreset(clean, payload)) throw new Error('settings could not be written')
                  onSaved({ ...presets, [clean]: payload })
                  setName('')
                })().catch((e: unknown) => onError(`Could not save preset: ${String(e)}`))}>
          Save
        </button>
      </div>
    </div>
  )
}
