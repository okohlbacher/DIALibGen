import { useCallback, useEffect, useMemo, useRef, useState, type JSX } from 'react'
import ParamField from './ParamField'
import { buildSpecs, inertBecause, invalidValue, refinementSpecs, tuningError, tuningSpecs, type ParamSpec } from './paramLayout'
import type { BinaryInfo, Mode, ModelHead, ModelStatus, RunResult, TuningOption } from './types'

type Values = Record<string, unknown>
type Configs = Partial<Record<Mode, Values>>
interface Outcome { ok: boolean; text: string; path?: string }
const MAX_LOG_LINES = 2000
const modes: Mode[] = ['generate', 'refine', 'tune']
const record = (v: unknown): v is Values => !!v && typeof v === 'object' && !Array.isArray(v)

function humanBytes(n: number): string {
  const units = ['B', 'kB', 'MB', 'GB', 'TB']
  let i = 0
  while (n >= 1000 && i < units.length - 1) { n /= 1000; i += 1 }
  return `${n < 10 && i > 0 ? n.toFixed(1) : Math.round(n)} ${units[i]}`
}
export function withExtension(path: string, ext: 'parquet' | 'tsv'): string {
  return path.replace(/\.(parquet|tsv)$/i, '') + '.' + ext
}
export function stripPrivate(v: Values): Values {
  return Object.fromEntries(Object.entries(v).filter(([key]) => !key.startsWith('__')))
}
function toolConfig(v: Values): Values {
  const config = stripPrivate(v)
  if (config.nce == null) delete config.nce
  return config
}

export default function App(): JSX.Element {
  const [info, setInfo] = useState<BinaryInfo | null>(null)
  const [initialized, setInitialized] = useState(false)
  const [mode, setMode] = useState<Mode>('generate')
  const [nativeDefaults, setNativeDefaults] = useState<Configs>({})
  const [configs, setConfigs] = useState<Configs>({})
  const [configError, setConfigError] = useState<string | null>(null)
  const [options, setOptions] = useState<TuningOption[] | null>(null)
  const [tuning, setTuning] = useState<Values>({})
  const [tune, setTune] = useState(false)
  const [tuningFailure, setTuningFailure] = useState<string | null>(null)
  const [fasta, setFasta] = useState('')
  const [library, setLibrary] = useState('')
  const [ids, setIds] = useState('')
  const [out, setOut] = useState('')
  const [outReport, setOutReport] = useState('')
  const [tuneOutModels, setTuneOutModels] = useState('')
  const [modelDir, setModelDir] = useState('')
  const [modelCheck, setModelCheck] = useState<{ key: string; status: ModelStatus } | null>(null)
  const [threads, setThreads] = useState<number | null>(1)
  const [showAdvanced, setShowAdvanced] = useState(false)
  const [showTrainingAdvanced, setShowTrainingAdvanced] = useState(false)
  const [log, setLog] = useState<string[]>([])
  const [running, setRunning] = useState(false)
  const [outcome, setOutcome] = useState<Outcome | null>(null)
  const [presets, setPresets] = useState<Record<string, Values>>({})
  const logRef = useRef<HTMLPreElement | null>(null)
  const outRef = useRef(out)
  // Invalidates a file dialog/default response when a newer mode/recipe wins.
  const revision = useRef(0)
  const defaults = nativeDefaults[mode] ?? null
  const values = configs[mode] ?? {}
  const usesTuning = mode === 'tune' || (mode === 'refine' && tune)
  const specs = useMemo(() => !defaults ? [] : mode === 'generate' ? buildSpecs(defaults) : refinementSpecs(defaults), [defaults, mode])
  const trainingSpecs = useMemo(() => tuningSpecs(options ?? []), [options])
  const heads: ModelHead[] = tuning.tune_heads === 'rt' ? ['rt'] : tuning.tune_heads === 'ccs' ? ['ccs'] : ['rt', 'ccs']
  const modelKey = mode === 'generate' ? `${modelDir}:all` : usesTuning ? `${modelDir}:${heads.join(',')}` : ''
  const models = modelCheck?.key === modelKey ? modelCheck.status : null
  const modelsReady = !modelKey || (!!models && models.missing.length === 0)

  const applySaved = useCallback((payload: Values, fallbackMode: Mode) => {
    const gui = payload.__gui
    if (gui != null && (!record(gui) || gui.version !== 1 || !modes.includes(gui.mode as Mode))) {
      throw new Error('unsupported desktop settings format')
    }
    revision.current += 1
    if (record(gui)) {
      const selected = gui.mode as Mode
      const restored: Configs = {}
      if (record(gui.configs)) {
        for (const name of modes) if (record(gui.configs[name])) restored[name] = stripPrivate(gui.configs[name] as Values)
      }
      // The visible native config is authoritative; the metadata preserves other modes.
      restored[selected] = { ...restored[selected], ...stripPrivate(payload) }
      setConfigs((old) => Object.fromEntries(modes.map((name) => [name, { ...old[name], ...restored[name] }])))
      setMode(selected)
      if (record(gui.tuning)) setTuning((old) => ({ ...old, ...gui.tuning as Values }))
      if (typeof gui.tune === 'boolean') setTune(gui.tune)
      for (const [key, setter] of [['modelDir', setModelDir], ['tuneOutModels', setTuneOutModels], ['outReport', setOutReport],
        ['fasta', setFasta], ['library', setLibrary], ['ids', setIds], ['out', setOut]] as const) {
        if (typeof gui[key] === 'string') setter(gui[key] as string)
      }
      if (typeof gui.threads === 'number') setThreads(gui.threads)
    } else {
      setMode(fallbackMode)
      setConfigs((old) => ({ ...old, [fallbackMode]: { ...old[fallbackMode], ...stripPrivate(payload),
        ...(fallbackMode === 'generate' ? { nce: payload.nce ?? null } : {}) } }))
      if (typeof payload.__modelDir === 'string') setModelDir(payload.__modelDir)
    }
  }, [])

  useEffect(() => {
    let alive = true
    void (async () => {
      try {
        const [binary, settings] = await Promise.all([window.dialibgen.probe(), window.dialibgen.loadSettings()])
        if (!alive) return
        setInfo(binary)
        setPresets(settings.presets ?? {})
        if (settings.lastUsed) applySaved(settings.lastUsed, 'generate')
      } catch (e) { if (alive) setOutcome({ ok: false, text: `Could not load settings: ${String(e)}` }) }
      if (alive) setInitialized(true)
    })()
    return () => { alive = false }
  }, [applySaved])

  useEffect(() => {
    if (!initialized || nativeDefaults[mode]) return
    let alive = true
    setConfigError(null)
    void window.dialibgen.defaultConfig(mode).then((d) => {
      if (!alive) return
      setNativeDefaults((old) => ({ ...old, [mode]: d }))
      setConfigs((old) => ({ ...old, [mode]: { ...d, ...(mode === 'generate' ? { nce: null } : {}), ...old[mode] } }))
    }).catch((e: unknown) => { if (alive) setConfigError(String(e)) })
    return () => { alive = false }
  }, [initialized, mode, nativeDefaults])

  useEffect(() => {
    if (!usesTuning || options) return
    let alive = true
    setTuningFailure(null)
    void window.dialibgen.tuningOptions().then((schema) => {
      if (!alive) return
      setOptions(schema)
      setTuning((old) => ({ ...Object.fromEntries(schema.map((s) => [s.name, s.value])), ...old }))
    }).catch((e: unknown) => { if (alive) setTuningFailure(String(e)) })
    return () => { alive = false }
  }, [usesTuning, options])

  useEffect(() => {
    if (!modelKey) return
    let alive = true
    const result = mode === 'generate' ? window.dialibgen.models(modelDir || undefined)
      : window.dialibgen.models(modelDir || undefined, heads)
    void result.then((status) => { if (alive) setModelCheck({ key: modelKey, status }) })
      .catch((e: unknown) => { if (alive) setModelCheck({ key: modelKey, status: { dir: null, missing: [String(e)] } }) })
    return () => { alive = false }
  // The key includes the directory and every required head; no array identity dependency.
  }, [modelKey])

  useEffect(() => {
    const offLog = window.dialibgen.onLog((line) => setLog((old) => [...old.slice(-(MAX_LOG_LINES - 1)), line]))
    const offDone = window.dialibgen.onDone((r: RunResult) => {
      setRunning(false)
      setOutcome(r.ok && r.bytes ? { ok: true, text: `Wrote ${humanBytes(r.bytes)}`, path: outRef.current }
        : r.ok ? { ok: false, text: 'The tool exited cleanly but no library is on disk.' }
          : { ok: false, text: `Failed (exit ${r.code ?? 'killed'}). See the log.` })
    })
    return () => { offLog(); offDone() }
  }, [])
  useEffect(() => { if (logRef.current) logRef.current.scrollTop = logRef.current.scrollHeight }, [log])

  const set = (name: string, value: unknown) => setConfigs((old) => ({ ...old, [mode]: { ...old[mode], [name]: value } }))
  const outExt = /\.parquet$/.test(out) ? 'parquet' : /\.tsv$/.test(out) ? 'tsv' : null
  const threadsValid = threads !== null && Number.isSafeInteger(threads) && threads >= 0 && threads <= 2147483647
  const configValidation = defaults && (Object.keys(toolConfig(values)).find((key) => !(key in defaults))
    ? 'The configuration contains an option from another mode or an unknown option.'
    : specs.map((spec) => invalidValue(spec, values[spec.name])).find(Boolean))
  const tuneValidation = usesTuning && (Object.keys(tuning).some((key) => !options?.some((s) => s.name === key))
    ? 'The fine-tuning settings contain an unknown option.'
    : trainingSpecs.map((spec) => invalidValue(spec, tuning[spec.name])).find(Boolean) || tuningError(tuning))
  const modeValidation = mode === 'tune' && defaults && Object.keys(defaults).some((key) => JSON.stringify(values[key]) !== JSON.stringify(defaults[key]))
    ? 'Fine-tune mode requires its native prediction-only configuration. Reset settings or select Refine to replace observed values.'
    : mode === 'refine' && usesTuning && values.rt_unit === 'minmax' ? 'Fine-tuning requires observed RT units, not minmax.'
      : mode === 'refine' && values.write_rt === true && values.filter === false && (!usesTuning || tuning.tune_heads === 'ccs')
        ? 'Observed RT with filtering off requires RT fine-tuning, or turn off write rt.' : null
  const validation = configValidation || tuneValidation || modeValidation
  const inputReady = mode === 'generate' ? !!fasta.trim() : /\.(parquet|tsv)$/.test(library) && /\.parquet$/.test(ids)
  const canRun = initialized && threadsValid && !running && inputReady && !!out.trim() && !!outExt && modelsReady
    && !!defaults && !!info?.ok && (!usesTuning || !!options?.length) && !validation && (mode === 'generate' || !outReport || /\.tsv$/.test(outReport))
  const snapshot = (): Values => ({ ...toolConfig(values), __gui: { version: 1, mode, configs, tuning, tune, modelDir,
    tuneOutModels, outReport, threads, fasta, library, ids, out } })

  async function start(): Promise<void> {
    if (!canRun) return
    setRunning(true); outRef.current = out; setLog([]); setOutcome(null)
    try {
      if (!await window.dialibgen.saveLast(snapshot())) throw new Error('could not save settings; check that the app config directory is writable')
      const result = await window.dialibgen.run({ mode, in: mode === 'generate' ? fasta : library, out,
        config: toolConfig(values), modelDir: mode === 'generate' || usesTuning ? modelDir || null : null, threads,
        ...(mode !== 'generate' ? { ids, ...(outReport ? { outReport } : {}) } : {}),
        ...(usesTuning ? { tuning, ...(tuneOutModels ? { tuneOutModels } : {}), ...(mode === 'refine' ? { tune: true } : {}) } : {}) })
      if (!result.started) { setRunning(false); setOutcome({ ok: false, text: result.reason ?? 'could not start' }) }
    } catch (e) { setRunning(false); setOutcome({ ok: false, text: `Could not start: ${String(e)}` }) }
  }
  async function load(): Promise<void> {
    const token = revision.current
    const path = await window.dialibgen.pickConfig()
    if (!path) return
    const parsed = await window.dialibgen.readConfig(path)
    if (token === revision.current) applySaved(parsed, mode)
  }
  async function save(): Promise<void> {
    const payload = snapshot()
    const path = await window.dialibgen.saveConfig('dialibgen-settings.json')
    if (path && !await window.dialibgen.writeConfig(path, payload)) throw new Error('settings could not be written')
  }
  function reset(): void {
    if (!defaults) return
    revision.current += 1
    setConfigs((old) => ({ ...old, [mode]: { ...defaults, ...(mode === 'generate' ? { nce: null } : {}) } }))
    if (usesTuning && options) setTuning(Object.fromEntries(options.map((s) => [s.name, s.value])))
  }
  const renderParam = (spec: ParamSpec, training = false) => <ParamField key={spec.name} spec={spec}
    value={(training ? tuning : values)[spec.name]} inert={training ? null : inertBecause(spec.name, values)}
    onChange={(value) => training ? setTuning((old) => ({ ...old, [spec.name]: value })) : set(spec.name, value)} />

  return <div className="app">
    <header><h1>DIALibGen</h1>
      <span className={`badge ${info?.ok ? 'ok' : 'bad'}`}>{info ? info.detail : 'probing…'}</span>
      <span className={`badge ${modelsReady ? 'ok' : 'bad'}`}>{!modelKey ? 'models not required' : !models ? 'checking models…'
        : modelsReady ? 'models ready' : `${models.missing.length} model${models.missing.length === 1 ? '' : 's'} missing`}</span>
    </header>
    <main><section className="controls">
      <div className="field"><label htmlFor="mode">workflow</label><select id="mode" value={mode} disabled={running || !initialized}
        onChange={(e) => { revision.current += 1; setMode(e.target.value as Mode); setConfigError(null); setShowAdvanced(false) }}>
        <option value="generate">Generate</option><option value="refine">Refine</option><option value="tune">Fine-tune</option>
      </select><p className="help">{mode === 'generate' ? 'Predict a spectral library from protein sequences.' : mode === 'refine'
        ? 'Refine a library with observations from a DIA-NN report.' : 'Train RT/CCS models and re-predict the whole library. Observed-value replacement belongs to Refine.'}</p></div>
      {mode === 'generate' ? <PathField id="fasta" label="protein FASTA" value={fasta} set={setFasta} pick={window.dialibgen.pickFasta} />
        : <><PathField id="library" label="input library" value={library} set={setLibrary} pick={window.dialibgen.pickLibrary} help="Existing library ending in .parquet or .tsv." />
          <PathField id="ids" label="DIA-NN report (Parquet)" value={ids} set={setIds} pick={window.dialibgen.pickReport} help="DIA-NN report ending in .parquet." /></>}
      <PathField id="out" label="output library" value={out} set={setOut} pick={() => window.dialibgen.pickOutput(out || undefined)} />
      <div className="row tight formats">{(['parquet', 'tsv'] as const).map((ext) => <button type="button" key={ext}
        className={`chip${outExt === ext ? ' on' : ''}`} disabled={!out} onClick={() => setOut(withExtension(out, ext))}>.{ext}</button>)}</div>
      <p className="help">{outExt === 'tsv' ? 'DIA-NN dialect. A TSV cannot carry the recipe, so this library is not reproducible from itself.'
        : outExt === 'parquet' ? 'Carries the config verbatim in its schema metadata, so the library states how it was built.' : 'The name must end in .parquet or .tsv.'}</p>
      {mode === 'refine' && <div className="field check"><input id="tune" type="checkbox" checked={tune} disabled={running} onChange={(e) => setTune(e.target.checked)} />
        <label htmlFor="tune">fine-tune before refinement</label><p className="help">Train selected models and re-predict before applying observed values.</p></div>}
      {modelKey && <><PathField id="models" label="model directory" value={modelDir} set={setModelDir} pick={window.dialibgen.pickDirectory} placeholder="(search the usual places)" />
        <p className={`help${models?.missing.length ? ' bad' : ''}`}>{models?.missing.length ? `Not found: ${models.missing.join(', ')}. Release builds include these models.` : models?.dir ?? ''}</p></>}
      {usesTuning && <PathField id="tuned-models" label="keep tuned models in (optional)" value={tuneOutModels} set={setTuneOutModels} pick={window.dialibgen.pickDirectory}
        help="Keep ONNX models, training summaries and trajectories. Empty uses temporary models; the library is always written." />}
      {mode !== 'generate' && <PathField id="out-report" label="residual report (optional TSV)" value={outReport} set={setOutReport} />}
      <div className="field"><label htmlFor="threads">threads</label><input id="threads" type="number" min={0} max={2147483647} step={1}
        value={threads ?? ''} aria-invalid={!threadsValid} onChange={(e) => setThreads(e.target.value === '' ? null : Number(e.target.value))} />
        <p className={`help${threadsValid ? '' : ' bad'}`}>{threadsValid ? '1 is the default. 0 uses all available cores, with more memory.' : 'Enter a whole number of threads (0 or more).'}</p></div>
      {validation && <p className="help bad" role="alert">{validation}</p>}
      {mode !== 'generate' && outReport && !/\.tsv$/.test(outReport) && <p className="help bad">The residual report name must end in .tsv.</p>}
      <div className="row"><button type="button" disabled={!canRun} onClick={() => void start()}>{running ? 'Running…' : mode === 'generate' ? 'Generate library' : mode === 'refine' ? 'Refine library' : 'Fine-tune library'}</button>
        <button type="button" className="secondary" disabled={!running} onClick={() => void window.dialibgen.cancel().catch((e: unknown) => setOutcome({ ok: false, text: `Could not cancel: ${String(e)}` }))}>Cancel</button></div>
      <Presets presets={presets} values={snapshot()} disabled={running || !defaults} onApply={(payload) => {
        try { applySaved(payload, 'generate') } catch (e) { setOutcome({ ok: false, text: `Could not load preset: ${String(e)}` }) }
      }} onSaved={setPresets} onError={(text) => setOutcome({ ok: false, text })} />
      {outcome && <div className={`outcome ${outcome.ok ? 'ok' : 'bad'}`}><span>{outcome.text}</span>
        {outcome.ok && outcome.path && <button type="button" className="secondary slim" onClick={() => void window.dialibgen.revealPath(outcome.path as string)}>Show in folder</button>}</div>}
    </section><section className="right"><div className="params">
      <div className="section-head"><h2>Desktop settings</h2><div className="row tight actions">
        <button type="button" className="secondary slim" disabled={!defaults || running} onClick={reset}>Reset</button>
        <button type="button" className="secondary slim" disabled={running || !initialized} onClick={() => void load().catch((e: unknown) => setOutcome({ ok: false, text: `Could not load config: ${String(e)}` }))}>Load…</button>
        <button type="button" className="secondary slim" disabled={running || !defaults} onClick={() => void save().catch((e: unknown) => setOutcome({ ok: false, text: `Could not save config: ${String(e)}` }))}>Save…</button>
      </div></div>
      <p className="help">Saved desktop settings include the workflow and fine-tuning controls. Reopen them here; use the CLI’s own configuration export for -config.</p>
      {configError && <p className="help bad">Could not read the tool’s own defaults: {configError}. The form is unavailable until the binary runs — check the path in the badge above.</p>}
      {mode !== 'tune' && <><h3>{mode === 'generate' ? 'Library configuration' : 'Refinement configuration'}</h3>
        {specs.filter((s) => !s.hidden && s.group === 'core').map((s) => renderParam(s))}
        {specs.some((s) => !s.hidden && s.group === 'advanced') && <><button type="button" className="disclose" aria-expanded={showAdvanced} onClick={() => setShowAdvanced((old) => !old)}>
          {showAdvanced ? '▾' : '▸'} Advanced ({specs.filter((s) => !s.hidden && s.group === 'advanced').length})</button>
          {showAdvanced && specs.filter((s) => !s.hidden && s.group === 'advanced').map((s) => renderParam(s))}</>}
      </>}
      {usesTuning && <><h3>Fine-tuning</h3><p className="help">CPU training is included in release builds. CUDA training and GPU re-prediction require a suitable source build and matching GPU libraries.</p>
        {tuningFailure && <p className="help bad">Could not read fine-tuning options: {tuningFailure}</p>}
        {!options && !tuningFailure && <p className="help">Loading the tool’s fine-tuning options…</p>}
        {trainingSpecs.filter((s) => s.group === 'core').map((s) => renderParam(s, true))}
        {trainingSpecs.some((s) => s.group === 'advanced') && <><button type="button" className="disclose" aria-expanded={showTrainingAdvanced} onClick={() => setShowTrainingAdvanced((old) => !old)}>
          {showTrainingAdvanced ? '▾' : '▸'} Advanced fine-tuning ({trainingSpecs.filter((s) => s.group === 'advanced').length})</button>
          {showTrainingAdvanced && trainingSpecs.filter((s) => s.group === 'advanced').map((s) => renderParam(s, true))}</>}
      </>}
    </div><div className="logpane"><div className="section-head"><h2>Log</h2></div><pre ref={logRef} className="log">{log.join('\n')}</pre></div></section></main>
  </div>
}

function PathField({ id, label, value, set, pick, placeholder, help }: { id: string; label: string; value: string;
  set: (value: string) => void; pick?: () => Promise<string | null>; placeholder?: string; help?: string }): JSX.Element {
  const [error, setError] = useState('')
  return <div className="field"><label htmlFor={id}>{label}</label><div className="row tight">
    <input id={id} type="text" value={value} placeholder={placeholder ?? '(not set)'} onChange={(e) => set(e.target.value)} />
    {pick && <button type="button" className="secondary slim" onClick={() => void pick().then((path) => { setError(''); if (path) set(path) }).catch((e: unknown) => setError(String(e)))}>Browse…</button>}
  </div>{help && <p className="help">{help}</p>}{error && <p className="help bad">{error}</p>}</div>
}

function Presets({
  presets,
  values,
  onApply,
  onSaved,
  onError,
  disabled
}: {
  presets: Record<string, Values>
  values: Values
  onApply: (v: Values) => void
  onSaved: (p: Record<string, Values>) => void
  onError: (message: string) => void
  disabled: boolean
}): JSX.Element {
  const [name, setName] = useState('')
  const [selected, setSelected] = useState('')
  const names = Object.keys(presets).sort()
  return (
    <div className="field presets">
      <label htmlFor="preset">presets</label>
      <div className="row tight">
        <select id="preset" disabled={disabled} value={selected} onChange={(e) => {
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
        <button type="button" className="secondary slim" disabled={disabled || !selected}
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
        <button type="button" className="secondary slim" disabled={disabled || !name.trim()}
                onClick={() => void (async () => {
                  const clean = name.trim()
                  const payload = values
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
