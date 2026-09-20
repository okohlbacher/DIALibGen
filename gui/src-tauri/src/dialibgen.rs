// DIALibGen CLI integration: locate the binary, ask it for its own
// effective config, run it, stream its stderr as events, and cancel by killing
// the child.
//
// The tool's `-write_config` output is the single source of truth for the
// parameter set: the frontend builds its form from it and build_config()
// validates against it, so neither this file nor the frontend carries a second,
// drifting copy of the defaults.
//
// Security: the binary is spawned with an explicit argv array (never a shell
// string), the config reaches it as a FILE rather than as command-line values,
// and only keys the tool itself declared are written into that file.

use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Arc, Mutex};

use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use tauri::{AppHandle, Emitter, Manager, State};

fn exe_name() -> &'static str {
    if cfg!(windows) {
        "DIALibGen.exe"
    } else {
        "DIALibGen"
    }
}

/// The three ONNX exports bundled with release builds.
pub const MODEL_FILES: [&str; 3] = [
    "peptdeep_rt_dynamic.onnx",
    "peptdeep_ms2_dynamic.onnx",
    "peptdeep_ccs_dynamic.onnx",
];

/// The single in-flight run. Library generation is long and memory-hungry; a
/// second `run` while one is live is refused rather than queued.
#[derive(Default)]
pub struct RunManager {
    pub current: Mutex<Option<CurrentRun>>,
    pub launch: Mutex<()>,
    pub seq: AtomicU64,
    cancel_seq: AtomicU64,
    shutting_down: AtomicBool,
}

pub struct CurrentRun {
    pub run_id: u64,
    pub child: Arc<Mutex<Child>>,
    pub config_dir: Option<tempfile::TempDir>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct BinaryInfo {
    bin: String,
    source: String,
    ok: bool,
    version: Option<String>,
    detail: String,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ModelStatus {
    /// Where the models were found, if all three are in one directory.
    dir: Option<String>,
    /// The file names that are NOT there. Empty means ready to run.
    missing: Vec<String>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct RunStarted {
    started: bool,
    reason: Option<String>,
}

#[derive(Clone, Copy, Default, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Mode {
    #[default]
    Generate,
    Refine,
    Tune,
}

impl Mode {
    fn as_str(self) -> &'static str {
        match self { Self::Generate => "generate", Self::Refine => "refine", Self::Tune => "tune" }
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RunParams {
    #[serde(default)]
    mode: Mode,
    #[serde(rename = "in")]
    input: String,
    out: String,
    #[serde(default)]
    ids: Option<String>,
    #[serde(default)]
    tune: bool,
    #[serde(default)]
    tuning: Map<String, Value>,
    #[serde(default)]
    tune_out_models: Option<String>,
    #[serde(default)]
    out_report: Option<String>,
    /// The library config, as the form produced it. Filtered against the tool's
    /// own effective config before it is written out.
    #[serde(default)]
    config: Map<String, Value>,
    /// Optional directory holding the three .onnx files; passed to the child as
    /// DIALIBGEN_MODEL_DIR. The tool searches its own locations when unset.
    #[serde(default)]
    model_dir: Option<String>,
    /// 1 matches the TOPP default; 0 requests all available cores.
    #[serde(default)]
    threads: Option<i64>,
}

pub struct Resolved {
    pub bin: PathBuf,
    pub data: Option<PathBuf>,
    /// The bundled share/DIALibGen, if present.
    pub share: Option<PathBuf>,
    pub source: &'static str,
}

/// Resolve the binary in priority order: DIALIBGEN_BIN env override, then
/// bundled beside the app, then the bare name on PATH.
pub fn resolve_binary<R: tauri::Runtime>(app: &AppHandle<R>) -> Resolved {
    if let Ok(env_bin) = std::env::var("DIALIBGEN_BIN") {
        if !env_bin.is_empty() && PathBuf::from(&env_bin).exists() {
            let data = std::env::var("OPENMS_DATA_PATH").ok().map(PathBuf::from);
            return Resolved { bin: PathBuf::from(env_bin), data, share: None, source: "env" };
        }
    }
    let mut roots: Vec<PathBuf> = Vec::new();
    if let Ok(rd) = app.path().resource_dir() {
        roots.push(rd.join("resources").join("dialibgen"));
    }
    // Dev fallback: `tauri dev` does not populate resource_dir the way a bundle
    // does, so the resources dir next to this crate is baked in at compile time.
    #[cfg(debug_assertions)]
    roots.push(PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/resources/dialibgen")));
    for root in roots {
        let bin = root.join("bin").join(exe_name());
        if bin.exists() {
            let data = root.join("share").join("OpenMS");
            let data = if data.exists() { Some(data) } else { None };
            let share = root.join("share").join("DIALibGen");
            let share = if share.exists() { Some(share) } else { None };
            return Resolved { bin, data, share, source: "bundled" };
        }
    }
    Resolved { bin: PathBuf::from(exe_name()), data: None, share: None, source: "path" }
}

/// Everything the child needs that this process knows and it does not.
fn apply_env(cmd: &mut Command, r: &Resolved) {
    #[cfg(target_os = "linux")]
    if r.source == "bundled" {
        // AppRun prepends the GUI's GTK libraries. The private CLI has its own
        // verified relative RUNPATHs and must not inherit that different SDK.
        cmd.env_remove("LD_LIBRARY_PATH");
    }
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        cmd.creation_flags(0x08000000); // CREATE_NO_WINDOW: the GUI owns the log.
    }
    if let Some(d) = &r.data {
        cmd.env("OPENMS_DATA_PATH", d);
    }
    if let Some(s) = &r.share {
        // Point the tool's data lookup at the bundled share explicitly. Its own
        // <bin>/../share arithmetic is correct for an install but a copied or
        // symlinked binary -- which is exactly what `tauri dev` sets up -- can
        // land somewhere that arithmetic does not reach.
        cmd.env("DIALIBGEN_DATA_DIR", s);
    }
}

/// TOPP --help prints "... Version: X.Y.Z ...". Capture the version token as
/// proof the binary actually runs and links its libraries on this machine.
fn parse_version(text: &str) -> Option<String> {
    let idx = text.find("Version:")?;
    let rest = text[idx + "Version:".len()..].trim_start();
    let tok: String = rest.chars().take_while(|c| !c.is_whitespace() && *c != ',').collect();
    if tok.chars().next().is_some_and(|c| c.is_ascii_digit()) {
        Some(tok)
    } else {
        None
    }
}

#[tauri::command(async)]
pub fn probe<R: tauri::Runtime>(app: AppHandle<R>) -> BinaryInfo {
    let r = resolve_binary(&app);
    probe_resolved(&r)
}

fn probe_resolved(r: &Resolved) -> BinaryInfo {
    let bin_s = r.bin.display().to_string();
    let mut cmd = Command::new(&r.bin);
    cmd.arg("--help");
    apply_env(&mut cmd, r);
    match cmd.output() {
        Ok(out) => {
            let text = format!(
                "{}\n{}",
                String::from_utf8_lossy(&out.stdout),
                String::from_utf8_lossy(&out.stderr)
            );
            let version = parse_version(&text);
            let ok = out.status.success() && version.is_some();
            let detail = match (&version, out.status.success()) {
                (Some(v), true) => format!("DIALibGen {v} ({})", r.source),
                (_, false) => format!("cannot execute successfully ({}): {}", out.status, text.trim()),
                (None, true) => format!("unrecognized DIALibGen version ({})", r.source),
            };
            BinaryInfo { bin: bin_s, source: r.source.into(), ok, version, detail }
        }
        Err(e) => BinaryInfo {
            bin: bin_s,
            source: r.source.into(),
            ok: false,
            version: None,
            detail: format!("cannot execute: {e}"),
        },
    }
}

/// Which of the three models are present in @p dir, or -- when @p dir is None --
/// in the places the tool itself searches that this process can also see.
#[tauri::command]
pub fn models<R: tauri::Runtime>(app: AppHandle<R>, dir: Option<String>, heads: Option<Vec<String>>) -> Result<ModelStatus, String> {
    let files: Vec<&str> = match heads {
        None => MODEL_FILES.to_vec(),
        Some(heads) if heads.is_empty() => return Err("select at least one model head".into()),
        Some(heads) => heads.iter().map(|head| match head.as_str() {
            "rt" => Ok(MODEL_FILES[0]), "ms2" => Ok(MODEL_FILES[1]), "ccs" => Ok(MODEL_FILES[2]),
            _ => Err(format!("unknown model head: {head}")),
        }).collect::<Result<_, _>>()?,
    };
    let mut candidates: Vec<PathBuf> = Vec::new();
    if let Some(d) = dir.as_deref().filter(|d| !d.is_empty()) {
        candidates.push(PathBuf::from(d));
    } else {
        if let Ok(env) = std::env::var("DIALIBGEN_MODEL_DIR") {
            if !env.is_empty() {
                candidates.push(PathBuf::from(env));
            }
        }
        let r = resolve_binary(&app);
        if let Some(s) = &r.share {
            candidates.push(s.join("models"));
        }
        if let Some(parent) = r.bin.parent() {
            candidates.push(parent.join("..").join("share").join("DIALibGen").join("models"));
        }
        if let Some(d) = &r.data {
            candidates.push(d.join("models"));
        }
    }
    // Report the first directory that has ALL requested heads: a split set is not usable,
    // and naming a directory that holds one of them would be misleading.
    let mut best: Option<(PathBuf, Vec<String>)> = None;
    for c in candidates {
        let missing: Vec<String> = files
            .iter()
            .filter(|f| !c.join(f).is_file())
            .map(|f| (*f).to_string())
            .collect();
        if missing.is_empty() {
            return Ok(ModelStatus { dir: Some(c.display().to_string()), missing });
        }
        if best.is_none() || missing.len() < best.as_ref().unwrap().1.len() {
            best = Some((c, missing));
        }
    }
    Ok(match best {
        Some((c, missing)) => ModelStatus { dir: Some(c.display().to_string()), missing },
        None => ModelStatus { dir: None, missing: files.iter().map(|f| f.to_string()).collect() },
    })
}

/// The tool's own effective config: every default materialised. This is what
/// the form is built from, so the GUI can never offer a key the tool does not
/// have, nor default one differently from the CLI.
#[tauri::command(async)]
pub fn default_config<R: tauri::Runtime>(app: AppHandle<R>, mode: Option<Mode>) -> Result<Value, String> {
    let r = resolve_binary(&app);
    default_config_resolved(&r, mode.unwrap_or_default())
}

fn private_tempdir() -> std::io::Result<tempfile::TempDir> {
    let mut builder = tempfile::Builder::new();
    builder.prefix("dialibgen-");
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        builder.permissions(std::fs::Permissions::from_mode(0o700));
    }
    builder.tempdir()
}

fn write_run_config(dir: &Path, config: &Map<String, Value>) -> std::io::Result<PathBuf> {
    let path = dir.join("config.json");
    let mut file = std::fs::OpenOptions::new().write(true).create_new(true).open(&path)?;
    let json = serde_json::to_vec_pretty(config).map_err(std::io::Error::other)?;
    file.write_all(&json)?;
    Ok(path)
}

fn default_config_resolved(r: &Resolved, mode: Mode) -> Result<Value, String> {
    let dir = private_tempdir().map_err(|e| e.to_string())?;
    let path = dir.path().join("effective.json");
    let mut cmd = Command::new(&r.bin);
    cmd.args(["-mode", mode.as_str(), "-write_config"]).arg(&path);
    apply_env(&mut cmd, r);
    let out = cmd.output().map_err(|e| format!("cannot run the tool: {e}"))?;
    if !out.status.success() || !path.is_file() {
        return Err(format!(
            "-write_config failed ({}): {}",
            out.status,
            String::from_utf8_lossy(&out.stderr).trim()
        ));
    }
    let text = std::fs::read_to_string(&path).map_err(|e| e.to_string());
    let value: Value = serde_json::from_str(&text?)
        .map_err(|e| format!("cannot parse the effective config: {e}"))?;
    if !value.is_object() {
        return Err("the tool's effective config is not a JSON object".into());
    }
    Ok(value)
}

#[derive(Debug, Serialize)]
pub struct TuningOption {
    name: String,
    kind: String,
    value: Value,
    description: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    min: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    max: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    choices: Option<Vec<String>>,
}

fn training_key(name: &str) -> bool {
    matches!(name, "tune_heads" | "tune_predict_gpu" | "tune_predict_sessions" | "tune_keep_free_cysteine_offset")
        || ["filter:", "cohort:", "train:", "stop:", "machine:"].iter().any(|prefix| name.starts_with(prefix))
}

fn parse_tuning_options(xml: &[u8]) -> Result<Vec<TuningOption>, String> {
    use quick_xml::{events::Event, Reader};
    let mut reader = Reader::from_reader(quick_xml::encoding::DecodingReader::new(xml));
    reader.config_mut().expand_empty_elements = true;
    let mut buffer = Vec::new();
    let mut nodes = Vec::<String>::new();
    let mut options = Vec::<TuningOption>::new();
    loop {
        match reader.read_event_into(&mut buffer).map_err(|e| format!("invalid TOPP INI: {e}"))? {
            Event::Start(ref element) | Event::Empty(ref element) => {
                let attributes = element.attributes().map(|attribute| {
                    let attribute = attribute.map_err(|e| e.to_string())?;
                    let key = attribute.key.as_ref().to_string();
                    let value = attribute.normalized_value(quick_xml::XmlVersion::Explicit1_0).map_err(|e| e.to_string())?.into_owned();
                    Ok((key, value))
                }).collect::<Result<std::collections::HashMap<_, _>, String>>()?;
                if element.name().as_ref() == "NODE" {
                    nodes.push(attributes.get("name").ok_or("TOPP NODE has no name")?.clone());
                } else if element.name().as_ref() == "ITEM" && nodes.len() >= 2
                    && nodes[0] == "DIALibGen" && nodes[1] == "1" {
                    let leaf = attributes.get("name").ok_or("TOPP ITEM has no name")?;
                    let name = nodes[2..].iter().map(String::as_str).chain([leaf.as_str()]).collect::<Vec<_>>().join(":");
                    if !training_key(&name) { continue; }
                    if !name.chars().all(|c| c.is_ascii_alphanumeric() || c == '_' || c == ':')
                        || options.iter().any(|option| option.name == name) {
                        return Err(format!("invalid or duplicate tuning option: {name}"));
                    }
                    let kind = attributes.get("type").ok_or("TOPP ITEM has no type")?.clone();
                    let raw = attributes.get("value").ok_or("TOPP ITEM has no value")?;
                    let invalid = || format!("invalid {kind} default for {name}");
                    let value = match kind.as_str() {
                        "bool" => Value::Bool(raw.parse::<bool>().map_err(|_| invalid())?),
                        "int" => Value::from(raw.parse::<i32>().map_err(|_| invalid())?),
                        "double" => serde_json::Number::from_f64(raw.parse::<f64>().map_err(|_| invalid())?)
                            .map(Value::Number).ok_or_else(invalid)?,
                        "string" => Value::String(raw.clone()),
                        _ => return Err(format!("unsupported tuning option type: {kind}")),
                    };
                    let mut option = TuningOption { name, kind, value,
                        description: attributes.get("description").cloned().unwrap_or_default(),
                        min: None, max: None, choices: None };
                    if let Some(restriction) = attributes.get("restrictions").filter(|s| !s.is_empty()) {
                        if matches!(option.kind.as_str(), "int" | "double") {
                            let (low, high) = restriction.split_once(':').ok_or("invalid TOPP numeric range")?;
                            let bound = |text: &str| -> Result<Option<f64>, String> {
                                if text.is_empty() { return Ok(None); }
                                let value = text.parse::<f64>().map_err(|_| "invalid TOPP numeric bound")?;
                                if !value.is_finite() { return Err("non-finite TOPP numeric bound".into()); }
                                Ok(Some(value))
                            };
                            option.min = bound(low)?;
                            option.max = bound(high)?;
                        } else if option.kind == "string" {
                            option.choices = Some(restriction.split(',').map(str::to_string).collect());
                        }
                    }
                    options.push(option);
                }
            },
            Event::End(element) if element.name().as_ref() == "NODE" => { nodes.pop(); },
            Event::Decl(declaration) => {
                let encoding = declaration.encoder().ok_or("unsupported TOPP INI encoding")?;
                reader.get_mut().set_encoding(encoding);
            },
            Event::DocType(_) => return Err("TOPP INI must not contain a document type".into()),
            Event::Eof => break,
            _ => {},
        }
    }
    if !nodes.is_empty() || !options.iter().any(|option| option.name == "train:epochs") {
        return Err("the tool's TOPP INI has no complete training parameter schema".into());
    }
    Ok(options)
}

#[tauri::command(async)]
pub fn tuning_options<R: tauri::Runtime>(app: AppHandle<R>) -> Result<Vec<TuningOption>, String> {
    tuning_options_resolved(&resolve_binary(&app))
}

fn tuning_options_resolved(r: &Resolved) -> Result<Vec<TuningOption>, String> {
    let dir = private_tempdir().map_err(|e| e.to_string())?;
    let path = dir.path().join("DIALibGen.ini");
    let mut cmd = Command::new(&r.bin);
    cmd.arg("-write_ini").arg(&path);
    apply_env(&mut cmd, r);
    let output = cmd.output().map_err(|e| format!("cannot read training options: {e}"))?;
    if !output.status.success() {
        return Err(format!("-write_ini failed ({}): {}", output.status, String::from_utf8_lossy(&output.stderr)));
    }
    parse_tuning_options(&std::fs::read(path).map_err(|e| e.to_string())?)
}

fn tuning_argv(schema: &[TuningOption], values: &Map<String, Value>) -> Result<Vec<String>, String> {
    let mut args = Vec::new();
    for (key, value) in values {
        let option = schema.iter().find(|option| &option.name == key)
            .ok_or_else(|| format!("unknown tuning option: {key}"))?;
        let invalid = || format!("invalid value for {key} ({})", option.kind);
        if option.kind == "bool" {
            if value.as_bool().ok_or_else(invalid)? { args.push(format!("-{key}")); }
            continue;
        }
        let text = match option.kind.as_str() {
            "int" => {
                let integer = value.as_i64().filter(|&v| i32::try_from(v).is_ok()).ok_or_else(invalid)?;
                integer.to_string()
            },
            "double" => value.as_f64().filter(|v| v.is_finite()).ok_or_else(invalid)?.to_string(),
            "string" => {
                let text = value.as_str().ok_or_else(invalid)?;
                if let Some(choices) = &option.choices {
                    if !choices.iter().any(|choice| choice == text) { return Err(invalid()); }
                } else if key != "machine:device" || !(text == "cpu" || text == "cuda"
                    || text.strip_prefix("cuda:").is_some_and(|index| !index.is_empty() && index.bytes().all(|c| c.is_ascii_digit()) && index.parse::<u32>().is_ok())) {
                    return Err(invalid());
                }
                text.to_string()
            },
            _ => return Err(invalid()),
        };
        if let Some(number) = value.as_f64() {
            if option.min.is_some_and(|min| number < min) || option.max.is_some_and(|max| number > max) {
                return Err(format!("{key} is outside the tool's permitted range"));
            }
        }
        args.extend([format!("-{key}"), text]);
    }
    Ok(args)
}

/// Filter the form's values against the tool's own key set, so a key the tool
/// would reject (it refuses unknown keys outright) can never reach the file.
fn build_config(reference: &Map<String, Value>, wanted: &Map<String, Value>) -> Map<String, Value> {
    let mut out = Map::new();
    for (k, v) in wanted {
        if reference.contains_key(k) && !v.is_null() {
            out.insert(k.clone(), v.clone());
        }
    }
    out
}

/// A value beginning with '-' would be re-parsed by the CLI as a fresh option
/// (OpenMS treats any token starting '-' + non-digit as one), letting a hostile
/// renderer smuggle registered flags through a path. Genuinely dash-named paths
/// stay usable via the ./ prefix.
fn defang_path(p: &str) -> String {
    let b = p.as_bytes();
    let flaggish = b.first() == Some(&b'-')
        && !matches!(b.get(1), Some(c) if c.is_ascii_digit() || *c == b'.');
    if flaggish {
        format!("./{p}")
    } else {
        p.to_string()
    }
}

fn read_stream<R: tauri::Runtime>(stream: impl Read, app: &AppHandle<R>) {
    let reader = BufReader::new(stream);
    for line in reader.lines() {
        let Ok(line) = line else { break };
        let _ = app.emit("dialibgen:log", line);
    }
}

fn unused_output(path: &Path) -> Result<(), String> {
    match std::fs::symlink_metadata(path) {
        Ok(_) => Err(format!("the output path already exists; choose an unused file name: {}", path.display())),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(format!("cannot check the output path: {e}")),
    }
}

fn validate_run_inputs(params: &RunParams) -> Result<(), String> {
    let input = params.input.trim();
    let output = params.out.trim();
    let ids = params.ids.as_deref().unwrap_or("").trim();
    let report = params.out_report.as_deref().unwrap_or("").trim();
    let models = params.tune_out_models.as_deref().unwrap_or("").trim();
    let training = params.mode == Mode::Tune || params.tune;
    if input.is_empty() || output.is_empty() {
        return Err("an input and an output file are both required".into());
    }
    if !output.ends_with(".parquet") && !output.ends_with(".tsv") {
        return Err("the output file must end in .parquet or .tsv".into());
    }
    if params.threads.is_some_and(|n| n < 0 || n > i32::MAX as i64) {
        return Err("threads must be an integer between 0 and 2147483647".into());
    }
    if params.mode == Mode::Generate {
        if !ids.is_empty() || !report.is_empty() || training || !params.tuning.is_empty() || !models.is_empty() {
            return Err("generation cannot receive refinement or tuning parameters".into());
        }
    } else {
        if !input.ends_with(".tsv") && !input.ends_with(".parquet") {
            return Err("the input library must end in .parquet or .tsv".into());
        }
        if ids.is_empty() { return Err("a DIA-NN report is required for refinement or tuning".into()); }
        if !ids.ends_with(".parquet") { return Err("the DIA-NN report must end in .parquet".into()); }
        if !report.is_empty() && !report.ends_with(".tsv") {
            return Err("the residual report must end in .tsv".into());
        }
        if params.mode == Mode::Tune && params.tune {
            return Err("the tune flag is only used with refine mode".into());
        }
        if !training && (!params.tuning.is_empty() || !models.is_empty()) {
            return Err("enable tuning before setting training parameters".into());
        }
        unused_output(Path::new(&format!("{output}.refine.json")))?;
    }
    unused_output(Path::new(output))?;
    if !report.is_empty() { unused_output(Path::new(report))?; }
    if !models.is_empty() {
        let heads = params.tuning.get("tune_heads").and_then(Value::as_str).unwrap_or("both");
        for head in ["rt", "ccs"].iter().filter(|&&head| heads == "both" || heads == head) {
            for suffix in [".onnx", ".onnx.tune.json", ".onnx.trajectory.tsv"] {
                unused_output(&Path::new(models).join(format!("peptdeep_{head}_dynamic{suffix}")))?;
            }
        }
    }
    Ok(())
}

#[tauri::command(async)]
pub fn run<R: tauri::Runtime>(app: AppHandle<R>, state: State<'_, RunManager>, params: RunParams) -> RunStarted {
    let cancel_seq = state.cancel_seq.load(Ordering::SeqCst);
    let refuse = |reason: String| RunStarted { started: false, reason: Some(reason) };
    let Ok(_launch) = state.launch.try_lock() else {
        return refuse("a run is already starting".into());
    };
    {
        // Cancellation and the final spawn share this short critical section;
        // the potentially slow CLI preflight below never holds it.
        let current = state.current.lock().unwrap();
        if current.is_some() {
            return refuse("a run is already in progress".into());
        }
        if state.shutting_down.load(Ordering::SeqCst) {
            return refuse("the app is closing".into());
        }
    }
    if let Err(e) = validate_run_inputs(&params) { return refuse(e); }

    let r = resolve_binary(&app);
    let training_args = if params.mode == Mode::Tune || params.tune {
        match tuning_options_resolved(&r).and_then(|schema| tuning_argv(&schema, &params.tuning)) {
            Ok(args) => args,
            Err(e) => return refuse(e),
        }
    } else { Vec::new() };
    let reference = match default_config_resolved(&r, params.mode) {
        Ok(Value::Object(m)) => m,
        Ok(_) => return refuse("the tool's effective config is not a JSON object".into()),
        Err(e) => return refuse(e),
    };
    let mut config = build_config(&reference, &params.config);
    // Bind the run to the same complete directory the readiness check reports.
    // An explicit incomplete directory must never fall back to bundled models.
    let model_dir = if params.mode == Mode::Generate || params.mode == Mode::Tune || params.tune {
        let heads = if params.mode == Mode::Generate { None } else {
            Some(match params.tuning.get("tune_heads").and_then(Value::as_str).unwrap_or("both") {
                "rt" => vec!["rt".into()], "ccs" => vec!["ccs".into()], _ => vec!["rt".into(), "ccs".into()],
            })
        };
        let selected = match models(app.clone(), params.model_dir.clone(), heads) {
            Ok(status) if status.missing.is_empty() => status.dir,
            Ok(status) => return refuse(format!("missing models in the selected directory: {}", status.missing.join(", "))),
            Err(error) => return refuse(error),
        };
        let Some(selected) = selected else { return refuse("no model directory found".into()); };
        let selected = match std::fs::canonicalize(&selected) {
            Ok(path) => path,
            Err(error) => return refuse(format!("cannot resolve model directory {selected}: {error}")),
        };
        if params.mode == Mode::Generate {
            // Effective generation defaults contain resolved model paths. Replace
            // them rather than letting stale/imported paths override the picker.
            // Explicit file paths also stop the CLI's per-head fallback search.
            for (key, file) in ["rt_model", "ms2_model", "ccs_model"].into_iter().zip(MODEL_FILES) {
                config.insert(key.into(), Value::String(selected.join(file).display().to_string()));
            }
        }
        Some(selected)
    } else { None };

    // The config goes to the child as a FILE. Two reasons: the tool takes it
    // that way and embeds it verbatim in the Parquet output as the recipe, and
    // a file keeps values with spaces, quotes and unicode off the command line.
    let dir = match private_tempdir() {
        Ok(dir) => dir,
        Err(e) => return refuse(format!("cannot create a temporary config directory: {e}")),
    };
    let config_file = match write_run_config(dir.path(), &config) {
        Ok(path) => path,
        Err(e) => return refuse(format!("cannot write the config file: {e}")),
    };

    let mut cmd = Command::new(&r.bin);
    cmd.args(["-mode", params.mode.as_str(), "-in"])
        .arg(defang_path(params.input.trim()))
        .arg("-config")
        .arg(&config_file)
        .arg("-out")
        .arg(defang_path(params.out.trim()))
        .arg("-threads")
        .arg(params.threads.unwrap_or(1).to_string())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    if params.tune { cmd.arg("-tune"); }
    for (flag, value) in [("-ids", &params.ids), ("-out_report", &params.out_report), ("-tune_out_models", &params.tune_out_models)] {
        if let Some(value) = value.as_deref().map(str::trim).filter(|v| !v.is_empty()) {
            cmd.arg(flag).arg(defang_path(value));
        }
    }
    cmd.args(training_args);
    apply_env(&mut cmd, &r);
    if let Some(directory) = model_dir {
        cmd.env("DIALIBGEN_MODEL_DIR", directory);
    }

    let mut current = state.current.lock().unwrap();
    if state.shutting_down.load(Ordering::SeqCst) || state.cancel_seq.load(Ordering::SeqCst) != cancel_seq {
        return refuse("run cancelled before launch".into());
    }
    let mut child = match cmd.spawn() {
        Ok(c) => c,
        // A spawn failure is known synchronously: report it as a non-start so
        // the frontend surfaces "could not start" and settles cleanly.
        Err(e) => {
            return refuse(format!("cannot start: {e}"));
        }
    };

    let run_id = state.seq.fetch_add(1, Ordering::SeqCst) + 1;
    let stdout = child.stdout.take();
    let stderr = child.stderr.take();
    let child_arc = Arc::new(Mutex::new(child));
    *current =
        Some(CurrentRun { run_id, child: child_arc.clone(), config_dir: Some(dir) });
    drop(current);

    // One reader thread per stream. stdout is drained as log even though the
    // tool writes nothing there, so a full pipe can never block the child.
    let a_err = app.clone();
    let t_err = std::thread::spawn(move || {
        if let Some(s) = stderr {
            read_stream(s, &a_err);
        }
    });
    let a_out = app.clone();
    let t_out = std::thread::spawn(move || {
        if let Some(s) = stdout {
            read_stream(s, &a_out);
        }
    });

    // Coordinator: once both streams hit EOF the process has ended; reap it,
    // clean up the config, emit exactly one terminal event, release the slot.
    let a_done = app.clone();
    let out_path = params.out.trim().to_string();
    std::thread::spawn(move || {
        let _ = t_err.join();
        let _ = t_out.join();
        let status = child_arc.lock().unwrap().wait();
        let (ok, code) = match status {
            Ok(s) => (s.success(), s.code()),
            Err(_) => (false, None),
        };
        {
            let st = a_done.state::<RunManager>();
            let mut cur = st.current.lock().unwrap();
            if cur.as_ref().is_some_and(|c| c.run_id == run_id) {
                *cur = None;
            }
        }
        let bytes = std::fs::metadata(Path::new(&out_path)).map(|m| m.len()).ok();
        let _ = a_done.emit(
            "dialibgen:done",
            serde_json::json!({ "ok": ok, "code": code, "bytes": bytes }),
        );
    });

    RunStarted { started: true, reason: None }
}

/// Read a JSON library config the user picked.
///
/// A command rather than the fs plugin: the plugin's scope is a compile-time
/// allowlist of paths, and the whole point here is a file the user chose at
/// run time. Narrowed to JSON -- the only thing the form can consume -- so this
/// is not a general "read any file" capability handed to the renderer.
#[tauri::command]
pub fn read_config(path: String) -> Result<Value, String> {
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{path}: {e}"))?;
    let v: Value = serde_json::from_str(&text).map_err(|e| format!("{path}: {e}"))?;
    if v.is_object() {
        Ok(v)
    } else {
        Err(format!("{path}: not a JSON object"))
    }
}

/// Write a JSON library config to a path the user picked. Same reasoning as
/// read_config, and the value is re-serialised here rather than passed through
/// as text, so only well-formed JSON can ever be written.
#[tauri::command]
pub fn write_config(path: String, value: Value) -> Result<bool, String> {
    let text = serde_json::to_string_pretty(&value).map_err(|e| e.to_string())?;
    std::fs::write(&path, text + "\n").map_err(|e| format!("{path}: {e}"))?;
    Ok(true)
}

#[tauri::command]
pub fn cancel(state: State<'_, RunManager>) -> Value {
    cancel_current(&state)
}

fn cancel_current(state: &RunManager) -> Value {
    let cur = state.current.lock().unwrap();
    state.cancel_seq.fetch_add(1, Ordering::SeqCst);
    if let Some(run) = cur.as_ref() {
        if let Ok(mut child) = run.child.lock() {
            let _ = child.kill(); // the coordinator reaps, cleans up and emits done
        }
        serde_json::json!({ "cancelled": true })
    } else {
        serde_json::json!({ "cancelled": state.launch.try_lock().is_err() })
    }
}

/// Finish cancellation and release temporary files before the app exits.
pub fn shutdown(state: &RunManager) {
    let run = {
        let mut current = state.current.lock().unwrap();
        state.shutting_down.store(true, Ordering::SeqCst);
        state.cancel_seq.fetch_add(1, Ordering::SeqCst);
        current.take()
    };
    if let Some(run) = run {
        if let Ok(mut child) = run.child.lock() {
            let _ = child.kill();
            let _ = child.wait();
        }
        if let Some(dir) = run.config_dir {
            let _ = dir.close();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Exported by the verified 0.11.0 CLI; runtime options always come from the
    // selected executable, while this fixture detects parser/schema drift.
    const TOPP_INI: &[u8] = include_bytes!("../testdata/tuning-options.ini");

    fn obj(s: &str) -> Map<String, Value> {
        match serde_json::from_str(s).unwrap() {
            Value::Object(m) => m,
            _ => panic!("not an object"),
        }
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn bundled_cli_does_not_inherit_appimage_gui_libraries() {
        for (source, removed) in [("bundled", true), ("env", false), ("path", false)] {
            let resolved = Resolved { bin: "DIALibGen".into(), data: None, share: None, source };
            let mut command = Command::new(&resolved.bin);
            apply_env(&mut command, &resolved);
            assert_eq!(command.get_envs().any(|(key, value)| key == "LD_LIBRARY_PATH" && value.is_none()), removed);
        }
    }

    #[test]
    fn native_topp_training_schema_preserves_types_ranges_entities_and_encoding() {
        let options = parse_tuning_options(TOPP_INI).unwrap();
        assert_eq!(options.len(), 28);
        assert!(options.iter().all(|option| training_key(&option.name)));
        assert!(!options.iter().any(|option| ["tune", "tune_models", "tune_out_models"].contains(&option.name.as_str())));
        let option = |name| options.iter().find(|option| option.name == name).unwrap();
        assert_eq!(option("filter:min_charge").min, Some(1.0));
        assert_eq!(option("filter:min_charge").max, Some(8.0));
        assert_eq!(option("train:epochs").min, Some(1.0));
        assert_eq!(option("train:epochs").max, None);
        assert_eq!(option("train:lr").value, serde_json::json!(0.0001));
        assert_eq!(option("cohort:full_fit").value, false);
        assert_eq!(option("tune_heads").choices, Some(vec!["rt".into(), "ccs".into(), "both".into()]));
        assert!(option("filter:rt_max_minutes").description.contains("run's maximum"));
        let minimal = b"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?><PARAMETERS><NODE name=\"DIALibGen\"><NODE name=\"1\"><NODE name=\"train\"><ITEM name=\"epochs\" type=\"int\" value=\"100\" description=\"caf\xe9 &amp; &lt;epochs&gt;\"/></NODE></NODE></NODE></PARAMETERS>";
        assert_eq!(parse_tuning_options(minimal).unwrap()[0].description, "café & <epochs>");
        for invalid in [b"<broken>".as_slice(), b"<!DOCTYPE x><PARAMETERS/>",
                        b"<NODE name=\"another-tool\"><NODE name=\"1\"><ITEM name=\"train:epochs\" type=\"int\" value=\"1\"/></NODE></NODE>"] {
            assert!(parse_tuning_options(invalid).is_err());
        }
        let duplicate = String::from_utf8(TOPP_INI.to_vec()).unwrap().replace(
            "<ITEM name=\"warmup\"", "<ITEM name=\"epochs\"");
        assert!(parse_tuning_options(duplicate.as_bytes()).unwrap_err().contains("duplicate"));
    }

    #[test]
    fn training_arguments_reject_wrong_types_ranges_and_option_injection() {
        let options = parse_tuning_options(TOPP_INI).unwrap();
        let args = tuning_argv(&options, &obj(r#"{"train:epochs":25,"train:lr":0.0002,"tune_heads":"rt","cohort:full_fit":true,"machine:no_cudnn":false,"machine:device":"cuda:1"}"#)).unwrap();
        assert!(args.contains(&"-cohort:full_fit".into()));
        assert!(!args.contains(&"true".into()));
        assert!(!args.contains(&"-machine:no_cudnn".into()));
        assert!(args.windows(2).any(|pair| pair == ["-train:epochs", "25"]));
        assert!(args.windows(2).any(|pair| pair == ["-machine:device", "cuda:1"]));
        for invalid in [r#"{"-force":true}"#, r#"{"tune_models":"-force"}"#, r#"{"train:epochs":2.5}"#,
            r#"{"train:epochs":"25"}"#, r#"{"train:epochs":0}"#, r#"{"train:epochs":2147483648}"#,
            r#"{"filter:min_charge":9}"#, r#"{"train:lr":"NaN"}"#, r#"{"train:lr":null}"#,
            r#"{"cohort:full_fit":"false"}"#, r#"{"tune_heads":"rt -force"}"#,
            r#"{"machine:device":"-force"}"#, r#"{"machine:device":"cuda:1 -force"}"#,
            r#"{"machine:device":"cuda:"}"#, r#"{"machine:device":"cuda:999999999999999"}"#] {
            assert!(tuning_argv(&options, &obj(invalid)).is_err(), "{invalid}");
        }
    }

    #[test]
    fn mode_inputs_and_all_output_artifacts_are_checked_before_launch() {
        let dir = private_tempdir().unwrap();
        let base = serde_json::json!({"mode":"refine", "in":"input.tsv", "ids":"ids.parquet", "out":dir.path().join("out.tsv")});
        let validate = |value: Value| validate_run_inputs(&serde_json::from_value::<RunParams>(value).unwrap());
        assert!(validate(base.clone()).is_ok());
        for (key, value) in [("mode", serde_json::json!("generate")), ("in", serde_json::json!("in.fasta")),
            ("ids", serde_json::json!("")), ("ids", serde_json::json!("in.tsv")),
            ("outReport", serde_json::json!("report.parquet")), ("threads", serde_json::json!(-1)),
            ("tuning", serde_json::json!({"train:epochs":1}))] {
            let mut invalid = base.clone(); invalid[key] = value;
            assert!(validate(invalid).is_err(), "{key}");
        }
        let mut training = base.clone(); training["tune"] = true.into();
        training["tuning"] = serde_json::json!({"tune_heads":"rt"});
        training["tuneOutModels"] = dir.path().display().to_string().into();
        assert!(validate(training.clone()).is_ok());
        for filename in ["out.tsv", "out.tsv.refine.json", "peptdeep_rt_dynamic.onnx",
            "peptdeep_rt_dynamic.onnx.tune.json", "peptdeep_rt_dynamic.onnx.trajectory.tsv"] {
            let path = dir.path().join(filename);
            std::fs::write(&path, "preserve").unwrap();
            assert!(validate(training.clone()).is_err(), "{filename}");
            assert_eq!(std::fs::read_to_string(&path).unwrap(), "preserve");
            std::fs::remove_file(path).unwrap();
        }
        let report = dir.path().join("residuals.tsv");
        training["outReport"] = report.display().to_string().into();
        #[cfg(unix)]
        {
            std::os::unix::fs::symlink(dir.path().join("missing"), &report).unwrap();
            assert!(validate(training).is_err());
        }
    }

    #[cfg(unix)]
    fn fake_cli(script: &str) -> (tempfile::TempDir, Resolved) {
        use std::os::unix::fs::PermissionsExt;
        let dir = private_tempdir().unwrap();
        let bin = dir.path().join("DIALibGen");
        std::fs::write(&bin, format!("#!/bin/sh\n{script}\n")).unwrap();
        std::fs::set_permissions(&bin, std::fs::Permissions::from_mode(0o700)).unwrap();
        (dir, Resolved { bin, source: "env", data: None, share: None })
    }

    #[cfg(unix)]
    #[test]
    fn health_probe_requires_success_and_a_version_banner() {
        for (script, ok) in [
            ("echo 'Version: 0.11.0'; exit 0", true),
            ("echo 'Version: 0.11.0'; exit 13", false),
            ("echo 'dyld: missing library' >&2; exit 1", false),
            ("echo 'some other program'; exit 0", false),
        ] {
            let (_dir, resolved) = fake_cli(script);
            assert_eq!(probe_resolved(&resolved).ok, ok, "{script}");
        }
        let missing = Resolved { bin: PathBuf::from("/nonexistent/dialibgen"), source: "env", data: None, share: None };
        assert!(!probe_resolved(&missing).ok);
    }

    #[cfg(unix)]
    #[test]
    fn effective_config_requires_success_and_valid_json() {
        let (_dir, mut resolved) = fake_cli("printf '{\"nce\":30,\"data\":\"%s\"}' \"$OPENMS_DATA_PATH\" > \"$4\"");
        resolved.data = Some(PathBuf::from("/data with spaces"));
        assert_eq!(default_config_resolved(&resolved, Mode::Generate).unwrap(), serde_json::json!({"nce":30,"data":"/data with spaces"}));
        for script in ["exit 0", "echo broken > \"$4\"", "echo '[]' > \"$4\"", "echo '{}' > \"$4\"; exit 13"] {
            let (_dir, resolved) = fake_cli(script);
            assert!(default_config_resolved(&resolved, Mode::Generate).is_err(), "{script}");
        }
    }

    #[test]
    fn config_directories_are_unique_private_and_removed_on_drop() {
        let first = private_tempdir().unwrap();
        let second = private_tempdir().unwrap();
        assert_ne!(first.path(), second.path());
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            assert_eq!(first.path().metadata().unwrap().permissions().mode() & 0o777, 0o700);
        }
        let path = first.path().to_path_buf();
        std::fs::write(path.join("config.json"), "{}").unwrap();
        drop(first);
        assert!(!path.exists());
        assert!(second.path().exists());
    }

    #[test]
    fn run_config_never_overwrites_an_existing_file() {
        let dir = private_tempdir().unwrap();
        let config = obj(r#"{"enzyme":"Lys-C","label":"quoted \\\" μ"}"#);
        let path = write_run_config(dir.path(), &config).unwrap();
        assert_eq!(read_config(path.display().to_string()).unwrap(), Value::Object(config));
        let before = std::fs::read(&path).unwrap();
        assert!(write_run_config(dir.path(), &Map::new()).is_err());
        assert_eq!(std::fs::read(&path).unwrap(), before);
    }

    #[cfg(unix)]
    #[test]
    fn run_config_refuses_a_symlink_without_touching_its_target() {
        let dir = private_tempdir().unwrap();
        let victim = dir.path().join("protected.json");
        std::fs::write(&victim, "keep this file").unwrap();
        let work = dir.path().join("run");
        std::fs::create_dir(&work).unwrap();
        std::os::unix::fs::symlink(&victim, work.join("config.json")).unwrap();
        assert!(write_run_config(&work, &Map::new()).is_err());
        assert_eq!(std::fs::read_to_string(&victim).unwrap(), "keep this file");
    }

    #[test]
    fn release_csp_allows_local_ipc_without_remote_or_inline_scripts() {
        let config: Value = serde_json::from_str(include_str!("../tauri.conf.json")).unwrap();
        let csp = config["app"]["security"]["csp"].as_str().unwrap();
        assert!(csp.contains("connect-src ipc: http://ipc.localhost"));
        assert!(!csp.contains("unsafe-inline"));
        assert!(!csp.contains("unsafe-eval"));
        assert!(!csp.contains('*'));
    }

    #[cfg(unix)]
    #[test]
    fn ipc_generation_run_checks_inputs_streams_results_and_cancels() {
        use tauri::{Listener, WebviewWindowBuilder};
        use tauri::test::{get_ipc_response, mock_builder, mock_context, noop_assets};
        use std::time::Duration;
        // One end-to-end command test owns this process-wide override. Other
        // tests use explicit Resolved paths and never change the environment.
        struct RestoreEnv(&'static str, Option<std::ffi::OsString>);
        impl Drop for RestoreEnv {
            fn drop(&mut self) {
                if let Some(value) = &self.1 { std::env::set_var(self.0, value); }
                else { std::env::remove_var(self.0); }
            }
        }
        let (dir, resolved) = fake_cli(r#"
work=${0%/*}
printf '%s\n' "$@" >> "$work/calls.txt"
if [ "$1" = --help ]; then echo 'DIALibGen Version: 0.11.0'; exit 0; fi
if [ "$1" = -write_ini ]; then cp "$work/schema.ini" "$2"; exit 0; fi
if [ "$3" = -write_config ]; then
  if [ -f "$work/fail-config" ]; then echo 'config failed' >&2; exit 12; fi
  if [ -f "$work/pause-config" ]; then
    touch "$work/config-entered"
    i=0
    while [ -f "$work/pause-config" ] && [ "$i" -lt 1000 ]; do sleep 0.01; i=$((i+1)); done
  fi
  printf '{"instrument":"QE","nce":30,"precursor_charges":[2,3],"requested_mode":"%s","rt_model":"/stale/rt.onnx","ms2_model":"/stale/ms2.onnx","ccs_model":"/stale/ccs.onnx"}' "$2" > "$4"
  exit 0
fi
printf '%s\n' "$@" > "$work/argv.txt"
while [ "$#" -gt 0 ]; do
  case "$1" in
    -config) cp "$2" "$work/received.json"; printf '%s' "$2" > "$work/config-path.txt"; shift;;
    -out) out=$2; shift;;
  esac
  shift
done
printf '%s' "$DIALIBGEN_MODEL_DIR" > "$work/model-env.txt"
echo 'stdout log'
echo 'stderr log' >&2
if [ -f "$work/wait" ]; then exec sleep 30; fi
if [ -f "$work/fail-run" ]; then exit 17; fi
printf 'library' > "$out"
"#);
        std::fs::write(dir.path().join("schema.ini"), TOPP_INI).unwrap();
        let _restore = RestoreEnv("DIALIBGEN_BIN", std::env::var_os("DIALIBGEN_BIN"));
        let _restore_models = RestoreEnv("DIALIBGEN_MODEL_DIR", std::env::var_os("DIALIBGEN_MODEL_DIR"));
        let _restore_data = RestoreEnv("OPENMS_DATA_PATH", std::env::var_os("OPENMS_DATA_PATH"));
        std::env::set_var("DIALIBGEN_BIN", &resolved.bin);
        let app = mock_builder().manage(RunManager::default())
            .invoke_handler(tauri::generate_handler![probe, models, default_config, tuning_options, run, cancel, read_config, write_config])
            .build(mock_context(noop_assets())).unwrap();
        let window = WebviewWindowBuilder::new(&app, "main", Default::default()).build().unwrap();
        let invoke = |command: &str, body: Value| {
            get_ipc_response(&window, tauri::webview::InvokeRequest {
                cmd: command.into(), callback: tauri::ipc::CallbackFn(0), error: tauri::ipc::CallbackFn(1),
                url: "tauri://localhost".parse().unwrap(), body: tauri::ipc::InvokeBody::Json(body),
                headers: Default::default(), invoke_key: tauri::test::INVOKE_KEY.into(),
            }).map(|response| response.deserialize::<Value>().unwrap())
        };
        assert_eq!(invoke("probe", serde_json::json!({})).unwrap()["ok"], true);
        assert_eq!(invoke("default_config", serde_json::json!({})).unwrap()["nce"], 30);
        for mode in ["generate", "refine", "tune"] {
            assert_eq!(invoke("default_config", serde_json::json!({"mode":mode})).unwrap()["requested_mode"], mode);
        }
        assert!(invoke("default_config", serde_json::json!({"mode":"-force"})).is_err());
        assert_eq!(invoke("tuning_options", serde_json::json!({})).unwrap().as_array().unwrap().len(), 28);
        let model_dir = dir.path().join("models with spaces");
        std::fs::create_dir(&model_dir).unwrap();
        let model_query = serde_json::json!({"dir": model_dir});
        assert_eq!(invoke("models", model_query.clone()).unwrap()["missing"].as_array().unwrap().len(), 3);
        for name in MODEL_FILES { std::fs::write(model_dir.join(name), "model").unwrap(); }
        assert_eq!(invoke("models", model_query).unwrap()["missing"], serde_json::json!([]));
        std::env::set_var("DIALIBGEN_MODEL_DIR", &model_dir);
        let automatic = invoke("models", serde_json::json!({})).unwrap();
        assert_eq!(automatic["dir"], model_dir.display().to_string());
        assert_eq!(automatic["missing"], serde_json::json!([]));
        let incomplete = dir.path().join("incomplete models");
        std::fs::create_dir(&incomplete).unwrap();
        std::fs::write(incomplete.join(MODEL_FILES[0]), "model").unwrap();
        // An explicit choice must report its missing files rather than silently
        // substituting a complete directory from the environment.
        let explicit = invoke("models", serde_json::json!({"dir":incomplete})).unwrap();
        assert_eq!(explicit["dir"], incomplete.display().to_string());
        assert_eq!(explicit["missing"], serde_json::json!([MODEL_FILES[1], MODEL_FILES[2]]));
        assert_eq!(invoke("models", serde_json::json!({"dir":incomplete,"heads":["rt"]})).unwrap()["missing"], serde_json::json!([]));
        assert_eq!(invoke("models", serde_json::json!({"dir":incomplete,"heads":["ccs"]})).unwrap()["missing"], serde_json::json!([MODEL_FILES[2]]));
        assert!(invoke("models", serde_json::json!({"heads":["-force"]})).is_err());
        assert!(invoke("models", serde_json::json!({"heads":[]})).is_err());
        let config_path = dir.path().join("picked config.json");
        let picked = serde_json::json!({"instrument":"Lumos"});
        assert_eq!(invoke("write_config", serde_json::json!({"path":config_path,"value":picked})).unwrap(), true);
        assert_eq!(invoke("read_config", serde_json::json!({"path":config_path})).unwrap(), picked);
        assert!(invoke("read_config", serde_json::json!({"path":dir.path().join("missing.json")})).is_err());
        let output = dir.path().join("output with spaces.tsv");
        let mut params = serde_json::json!({"in":"-fasta with spaces","out":output,
            "threads":3,"modelDir":model_dir,"config":{"instrument":"Lumos","nce":null,"unknown":"drop","rt_model":"/imported/rt.onnx","ms2_model":"/imported/ms2.onnx","ccs_model":"/imported/ccs.onnx"}});
        for (key, value, reason) in [("in", "", "required"), ("out", "output.TSV", "must end")] {
            let mut invalid = params.clone(); invalid[key] = value.into();
            let result = invoke("run", serde_json::json!({"params":invalid})).unwrap();
            assert_eq!(result["started"], false);
            assert!(result["reason"].as_str().unwrap().contains(reason));
        }
        let mut missing_model = params.clone();
        missing_model["modelDir"] = incomplete.display().to_string().into();
        assert_eq!(invoke("run", serde_json::json!({"params":missing_model})).unwrap()["started"], false);
        std::fs::remove_file(model_dir.join(MODEL_FILES[2])).unwrap();
        assert_eq!(invoke("run", serde_json::json!({"params":params})).unwrap()["started"], false);
        std::fs::write(model_dir.join(MODEL_FILES[2]), "model").unwrap();
        let calls_before = std::fs::read(dir.path().join("calls.txt")).unwrap();
        std::fs::write(&output, "keep existing library").unwrap();
        let result = invoke("run", serde_json::json!({"params":params})).unwrap();
        assert_eq!(result["started"], false);
        assert!(result["reason"].as_str().unwrap().contains("choose an unused file name"));
        assert_eq!(std::fs::read_to_string(&output).unwrap(), "keep existing library");
        assert_eq!(std::fs::read(dir.path().join("calls.txt")).unwrap(), calls_before);
        assert!(app.state::<RunManager>().current.lock().unwrap().is_none());
        std::fs::remove_file(&output).unwrap();
        // A dangling symlink is also an occupied output name.
        std::os::unix::fs::symlink(dir.path().join("missing-target"), &output).unwrap();
        assert_eq!(invoke("run", serde_json::json!({"params":params})).unwrap()["started"], false);
        assert_eq!(std::fs::read(dir.path().join("calls.txt")).unwrap(), calls_before);
        std::fs::remove_file(&output).unwrap();
        let mut fractional_threads = params.clone();
        fractional_threads["threads"] = serde_json::json!(2.5);
        assert!(invoke("run", serde_json::json!({"params":fractional_threads})).is_err());
        assert_eq!(std::fs::read(dir.path().join("calls.txt")).unwrap(), calls_before);
        {
            let state = app.state::<RunManager>();
            let _launch = state.launch.lock().unwrap();
            assert_eq!(invoke("run", serde_json::json!({"params":params})).unwrap()["reason"], "a run is already starting");
        }
        std::fs::write(dir.path().join("fail-config"), "").unwrap();
        assert_eq!(invoke("run", serde_json::json!({"params":params})).unwrap()["started"], false);
        std::fs::remove_file(dir.path().join("fail-config")).unwrap();
        let (tx, rx) = std::sync::mpsc::channel();
        app.listen("dialibgen:done", move |event| { let _ = tx.send(serde_json::from_str::<Value>(event.payload()).unwrap()); });
        let (log_tx, logs) = std::sync::mpsc::channel();
        app.listen("dialibgen:log", move |event| { let _ = log_tx.send(serde_json::from_str::<String>(event.payload()).unwrap()); });
        assert_eq!(invoke("run", serde_json::json!({"params":params})).unwrap()["started"], true);
        let done = rx.recv_timeout(Duration::from_secs(10)).unwrap();
        assert_eq!(done, serde_json::json!({"ok":true,"code":0,"bytes":7}));
        assert_eq!(std::fs::read_to_string(&output).unwrap(), "library");
        let received = read_config(dir.path().join("received.json").display().to_string()).unwrap();
        assert_eq!(received["instrument"], "Lumos");
        assert!(received.get("nce").is_none());
        assert!(received.get("unknown").is_none());
        for (key, file) in ["rt_model", "ms2_model", "ccs_model"].into_iter().zip(MODEL_FILES) {
            assert_eq!(received[key], model_dir.join(file).canonicalize().unwrap().display().to_string());
        }
        let argv = std::fs::read_to_string(dir.path().join("argv.txt")).unwrap();
        let args: Vec<_> = argv.lines().collect();
        assert_eq!(&args[..4], &["-mode", "generate", "-in", "./-fasta with spaces"]);
        assert_eq!(&args[args.len()-2..], &["-threads", "3"]);
        assert_eq!(std::fs::read_to_string(dir.path().join("model-env.txt")).unwrap(), model_dir.canonicalize().unwrap().display().to_string());
        let temporary = std::fs::read_to_string(dir.path().join("config-path.txt")).unwrap();
        assert!(!Path::new(&temporary).exists());
        let mut lines: Vec<_> = logs.try_iter().collect(); lines.sort();
        assert_eq!(lines, ["stderr log", "stdout log"]);
        assert!(app.state::<RunManager>().current.lock().unwrap().is_none());
        assert_eq!(invoke("cancel", serde_json::json!({})).unwrap()["cancelled"], false);

        // Readiness may choose a complete data-directory fallback after an
        // incomplete environment override. Launch must use that exact choice.
        let data = dir.path().join("fallback data");
        let fallback_models = data.join("models");
        std::fs::create_dir_all(&fallback_models).unwrap();
        for file in MODEL_FILES { std::fs::write(fallback_models.join(file), "fallback model").unwrap(); }
        std::env::set_var("DIALIBGEN_MODEL_DIR", &incomplete);
        std::env::set_var("OPENMS_DATA_PATH", &data);
        assert_eq!(invoke("models", serde_json::json!({})).unwrap()["dir"], fallback_models.display().to_string());
        for mode in ["generate", "tune"] {
            let mut fallback = params.clone();
            fallback["mode"] = mode.into();
            fallback["modelDir"] = Value::Null;
            fallback["out"] = dir.path().join(format!("fallback-{mode}.tsv")).display().to_string().into();
            if mode == "tune" {
                fallback["in"] = "library.tsv".into();
                fallback["ids"] = "report.parquet".into();
                fallback["config"] = serde_json::json!({});
            }
            assert_eq!(invoke("run", serde_json::json!({"params":fallback})).unwrap()["started"], true);
            assert_eq!(rx.recv_timeout(Duration::from_secs(10)).unwrap()["ok"], true);
            assert_eq!(std::fs::read_to_string(dir.path().join("model-env.txt")).unwrap(), fallback_models.canonicalize().unwrap().display().to_string());
            if mode == "generate" {
                let config = read_config(dir.path().join("received.json").display().to_string()).unwrap();
                for (key, file) in ["rt_model", "ms2_model", "ccs_model"].into_iter().zip(MODEL_FILES) {
                    assert_eq!(config[key], fallback_models.join(file).canonicalize().unwrap().display().to_string());
                }
            }
        }
        std::env::set_var("DIALIBGEN_MODEL_DIR", &model_dir);

        for mode in ["refine", "tune"] {
            let trained = serde_json::json!({"mode":mode,"in":"-library with spaces.tsv", "ids":"-force.parquet",
                "out":dir.path().join(format!("{mode}.tsv")), "outReport":"-residual report.tsv",
                "tune":mode == "refine", "tuneOutModels":"-models with spaces",
                "tuning":{"tune_heads":"rt","train:epochs":25,"cohort:full_fit":true,"machine:no_cudnn":false}});
            assert_eq!(invoke("run", serde_json::json!({"params":trained})).unwrap()["started"], true);
            assert_eq!(rx.recv_timeout(Duration::from_secs(10)).unwrap()["ok"], true);
            let arguments = std::fs::read_to_string(dir.path().join("argv.txt")).unwrap();
            let arguments: Vec<_> = arguments.lines().collect();
            assert_eq!(&arguments[..4], &["-mode", mode, "-in", "./-library with spaces.tsv"]);
            for pair in [["-ids", "./-force.parquet"], ["-out_report", "./-residual report.tsv"],
                         ["-tune_out_models", "./-models with spaces"], ["-train:epochs", "25"], ["-tune_heads", "rt"]] {
                assert!(arguments.windows(2).any(|args| args == pair));
            }
            assert_eq!(arguments.contains(&"-tune"), mode == "refine");
            assert!(arguments.contains(&"-cohort:full_fit"));
            assert!(!arguments.contains(&"-machine:no_cudnn"));
            let mut invalid = trained;
            invalid["out"] = dir.path().join("invalid-tuning.tsv").display().to_string().into();
            invalid["tuning"]["machine:device"] = "-force".into();
            let before = std::fs::read(dir.path().join("argv.txt")).unwrap();
            assert_eq!(invoke("run", serde_json::json!({"params":invalid})).unwrap()["started"], false);
            assert_eq!(std::fs::read(dir.path().join("argv.txt")).unwrap(), before);
        }

        std::fs::write(dir.path().join("fail-run"), "").unwrap();
        params["out"] = dir.path().join("failed.tsv").display().to_string().into();
        assert_eq!(invoke("run", serde_json::json!({"params":params})).unwrap()["started"], true);
        assert_eq!(rx.recv_timeout(Duration::from_secs(10)).unwrap(), serde_json::json!({"ok":false,"code":17,"bytes":null}));
        std::fs::remove_file(dir.path().join("fail-run")).unwrap();
        std::fs::write(dir.path().join("wait"), "").unwrap();
        assert_eq!(invoke("run", serde_json::json!({"params":params})).unwrap()["started"], true);
        assert_eq!(invoke("run", serde_json::json!({"params":params})).unwrap()["reason"], "a run is already in progress");
        assert_eq!(invoke("cancel", serde_json::json!({})).unwrap()["cancelled"], true);
        assert_eq!(rx.recv_timeout(Duration::from_secs(10)).unwrap()["ok"], false);
        assert!(app.state::<RunManager>().current.lock().unwrap().is_none());
        std::fs::remove_file(dir.path().join("wait")).unwrap();

        // Block the real native preflight helper, then cancel/close while no
        // generation child exists. Neither path may launch after it returns.
        for closing in [false, true] {
            let pause = dir.path().join("pause-config");
            let entered = dir.path().join("config-entered");
            let _ = std::fs::remove_file(&entered);
            std::fs::write(&pause, "").unwrap();
            let before = std::fs::read(dir.path().join("argv.txt")).unwrap();
            let handle = app.handle().clone();
            let pending_params = serde_json::from_value(params.clone()).unwrap();
            let pending = std::thread::spawn(move || run(handle.clone(), handle.state::<RunManager>(), pending_params));
            let deadline = std::time::Instant::now() + Duration::from_secs(5);
            while !entered.exists() && std::time::Instant::now() < deadline {
                std::thread::sleep(Duration::from_millis(10));
            }
            let preflight_entered = entered.exists();
            if closing { shutdown(&app.state::<RunManager>()); }
            else { assert_eq!(invoke("cancel", serde_json::json!({})).unwrap()["cancelled"], true); }
            std::fs::remove_file(&pause).unwrap();
            let result = pending.join().unwrap();
            assert!(preflight_entered, "the CLI helper must enter preflight before cancellation");
            assert!(!result.started);
            assert_eq!(result.reason.as_deref(), Some("run cancelled before launch"));
            assert_eq!(std::fs::read(dir.path().join("argv.txt")).unwrap(), before);
            assert!(!Path::new(params["out"].as_str().unwrap()).exists());
            assert!(app.state::<RunManager>().current.lock().unwrap().is_none());
            assert!(rx.try_recv().is_err());
            if !closing {
                let mut retry = params.clone();
                retry["out"] = dir.path().join("after-cancel.tsv").display().to_string().into();
                assert_eq!(invoke("run", serde_json::json!({"params":retry})).unwrap()["started"], true);
                assert_eq!(rx.recv_timeout(Duration::from_secs(10)).unwrap()["ok"], true);
            }
        }
        assert_eq!(invoke("run", serde_json::json!({"params":params})).unwrap()["reason"], "the app is closing");
    }

    // build_config is the trust boundary: the tool REFUSES an unknown config
    // key outright, so a stray one from the renderer would fail every run.
    #[test]
    fn config_keeps_only_declared_keys() {
        let reference = obj(r#"{"enzyme":"Trypsin/P","missed_cleavages":2,"decoys":"none"}"#);
        let wanted = obj(r#"{"enzyme":"Trypsin","no_such_key":1,"decoys":"reverse"}"#);
        let got = build_config(&reference, &wanted);
        assert_eq!(got.get("enzyme").unwrap(), "Trypsin");
        assert_eq!(got.get("decoys").unwrap(), "reverse");
        assert!(!got.contains_key("no_such_key"));
        assert_eq!(got.len(), 2);
    }

    #[test]
    fn config_drops_nulls() {
        // A cleared form field arrives as null; writing it would mean "set this
        // key to null", which the tool's typed reader rejects.
        let reference = obj(r#"{"nce":30.0,"instrument":"QE"}"#);
        let got = build_config(&reference, &obj(r#"{"nce":null,"instrument":"timsTOF"}"#));
        assert!(!got.contains_key("nce"));
        assert_eq!(got.len(), 1);
    }

    #[test]
    fn config_preserves_arrays_and_nesting() {
        let reference = obj(r#"{"precursor_charges":[1,2],"peptide_length":[7,30]}"#);
        let got = build_config(&reference, &obj(r#"{"precursor_charges":[2,3,4],"peptide_length":[8,25]}"#));
        assert_eq!(got.get("precursor_charges").unwrap().as_array().unwrap().len(), 3);
        assert_eq!(got.get("peptide_length").unwrap()[0], 8);
    }

    #[test]
    fn dash_leading_paths_are_defanged_not_dropped() {
        assert_eq!(defang_path("-force"), "./-force");
        assert_eq!(defang_path("./-force"), "./-force");
        assert_eq!(defang_path("-5.fasta"), "-5.fasta"); // digit after the dash
        assert_eq!(defang_path("/tmp/a.fasta"), "/tmp/a.fasta");
    }

    #[test]
    fn version_parses_from_help_banner() {
        assert_eq!(
            parse_version("DIALibGen -- blah\nVersion: 0.2.0 Sep 11 2026"),
            Some("0.2.0".to_string())
        );
        assert_eq!(parse_version("no banner here"), None);
    }

    #[test]
    fn config_file_roundtrip_and_invalid_inputs() {
        let dir = std::env::temp_dir().join(format!("dialibgen-config-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("library config.json").display().to_string();
        let value = serde_json::json!({"enzyme": "Lys-C", "charges": [2, 3], "label": "μ sample"});
        assert!(write_config(path.clone(), value.clone()).unwrap());
        assert_eq!(read_config(path.clone()).unwrap(), value);
        std::fs::write(&path, "[1,2]").unwrap();
        assert!(read_config(path.clone()).unwrap_err().contains("not a JSON object"));
        std::fs::write(&path, "{broken json").unwrap();
        assert!(read_config(path.clone()).is_err());
        std::fs::remove_file(&path).unwrap();
        assert!(read_config(path).is_err());
        assert!(write_config(dir.display().to_string(), value).is_err());
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn cancellation_kills_child_and_leaves_reaping_to_coordinator() {
        let state = RunManager::default();
        assert_eq!(cancel_current(&state), serde_json::json!({"cancelled": false}));
        #[cfg(windows)]
        let child = Command::new("powershell.exe")
            .args(["-NoProfile", "-Command", "Start-Sleep -Seconds 30"])
            .stdout(Stdio::null()).stderr(Stdio::null()).spawn().unwrap();
        #[cfg(not(windows))]
        let child = Command::new("sleep").arg("30").spawn().unwrap();
        let child = Arc::new(Mutex::new(child));
        let dir = private_tempdir().unwrap();
        let path = dir.path().to_path_buf();
        *state.current.lock().unwrap() = Some(CurrentRun { run_id: 1, child: child.clone(), config_dir: Some(dir) });
        assert_eq!(cancel_current(&state), serde_json::json!({"cancelled": true}));
        assert!(!child.lock().unwrap().wait().unwrap().success());
        assert!(state.current.lock().unwrap().is_some());
        shutdown(&state);
        assert!(state.current.lock().unwrap().is_none());
        assert!(!path.exists());
    }
}
