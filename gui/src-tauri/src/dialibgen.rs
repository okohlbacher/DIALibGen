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
use std::sync::atomic::{AtomicU64, Ordering};
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

/// The three ONNX exports the tool needs. No tagged OpenMS release ships them,
/// so the GUI has to be able to say which one is missing and let the user point
/// at a directory holding all three.
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
    pub seq: AtomicU64,
}

pub struct CurrentRun {
    pub run_id: u64,
    pub child: Arc<Mutex<Child>>,
    /// Deleted when the run ends, whatever the outcome.
    pub config_file: Option<PathBuf>,
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

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RunParams {
    #[serde(rename = "in")]
    input: String,
    out: String,
    /// The library config, as the form produced it. Filtered against the tool's
    /// own effective config before it is written out.
    #[serde(default)]
    config: Map<String, Value>,
    /// Optional directory holding the three .onnx files; passed to the child as
    /// DIALIBGEN_MODEL_DIR. The tool searches its own locations when unset.
    #[serde(default)]
    model_dir: Option<String>,
    /// 0 = all cores, which is the tool's own default.
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
pub fn resolve_binary(app: &AppHandle) -> Resolved {
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

#[tauri::command]
pub fn probe(app: AppHandle) -> BinaryInfo {
    let r = resolve_binary(&app);
    let bin_s = r.bin.display().to_string();
    let mut cmd = Command::new(&r.bin);
    cmd.arg("--help");
    apply_env(&mut cmd, &r);
    match cmd.output() {
        Ok(out) => {
            let text = format!(
                "{}\n{}",
                String::from_utf8_lossy(&out.stdout),
                String::from_utf8_lossy(&out.stderr)
            );
            let version = parse_version(&text);
            let detail = match &version {
                Some(v) => format!("DIALibGen {v} ({})", r.source),
                None => format!("runs, version unknown ({})", r.source),
            };
            BinaryInfo { bin: bin_s, source: r.source.into(), ok: true, version, detail }
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
pub fn models(app: AppHandle, dir: Option<String>) -> ModelStatus {
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
    // Report the first directory that has ALL three: a split set is not usable,
    // and naming a directory that holds one of them would be misleading.
    let mut best: Option<(PathBuf, Vec<String>)> = None;
    for c in candidates {
        let missing: Vec<String> = MODEL_FILES
            .iter()
            .filter(|f| !c.join(f).exists())
            .map(|f| (*f).to_string())
            .collect();
        if missing.is_empty() {
            return ModelStatus { dir: Some(c.display().to_string()), missing };
        }
        if best.is_none() || missing.len() < best.as_ref().unwrap().1.len() {
            best = Some((c, missing));
        }
    }
    match best {
        Some((c, missing)) => ModelStatus { dir: Some(c.display().to_string()), missing },
        None => ModelStatus { dir: None, missing: MODEL_FILES.iter().map(|f| f.to_string()).collect() },
    }
}

/// The tool's own effective config: every default materialised. This is what
/// the form is built from, so the GUI can never offer a key the tool does not
/// have, nor default one differently from the CLI.
#[tauri::command]
pub fn default_config(app: AppHandle) -> Result<Value, String> {
    let r = resolve_binary(&app);
    let dir = std::env::temp_dir().join(format!("dialibgen-cfg-{}", std::process::id()));
    std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
    let path = dir.join("effective.json");
    let mut cmd = Command::new(&r.bin);
    cmd.arg("-write_config").arg(&path);
    apply_env(&mut cmd, &r);
    let out = cmd.output().map_err(|e| format!("cannot run the tool: {e}"))?;
    if !path.exists() {
        let _ = std::fs::remove_dir_all(&dir);
        return Err(format!(
            "-write_config produced nothing: {}",
            String::from_utf8_lossy(&out.stderr).trim()
        ));
    }
    let text = std::fs::read_to_string(&path).map_err(|e| e.to_string());
    let _ = std::fs::remove_dir_all(&dir);
    serde_json::from_str(&text?).map_err(|e| format!("cannot parse the effective config: {e}"))
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

/// OpenMS's ProgressLogger writes "Progress of 'xyz': 42.13 %" to stderr. There
/// is no machine-readable progress contract to parse, and inventing one would
/// mean a CLI change that exists only for this GUI.
fn parse_progress_line(line: &str) -> Option<(String, f64)> {
    let rest = line.trim().strip_prefix("Progress of ")?;
    let (label, tail) = rest.split_once(':')?;
    let pct = tail.trim().strip_suffix('%')?.trim().parse::<f64>().ok()?;
    Some((label.trim().trim_matches('\'').to_string(), pct))
}

fn read_stream<R: Read>(stream: R, app: &AppHandle, parse_progress: bool) {
    let reader = BufReader::new(stream);
    for line in reader.lines() {
        let Ok(line) = line else { break };
        if parse_progress {
            if let Some((label, pct)) = parse_progress_line(&line) {
                let _ = app.emit(
                    "dialibgen:progress",
                    serde_json::json!({ "label": label, "percent": pct }),
                );
                continue;
            }
        }
        let _ = app.emit("dialibgen:log", line);
    }
}

#[tauri::command]
pub fn run(app: AppHandle, state: State<'_, RunManager>, params: RunParams) -> RunStarted {
    let refuse = |reason: String| RunStarted { started: false, reason: Some(reason) };
    {
        if state.current.lock().unwrap().is_some() {
            return refuse("a run is already in progress".into());
        }
    }
    if params.input.trim().is_empty() || params.out.trim().is_empty() {
        return refuse("a FASTA and an output file are both required".into());
    }
    // The CLI refuses anything else, but saying so here costs one line and
    // saves the user a process start and a stack of log lines.
    let out_ok = params.out.ends_with(".parquet") || params.out.ends_with(".tsv");
    if !out_ok {
        return refuse("the output file must end in .parquet or .tsv".into());
    }

    let reference = match default_config(app.clone()) {
        Ok(Value::Object(m)) => m,
        Ok(_) => return refuse("the tool's effective config is not a JSON object".into()),
        Err(e) => return refuse(e),
    };
    let config = build_config(&reference, &params.config);

    // The config goes to the child as a FILE. Two reasons: the tool takes it
    // that way and embeds it verbatim in the Parquet output as the recipe, and
    // a file keeps values with spaces, quotes and unicode off the command line.
    let dir = std::env::temp_dir().join(format!(
        "dialibgen-run-{}-{}",
        std::process::id(),
        state.seq.load(Ordering::SeqCst)
    ));
    if std::fs::create_dir_all(&dir).is_err() {
        return refuse("cannot create a temporary directory for the config".into());
    }
    let config_file = dir.join("config.json");
    let written = std::fs::File::create(&config_file).and_then(|mut f| {
        f.write_all(serde_json::to_string_pretty(&Value::Object(config)).unwrap_or_default().as_bytes())
    });
    if written.is_err() {
        let _ = std::fs::remove_dir_all(&dir);
        return refuse("cannot write the config file".into());
    }

    let r = resolve_binary(&app);
    let mut cmd = Command::new(&r.bin);
    cmd.arg("-in")
        .arg(defang_path(params.input.trim()))
        .arg("-config")
        .arg(&config_file)
        .arg("-out")
        .arg(defang_path(params.out.trim()))
        .arg("-threads")
        .arg(params.threads.unwrap_or(0).max(0).to_string())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    apply_env(&mut cmd, &r);
    if let Some(d) = params.model_dir.as_deref().filter(|d| !d.trim().is_empty()) {
        cmd.env("DIALIBGEN_MODEL_DIR", d.trim());
    }

    let mut child = match cmd.spawn() {
        Ok(c) => c,
        // A spawn failure is known synchronously: report it as a non-start so
        // the frontend surfaces "could not start" and settles cleanly.
        Err(e) => {
            let _ = std::fs::remove_dir_all(&dir);
            return refuse(format!("cannot start: {e}"));
        }
    };

    let run_id = state.seq.fetch_add(1, Ordering::SeqCst) + 1;
    let stdout = child.stdout.take();
    let stderr = child.stderr.take();
    let child_arc = Arc::new(Mutex::new(child));
    *state.current.lock().unwrap() =
        Some(CurrentRun { run_id, child: child_arc.clone(), config_file: Some(config_file) });

    // One reader thread per stream. stdout is drained as log even though the
    // tool writes nothing there, so a full pipe can never block the child.
    let a_err = app.clone();
    let t_err = std::thread::spawn(move || {
        if let Some(s) = stderr {
            read_stream(s, &a_err, true);
        }
    });
    let a_out = app.clone();
    let t_out = std::thread::spawn(move || {
        if let Some(s) = stdout {
            read_stream(s, &a_out, false);
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
                if let Some(f) = cur.as_ref().and_then(|c| c.config_file.clone()) {
                    if let Some(d) = f.parent() {
                        let _ = std::fs::remove_dir_all(d);
                    }
                }
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
    let cur = state.current.lock().unwrap();
    if let Some(run) = cur.as_ref() {
        if let Ok(mut child) = run.child.lock() {
            let _ = child.kill(); // the coordinator reaps, cleans up and emits done
        }
        serde_json::json!({ "cancelled": true })
    } else {
        serde_json::json!({ "cancelled": false })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn obj(s: &str) -> Map<String, Value> {
        match serde_json::from_str(s).unwrap() {
            Value::Object(m) => m,
            _ => panic!("not an object"),
        }
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
    fn progress_line_parses() {
        assert_eq!(
            parse_progress_line("Progress of 'predicting MS2': 42.13 %"),
            Some(("predicting MS2".to_string(), 42.13))
        );
        assert_eq!(parse_progress_line("digest: 12 proteins"), None);
        assert_eq!(parse_progress_line("Progress of 'x': nope %"), None);
    }

    #[test]
    fn version_parses_from_help_banner() {
        assert_eq!(
            parse_version("DIALibGen -- blah\nVersion: 0.2.0 Sep 11 2026"),
            Some("0.2.0".to_string())
        );
        assert_eq!(parse_version("no banner here"), None);
    }
}
