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
pub fn models<R: tauri::Runtime>(app: AppHandle<R>, dir: Option<String>) -> ModelStatus {
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
#[tauri::command(async)]
pub fn default_config<R: tauri::Runtime>(app: AppHandle<R>) -> Result<Value, String> {
    let r = resolve_binary(&app);
    default_config_resolved(&r)
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

fn default_config_resolved(r: &Resolved) -> Result<Value, String> {
    let dir = private_tempdir().map_err(|e| e.to_string())?;
    let path = dir.path().join("effective.json");
    let mut cmd = Command::new(&r.bin);
    cmd.args(["-mode", "generate", "-write_config"]).arg(&path);
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

#[tauri::command(async)]
pub fn run<R: tauri::Runtime>(app: AppHandle<R>, state: State<'_, RunManager>, params: RunParams) -> RunStarted {
    let refuse = |reason: String| RunStarted { started: false, reason: Some(reason) };
    let Ok(_launch) = state.launch.try_lock() else {
        return refuse("a run is already starting".into());
    };
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
    let dir = match private_tempdir() {
        Ok(dir) => dir,
        Err(e) => return refuse(format!("cannot create a temporary config directory: {e}")),
    };
    let config_file = match write_run_config(dir.path(), &config) {
        Ok(path) => path,
        Err(e) => return refuse(format!("cannot write the config file: {e}")),
    };

    let r = resolve_binary(&app);
    let mut cmd = Command::new(&r.bin);
    cmd.args(["-mode", "generate", "-in"])
        .arg(defang_path(params.input.trim()))
        .arg("-config")
        .arg(&config_file)
        .arg("-out")
        .arg(defang_path(params.out.trim()))
        .arg("-threads")
        .arg(params.threads.unwrap_or(1).max(0).to_string())
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
            return refuse(format!("cannot start: {e}"));
        }
    };

    let run_id = state.seq.fetch_add(1, Ordering::SeqCst) + 1;
    let stdout = child.stdout.take();
    let stderr = child.stderr.take();
    let child_arc = Arc::new(Mutex::new(child));
    *state.current.lock().unwrap() =
        Some(CurrentRun { run_id, child: child_arc.clone(), config_dir: Some(dir) });

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
    if let Some(run) = cur.as_ref() {
        if let Ok(mut child) = run.child.lock() {
            let _ = child.kill(); // the coordinator reaps, cleans up and emits done
        }
        serde_json::json!({ "cancelled": true })
    } else {
        serde_json::json!({ "cancelled": false })
    }
}

/// Finish cancellation and release temporary files before the app exits.
pub fn shutdown(state: &RunManager) {
    if let Some(run) = state.current.lock().unwrap().take() {
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

    fn obj(s: &str) -> Map<String, Value> {
        match serde_json::from_str(s).unwrap() {
            Value::Object(m) => m,
            _ => panic!("not an object"),
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
        assert_eq!(default_config_resolved(&resolved).unwrap(), serde_json::json!({"nce":30,"data":"/data with spaces"}));
        for script in ["exit 0", "echo broken > \"$4\"", "echo '[]' > \"$4\"", "echo '{}' > \"$4\"; exit 13"] {
            let (_dir, resolved) = fake_cli(script);
            assert!(default_config_resolved(&resolved).is_err(), "{script}");
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
        struct RestoreBin(Option<std::ffi::OsString>);
        impl Drop for RestoreBin {
            fn drop(&mut self) {
                if let Some(value) = &self.0 { std::env::set_var("DIALIBGEN_BIN", value); }
                else { std::env::remove_var("DIALIBGEN_BIN"); }
            }
        }
        let (dir, resolved) = fake_cli(r#"
work=${0%/*}
if [ "$1" = --help ]; then echo 'DIALibGen Version: 0.11.0'; exit 0; fi
if [ "$3" = -write_config ]; then
  if [ -f "$work/fail-config" ]; then echo 'config failed' >&2; exit 12; fi
  printf '{"instrument":"QE","nce":30,"precursor_charges":[2,3]}' > "$4"
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
        let _restore = RestoreBin(std::env::var_os("DIALIBGEN_BIN"));
        std::env::set_var("DIALIBGEN_BIN", &resolved.bin);
        let app = mock_builder().manage(RunManager::default())
            .invoke_handler(tauri::generate_handler![probe, models, default_config, run, cancel, read_config, write_config])
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
        let model_dir = dir.path().join("models with spaces");
        std::fs::create_dir(&model_dir).unwrap();
        let model_query = serde_json::json!({"dir": model_dir});
        assert_eq!(invoke("models", model_query.clone()).unwrap()["missing"].as_array().unwrap().len(), 3);
        for name in MODEL_FILES { std::fs::write(model_dir.join(name), "model").unwrap(); }
        assert_eq!(invoke("models", model_query).unwrap()["missing"], serde_json::json!([]));
        let config_path = dir.path().join("picked config.json");
        let picked = serde_json::json!({"instrument":"Lumos"});
        assert_eq!(invoke("write_config", serde_json::json!({"path":config_path,"value":picked})).unwrap(), true);
        assert_eq!(invoke("read_config", serde_json::json!({"path":config_path})).unwrap(), picked);
        assert!(invoke("read_config", serde_json::json!({"path":dir.path().join("missing.json")})).is_err());
        let output = dir.path().join("output with spaces.tsv");
        let mut params = serde_json::json!({"in":"-fasta with spaces","out":output,
            "threads":3,"modelDir":model_dir,"config":{"instrument":"Lumos","nce":null,"unknown":"drop"}});
        for (key, value, reason) in [("in", "", "required"), ("out", "output.TSV", "must end")] {
            let mut invalid = params.clone(); invalid[key] = value.into();
            let result = invoke("run", serde_json::json!({"params":invalid})).unwrap();
            assert_eq!(result["started"], false);
            assert!(result["reason"].as_str().unwrap().contains(reason));
        }
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
        assert_eq!(read_config(dir.path().join("received.json").display().to_string()).unwrap(), picked);
        let argv = std::fs::read_to_string(dir.path().join("argv.txt")).unwrap();
        let args: Vec<_> = argv.lines().collect();
        assert_eq!(&args[..4], &["-mode", "generate", "-in", "./-fasta with spaces"]);
        assert_eq!(&args[args.len()-2..], &["-threads", "3"]);
        assert_eq!(std::fs::read_to_string(dir.path().join("model-env.txt")).unwrap(), model_dir.display().to_string());
        let temporary = std::fs::read_to_string(dir.path().join("config-path.txt")).unwrap();
        assert!(!Path::new(&temporary).exists());
        let mut lines: Vec<_> = logs.try_iter().collect(); lines.sort();
        assert_eq!(lines, ["stderr log", "stdout log"]);
        assert!(app.state::<RunManager>().current.lock().unwrap().is_none());
        assert_eq!(invoke("cancel", serde_json::json!({})).unwrap()["cancelled"], false);

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
