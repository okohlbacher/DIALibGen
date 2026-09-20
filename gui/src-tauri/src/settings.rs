// Named config presets + last-used state, persisted as one JSON file in the app
// config dir. Writes are atomic (temp + rename). A corrupt file is moved aside
// once rather than silently overwritten, so recoverable presets survive.

use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::path::{Path, PathBuf};
use tauri::{AppHandle, Manager};

#[derive(Serialize, Deserialize, Clone)]
#[serde(rename_all = "camelCase")]
pub struct Settings {
    schema_version: u32,
    last_used: Option<Value>,
    presets: Map<String, Value>,
}

impl Default for Settings {
    fn default() -> Self {
        Settings { schema_version: 1, last_used: None, presets: Map::new() }
    }
}

fn is_reserved(k: &str) -> bool {
    matches!(k, "__proto__" | "constructor" | "prototype")
}

// Drop keys that would poison Object.prototype once merged in the renderer.
fn sanitize_map(m: &Map<String, Value>) -> Map<String, Value> {
    m.iter().filter(|(k, _)| !is_reserved(k)).map(|(k, v)| (k.clone(), v.clone())).collect()
}

fn sanitize_value(v: Value) -> Value {
    match v {
        Value::Object(m) => Value::Object(sanitize_map(&m)),
        other => other,
    }
}

fn settings_file(app: &AppHandle) -> Option<PathBuf> {
    app.path().app_config_dir().ok().map(|d| d.join("dialibgen-settings.json"))
}

fn supported_settings(v: &Value) -> bool {
    v.is_object()
        && v.get("schemaVersion").is_none_or(|x| x == 1)
        && v.get("lastUsed").is_none_or(|x| x.is_null() || x.is_object())
        && v.get("presets").is_none_or(Value::is_object)
}

pub fn load(app: &AppHandle) -> Settings {
    let Some(path) = settings_file(app) else { return Settings::default() };
    load_path(&path)
}

fn load_path(path: &Path) -> Settings {
    if !path.exists() {
        return Settings::default();
    }
    let Ok(raw) = std::fs::read_to_string(path) else { return Settings::default() };
    match serde_json::from_str::<Value>(&raw) {
        Ok(v) if supported_settings(&v) => {
            let last_used =
                v.get("lastUsed").and_then(|x| x.as_object()).map(|m| Value::Object(sanitize_map(m)));
            let presets = v
                .get("presets")
                .and_then(|x| x.as_object())
                .map(|m| {
                    sanitize_map(m)
                        .into_iter()
                        .filter_map(|(k, val)| {
                            val.as_object().map(|vm| (k, Value::Object(sanitize_map(vm))))
                        })
                        .collect::<Map<String, Value>>()
                })
                .unwrap_or_default();
            Settings { schema_version: 1, last_used, presets }
        }
        _ => {
            // Preserve malformed or incompatible settings before starting
            // fresh. Never replace an existing recovery copy.
            let mut corrupt = path.to_path_buf();
            if let Some(name) = path.file_name().map(|n| n.to_string_lossy().to_string()) {
                corrupt.set_file_name(format!("{name}.corrupt"));
                if !corrupt.exists() {
                    let _ = std::fs::rename(path, &corrupt);
                }
            }
            Settings::default()
        }
    }
}

fn save(app: &AppHandle, s: &Settings) -> bool {
    let Some(path) = settings_file(app) else { return false };
    save_path(&path, s)
}

fn save_path(path: &Path, s: &Settings) -> bool {
    // If load could not preserve a future-format file because a recovery copy
    // already exists, do not overwrite it with this version's schema.
    if std::fs::read_to_string(path).ok()
        .and_then(|text| serde_json::from_str::<Value>(&text).ok())
        .and_then(|v| v.get("schemaVersion").cloned())
        .is_some_and(|version| version != 1)
    {
        return false;
    }
    if let Some(dir) = path.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let Ok(json) = serde_json::to_string_pretty(s) else { return false };
    let mut tmp = path.to_path_buf();
    tmp.set_extension("json.tmp");
    if std::fs::write(&tmp, json).is_err() {
        return false;
    }
    std::fs::rename(&tmp, path).is_ok()
}

#[tauri::command]
pub fn load_settings(app: AppHandle) -> Settings {
    load(&app)
}

#[tauri::command]
pub fn save_last(app: AppHandle, values: Value) -> bool {
    let mut s = load(&app);
    s.last_used = Some(sanitize_value(values));
    save(&app, &s)
}

#[tauri::command]
pub fn save_preset(app: AppHandle, name: String, values: Value) -> bool {
    let clean = name.trim();
    if clean.is_empty() || is_reserved(clean) {
        return false;
    }
    let mut s = load(&app);
    s.presets.insert(clean.to_string(), sanitize_value(values));
    save(&app, &s)
}

#[tauri::command]
pub fn delete_preset(app: AppHandle, name: String) -> bool {
    let mut s = load(&app);
    if s.presets.remove(&name).is_none() {
        return false;
    }
    save(&app, &s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicU64, Ordering};

    fn test_dir() -> PathBuf {
        static SEQUENCE: AtomicU64 = AtomicU64::new(0);
        let dir = std::env::temp_dir().join(format!(
            "dialibgen-settings-test-{}-{}", std::process::id(), SEQUENCE.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn reserved_keys_are_stripped() {
        let m: Map<String, Value> = serde_json::from_str(
            r#"{"__proto__":{"x":1},"constructor":2,"prototype":3,"enzyme":"Trypsin/P"}"#,
        )
        .unwrap();
        let s = sanitize_map(&m);
        assert_eq!(s.len(), 1);
        assert!(s.contains_key("enzyme"));
    }

    #[test]
    fn sanitize_value_recurses_one_level() {
        let v: Value = serde_json::from_str(r#"{"__proto__":1,"ok":true}"#).unwrap();
        let s = sanitize_value(v);
        assert!(s.get("__proto__").is_none());
        assert_eq!(s.get("ok"), Some(&Value::Bool(true)));
    }

    #[test]
    fn settings_roundtrip_replaces_existing_file_and_cleans_temporary_file() {
        let dir = test_dir();
        let path = dir.join("nested/settings.json");
        assert!(load_path(&path).presets.is_empty());
        let mut settings = Settings {
            last_used: Some(serde_json::json!({"enzyme": "Lys-C", "__modelDir": "a path/μ"})),
            ..Settings::default()
        };
        settings.presets.insert("my preset".into(), serde_json::json!({"charges": [2, 3]}));
        assert!(save_path(&path, &settings));
        let loaded = load_path(&path);
        assert_eq!(loaded.last_used, settings.last_used);
        assert_eq!(loaded.presets, settings.presets);
        assert!(save_path(&path, &Settings::default()));
        assert!(load_path(&path).last_used.is_none());
        assert!(!path.with_extension("json.tmp").exists());
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn corrupt_settings_are_preserved_once_before_recovery() {
        let dir = test_dir();
        let path = dir.join("settings.json");
        let backup = dir.join("settings.json.corrupt");
        std::fs::write(&path, "recoverable original").unwrap();
        assert!(load_path(&path).presets.is_empty());
        assert_eq!(std::fs::read_to_string(&backup).unwrap(), "recoverable original");
        std::fs::write(&path, "second corrupt file").unwrap();
        assert!(load_path(&path).presets.is_empty());
        assert_eq!(std::fs::read_to_string(&backup).unwrap(), "recoverable original");
        assert!(save_path(&path, &Settings::default()));
        assert!(load_path(&path).presets.is_empty());
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn loaded_settings_filter_reserved_keys_and_invalid_presets() {
        let dir = test_dir();
        let path = dir.join("settings.json");
        std::fs::write(&path, r#"{"lastUsed":{"__proto__":1,"enzyme":"Lys-C"},"presets":{"bad":7,"constructor":{},"good":{"prototype":8,"nce":30}}}"#).unwrap();
        let settings = load_path(&path);
        assert_eq!(settings.last_used.unwrap(), serde_json::json!({"enzyme":"Lys-C"}));
        assert_eq!(settings.presets.len(), 1);
        assert_eq!(settings.presets["good"], serde_json::json!({"nce":30}));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn incompatible_settings_are_preserved_and_never_overwrite_a_future_schema() {
        for text in ["[]", r#"{"presets":7}"#, r#"{"lastUsed":false}"#, r#"{"schemaVersion":2,"presets":{"future":{}}}"#] {
            let dir = test_dir();
            let path = dir.join("settings.json");
            std::fs::write(&path, text).unwrap();
            assert!(load_path(&path).presets.is_empty());
            assert_eq!(std::fs::read_to_string(dir.join("settings.json.corrupt")).unwrap(), text);
            let future = r#"{"schemaVersion":2,"presets":{"future":{}}}"#;
            std::fs::write(&path, future).unwrap();
            assert!(load_path(&path).presets.is_empty());
            assert!(!save_path(&path, &Settings::default()));
            assert_eq!(std::fs::read_to_string(&path).unwrap(), future);
            std::fs::remove_dir_all(dir).unwrap();
        }
    }
}
