// Named config presets + last-used state, persisted as one JSON file in the app
// config dir. Writes are atomic (temp + rename). A corrupt file is moved aside
// once rather than silently overwritten, so recoverable presets survive.

use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use std::path::PathBuf;
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

pub fn load(app: &AppHandle) -> Settings {
    let Some(path) = settings_file(app) else { return Settings::default() };
    if !path.exists() {
        return Settings::default();
    }
    let Ok(raw) = std::fs::read_to_string(&path) else { return Settings::default() };
    match serde_json::from_str::<Value>(&raw) {
        Ok(v) => {
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
        Err(_) => {
            // Move the unparseable file aside ONCE so the next save() does not
            // destroy recoverable presets; then start from empty.
            let mut corrupt = path.clone();
            if let Some(name) = path.file_name().map(|n| n.to_string_lossy().to_string()) {
                corrupt.set_file_name(format!("{name}.corrupt"));
                if !corrupt.exists() {
                    let _ = std::fs::rename(&path, &corrupt);
                }
            }
            Settings::default()
        }
    }
}

fn save(app: &AppHandle, s: &Settings) -> bool {
    let Some(path) = settings_file(app) else { return false };
    if let Some(dir) = path.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    let Ok(json) = serde_json::to_string_pretty(s) else { return false };
    let mut tmp = path.clone();
    tmp.set_extension("json.tmp");
    if std::fs::write(&tmp, json).is_err() {
        return false;
    }
    std::fs::rename(&tmp, &path).is_ok()
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
}
