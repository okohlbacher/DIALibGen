mod dialibgen;
mod settings;

use dialibgen::RunManager;
use tauri::{Manager, WindowEvent};

fn on_run_event<R: tauri::Runtime>(app: &tauri::AppHandle<R>, event: tauri::RunEvent) {
    // Application Quit can bypass window destruction; finish cleanup before
    // Tauri exits the process and drops the worker threads.
    if let tauri::RunEvent::Exit = event {
        dialibgen::shutdown(&app.state::<RunManager>());
    }
}

#[cfg_attr(mobile, tauri::mobile_entry_point)]
pub fn run() {
    tauri::Builder::default()
        // single-instance must be registered first: a second launch focuses the
        // existing window instead of racing over the settings file and run state.
        .plugin(tauri_plugin_single_instance::init(|app, _args, _cwd| {
            if let Some(w) = app.get_webview_window("main") {
                let _ = w.unminimize();
                let _ = w.set_focus();
            }
        }))
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_opener::init())
        .manage(RunManager::default())
        .invoke_handler(tauri::generate_handler![
            dialibgen::probe,
            dialibgen::models,
            dialibgen::default_config,
            dialibgen::tuning_options,
            dialibgen::run,
            dialibgen::cancel,
            dialibgen::read_config,
            dialibgen::write_config,
            settings::load_settings,
            settings::save_last,
            settings::save_preset,
            settings::delete_preset,
        ])
        .on_window_event(|window, event| {
            // Closing the window must not orphan a running CLI process -- a
            // library build can hold gigabytes and run for an hour.
            if let WindowEvent::Destroyed = event {
                let app = window.app_handle().clone();
                dialibgen::shutdown(&app.state::<RunManager>());
            }
        })
        .build(tauri::generate_context!())
        .expect("error while building tauri application")
        .run(on_run_event);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::process::{Command, Stdio};
    use std::sync::{Arc, Mutex};
    use tauri::test::{mock_builder, mock_context, noop_assets};

    #[test]
    fn application_exit_reaps_the_child_and_removes_its_config_without_a_window_event() {
        let app = mock_builder().manage(RunManager::default()).build(mock_context(noop_assets())).unwrap();
        #[cfg(windows)]
        let child = Command::new("powershell.exe")
            .args(["-NoProfile", "-Command", "Start-Sleep -Seconds 30"])
            .stdout(Stdio::null()).stderr(Stdio::null()).spawn().unwrap();
        #[cfg(not(windows))]
        let child = Command::new("sleep").arg("30")
            .stdout(Stdio::null()).stderr(Stdio::null()).spawn().unwrap();
        let child = Arc::new(Mutex::new(child));
        let dir = tempfile::tempdir().unwrap();
        let config = dir.path().join("config.json");
        std::fs::write(&config, "{}").unwrap();
        *app.state::<RunManager>().current.lock().unwrap() = Some(dialibgen::CurrentRun {
            run_id: 1, child: child.clone(), config_dir: Some(dir),
        });

        on_run_event(app.handle(), tauri::RunEvent::Ready);
        assert!(child.lock().unwrap().try_wait().unwrap().is_none());
        assert!(config.exists());
        on_run_event(app.handle(), tauri::RunEvent::Exit);
        assert!(!child.lock().unwrap().try_wait().unwrap().unwrap().success());
        assert!(app.state::<RunManager>().current.lock().unwrap().is_none());
        assert!(!config.parent().unwrap().exists());
        // A window-destruction callback and app-exit callback may both fire.
        on_run_event(app.handle(), tauri::RunEvent::Exit);
        assert!(app.state::<RunManager>().current.lock().unwrap().is_none());
    }
}
