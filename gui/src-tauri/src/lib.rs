mod dialibgen;
mod settings;

use dialibgen::RunManager;
use tauri::{Manager, WindowEvent};

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
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
