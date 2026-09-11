import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// 1420 is the port tauri.conf.json's devUrl points at; `strictPort` makes a
// clash fail loudly instead of silently serving the app somewhere Tauri is not
// looking.
export default defineConfig({
  plugins: [react()],
  clearScreen: false,
  server: { port: 1420, strictPort: true }
})
