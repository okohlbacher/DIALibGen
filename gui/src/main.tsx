import React from 'react'
import { createRoot } from 'react-dom/client'
import './api'
import App from './App'
import './index.css'

createRoot(document.getElementById('root') as HTMLElement).render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
)
