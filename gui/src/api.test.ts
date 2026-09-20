import { beforeEach, describe, expect, it, vi } from 'vitest'
import { invoke } from '@tauri-apps/api/core'
import { listen } from '@tauri-apps/api/event'
import { open, save } from '@tauri-apps/plugin-dialog'
import { revealItemInDir } from '@tauri-apps/plugin-opener'
import './api'

vi.mock('@tauri-apps/api/core', () => ({ invoke: vi.fn().mockResolvedValue(true) }))
vi.mock('@tauri-apps/api/event', () => ({ listen: vi.fn() }))
vi.mock('@tauri-apps/plugin-dialog', () => ({ open: vi.fn(), save: vi.fn() }))
vi.mock('@tauri-apps/plugin-opener', () => ({ revealItemInDir: vi.fn().mockResolvedValue(undefined) }))

beforeEach(() => vi.clearAllMocks())

describe('native bridge', () => {
  it('passes structured config and paths through named commands without shell interpolation', async () => {
    const api = window.dialibgen
    const config = { enzyme: 'Lys-C', fixed_modifications: ['Acetyl (Protein N-term)'] }
    const params = { in: 'a path/$(unsafe).fasta', out: 'out.tsv', config, threads: 1, modelDir: null }
    await api.probe(); await api.models(); await api.models('/custom models'); await api.defaultConfig()
    await api.readConfig('some config.json'); await api.writeConfig('saved.json', config)
    await api.run(params); await api.cancel(); await api.loadSettings()
    await api.saveLast(config); await api.savePreset('preset', config); await api.deletePreset('preset')
    expect(vi.mocked(invoke).mock.calls).toEqual([
      ['probe'], ['models', { dir: null }], ['models', { dir: '/custom models' }], ['default_config', { mode: 'generate' }],
      ['read_config', { path: 'some config.json' }], ['write_config', { path: 'saved.json', value: config }],
      ['run', { params }], ['cancel'], ['load_settings'], ['save_last', { values: config }],
      ['save_preset', { name: 'preset', values: config }], ['delete_preset', { name: 'preset' }]
    ])
  })

  it('passes selected mode, head requirements and tuning metadata through the native bridge', async () => {
    await window.dialibgen.defaultConfig('refine')
    await window.dialibgen.defaultConfig('tune')
    await window.dialibgen.models('/rt models', ['rt'])
    await window.dialibgen.tuningOptions()
    expect(vi.mocked(invoke).mock.calls).toEqual([
      ['default_config', { mode: 'refine' }], ['default_config', { mode: 'tune' }],
      ['models', { dir: '/rt models', heads: ['rt'] }], ['tuning_options']
    ])
  })

  it('returns only selected single paths and treats dialog cancellation as null', async () => {
    for (const picker of [window.dialibgen.pickFasta, window.dialibgen.pickLibrary, window.dialibgen.pickReport, window.dialibgen.pickDirectory, window.dialibgen.pickConfig]) {
      vi.mocked(open).mockResolvedValueOnce('/selected').mockResolvedValueOnce(null).mockResolvedValueOnce(['/multiple'])
      expect(await picker()).toBe('/selected')
      expect(await picker()).toBeNull()
      expect(await picker()).toBeNull()
    }
    for (const picker of [window.dialibgen.pickOutput, window.dialibgen.saveConfig]) {
      vi.mocked(save).mockResolvedValueOnce('/selected').mockResolvedValueOnce(null)
      expect(await picker('/suggested')).toBe('/selected')
      expect(await picker()).toBeNull()
    }
    expect(vi.mocked(open).mock.calls[0][0]).toMatchObject({ multiple: false, filters: [{ extensions: ['fasta', 'fa', 'faa'] }] })
    expect(vi.mocked(save).mock.calls[0][0]).toMatchObject({ defaultPath: '/suggested', filters: [{ extensions: ['parquet', 'tsv'] }] })
    expect(await window.dialibgen.revealPath('/library.parquet')).toBe(true)
    expect(revealItemInDir).toHaveBeenCalledWith('/library.parquet')
  })

  it('forwards event payloads and removes resolved subscriptions', async () => {
    const unlisten = vi.fn()
    vi.mocked(listen).mockResolvedValue(unlisten)
    const received = vi.fn()
    const off = window.dialibgen.onLog(received)
    await Promise.resolve()
    const handler = vi.mocked(listen).mock.calls[0][1]
    handler({ event: 'dialibgen:log', id: 1, payload: 'finished digest' })
    expect(received).toHaveBeenCalledWith('finished digest')
    off()
    expect(unlisten).toHaveBeenCalledOnce()
  })

  it('removes a pending subscription when the component unmounts before listen resolves', async () => {
    let complete!: (unlisten: () => void) => void
    vi.mocked(listen).mockImplementation(() => new Promise((resolve) => { complete = resolve }))
    const unlisten = vi.fn()
    const off = window.dialibgen.onLog(vi.fn())
    off()
    complete(unlisten)
    await Promise.resolve()
    expect(unlisten).toHaveBeenCalledOnce()
    const doneOff = window.dialibgen.onDone(vi.fn())
    expect(listen).toHaveBeenLastCalledWith('dialibgen:done', expect.any(Function))
    doneOff(); complete(vi.fn()); await Promise.resolve()
  })
})
