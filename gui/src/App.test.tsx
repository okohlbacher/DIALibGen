import { act, fireEvent, render, screen, waitFor } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { beforeEach, describe, expect, it, vi } from 'vitest'
import App, { stripPrivate, withExtension } from './App'
import { installMockBridge, SAMPLE_CONFIG, type MockBridge } from './testing/mockBridge'

let bridge: MockBridge

async function renderApp(opts: Parameters<typeof installMockBridge>[0] = {}, configure?: (b: MockBridge) => void): Promise<void> {
  bridge = installMockBridge(opts)
  configure?.(bridge)
  // The startup effect resolves three promises before the first paint; without
  // act() every one of them logs a "not wrapped in act(...)" warning that
  // buries the real output.
  await act(async () => {
    render(<App />)
  })
  await waitFor(() => expect(screen.getByText(/DIALibGen 0\.11\.0|cannot execute/)).toBeTruthy())
}

describe('withExtension', () => {
  it('swaps the format without making the user re-pick the file', () => {
    expect(withExtension('/a/lib.parquet', 'tsv')).toBe('/a/lib.tsv')
    expect(withExtension('/a/lib.tsv', 'parquet')).toBe('/a/lib.parquet')
    // A name with no known extension gains one rather than losing its tail.
    expect(withExtension('/a/lib', 'parquet')).toBe('/a/lib.parquet')
    expect(withExtension('/a/my.v2.lib', 'tsv')).toBe('/a/my.v2.lib.tsv')
  })
})

describe('stripPrivate', () => {
  it('removes the GUI’s own keys, which the tool would refuse', () => {
    // The CLI rejects an unknown config key outright, so a leaked __modelDir
    // would fail every run rather than be ignored.
    expect(stripPrivate({ enzyme: 'Trypsin', __modelDir: '/m' })).toEqual({ enzyme: 'Trypsin' })
  })
})

describe('App', () => {
  beforeEach(() => {
    document.body.innerHTML = ''
  })

  it('builds the form from the tool’s own effective config', async () => {
    await renderApp()
    await waitFor(() => expect(screen.getByLabelText('enzyme')).toBeTruthy())
    expect((screen.getByLabelText('enzyme') as HTMLSelectElement).value).toBe('Trypsin/P')
    expect((screen.getByLabelText('missed cleavages') as HTMLInputElement).value).toBe('1')
    // The range widget is two inputs over one [min, max] key.
    expect((screen.getByLabelText('peptide length minimum') as HTMLInputElement).value).toBe('7')
    expect((screen.getByLabelText('peptide length maximum') as HTMLInputElement).value).toBe('30')
  })

  it('keeps advanced parameters behind a disclosure', async () => {
    await renderApp()
    await waitFor(() => expect(screen.getByLabelText('enzyme')).toBeTruthy())
    expect(screen.queryByLabelText('min relative intensity')).toBeNull()
    await userEvent.click(screen.getByRole('button', { name: /advanced/i }))
    expect(screen.getByLabelText('min relative intensity')).toBeTruthy()
  })

  it('refuses to run until a FASTA, an output and the models are all there', async () => {
    await renderApp({ models: { dir: null, missing: ['peptdeep_rt_dynamic.onnx'] } })
    const run = await screen.findByRole('button', { name: /generate library/i })
    expect((run as HTMLButtonElement).disabled).toBe(true)
    expect(screen.getByText(/1 model missing/)).toBeTruthy()
  })

  it('refuses an output name the CLI would refuse', async () => {
    await renderApp()
    await userEvent.type(screen.getByLabelText('protein FASTA'), '/d/p.fasta')
    await userEvent.type(screen.getByLabelText('output library'), '/d/lib.parqet')
    const run = screen.getByRole('button', { name: /generate library/i }) as HTMLButtonElement
    expect(run.disabled).toBe(true)
    expect(screen.getByText(/must end in \.parquet or \.tsv/)).toBeTruthy()
  })

  it('sends the config as an object with no GUI-private keys', async () => {
    await renderApp()
    await waitFor(() => expect(screen.getByLabelText('enzyme')).toBeTruthy())
    await userEvent.type(screen.getByLabelText('protein FASTA'), '/d/p.fasta')
    await userEvent.type(screen.getByLabelText('output library'), '/d/lib.parquet')
    await userEvent.click(screen.getByRole('button', { name: /generate library/i }))
    await waitFor(() => expect(bridge.runs.length).toBe(1))
    const sent = bridge.runs[0] as { in: string; out: string; threads: number; config: Record<string, unknown> }
    expect(sent.in).toBe('/d/p.fasta')
    expect(sent.out).toBe('/d/lib.parquet')
    expect(sent.threads).toBe(1)
    expect(sent.config.enzyme).toBe('Trypsin/P')
    expect(Object.keys(sent.config).some((k) => k.startsWith('__'))).toBe(false)
  })

  it('accepts fractional collision energy and m/z despite whole-number native defaults', async () => {
    await renderApp()
    await userEvent.click(screen.getByRole('button', { name: /advanced/i }))
    for (const [label, value] of [['nce', '27.5'], ['precursor mz minimum', '350.25'], ['fragment mz maximum', '1750.75']]) {
      const input = screen.getByLabelText(label) as HTMLInputElement
      fireEvent.change(input, { target: { value } })
      expect(input.step).toBe('any')
      expect(input.validity.valid).toBe(true)
    }
    await startRun()
    expect(bridge.runs[0]).toMatchObject({ config: { nce: 27.5, precursor_mz: [350.25, 1200], fragment_mz: [200, 1750.75] } })
  })

  it('rejects fractional and negative thread counts before invoking the native command', async () => {
    await renderApp()
    fireEvent.change(screen.getByLabelText('protein FASTA'), { target: { value: '/d/p.fasta' } })
    fireEvent.change(screen.getByLabelText('output library'), { target: { value: '/d/lib.tsv' } })
    const threads = screen.getByLabelText('threads')
    const run = screen.getByRole('button', { name: /generate library/i }) as HTMLButtonElement
    for (const value of ['2.5', '-1']) {
      fireEvent.change(threads, { target: { value } })
      expect(run.disabled).toBe(true)
      expect(screen.getByText('Enter a whole number of threads (0 or more).')).toBeTruthy()
      await userEvent.click(run)
      expect(bridge.runs).toHaveLength(0)
    }
    fireEvent.change(threads, { target: { value: '2' } })
    expect(run.disabled).toBe(false)
    await userEvent.click(run)
    expect(bridge.runs[0]).toMatchObject({ threads: 2 })
  })

  it('streams the log and settles on the terminal event', async () => {
    await renderApp()
    await userEvent.type(screen.getByLabelText('protein FASTA'), '/d/p.fasta')
    await userEvent.type(screen.getByLabelText('output library'), '/d/lib.parquet')
    await userEvent.click(screen.getByRole('button', { name: /generate library/i }))
    await waitFor(() => expect(bridge.runs.length).toBe(1))

    act(() => bridge.emitLog('digest: 3 proteins, 120 peptides'))
    await waitFor(() => expect(screen.getByText(/digest: 3 proteins/)).toBeTruthy())
    act(() => bridge.emitDone({ ok: true, code: 0, bytes: 2_500_000 }))
    await waitFor(() => expect(screen.getByText(/Wrote 2\.5 MB/)).toBeTruthy())
  })

  it('calls exit-0-with-no-file a failure', async () => {
    // A tool that exits cleanly but writes nothing has not produced a library,
    // and reporting success there is how an empty output goes unnoticed.
    await renderApp()
    await userEvent.type(screen.getByLabelText('protein FASTA'), '/d/p.fasta')
    await userEvent.type(screen.getByLabelText('output library'), '/d/lib.parquet')
    await userEvent.click(screen.getByRole('button', { name: /generate library/i }))
    await waitFor(() => expect(bridge.runs.length).toBe(1))
    act(() => bridge.emitDone({ ok: true, code: 0, bytes: null }))
    await waitFor(() => expect(screen.getByText(/no library is on disk/)).toBeTruthy())
  })

  it('says so when the tool cannot be asked for its defaults', async () => {
    await renderApp({ configError: 'cannot run the tool: No such file' })
    await waitFor(() => expect(screen.getByText(/Could not read the tool/)).toBeTruthy())
    expect(screen.queryByLabelText('enzyme')).toBeNull()
  })

  it('disables a parameter that would do nothing, and says why', async () => {
    await renderApp()
    await userEvent.click(await screen.findByRole('button', { name: /advanced/i }))
    const cap = screen.getByLabelText('max variable modifications') as HTMLInputElement
    expect(cap.disabled).toBe(true)
    expect(screen.getByText(/no variable modifications are configured/i)).toBeTruthy()
  })

  it('offers a value the tool has even when the choice list lacks it', async () => {
    // A tool that grows an enzyme this file has not listed must not have it
    // silently rewritten to the first option the moment the form renders.
    await renderApp({ config: { ...SAMPLE_CONFIG, enzyme: 'Some-New-Protease' } })
    const sel = (await screen.findByLabelText('enzyme')) as HTMLSelectElement
    expect(sel.value).toBe('Some-New-Protease')
  })

  async function startRun(): Promise<void> {
    fireEvent.change(screen.getByLabelText('protein FASTA'), { target: { value: '/d/p.fasta' } })
    fireEvent.change(screen.getByLabelText('output library'), { target: { value: '/d/original.parquet' } })
    await userEvent.click(screen.getByRole('button', { name: /generate library/i }))
  }

  it('locks submission while launch is pending and handles completion before launch returns', async () => {
    let started!: () => void
    await renderApp({}, (b) => {
      b.api.run = vi.fn(() => new Promise<{ started: boolean }>((resolve) => { started = () => resolve({ started: true }) }))
    })
    await startRun()
    expect((screen.getByRole('button', { name: 'Running…' }) as HTMLButtonElement).disabled).toBe(true)
    act(() => bridge.emitDone({ ok: true, code: 0, bytes: 20 }))
    await act(async () => started())
    expect(screen.getByText('Wrote 20 B')).toBeTruthy()
    expect((screen.getByRole('button', { name: /generate library/i }) as HTMLButtonElement).disabled).toBe(false)
    expect(bridge.api.run).toHaveBeenCalledTimes(1)
  })

  it('reports launch rejection and permits a retry', async () => {
    await renderApp({}, (b) => { b.api.run = vi.fn().mockRejectedValue(new Error('spawn denied')) })
    await startRun()
    expect(await screen.findByText(/Could not start: Error: spawn denied/)).toBeTruthy()
    expect((screen.getByRole('button', { name: /generate library/i }) as HTMLButtonElement).disabled).toBe(false)
    bridge.api.run = vi.fn().mockResolvedValue({ started: false, reason: 'already running' })
    await userEvent.click(screen.getByRole('button', { name: /generate library/i }))
    expect(await screen.findByText('already running')).toBeTruthy()
  })

  it.each(['false', 'rejected'])('reports a %s settings save before starting and permits retry', async (failure) => {
    await renderApp({}, (b) => {
      b.api.saveLast = failure === 'false' ? vi.fn().mockResolvedValue(false) : vi.fn().mockRejectedValue('settings read-only')
    })
    await startRun()
    expect(await screen.findByText(/Could not start:.*settings/)).toBeTruthy()
    expect(bridge.runs).toHaveLength(0)
    const run = screen.getByRole('button', { name: /generate library/i }) as HTMLButtonElement
    expect(run.disabled).toBe(false)
    bridge.api.saveLast = vi.fn().mockResolvedValue(true)
    await userEvent.click(run)
    expect(bridge.runs).toHaveLength(1)
  })

  it.each(['false', 'rejected'])('reports a %s preset write and preserves unsaved input and existing presets', async (failure) => {
    await renderApp({}, (b) => {
      b.api.loadSettings = async () => ({ schemaVersion: 1, lastUsed: null, presets: { existing: { enzyme: 'Lys-C' } } })
      b.api.savePreset = failure === 'false' ? vi.fn().mockResolvedValue(false) : vi.fn().mockRejectedValue('settings read-only')
      b.api.deletePreset = failure === 'false' ? vi.fn().mockResolvedValue(false) : vi.fn().mockRejectedValue('settings read-only')
    })
    const name = screen.getByLabelText('preset name') as HTMLInputElement
    await userEvent.type(name, 'new preset')
    await userEvent.click(screen.getByRole('button', { name: 'Save' }))
    expect(await screen.findByText(/Could not save preset:.*settings/)).toBeTruthy()
    expect(name.value).toBe('new preset')
    expect(screen.queryByRole('option', { name: 'new preset' })).toBeNull()
    await userEvent.selectOptions(screen.getByLabelText('presets'), 'existing')
    await userEvent.click(screen.getByRole('button', { name: 'Delete' }))
    expect(await screen.findByText(/Could not delete preset:.*settings/)).toBeTruthy()
    expect((screen.getByLabelText('presets') as HTMLSelectElement).value).toBe('existing')
    expect(screen.getByRole('option', { name: 'existing' })).toBeTruthy()
    bridge.api.deletePreset = vi.fn().mockResolvedValue(true)
    await userEvent.click(screen.getByRole('button', { name: 'Delete' }))
    expect(screen.queryByRole('option', { name: 'existing' })).toBeNull()
  })

  it('cancels the live job and waits for its terminal event', async () => {
    await renderApp({}, (b) => { b.api.cancel = vi.fn().mockResolvedValue({ cancelled: true }) })
    await startRun()
    act(() => bridge.emitLog('Predicting MS2'))
    expect(screen.getByText('Predicting MS2')).toBeTruthy()
    await userEvent.click(screen.getByRole('button', { name: 'Cancel' }))
    expect(bridge.api.cancel).toHaveBeenCalledOnce()
    expect(screen.getByRole('button', { name: 'Running…' })).toBeTruthy()
    act(() => bridge.emitDone({ ok: false, code: null, bytes: null }))
    expect(await screen.findByText(/Failed \(exit killed\)/)).toBeTruthy()
    expect(screen.getByText('Predicting MS2')).toBeTruthy()
  })

  it('reports a cancellation error without claiming the job stopped', async () => {
    await renderApp({}, (b) => { b.api.cancel = vi.fn().mockRejectedValue('permission denied') })
    await startRun()
    await userEvent.click(screen.getByRole('button', { name: 'Cancel' }))
    expect(await screen.findByText(/Could not cancel: permission denied/)).toBeTruthy()
    expect(screen.getByRole('button', { name: 'Running…' })).toBeTruthy()
  })

  it('reveals the completed job’s output even if the form was edited during the run', async () => {
    await renderApp({}, (b) => { b.api.revealPath = vi.fn().mockResolvedValue(true) })
    await startRun()
    fireEvent.change(screen.getByLabelText('output library'), { target: { value: '/d/next.tsv' } })
    act(() => bridge.emitDone({ ok: true, code: 0, bytes: 12000 }))
    await userEvent.click(await screen.findByRole('button', { name: 'Show in folder' }))
    expect(bridge.api.revealPath).toHaveBeenCalledWith('/d/original.parquet')
  })

  it('merges a partial loaded config, reports malformed files, and resets to native defaults', async () => {
    await renderApp()
    await userEvent.click(screen.getByRole('button', { name: 'Load…' }))
    expect((screen.getByLabelText('enzyme') as HTMLSelectElement).value).toBe('Lys-C')
    expect((screen.getByLabelText('peptide length minimum') as HTMLInputElement).value).toBe('7')
    bridge.api.readConfig = vi.fn().mockRejectedValue('invalid JSON')
    await userEvent.click(screen.getByRole('button', { name: 'Load…' }))
    expect(await screen.findByText(/Could not load config: invalid JSON/)).toBeTruthy()
    expect((screen.getByLabelText('enzyme') as HTMLSelectElement).value).toBe('Lys-C')
    await userEvent.click(screen.getByRole('button', { name: 'Reset' }))
    expect((screen.getByLabelText('enzyme') as HTMLSelectElement).value).toBe('Trypsin/P')
  })

  it('saves edited config and reports an unwritable destination', async () => {
    await renderApp({}, (b) => { b.api.writeConfig = vi.fn().mockResolvedValue(true) })
    await userEvent.selectOptions(screen.getByLabelText('enzyme'), 'Lys-C')
    await userEvent.click(screen.getByRole('button', { name: 'Save…' }))
    expect(bridge.api.writeConfig).toHaveBeenCalledWith('/data/config.json', expect.objectContaining({ enzyme: 'Lys-C' }))
    bridge.api.writeConfig = vi.fn().mockRejectedValue('read-only directory')
    await userEvent.click(screen.getByRole('button', { name: 'Save…' }))
    expect(await screen.findByText(/Could not save config: read-only directory/)).toBeTruthy()
  })

  it('applies, saves, and deletes presets while preserving defaults and private model state', async () => {
    await renderApp({}, (b) => {
      b.api.loadSettings = async () => ({ schemaVersion: 1, lastUsed: { enzyme: 'Lys-C', __modelDir: '/saved/models' }, presets: {} })
    })
    expect((screen.getByLabelText('enzyme') as HTMLSelectElement).value).toBe('Lys-C')
    expect((screen.getByLabelText('model directory') as HTMLInputElement).value).toBe('/saved/models')
    await userEvent.type(screen.getByLabelText('preset name'), '  my preset  ')
    await userEvent.click(screen.getByRole('button', { name: 'Save' }))
    await screen.findByRole('option', { name: 'my preset' })
    await userEvent.click(screen.getByRole('button', { name: 'Reset' }))
    await userEvent.selectOptions(screen.getByLabelText('presets'), 'my preset')
    expect((screen.getByLabelText('enzyme') as HTMLSelectElement).value).toBe('Lys-C')
    await userEvent.click(screen.getByRole('button', { name: 'Delete' }))
    expect(screen.queryByRole('option', { name: 'my preset' })).toBeNull()
  })

  it('uses native file pickers and format toggles without discarding the output name', async () => {
    await renderApp()
    const browse = screen.getAllByRole('button', { name: 'Browse…' })
    for (const button of browse) await userEvent.click(button)
    expect((screen.getByLabelText('protein FASTA') as HTMLInputElement).value).toBe('/data/proteins.fasta')
    await userEvent.click(screen.getByRole('button', { name: '.tsv' }))
    expect((screen.getByLabelText('output library') as HTMLInputElement).value).toBe('/data/library.tsv')
    await userEvent.click(screen.getByRole('button', { name: '.parquet' }))
    expect((screen.getByLabelText('output library') as HTMLInputElement).value).toBe('/data/library.parquet')
    fireEvent.change(screen.getByLabelText('threads'), { target: { value: '0' } })
    await userEvent.click(screen.getByRole('button', { name: /generate library/i }))
    expect(bridge.runs[0]).toMatchObject({ threads: 0, modelDir: '/models' })
  })

  it('sends edited advanced values and rechecks a manually selected model directory', async () => {
    await renderApp({}, (b) => { b.api.models = vi.fn().mockResolvedValue({ dir: '/models', missing: [] }) })
    fireEvent.change(screen.getByLabelText('model directory'), { target: { value: '/other models' } })
    await waitFor(() => expect(bridge.api.models).toHaveBeenCalledWith('/other models'))
    await userEvent.click(screen.getByRole('button', { name: /advanced/i }))
    fireEvent.change(screen.getByLabelText('min relative intensity'), { target: { value: '0.05' } })
    fireEvent.change(screen.getByLabelText('variable modifications'), { target: { value: 'Oxidation (M)' } })
    expect((screen.getByLabelText('max variable modifications') as HTMLInputElement).disabled).toBe(false)
    await startRun()
    expect(bridge.runs[0]).toMatchObject({ modelDir: '/other models', config: { min_relative_intensity: 0.05, variable_modifications: ['Oxidation (M)'] } })
  })

  it('leaves NCE automatic when selecting Lumos and omits it from saved and submitted configs', async () => {
    await renderApp({}, (b) => { b.api.writeConfig = vi.fn().mockResolvedValue(true) })
    expect((screen.getByLabelText('nce') as HTMLInputElement).value).toBe('')
    expect(screen.getByPlaceholderText('automatic')).toBeTruthy()
    await userEvent.selectOptions(screen.getByLabelText('instrument'), 'Lumos')
    await userEvent.click(screen.getByRole('button', { name: 'Save…' }))
    const config = vi.mocked(bridge.api.writeConfig).mock.calls[0][1]
    expect(config.instrument).toBe('Lumos')
    expect(config).not.toHaveProperty('nce')
    await startRun()
    expect((bridge.runs[0] as { config: object }).config).not.toHaveProperty('nce')
  })

  it('preserves explicit NCE on instrument changes and returns to automatic when cleared', async () => {
    await renderApp()
    fireEvent.change(screen.getByLabelText('nce'), { target: { value: '32' } })
    await userEvent.selectOptions(screen.getByLabelText('instrument'), 'Lumos')
    await startRun()
    expect(bridge.runs[0]).toMatchObject({ config: { instrument: 'Lumos', nce: 32 } })
    act(() => bridge.emitDone({ ok: true, code: 0, bytes: 50 }))
    fireEvent.change(screen.getByLabelText('nce'), { target: { value: '' } })
    await userEvent.click(screen.getByRole('button', { name: /generate library/i }))
    expect((bridge.runs[1] as { config: object }).config).not.toHaveProperty('nce')
  })

  it('treats NCE in an imported full config as explicit until Reset', async () => {
    await renderApp({}, (b) => { b.api.readConfig = async () => ({ ...SAMPLE_CONFIG, nce: 30 }) })
    await userEvent.click(screen.getByRole('button', { name: 'Load…' }))
    await userEvent.selectOptions(screen.getByLabelText('instrument'), 'Lumos')
    expect((screen.getByLabelText('nce') as HTMLInputElement).value).toBe('30')
    await userEvent.click(screen.getByRole('button', { name: 'Reset' }))
    expect((screen.getByLabelText('nce') as HTMLInputElement).value).toBe('')
  })

  it('keeps generation disabled for a failed health probe', async () => {
    await renderApp({ probeFails: true })
    await startRun()
    expect(bridge.runs).toHaveLength(0)
    expect((screen.getByRole('button', { name: /generate library/i }) as HTMLButtonElement).disabled).toBe(true)
  })

  it('requires the same output extension spelling as the CLI', async () => {
    await renderApp()
    fireEvent.change(screen.getByLabelText('protein FASTA'), { target: { value: '/d/p.fasta' } })
    fireEvent.change(screen.getByLabelText('output library'), { target: { value: '/d/output.PARQUET' } })
    expect((screen.getByRole('button', { name: /generate library/i }) as HTMLButtonElement).disabled).toBe(true)
    await userEvent.click(screen.getByRole('button', { name: '.parquet' }))
    expect((screen.getByLabelText('output library') as HTMLInputElement).value).toBe('/d/output.parquet')
  })
})
