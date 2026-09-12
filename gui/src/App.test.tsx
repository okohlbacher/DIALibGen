import { act, render, screen, waitFor } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { beforeEach, describe, expect, it } from 'vitest'
import App, { stripPrivate, withExtension } from './App'
import { installMockBridge, SAMPLE_CONFIG, type MockBridge } from './testing/mockBridge'

let bridge: MockBridge

async function renderApp(opts: Parameters<typeof installMockBridge>[0] = {}): Promise<void> {
  bridge = installMockBridge(opts)
  // The startup effect resolves three promises before the first paint; without
  // act() every one of them logs a "not wrapped in act(...)" warning that
  // buries the real output.
  await act(async () => {
    render(<App />)
  })
  await waitFor(() => expect(screen.getByText(/DIALibGen 0\.2\.0|cannot execute/)).toBeTruthy())
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
    const sent = bridge.runs[0] as { in: string; out: string; config: Record<string, unknown> }
    expect(sent.in).toBe('/d/p.fasta')
    expect(sent.out).toBe('/d/lib.parquet')
    expect(sent.config.enzyme).toBe('Trypsin/P')
    expect(Object.keys(sent.config).some((k) => k.startsWith('__'))).toBe(false)
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
})
