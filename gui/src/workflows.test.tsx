import { act, fireEvent, render, screen, waitFor } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import { describe, expect, it, vi } from 'vitest'
import App from './App'
import { installMockBridge, SAMPLE_CONFIG, type MockBridge } from './testing/mockBridge'
import refinementDefaults from './testing/refinement-defaults.json'
import nativeTuning from './testing/tuning-options.json'
import type { Mode, TuningOption } from './types'

let bridge: MockBridge
async function setup(configure?: (b: MockBridge) => void): Promise<void> {
  bridge = installMockBridge()
  configure?.(bridge)
  await act(async () => { render(<App />) })
  await screen.findByLabelText('enzyme')
}
async function mode(value: Mode): Promise<void> {
  await userEvent.selectOptions(screen.getByLabelText('workflow'), value)
  if (value === 'tune') await screen.findByLabelText('training heads')
  else await screen.findByLabelText(value === 'generate' ? 'enzyme' : 'q precursor')
}
function fillRefinement(): void {
  fireEvent.change(screen.getByLabelText('input library'), { target: { value: '/data/library.parquet' } })
  fireEvent.change(screen.getByLabelText('DIA-NN report (Parquet)'), { target: { value: '/data/report.parquet' } })
  fireEvent.change(screen.getByLabelText('output library'), { target: { value: '/data/tuned.parquet' } })
}
const button = (name: 'Refine library' | 'Fine-tune library') => screen.getByRole('button', { name }) as HTMLButtonElement

describe('desktop refinement and fine-tuning', () => {
  it('refines a library/report without requiring any prediction model or FASTA', async () => {
    await setup((b) => { b.api.models = vi.fn().mockResolvedValue({ dir: null, missing: ['rt', 'ccs', 'ms2'] }) })
    await mode('refine')
    expect(screen.queryByLabelText('protein FASTA')).toBeNull()
    expect(screen.queryByLabelText('model directory')).toBeNull()
    expect(screen.getByText('models not required')).toBeTruthy()
    expect(button('Refine library').disabled).toBe(true)
    fillRefinement()
    fireEvent.change(screen.getByLabelText('residual report (optional TSV)'), { target: { value: '/data/residuals.tsv' } })
    expect(button('Refine library').disabled).toBe(false)
    await userEvent.click(button('Refine library'))
    expect(bridge.runs[0]).toEqual({ mode: 'refine', in: '/data/library.parquet', ids: '/data/report.parquet', out: '/data/tuned.parquet',
      outReport: '/data/residuals.tsv', config: refinementDefaults, modelDir: null, threads: 1 })
  })

  it.each(['rt', 'ccs', 'both'])('fine-tunes the %s head selection with only required model files', async (head) => {
    await setup((b) => {
      b.api.models = vi.fn(async (_dir, heads) => ({ dir: '/models', missing: heads ? [] : ['ms2.onnx'] }))
    })
    await mode('tune')
    fillRefinement()
    await userEvent.selectOptions(screen.getByLabelText('training heads'), head)
    const wanted = head === 'both' ? ['rt', 'ccs'] : [head]
    await waitFor(() => expect(bridge.api.models).toHaveBeenLastCalledWith(undefined, wanted))
    expect(screen.queryByLabelText('write rt')).toBeNull()
    expect(screen.queryByLabelText('q precursor')).toBeNull()
    expect((screen.getByLabelText('machine device') as HTMLInputElement).value).toBe('cpu')
    expect(screen.getByText(/CUDA training.*require a suitable source build/)).toBeTruthy()
    fireEvent.change(screen.getByLabelText('keep tuned models in (optional)'), { target: { value: '/data/tuned models' } })
    expect(button('Fine-tune library').disabled).toBe(false)
    await userEvent.click(button('Fine-tune library'))
    expect(bridge.runs[0]).toMatchObject({ mode: 'tune', in: '/data/library.parquet', ids: '/data/report.parquet',
      tuneOutModels: '/data/tuned models', config: { filter: false, write_rt: false }, tuning: { tune_heads: head, 'machine:device': 'cpu' } })
    expect(bridge.runs[0]).not.toHaveProperty('tune')
    expect((bridge.runs[0] as { config: object }).config).not.toHaveProperty('__gui')
  })

  it('requires Parquet report input and a supported library extension', async () => {
    await setup(); await mode('tune'); fillRefinement()
    fireEvent.change(screen.getByLabelText('DIA-NN report (Parquet)'), { target: { value: '/data/report.tsv' } })
    expect(button('Fine-tune library').disabled).toBe(true)
    fireEvent.change(screen.getByLabelText('DIA-NN report (Parquet)'), { target: { value: '/data/report.parquet' } })
    fireEvent.change(screen.getByLabelText('input library'), { target: { value: '/data/input.csv' } })
    expect(button('Fine-tune library').disabled).toBe(true)
  })

  it('exposes every native tuning option and preserves scalar kinds and restrictions', async () => {
    await setup(); await mode('tune')
    expect(screen.queryByLabelText('filter q value')).toBeNull()
    await userEvent.click(screen.getByRole('button', { name: /Advanced fine-tuning/ }))
    for (const option of nativeTuning) {
      const field = document.getElementById(`p-${option.name}`) as HTMLInputElement
      expect(field, option.name).toBeTruthy()
      if (option.kind === 'int') expect(field.step).toBe('1')
      if (option.kind === 'double') expect(field.step).toBe('any')
    }
    fireEvent.change(screen.getByLabelText('filter rt max minutes'), { target: { value: '32.5' } })
    fillRefinement(); await userEvent.click(button('Fine-tune library'))
    expect(bridge.runs[0]).toMatchObject({ tuning: { 'filter:rt_max_minutes': 32.5 } })
  })

  it('supports fine-tuning before refinement and keeps per-mode configurations separate', async () => {
    await setup()
    await userEvent.selectOptions(screen.getByLabelText('enzyme'), 'Lys-C')
    await mode('refine')
    fireEvent.change(screen.getByLabelText('q precursor'), { target: { value: '0.005' } })
    await userEvent.click(screen.getByLabelText('fine-tune before refinement'))
    await screen.findByLabelText('training heads')
    await userEvent.selectOptions(screen.getByLabelText('training heads'), 'rt')
    await mode('generate')
    expect((screen.getByLabelText('enzyme') as HTMLSelectElement).value).toBe('Lys-C')
    expect(screen.queryByLabelText('training heads')).toBeNull()
    await mode('refine')
    expect((screen.getByLabelText('q precursor') as HTMLInputElement).value).toBe('0.005')
    fillRefinement(); await userEvent.click(button('Refine library'))
    expect(bridge.runs[0]).toMatchObject({ mode: 'refine', tune: true, config: { q_precursor: 0.005 }, tuning: { tune_heads: 'rt' } })
    expect((bridge.runs[0] as { config: object }).config).not.toHaveProperty('enzyme')
  })

  it('keeps refinement float defaults editable as decimals', async () => {
    await setup(); await mode('refine')
    await userEvent.click(screen.getByRole('button', { name: /Advanced/ }))
    for (const label of ['im ramp top', 'min match fraction', 'intensity min correlation', 'intensity mz tol ppm']) {
      const field = screen.getByLabelText(label) as HTMLInputElement
      expect(field.step).toBe('any')
      fireEvent.change(field, { target: { value: '0.25' } })
    }
    fillRefinement(); await userEvent.click(button('Refine library'))
    expect(bridge.runs[0]).toMatchObject({ config: { im_ramp_top: 0.25, min_match_fraction: 0.25 } })
  })

  it('rejects fractional, empty and out-of-range training integers and nonpositive learning rates', async () => {
    await setup(); await mode('tune'); fillRefinement()
    for (const value of ['1.5', '0', '2147483648', '']) {
      fireEvent.change(screen.getByLabelText('training epochs'), { target: { value } })
      expect(button('Fine-tune library').disabled).toBe(true)
      await userEvent.click(button('Fine-tune library'))
    }
    fireEvent.change(screen.getByLabelText('training epochs'), { target: { value: '100' } })
    fireEvent.change(screen.getByLabelText('training lr'), { target: { value: '0' } })
    expect(screen.getByRole('alert').textContent).toContain('learning rate must be positive')
    expect(bridge.runs).toHaveLength(0)
    fireEvent.change(screen.getByLabelText('training lr'), { target: { value: '0.0002' } })
    expect(button('Fine-tune library').disabled).toBe(false)
  })

  it('rejects invalid ranges and incompatible training settings before starting', async () => {
    await setup(); await mode('tune'); fillRefinement()
    await userEvent.click(screen.getByRole('button', { name: /Advanced fine-tuning/ }))
    fireEvent.change(screen.getByLabelText('filter q value'), { target: { value: '1.1' } })
    expect(button('Fine-tune library').disabled).toBe(true)
    fireEvent.change(screen.getByLabelText('filter q value'), { target: { value: '0.01' } })
    fireEvent.change(screen.getByLabelText('training warmup'), { target: { value: '101' } })
    expect(screen.getByRole('alert').textContent).toContain('warmup cannot exceed')
    fireEvent.change(screen.getByLabelText('training warmup'), { target: { value: '10' } })
    fireEvent.change(screen.getByLabelText('cohort train size'), { target: { value: '200' } })
    fireEvent.change(screen.getByLabelText('cohort train frac'), { target: { value: '0.5' } })
    expect(screen.getByRole('alert').textContent).toContain('size or fraction')
    fireEvent.change(screen.getByLabelText('cohort train frac'), { target: { value: '0' } })
    await userEvent.click(screen.getByLabelText('cohort full fit'))
    expect(screen.getByRole('alert').textContent).toContain('Full fit')
    await userEvent.click(screen.getByLabelText('cohort full fit'))
    fireEvent.change(screen.getByLabelText('filter min charge'), { target: { value: '1' } })
    expect(screen.getByRole('alert').textContent).toContain('charge 1')
    await userEvent.click(screen.getByLabelText('filter allow z1'))
    fireEvent.change(screen.getByLabelText('machine device'), { target: { value: 'mps' } })
    expect(screen.getByRole('alert').textContent).toContain('cpu, cuda or cuda:N')
    expect(bridge.runs).toHaveLength(0)
  })

  it('blocks incompatible observed RT and minmax choices when tuning during refinement', async () => {
    await setup(); await mode('refine'); fillRefinement()
    await userEvent.click(screen.getByLabelText('filter'))
    expect(button('Refine library').disabled).toBe(true)
    await userEvent.click(screen.getByLabelText('fine-tune before refinement'))
    await screen.findByLabelText('training heads')
    expect(button('Refine library').disabled).toBe(false)
    await userEvent.selectOptions(screen.getByLabelText('training heads'), 'ccs')
    expect(button('Refine library').disabled).toBe(true)
    await userEvent.click(screen.getByLabelText('filter'))
    await userEvent.click(screen.getByRole('button', { name: /^.*Advanced \(/ }))
    await userEvent.selectOptions(screen.getByLabelText('rt unit'), 'minmax')
    expect(button('Refine library').disabled).toBe(true)
  })

  it('ignores late defaults for a mode that was left before they arrived', async () => {
    let resolve!: (value: Record<string, unknown>) => void
    await setup((b) => {
      const original = b.api.defaultConfig
      b.api.defaultConfig = vi.fn((selected) => selected === 'refine' ? new Promise<Record<string, unknown>>((done) => { resolve = done }) : original(selected))
    })
    await userEvent.selectOptions(screen.getByLabelText('workflow'), 'refine')
    await mode('tune')
    await act(async () => resolve({ ...refinementDefaults, q_precursor: 0.9 }))
    fillRefinement(); await userEvent.click(button('Fine-tune library'))
    expect(bridge.runs[0]).toMatchObject({ mode: 'tune', config: { q_precursor: 0.01, filter: false } })
  })

  it('does not reuse model readiness while a different head is being checked', async () => {
    let resolve!: (value: { dir: string; missing: string[] }) => void
    await setup((b) => { b.api.models = vi.fn(async (_dir, heads) => heads?.join() === 'ccs'
      ? new Promise<{ dir: string; missing: string[] }>((done) => { resolve = done }) : { dir: '/models', missing: [] }) })
    await mode('tune'); fillRefinement()
    expect(button('Fine-tune library').disabled).toBe(false)
    await userEvent.selectOptions(screen.getByLabelText('training heads'), 'ccs')
    expect(button('Fine-tune library').disabled).toBe(true)
    await act(async () => resolve({ dir: '/models', missing: ['ccs.onnx'] }))
    expect(button('Fine-tune library').disabled).toBe(true)
    expect(screen.getByText('1 model missing')).toBeTruthy()
  })

  it('reports tuning metadata and model-check failures without enabling submission', async () => {
    await setup((b) => {
      b.api.tuningOptions = vi.fn().mockRejectedValue('unavailable INI')
      b.api.models = vi.fn().mockRejectedValue('unreadable models')
    })
    await userEvent.selectOptions(screen.getByLabelText('workflow'), 'tune')
    await screen.findByText(/Could not read fine-tuning options: unavailable INI/)
    fillRefinement()
    expect(button('Fine-tune library').disabled).toBe(true)
    expect(screen.getByText(/Not found: unreadable models/)).toBeTruthy()
  })

  it('round-trips desktop files, named presets and last-used settings without leaking metadata to native config', async () => {
    let saved: Record<string, unknown> = {}
    await setup((b) => {
      b.api.writeConfig = vi.fn(async (_path, value) => { saved = value; return true })
      b.api.readConfig = async () => saved
    })
    await mode('tune'); fillRefinement()
    await userEvent.selectOptions(screen.getByLabelText('training heads'), 'ccs')
    fireEvent.change(screen.getByLabelText('training epochs'), { target: { value: '33' } })
    await userEvent.click(screen.getByRole('button', { name: 'Save…' }))
    expect(saved).toMatchObject({ __gui: { mode: 'tune', tuning: { tune_heads: 'ccs', 'train:epochs': 33 } } })
    await userEvent.type(screen.getByLabelText('preset name'), 'CCS recipe')
    await userEvent.click(screen.getByRole('button', { name: 'Save' }))
    await mode('generate')
    await userEvent.click(screen.getByRole('button', { name: 'Load…' }))
    await screen.findByLabelText('training heads')
    expect((screen.getByLabelText('workflow') as HTMLSelectElement).value).toBe('tune')
    expect((screen.getByLabelText('training epochs') as HTMLInputElement).value).toBe('33')
    await userEvent.click(screen.getByRole('button', { name: 'Reset' }))
    await userEvent.selectOptions(screen.getByLabelText('presets'), 'CCS recipe')
    expect((screen.getByLabelText('training heads') as HTMLSelectElement).value).toBe('ccs')
    await userEvent.click(button('Fine-tune library'))
    expect(bridge.saved[0]).toMatchObject({ __gui: { mode: 'tune', tuning: { 'train:epochs': 33 } } })
    expect((bridge.runs[0] as { config: object }).config).not.toHaveProperty('__gui')
  })

  it('restores the selected mode and training controls from saved settings at startup', async () => {
    bridge = installMockBridge()
    bridge.api.loadSettings = async () => ({ schemaVersion: 1, presets: {}, lastUsed: { __gui: { version: 1, mode: 'tune',
      tuning: { tune_heads: 'rt', 'train:epochs': 42 }, modelDir: '/saved/rt', threads: 2 } } })
    await act(async () => { render(<App />) })
    await screen.findByLabelText('training heads')
    expect((screen.getByLabelText('workflow') as HTMLSelectElement).value).toBe('tune')
    expect((screen.getByLabelText('training epochs') as HTMLInputElement).value).toBe('42')
    expect((screen.getByLabelText('model directory') as HTMLInputElement).value).toBe('/saved/rt')
    fillRefinement(); await userEvent.click(button('Fine-tune library'))
    expect(bridge.runs[0]).toMatchObject({ mode: 'tune', threads: 2, tuning: { tune_heads: 'rt', 'train:epochs': 42 } })
  })

  it('loads plain native JSON into the current mode, but applies legacy presets to generation', async () => {
    await setup((b) => {
      b.api.readConfig = async () => ({ q_precursor: 0.002 })
      b.api.loadSettings = async () => ({ schemaVersion: 1, lastUsed: null, presets: { legacy: { enzyme: 'Lys-C' } } })
    })
    await mode('refine')
    await userEvent.click(screen.getByRole('button', { name: 'Load…' }))
    expect((screen.getByLabelText('q precursor') as HTMLInputElement).value).toBe('0.002')
    await userEvent.selectOptions(screen.getByLabelText('presets'), 'legacy')
    expect((screen.getByLabelText('workflow') as HTMLSelectElement).value).toBe('generate')
    expect((screen.getByLabelText('enzyme') as HTMLSelectElement).value).toBe('Lys-C')
  })

  it('ignores a late file read after changing mode', async () => {
    let resolve!: (value: Record<string, unknown>) => void
    await setup((b) => { b.api.readConfig = vi.fn(() => new Promise<Record<string, unknown>>((done) => { resolve = done })) })
    await userEvent.click(screen.getByRole('button', { name: 'Load…' }))
    await mode('tune')
    await act(async () => resolve({ enzyme: 'Lys-C' }))
    fillRefinement(); await userEvent.click(button('Fine-tune library'))
    expect((bridge.runs[0] as { config: object }).config).not.toHaveProperty('enzyme')
  })

  it.each(['refine', 'tune'] as const)('locks %s mode during the run, preserves logging and waits for cancel completion', async (selected) => {
    await setup((b) => { b.api.cancel = vi.fn().mockResolvedValue({ cancelled: true }) })
    await mode(selected); fillRefinement()
    await userEvent.click(button(selected === 'refine' ? 'Refine library' : 'Fine-tune library'))
    expect((screen.getByLabelText('workflow') as HTMLSelectElement).disabled).toBe(true)
    expect((screen.getByLabelText('presets') as HTMLSelectElement).disabled).toBe(true)
    act(() => bridge.emitLog('training/refinement progress'))
    await userEvent.click(screen.getByRole('button', { name: 'Cancel' }))
    expect(bridge.api.cancel).toHaveBeenCalledOnce()
    expect(screen.getByRole('button', { name: 'Running…' })).toBeTruthy()
    act(() => bridge.emitDone({ ok: false, code: null }))
    expect(screen.getByText(/Failed \(exit killed\)/)).toBeTruthy()
    expect(screen.getByText('training/refinement progress')).toBeTruthy()
    expect((screen.getByLabelText('workflow') as HTMLSelectElement).disabled).toBe(false)
  })

  it('rejects unsupported desktop settings and reports startup settings failures', async () => {
    await setup((b) => { b.api.readConfig = async () => ({ __gui: { version: 99, mode: 'tune' } }) })
    await userEvent.click(screen.getByRole('button', { name: 'Load…' }))
    expect(await screen.findByText(/unsupported desktop settings format/)).toBeTruthy()
    expect((screen.getByLabelText('workflow') as HTMLSelectElement).value).toBe('generate')
  })

  it('reports an unavailable settings store during startup', async () => {
    bridge = installMockBridge()
    bridge.api.loadSettings = vi.fn().mockRejectedValue('settings unavailable')
    await act(async () => { render(<App />) })
    expect(await screen.findByText(/Could not load settings: settings unavailable/)).toBeTruthy()
  })

  it('rejects native refinement changes imported into prediction-only fine-tuning mode', async () => {
    await setup((b) => { b.api.readConfig = async () => ({ write_rt: true }) })
    await mode('tune'); fillRefinement()
    await userEvent.click(screen.getByRole('button', { name: 'Load…' }))
    expect(button('Fine-tune library').disabled).toBe(true)
    expect(screen.getByRole('alert').textContent).toContain('prediction-only configuration')
    await userEvent.click(screen.getByRole('button', { name: 'Reset' }))
    expect(button('Fine-tune library').disabled).toBe(false)
  })

  it('uses library/report/model-output pickers and reports a rejected native dialog', async () => {
    await setup((b) => { b.api.pickReport = vi.fn().mockRejectedValueOnce('dialog unavailable').mockResolvedValue('/data/observations.parquet') })
    await mode('tune')
    const browseFor = (id: string) => document.getElementById(id)!.parentElement!.querySelector('button')!
    await userEvent.click(browseFor('library'))
    await userEvent.click(browseFor('tuned-models'))
    expect((screen.getByLabelText('input library') as HTMLInputElement).value).toBe('/data/library.parquet')
    expect((screen.getByLabelText('keep tuned models in (optional)') as HTMLInputElement).value).toBe('/models')
    await userEvent.click(browseFor('ids'))
    expect(await screen.findByText('dialog unavailable')).toBeTruthy()
    await userEvent.click(browseFor('ids'))
    expect((screen.getByLabelText('DIA-NN report (Parquet)') as HTMLInputElement).value).toBe('/data/observations.parquet')
    expect(screen.queryByText('dialog unavailable')).toBeNull()
  })

  it('does not let nonempty native/imported model paths override the desktop directory', async () => {
    await setup((b) => {
      const original = b.api.defaultConfig
      b.api.defaultConfig = async (selected) => selected === 'generate'
        ? { ...SAMPLE_CONFIG, rt_model: '/startup/rt.onnx', ms2_model: '/startup/ms2.onnx', ccs_model: '/startup/ccs.onnx' }
        : original(selected)
      b.api.models = vi.fn(async (directory) => ({ dir: directory || '/checked/bundled/models', missing: [] }))
      b.api.readConfig = async () => ({ rt_model: '/imported/rt.onnx', ms2_model: '/imported/ms2.onnx', ccs_model: '/imported/ccs.onnx' })
    })
    await userEvent.click(screen.getByRole('button', { name: 'Load…' }))
    fireEvent.change(screen.getByLabelText('model directory'), { target: { value: '/selected/models' } })
    fireEvent.change(screen.getByLabelText('protein FASTA'), { target: { value: '/data/proteins.fasta' } })
    fireEvent.change(screen.getByLabelText('output library'), { target: { value: '/data/generation.parquet' } })
    await userEvent.click(screen.getByRole('button', { name: 'Generate library' }))
    expect(bridge.runs[0]).toMatchObject({ modelDir: '/selected/models' })
    for (const key of ['rt_model', 'ms2_model', 'ccs_model']) expect((bridge.runs[0] as { config: object }).config).not.toHaveProperty(key)
  })

  it('binds an automatic tuning run to the checked fallback directory', async () => {
    await setup((b) => { b.api.models = vi.fn().mockResolvedValue({ dir: '/complete/fallback/models', missing: [] }) })
    await mode('tune'); fillRefinement()
    expect((screen.getByLabelText('model directory') as HTMLInputElement).value).toBe('')
    await userEvent.click(button('Fine-tune library'))
    expect(bridge.runs[0]).toMatchObject({ modelDir: '/complete/fallback/models' })
  })

  it('shows future native tuning options without hard-coding them in the frontend', async () => {
    await setup((b) => { b.api.tuningOptions = async () => [...nativeTuning as TuningOption[],
      { name: 'train:future', kind: 'double', value: 1.5, description: 'New native control', min: 0 }] })
    await mode('tune')
    await userEvent.click(screen.getByRole('button', { name: /Advanced fine-tuning/ }))
    expect((screen.getByLabelText('training future') as HTMLInputElement).value).toBe('1.5')
  })
})
