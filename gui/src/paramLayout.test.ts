import { describe, expect, it } from 'vitest'
import { buildSpecs, inferKind, inertBecause, invalidValue, tuningSpecs } from './paramLayout'
import nativeTuning from './testing/tuning-options.json'
import type { TuningOption } from './types'
import { SAMPLE_CONFIG } from './testing/mockBridge'

describe('inferKind', () => {
  it('preserves decimal semantics when JSON defaults are whole numbers', () => {
    // nce defaults to 30.0 but arrives as JSON 30, so an int/double split that
    // trusted the wire type alone would give it a step of 1.
    expect(inferKind('missed_cleavages', 1)).toBe('int')
    expect(inferKind('nce', 30)).toBe('double')
    expect(inferKind('min_relative_intensity', 0)).toBe('double')
    expect(inferKind('min_relative_intensity', 0.0001)).toBe('double')
  })

  it('treats a two-number array as a range only for the documented keys', () => {
    expect(inferKind('peptide_length', [7, 30])).toBe('int-range')
    expect(inferKind('precursor_mz', [350, 1200])).toBe('double-range')
    expect(inferKind('fragment_mz', [200, 1800])).toBe('double-range')
    // A charge list that happened to hold two entries must NOT become a range,
    // or a [2,3] default would render as "2 to 3" and write back a range.
    expect(inferKind('precursor_charges', [2, 3])).toBe('int-list')
  })

  it('calls an empty list a string list, not unsupported', () => {
    expect(inferKind('variable_modifications', [])).toBe('string-list')
  })

  it('does not drop a value it cannot edit', () => {
    expect(inferKind('whatever', { a: 1 })).toBe('unsupported')
    expect(inferKind('whatever', null)).toBe('unsupported')
    expect(inferKind('whatever', [1, 'mixed'])).toBe('unsupported')
  })
})

describe('buildSpecs', () => {
  it('covers every key the tool emitted', () => {
    const specs = buildSpecs(SAMPLE_CONFIG)
    expect(specs.map((s) => s.name).sort()).toEqual(Object.keys(SAMPLE_CONFIG).sort())
  })

  it('surfaces a key this file has never heard of instead of hiding it', () => {
    // The property that keeps the GUI from drifting behind the CLI: a new
    // config key must appear, described as undescribed, rather than vanish.
    const specs = buildSpecs({ ...SAMPLE_CONFIG, brand_new_knob: 5 })
    const s = specs.find((x) => x.name === 'brand_new_knob')
    expect(s).toBeDefined()
    expect(s?.hidden).toBe(false)
    expect(s?.group).toBe('advanced')
    expect(s?.kind).toBe('int')
    expect(s?.description).toMatch(/not yet described/i)
  })

  it('hides the three model paths, which the model picker owns', () => {
    const specs = buildSpecs({ ...SAMPLE_CONFIG, nce_source: 'instrument-default:QE', instrument_named: 'QE' })
    for (const n of ['rt_model', 'ms2_model', 'ccs_model', 'schema_version', 'nce_source', 'instrument_named']) {
      expect(specs.find((s) => s.name === n)?.hidden).toBe(true)
    }
  })

  it('hides the recorded sources and the echoed training settings, and describes the iRT standards', () => {
    const specs = buildSpecs({ ...SAMPLE_CONFIG, tuning: { 'train:epochs': 1 } })
    for (const n of ['precursor_mz_source', 'fragment_mz_source', 'tuning']) expect(specs.find((s) => s.name === n)?.hidden).toBe(true)
    expect(specs.find((s) => s.name === 'irt_standards')).toMatchObject({ group: 'advanced', hidden: false, kind: 'string' })
    expect(specs.find((s) => s.name === 'irt_standards')?.description).not.toMatch(/not yet described/i)
  })

  it('lets the m/z windows and the collision energy be blank (the instrument decides)', () => {
    const specs = buildSpecs({ ...SAMPLE_CONFIG, precursor_mz: null, fragment_mz: null })
    for (const n of ['precursor_mz', 'fragment_mz', 'nce']) {
      const spec = specs.find((s) => s.name === n)!
      expect(spec.auto).toBe(true)
      expect(invalidValue(spec, null)).toBeNull()
    }
    expect(specs.find((s) => s.name === 'precursor_mz')).toMatchObject({ kind: 'double-range', group: 'core' })
    expect(specs.find((s) => s.name === 'fragment_mz')?.kind).toBe('double-range')
    expect(specs.find((s) => s.name === 'precursor_mz')?.description).toMatch(/Astral 380-980; others 350-1200/)
    expect(specs.find((s) => s.name === 'fragment_mz')?.description).toMatch(/Astral 150-2000; others 200-1800/)
    expect(invalidValue(specs.find((s) => s.name === 'precursor_mz')!, [900, 400])).not.toBeNull()
    expect(invalidValue(specs.find((s) => s.name === 'enzyme')!, null)).not.toBeNull()
  })

  it('puts the core settings first and in pipeline order', () => {
    const core = buildSpecs(SAMPLE_CONFIG).filter((s) => s.group === 'core' && !s.hidden)
    expect(core[0].name).toBe('enzyme')
    expect(core.map((s) => s.name)).toContain('fixed_modifications')
  })

  it('sorts new keys after known keys without dropping either kind', () => {
    expect(buildSpecs({ z_future: 1, enzyme: 'Trypsin', a_future: true }).map((s) => s.name))
      .toEqual(['enzyme', 'a_future', 'z_future'])
  })

  it('limits fragment charges to the two channels the model predicts', () => {
    expect(buildSpecs({ max_fragment_charge: 2 })[0]).toMatchObject({ min: 1, max: 2 })
  })
})

describe('tuningSpecs', () => {
  it('keeps every native training option, advanced or not, and shows the train threads as core', () => {
    const specs = tuningSpecs(nativeTuning as TuningOption[])
    expect(specs).toHaveLength((nativeTuning as TuningOption[]).length)
    expect(specs.find((s) => s.name === 'machine:train_threads')?.group).toBe('core')
    expect(specs.find((s) => s.name === 'machine:threads')).toBeUndefined()
    expect(specs.find((s) => s.name === 'tune_heads')?.choices).toEqual(['auto', 'rt', 'ccs', 'both'])
  })
})

describe('inertBecause', () => {
  it('disables the iRT standards unless iRT rescaling is on', () => {
    expect(inertBecause('irt_standards', { irt_rescale: false })).toMatch(/no effect/i)
    expect(inertBecause('irt_standards', { irt_rescale: true })).toBeNull()
  })

  it('disables the variable-mod cap when there are no variable mods', () => {
    expect(inertBecause('max_variable_modifications', { variable_modifications: [] })).toMatch(/no effect/i)
    expect(inertBecause('max_variable_modifications', { variable_modifications: ['Oxidation (M)'] })).toBeNull()
  })

  it('disables the decoy m/z switch when no decoys are generated', () => {
    expect(inertBecause('recompute_decoy_mz', { decoys: 'none' })).toMatch(/no effect/i)
    expect(inertBecause('recompute_decoy_mz', { decoys: 'reverse' })).toBeNull()
  })
})

// Values can also arrive from a JSON/settings file, so browser input constraints
// alone do not establish the CLI contract.
describe('parameter value validation', () => {
  it('rejects non-finite, fractional integer, wrong-type and reversed-range values', async () => {
    const { invalidValue } = await import('./paramLayout')
    const base = { name: 'value', label: 'value', description: '', group: 'core' as const }
    for (const value of [NaN, Infinity, -Infinity, null, '4', 1.5]) {
      expect(invalidValue({ ...base, kind: 'int', min: 0 }, value)).not.toBeNull()
    }
    expect(invalidValue({ ...base, kind: 'double', min: 0, max: 1 }, 0.25)).toBeNull()
    expect(invalidValue({ ...base, kind: 'double-range' }, [4, 2])).not.toBeNull()
    expect(invalidValue({ ...base, kind: 'int-range' }, [1.5, 2])).not.toBeNull()
    expect(invalidValue({ ...base, kind: 'int-list' }, [1, Infinity])).not.toBeNull()
    expect(invalidValue({ ...base, name: 'precursor_charges', kind: 'int-list' }, [2, 2])).not.toBeNull()
    expect(invalidValue({ ...base, name: 'precursor_charges', kind: 'int-list' }, [9])).not.toBeNull()
    expect(invalidValue({ ...base, kind: 'bool' }, 'true')).not.toBeNull()
    expect(invalidValue({ ...base, kind: 'string' }, 2)).not.toBeNull()
    expect(invalidValue({ ...base, kind: 'string-list' }, [2])).not.toBeNull()
    expect(invalidValue({ ...base, kind: 'string', choices: ['rt', 'ccs'] }, 'both')).not.toBeNull()
  })
})
