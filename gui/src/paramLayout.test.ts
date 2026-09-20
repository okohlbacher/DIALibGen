import { describe, expect, it } from 'vitest'
import { buildSpecs, inferKind, inertBecause } from './paramLayout'
import { SAMPLE_CONFIG } from './testing/mockBridge'

describe('inferKind', () => {
  it('separates int from double by the tool’s own default', () => {
    // nce defaults to 30.0 but arrives as JSON 30, so an int/double split that
    // trusted the wire type alone would give it a step of 1.
    expect(inferKind('missed_cleavages', 1)).toBe('int')
    expect(inferKind('min_relative_intensity', 0.0001)).toBe('double')
  })

  it('treats a two-number array as a range only for the documented keys', () => {
    expect(inferKind('peptide_length', [7, 30])).toBe('int-range')
    expect(inferKind('precursor_mz', [350.0, 1200.5])).toBe('double-range')
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

describe('inertBecause', () => {
  it('disables the variable-mod cap when there are no variable mods', () => {
    expect(inertBecause('max_variable_modifications', { variable_modifications: [] })).toMatch(/no effect/i)
    expect(inertBecause('max_variable_modifications', { variable_modifications: ['Oxidation (M)'] })).toBeNull()
  })

  it('disables the decoy m/z switch when no decoys are generated', () => {
    expect(inertBecause('recompute_decoy_mz', { decoys: 'none' })).toMatch(/no effect/i)
    expect(inertBecause('recompute_decoy_mz', { decoys: 'reverse' })).toBeNull()
  })
})
