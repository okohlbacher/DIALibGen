import type { JSX } from 'react'
import type { ParamSpec } from './paramLayout'

interface Props {
  spec: ParamSpec
  value: unknown
  onChange: (v: unknown) => void
  /// Why this parameter currently does nothing, from inertBecause(). The field
  /// is still shown -- hiding it would make the option undiscoverable -- but it
  /// is disabled and says why.
  inert?: string | null
}

const num = (s: string): number => s.trim() === '' ? Number.NaN : Number(s)

// One widget per inferred kind. The tool's own defaults drive the control, so
// the UI cannot offer a shape the CLI would reject.
export default function ParamField({ spec, value, onChange, inert }: Props): JSX.Element {
  const off = !!inert
  const id = `p-${spec.name}`
  const label = (
    <label htmlFor={id} title={spec.description}>
      {spec.label}
    </label>
  )
  // The reason replaces the description rather than joining it: a disabled
  // field's one useful piece of information is what would switch it on.
  const help = <p className={off ? 'help inert' : 'help'}>{inert ?? spec.description}</p>
  const cls = `field${off ? ' off' : ''}`

  if (spec.kind === 'bool') {
    return (
      <div className={`${cls} check`}>
        <input id={id} type="checkbox" disabled={off} checked={value === true}
               onChange={(e) => onChange(e.target.checked)} />
        {label}
        {help}
      </div>
    )
  }

  if (spec.kind === 'int-range' || spec.kind === 'double-range') {
    const step = spec.kind === 'int-range' ? 1 : 'any'
    const pair = Array.isArray(value) ? (value as number[]) : [0, 0]
    const set = (i: 0 | 1, s: string) => {
      const next = [pair[0], pair[1]]
      next[i] = num(s)
      onChange(next)
    }
    return (
      <div className={cls}>
        {label}
        <div className="row tight">
          <input id={id} type="number" step={step} min={spec.min} max={spec.max} disabled={off}
                 aria-label={`${spec.label} minimum`} value={Number.isFinite(pair[0]) ? String(pair[0]) : ''}
                 onChange={(e) => set(0, e.target.value)} />
          <span className="between">to</span>
          <input type="number" step={step} min={spec.min} max={spec.max} disabled={off}
                 aria-label={`${spec.label} maximum`} value={Number.isFinite(pair[1]) ? String(pair[1]) : ''}
                 onChange={(e) => set(1, e.target.value)} />
        </div>
        {help}
      </div>
    )
  }

  if (spec.kind === 'int-list') {
    const items = Array.isArray(value) ? (value as number[]) : []
    return (
      <div className={cls}>
        {label}
        <input id={id} type="text" disabled={off} value={items.join(', ')} placeholder="e.g. 1, 2, 3, 4"
               onChange={(e) =>
                 onChange(
                   e.target.value
                     .split(',')
                     .map((s) => s.trim())
                     .filter((s) => s !== '')
                     .map(Number)
                 )
               } />
        {help}
      </div>
    )
  }

  if (spec.kind === 'string-list') {
    const items = Array.isArray(value) ? (value as string[]) : []
    return (
      <div className={cls}>
        {label}
        <textarea id={id} rows={Math.max(2, items.length + 1)} disabled={off} value={items.join('\n')}
                  placeholder="one per line (empty = none)"
                  onChange={(e) =>
                    onChange(e.target.value.split('\n').map((s) => s.trim()).filter((s) => s !== ''))
                  } />
        {help}
      </div>
    )
  }

  if (spec.choices?.length) {
    // A value the tool has but this list does not must not be silently replaced
    // by the first option, which is what a plain <select> would do on render.
    const current = String(value ?? '')
    const choices = spec.choices.includes(current) ? spec.choices : [current, ...spec.choices]
    return (
      <div className={cls}>
        {label}
        <select id={id} disabled={off} value={current} onChange={(e) => onChange(e.target.value)}>
          {choices.map((c) => (
            <option key={c} value={c}>{c}</option>
          ))}
        </select>
        {help}
      </div>
    )
  }

  if (spec.kind === 'int' || spec.kind === 'double') {
    return (
      <div className={cls}>
        {label}
        <input id={id} type="number" step={spec.kind === 'double' ? 'any' : 1} min={spec.min} max={spec.max}
               placeholder={spec.name === 'nce' ? 'automatic' : undefined}
               disabled={off} value={value === null || value === undefined ? '' : String(value)}
               onChange={(e) => onChange(e.target.value === '' ? null : num(e.target.value))} />
        {help}
      </div>
    )
  }

  if (spec.kind === 'unsupported') {
    return (
      <div className={`${cls} off`}>
        {label}
        <input id={id} type="text" readOnly value={JSON.stringify(value)} />
        <p className="help inert">This form cannot edit this value — edit the config file directly.</p>
      </div>
    )
  }

  return (
    <div className={cls}>
      {label}
      <input id={id} type="text" disabled={off} value={String(value ?? '')}
             onChange={(e) => onChange(e.target.value)} />
      {help}
    </div>
  )
}
