import { fireEvent, render, screen } from '@testing-library/react'
import { describe, expect, it, vi } from 'vitest'
import ParamField from './ParamField'
import type { ParamKind } from './paramLayout'

function field(kind: ParamKind, value: unknown) {
  const onChange = vi.fn()
  render(<ParamField spec={{ name: 'value', label: 'setting', description: 'help', group: 'core', kind }}
                     value={value} onChange={onChange} />)
  return onChange
}

describe('parameter editing', () => {
  it.each<['int' | 'double', number]>([['int', 2], ['double', 1.25]])('edits and clears a %s', (kind, value) => {
    const change = field(kind, value)
    fireEvent.change(screen.getByLabelText('setting'), { target: { value: '3.5' } })
    expect(change).toHaveBeenLastCalledWith(3.5)
    fireEvent.change(screen.getByLabelText('setting'), { target: { value: '' } })
    expect(change).toHaveBeenLastCalledWith(null)
  })

  it.each<['int-range' | 'double-range', number[]]>([['int-range', [7, 30]], ['double-range', [350.5, 1200.5]]])(
    'edits each end of a %s without changing the other', (kind, value) => {
      const change = field(kind, value)
      fireEvent.change(screen.getByLabelText('setting minimum'), { target: { value: '8' } })
      expect(change).toHaveBeenLastCalledWith([8, value[1]])
      fireEvent.change(screen.getByLabelText('setting maximum'), { target: { value: '40' } })
      expect(change).toHaveBeenLastCalledWith([value[0], 40])
    })

  it('passes a checkbox as a boolean', () => {
    const change = field('bool', false)
    fireEvent.click(screen.getByLabelText('setting'))
    expect(change).toHaveBeenCalledWith(true)
  })

  it('converts a charge list to numbers and ignores empty entries', () => {
    const change = field('int-list', [1, 2])
    fireEvent.change(screen.getByLabelText('setting'), { target: { value: '2, 3, , 4' } })
    expect(change).toHaveBeenCalledWith([2, 3, 4])
  })

  it('trims modification names while keeping spaces inside a name', () => {
    const change = field('string-list', [])
    fireEvent.change(screen.getByLabelText('setting'), { target: { value: ' Oxidation (M)\n\n Acetyl (Protein N-term) ' } })
    expect(change).toHaveBeenCalledWith(['Oxidation (M)', 'Acetyl (Protein N-term)'])
  })

  it('edits an unlisted future text parameter', () => {
    const change = field('string', 'old')
    fireEvent.change(screen.getByLabelText('setting'), { target: { value: 'new value' } })
    expect(change).toHaveBeenCalledWith('new value')
  })

  it('shows unsupported structured values read-only', () => {
    field('unsupported', { nested: true })
    const input = screen.getByLabelText('setting') as HTMLInputElement
    expect(input.readOnly).toBe(true)
    expect(input.value).toBe('{"nested":true}')
    expect(screen.getByText(/edit the config file directly/)).toBeTruthy()
  })
})
