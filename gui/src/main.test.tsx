import { expect, it, vi } from 'vitest'
import { createRoot } from 'react-dom/client'

vi.mock('./api', () => ({}))
vi.mock('react-dom/client', () => ({ createRoot: vi.fn(() => ({ render: vi.fn() })) }))

it('mounts the desktop entry point at its document root', async () => {
  document.body.innerHTML = '<div id="root"></div>'
  await import('./main')
  expect(createRoot).toHaveBeenCalledWith(document.getElementById('root'))
  expect(vi.mocked(createRoot).mock.results[0].value.render).toHaveBeenCalledOnce()
})
