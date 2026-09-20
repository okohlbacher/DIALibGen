// Bundle the installed production npm closure and the selected target's
// non-dev Cargo closure (including build/proc-macro dependencies, conservatively).
// Cargo populates its checksum-verified registry cache; license supplements are local.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { execFileSync } from 'node:child_process';
import { createHash } from 'node:crypto';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const cargoRoot = path.join(root, 'src-tauri');
const supplementsDir = path.join(root, 'scripts/license-supplements');
const supplements = JSON.parse(fs.readFileSync(path.join(supplementsDir, 'index.json')));
const target = process.argv[2] || process.env.TAURI_ENV_TARGET_TRIPLE ||
  execFileSync('rustc', ['-vV'], { encoding: 'utf8' }).match(/^host: (.+)$/m)[1];
const metadata = JSON.parse(execFileSync('cargo', ['metadata', '--locked',
  '--format-version', '1', '--filter-platform', target], { cwd: cargoRoot, encoding: 'utf8', maxBuffer: 32 * 1024 * 1024 }));
const sha256 = value => createHash('sha256').update(value).digest('hex');
// These Cargo.lock fields are plain quoted strings for registry packages.
const cargoChecksums = new Map([...fs.readFileSync(path.join(cargoRoot, 'Cargo.lock'), 'utf8')
  .matchAll(/\[\[package\]\]([\s\S]*?)(?=\[\[package\]\]|$)/g)].map(([, block]) => {
    const field = key => block.match(new RegExp(`^${key} = "([^"\\n]+)"$`, 'm'))?.[1];
    return [`${field('name')}@${field('version')}`, field('checksum')];
  }));
const notices = [];
const entries = [];

function licenseFiles(directory) {
  const result = [];
  for (const entry of fs.readdirSync(directory, { withFileTypes: true })) {
    if (!/^(licen[cs]e|copying|notice|unlicense|copyright)/i.test(entry.name)) continue;
    const filename = path.join(directory, entry.name);
    if (entry.isFile()) result.push(filename);
    else if (entry.isDirectory()) {
      for (const child of fs.readdirSync(filename, { withFileTypes: true })) {
        if (child.isFile()) result.push(path.join(filename, child.name));
      }
    }
  }
  return result.sort();
}

function add(entry, directory, declaredFile) {
  const id = `${entry.ecosystem}:${entry.name}@${entry.version}`;
  const files = licenseFiles(directory);
  if (declaredFile) files.push(path.resolve(directory, declaredFile));
  const texts = [...new Set(files)].map(filename => ({
    file: path.relative(directory, filename), text: fs.readFileSync(filename, 'utf8'),
  }));
  for (const item of supplements[id] || []) {
    const text = fs.readFileSync(path.join(supplementsDir, item.file), 'utf8');
    if (sha256(text) !== item.sha256) throw new Error(`changed supplemental license: ${item.file}`);
    texts.push({ file: item.file, text, url: item.url });
  }
  if (!entry.license || texts.length === 0 || texts.every(item => item.file.endsWith('.spdx'))) {
    throw new Error(`missing license terms for ${id}; add a version-pinned supplement`);
  }
  entry.licenseFiles = texts.map(({ file, text, url }) => ({ file, sha256: sha256(text), ...(url && { url }) }));
  entries.push(entry);
  notices.push(`\n${'='.repeat(78)}\n${id}\nLicense: ${entry.license}\nSource: ${entry.source}\n` +
    (entry.authors?.length ? `Authors: ${entry.authors.join('; ')}\n` : '') +
    texts.map(item => `\n--- ${item.file}${item.url ? ` (${item.url})` : ''} ---\n${item.text}\n`).join(''));
}

const lock = JSON.parse(fs.readFileSync(path.join(root, 'package-lock.json')));
for (const [location, info] of Object.entries(lock.packages).sort()) {
  if (!location || info.dev) continue;
  const directory = path.join(root, location);
  const pkg = JSON.parse(fs.readFileSync(path.join(directory, 'package.json')));
  add({ ecosystem: 'npm', name: pkg.name, version: pkg.version, license: pkg.license,
    source: info.resolved, integrity: info.integrity,
    repository: typeof pkg.repository === 'string' ? pkg.repository : pkg.repository?.url }, directory);
}

const nodes = new Map(metadata.resolve.nodes.map(node => [node.id, node]));
const selected = new Set();
const pending = [metadata.resolve.root];
while (pending.length) {
  const id = pending.pop();
  if (selected.has(id)) continue;
  selected.add(id);
  for (const dep of nodes.get(id).deps) {
    if (dep.dep_kinds.some(kind => kind.kind !== 'dev')) pending.push(dep.pkg);
  }
}
for (const pkg of metadata.packages.sort((a, b) => a.id.localeCompare(b.id))) {
  if (!selected.has(pkg.id) || pkg.id === metadata.resolve.root) continue;
  const directory = path.dirname(pkg.manifest_path);
  const checksum = cargoChecksums.get(`${pkg.name}@${pkg.version}`);
  if (!pkg.source?.startsWith('registry+') || !checksum) {
    throw new Error(`untracked Cargo source/checksum: ${pkg.id}`);
  }
  add({ ecosystem: 'cargo', name: pkg.name, version: pkg.version, license: pkg.license,
    source: `https://crates.io/api/v1/crates/${pkg.name}/${pkg.version}/download`,
    repository: pkg.repository, authors: pkg.authors, checksum }, directory, pkg.license_file);
}
const output = path.join(cargoRoot, 'resources/third-party-licenses');
fs.mkdirSync(output, { recursive: true });
fs.writeFileSync(path.join(output, 'inventory.json'), JSON.stringify({ target,
  scope: 'Production npm and selected-target non-dev Cargo closure, including Cargo build dependencies. System libraries and bundled DIALibGen runtime have separate notices.',
  packages: entries }, null, 2) + '\n');
fs.writeFileSync(path.join(output, 'THIRD_PARTY_NOTICES.txt'),
  'DIALibGen desktop GUI: third-party dependency notices\n' +
  'Dependency sources are unmodified. Exact-version source downloads and integrity checksums are in inventory.json.\n' +
  'This conservative Cargo inventory also includes build-time dependencies.\n' + notices.join(''));
console.log(`Bundled licenses for ${entries.filter(x => x.ecosystem === 'npm').length} npm and ${entries.filter(x => x.ecosystem === 'cargo').length} Cargo packages (${target}).`);
