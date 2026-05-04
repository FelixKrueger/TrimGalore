// Rasterise the brand SVG logos to PNG at sensible high-resolution sizes.
//
// Run via `npm run logos` from the Docs/ workspace. Idempotent — overwrites
// existing PNGs in public/logos/. Re-run whenever the source SVGs change.
//
// `sharp` is already a docs dependency (used elsewhere by Astro for image
// optimisation), so no new package is required.

import sharp from 'sharp';
import { readFile } from 'node:fs/promises';
import { join, dirname } from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const LOGO_DIR = join(__dirname, '..', 'public', 'logos');

// One PNG per SVG variant. Width-only resize preserves the SVG's native
// aspect ratio. Square logos pin both dimensions for explicit clarity.
const VARIANTS = [
  { name: 'hero-light',     width: 1920 },
  { name: 'hero-dark',      width: 1920 },
  { name: 'wordmark-light', width: 1200 },
  { name: 'wordmark-dark',  width: 1200 },
  { name: 'square-light',   width: 1024, height: 1024 },
  { name: 'square-dark',    width: 1024, height: 1024 },
];

console.log(`Rasterising ${VARIANTS.length} SVG → PNG into ${LOGO_DIR}`);

for (const { name, width, height } of VARIANTS) {
  const svgPath = join(LOGO_DIR, `${name}.svg`);
  const pngPath = join(LOGO_DIR, `${name}.png`);
  const svg = await readFile(svgPath);
  // density 384 = 4× the SVG's default 96 DPI, giving sharp's libvips
  // rasteriser enough headroom to produce crisp edges at the target width.
  const pipeline = sharp(svg, { density: 384 }).resize({ width, height });
  await pipeline.png({ compressionLevel: 9 }).toFile(pngPath);
  const meta = await sharp(pngPath).metadata();
  console.log(`  ✓ ${name}.png  ${meta.width}×${meta.height}`);
}
