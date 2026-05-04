import { readFileSync } from 'node:fs';
import { getCollection } from 'astro:content';
import satori from 'satori';
import { Resvg } from '@resvg/resvg-js';
import type { APIRoute } from 'astro';

// V4 brand-aligned OG card. Cream paper, the trim/galore wordmark top-left,
// a diagonal lime band cutting through the middle, and the page title +
// description bottom-left. Rendered with satori (HTML/JSX → SVG) and
// resvg-js (SVG → PNG) at build time. Output: /TrimGalore/og/<slug>.png.

// Fonts and the entries collection are loaded at module scope so they're read
// once per build, not once per page route.
const interRegular = readFileSync('./src/assets/fonts/Inter-Regular.ttf');
const interBold = readFileSync('./src/assets/fonts/Inter-Bold.ttf');
const bricolageBold = readFileSync('./src/assets/fonts/Bricolage-Bold.ttf');

const BG = '#F7F2E5';
const INK = '#23201B';
const INK_SOFT = '#5C5042';
const LIME = '#C2EE49';

// Tiny helper for satori's React-element-shaped tree.
type El = {
  type: string;
  props: { style?: Record<string, any>; children?: El | string | (El | string)[] } & Record<string, any>;
};
const e = (
  type: string,
  style: Record<string, any>,
  children?: El | string | (El | string)[],
): El => ({
  type,
  props: { style, children },
});

function card(title: string, description: string): El {
  return e(
    'div',
    {
      width: 1200,
      height: 630,
      background: BG,
      position: 'relative',
      display: 'flex',
      flexDirection: 'column',
      padding: '64px 80px',
      fontFamily: 'Inter',
      color: INK,
      overflow: 'hidden',
    },
    [
      // Lime diagonal band, full bleed, slightly above centre so the title
      // below has breathing room.
      e('div', {
        position: 'absolute',
        left: -200,
        right: -200,
        top: 230,
        height: 96,
        background: LIME,
        transform: 'rotate(-14deg)',
        transformOrigin: 'center',
      }),

      // Wordmark: trim · slash · galore
      e(
        'div',
        {
          display: 'flex',
          alignItems: 'center',
          fontFamily: 'Bricolage',
          fontWeight: 700,
          fontSize: 32,
          color: INK,
          letterSpacing: '-0.025em',
          gap: 10,
        },
        [
          e('span', { display: 'flex' }, 'trim'),
          e('span', {
            display: 'flex',
            width: 8,
            height: 28,
            background: LIME,
            transform: 'skewX(-14deg)',
          }),
          e('span', { display: 'flex' }, 'galore'),
        ],
      ),

      // Spacer pushes title block to the bottom.
      e('div', { display: 'flex', flex: 1 }),

      // Page title (large display). Drops to 64px when the title is long
      // enough that 84px would wrap awkwardly past the 980px max-width.
      e(
        'div',
        {
          display: 'flex',
          fontFamily: 'Bricolage',
          fontWeight: 700,
          fontSize: title.length > 28 ? 64 : 84,
          lineHeight: 1.02,
          letterSpacing: '-0.035em',
          color: INK,
          maxWidth: 980,
        },
        title,
      ),

      // Description
      e(
        'div',
        {
          display: 'flex',
          marginTop: 22,
          fontFamily: 'Inter',
          fontWeight: 400,
          fontSize: 26,
          lineHeight: 1.4,
          color: INK_SOFT,
          maxWidth: 920,
        },
        description,
      ),
    ],
  );
}

const entries = await getCollection('docs');
const toSlug = (id: string) => id.replace(/\.(mdx?|markdown)$/, '') || 'index';

const pages = entries.map((entry) => ({
  slug: `${toSlug(entry.id)}.png`,
  title: entry.data.title || 'Trim Galore',
  description:
    entry.data.description ||
    'Adapter and quality trimming for NGS data. A single Rust binary, drop-in compatible with v0.6.x scripts and pipelines.',
}));

export async function getStaticPaths() {
  return pages.map((p) => ({
    params: { route: p.slug },
    props: { title: p.title, description: p.description },
  }));
}

export const GET: APIRoute = async ({ props }) => {
  const { title, description } = props as { title: string; description: string };

  const svg = await satori(card(title, description) as any, {
    width: 1200,
    height: 630,
    fonts: [
      { name: 'Inter', data: interRegular, weight: 400, style: 'normal' },
      { name: 'Inter', data: interBold, weight: 700, style: 'normal' },
      { name: 'Bricolage', data: bricolageBold, weight: 700, style: 'normal' },
    ],
  });

  const png = new Resvg(svg, { fitTo: { mode: 'width', value: 1200 } }).render().asPng();

  return new Response(png, {
    headers: { 'Content-Type': 'image/png', 'Cache-Control': 'public, max-age=31536000, immutable' },
  });
};
