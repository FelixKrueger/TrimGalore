# Trim Galore documentation

This directory holds the Astro Starlight site that publishes to
<https://felixkrueger.github.io/TrimGalore/>.

## Local development

```bash
cd docs
npm install
npm run dev      # dev server on http://localhost:4321/TrimGalore/
npm run build    # static build into docs/dist/
npm run preview  # serve docs/dist/ locally
```

## Deploying

The site auto-deploys from the `master` branch via
`.github/workflows/docs.yml`. To enable GitHub Pages for this repo:

1. Settings → Pages → Build and deployment. Set Source to GitHub Actions.
2. Push a change under `docs/`, `CHANGELOG.md`, or
   `.github/workflows/docs.yml`, or trigger the workflow manually.

The workflow runs `npm ci && npm run build` in `docs/` and uploads
`docs/dist/` as the Pages artifact.

## Content layout

| Folder | Purpose |
|--------|---------|
| `src/content/docs/` | Markdown and MDX pages, one file per route. |
| `src/assets/` | Images embedded by markdown pages. Astro optimises these to webp. |
| `public/` | Static files served at the site root (favicon, robots.txt, etc.). |
| `astro.config.mjs` | Site config: sidebar, base path, integrations. |

## Sync with the repo `CHANGELOG.md`

`docs/src/content/docs/reference/changelog.md` is a copy of the
top-level `CHANGELOG.md` with a Starlight frontmatter block prepended
and one inline image path rewritten. Keep them in sync when releasing.
