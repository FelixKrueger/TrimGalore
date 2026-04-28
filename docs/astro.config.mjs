// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

export default defineConfig({
  site: 'https://felixkrueger.github.io',
  base: '/TrimGalore',
  trailingSlash: 'ignore',
  // Disable smartypants. It rewrites `--` to em-dashes (—) even inside
  // JSX <code> elements in MDX, which silently mangles CLI flag names
  // like --fastqc → —fastqc on the homepage. Em-dashes we actually want
  // are typed literally throughout, so nothing is lost.
  markdown: { smartypants: false },
  integrations: [
    starlight({
      title: 'Trim Galore',
      description:
        'Consistent quality and adapter trimming for next-generation sequencing data, with special handling for RRBS libraries.',
      // No logo image: project is maintained solo outside Babraham; the
      // historical logo on the user guide was removed for the same reason.
      favicon: '/favicon.svg',
      social: [
        {
          icon: 'github',
          label: 'GitHub',
          href: 'https://github.com/FelixKrueger/TrimGalore',
        },
      ],
      editLink: {
        baseUrl:
          'https://github.com/FelixKrueger/TrimGalore/edit/master/docs/',
      },
      components: {
        Head: './src/components/Head.astro',
        Header: './src/components/Header.astro',
        Hero: './src/components/Hero.astro',
        SiteTitle: './src/components/SiteTitle.astro',
      },
      sidebar: [
        {
          label: 'Start here',
          items: [
            { label: 'Introduction', link: '/' },
            { label: 'Installation', slug: 'install' },
            { label: 'Quick start', slug: 'quickstart' },
          ],
        },
        {
          label: 'User guide',
          items: [
            { label: 'How Trim Galore works', slug: 'guide/overview' },
            { label: 'Quality trimming', slug: 'guide/quality' },
            { label: 'Adapter trimming', slug: 'guide/adapters' },
            { label: 'Length filtering', slug: 'guide/length' },
            { label: 'Paired-end data', slug: 'guide/paired-end' },
            { label: 'Output files', slug: 'guide/outputs' },
            { label: 'Trimming reports', slug: 'guide/reports' },
            { label: 'Flag reference notes', slug: 'guide/flags' },
          ],
        },
        {
          label: 'Specialty modes',
          items: [
            { label: 'RRBS mode', slug: 'modes/rrbs' },
            { label: 'Hard-trimming', slug: 'modes/hardtrim' },
            { label: 'Mouse Epigenetic Clock', slug: 'modes/clock' },
            { label: 'IMPLICON UMI transfer', slug: 'modes/implicon' },
            { label: 'Demultiplexing', slug: 'modes/demux' },
            { label: 'Poly-G / Poly-A', slug: 'modes/poly' },
          ],
        },
        {
          label: 'Bisulfite & RRBS',
          items: [
            { label: 'A brief guide to RRBS', slug: 'rrbs/guide' },
            { label: 'Directional libraries', slug: 'rrbs/directional' },
            {
              label: 'Non-directional & paired-end',
              slug: 'rrbs/non-directional',
            },
            { label: 'QC measures for RRBS', slug: 'rrbs/qc' },
          ],
        },
        {
          label: 'Performance',
          items: [
            { label: 'Benchmarks', slug: 'performance/benchmarks' },
            { label: 'Threading model', slug: 'performance/threading' },
          ],
        },
        {
          label: 'Reference',
          items: [
            { label: 'Migrating from v0.6.x', slug: 'reference/migration' },
            { label: 'Changelog', slug: 'reference/changelog' },
            { label: 'Logo', slug: 'reference/logo' },
            { label: 'Credits & license', slug: 'reference/credits' },
          ],
        },
      ],
      customCss: ['./src/styles/custom.css'],
    }),
  ],
});
