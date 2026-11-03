import starlight from '@astrojs/starlight';
import { satteri } from '@astrojs/markdown-satteri';
import AutoImport from 'astro-auto-import';
import { defineConfig } from 'astro/config';
import { siNextflow } from 'simple-icons';
import codeblocks, { linksValidatorExclude } from 'starlight-codeblocks';
import { nextflow } from 'starlight-codeblocks/adapters/nextflow';
import starlightLinksValidator from 'starlight-links-validator';
import starlightLlmsTxt from 'starlight-llms-txt';
import starlightPageActions from 'starlight-page-actions';
import starlightQuiz from 'starlight-quiz';
import links from './src/markdown/links.mjs';
import { redirects } from './src/redirects.mjs';
import { sidebar } from './src/sidebar.mjs';

// One build serves one prefix: /development/, /latest/ or /<version>/.
const base = process.env.BASE ?? '/';

const nfCoreModule = ({ name, path }) =>
  path?.includes('nf-core/')
    ? { href: `https://nf-co.re/modules/${name.toLowerCase()}`, kind: 'process', source: 'nf-core modules' }
    : undefined;

export default defineConfig({
  site: 'https://training.nextflow.io',
  base,
  trailingSlash: 'always',
  redirects,
  markdown: {
    processor: satteri({ features: { headingAttributes: true }, hastPlugins: [links({ base })] }),
  },
  integrations: [
    starlight({
      title: 'Nextflow Training',
      logo: { light: './src/assets/nextflow-logo.png', dark: './src/assets/nextflow-logo-white.png', replacesTitle: true },
      favicon: '/favicon.png',
      defaultLocale: 'root',
      locales: {
        root: { label: 'English', lang: 'en' },
        ca: { label: 'Català' },
        de: { label: 'Deutsch' },
        es: { label: 'Español' },
        fr: { label: 'Français' },
        hi: { label: 'हिन्दी' },
        it: { label: 'Italiano' },
        ko: { label: '한국어' },
        pl: { label: 'Polski' },
        pt: { label: 'Português' },
        tr: { label: 'Türkçe' },
      },
      sidebar,
      social: [{ icon: 'github', label: 'GitHub', href: 'https://github.com/nextflow-io/training' }],
      editLink: { baseUrl: 'https://github.com/nextflow-io/training/edit/master/docs/' },
      customCss: ['./src/styles/custom.css'],
      components: {
        Head: './src/components/Head.astro',
        Banner: './src/components/Banner.astro',
        EditLink: './src/components/EditLink.astro',
        SiteTitle: './src/components/VersionMenu.astro',
        PageSidebar: './src/components/TocResources.astro',
        Footer: './src/components/Consent.astro',
      },
      plugins: [
        codeblocks({
          fileIcons: {
            icons: { nextflow: siNextflow.svg },
            files: { '.nf': 'nextflow', 'nextflow.config': 'nextflow' },
            languages: { nextflow: { icon: 'nextflow', colour: '#31C9AC' } },
          },
          runnable: { outputDelay: 60 },
          apiLinks: { adapters: [nextflow({ modules: nfCoreModule })] },
        }),
        starlightQuiz({ progressPosition: 'bottom' }),
        starlightLlmsTxt(),
        // No baseUrl: page actions then leave llms.txt to starlight-llms-txt.
        starlightPageActions(),
        // Translated pages link to English fallback pages on purpose.
        starlightLinksValidator({ exclude: linksValidatorExclude, errorOnFallbackPages: false }),
      ],
    }),
    // Must come after starlight(): it only injects imports when MDX is already registered.
    AutoImport({
      imports: [
        { '@astrojs/starlight/components': ['Aside', 'Tabs', 'TabItem', 'Code', 'CardGrid', 'LinkCard', 'LinkButton'] },
        { 'starlight-codeblocks/components': ['CodeWalkthrough', 'Scrollycoding', 'Step'] },
        { 'starlight-quiz/components': ['Quiz', 'QuizResults'] },
        './src/components/CourseInfo.astro',
      ],
    }),
  ],
});
