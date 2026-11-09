import { createRequire } from 'node:module';

// A path, not the package name: Prettier imports plugin names from the working directory,
// which misses a prek hook's packages. require.resolve also searches NODE_PATH, which prek sets.
const astro = createRequire(import.meta.url).resolve('prettier-plugin-astro');

export default {
  singleQuote: true,
  printWidth: 120,
  plugins: [astro],
  overrides: [{ files: '*.astro', options: { parser: 'astro' } }],
};
