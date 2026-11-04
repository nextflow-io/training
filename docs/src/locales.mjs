// Every site language. Translations live in src/content/docs/<key>/.
export const allLocales = {
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
};

// TRANSLATIONS=false builds English only: about a tenth of the pages, for faster previews.
export const translations = process.env.TRANSLATIONS !== 'false';

export const translatedLocales = Object.keys(allLocales).filter((key) => key !== 'root');

export const locales = translations ? allLocales : { root: allLocales.root };
