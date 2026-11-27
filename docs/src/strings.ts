import en from './content/i18n/en.json';

/** UI strings under `prefix` for client scripts, keyed without it: this language's, or English where missing. */
export function uiStrings(all: Record<string, string>, prefix: string) {
  const pick = (o: Record<string, string>) =>
    Object.entries(o)
      .filter(([k]) => k.startsWith(prefix))
      .map(([k, v]) => [k.slice(prefix.length), v]);
  return Object.fromEntries([...pick(en), ...pick(all)]);
}
