import { pluginLineNumbers } from '@expressive-code/plugin-line-numbers';
import { defineEcConfig } from '@astrojs/starlight/expressive-code';
import { pluginCodeblocks } from 'starlight-codeblocks/expressive-code';

/** A dark-theme value, then a light-theme one. */
const pick =
  (dark, light) =>
  ({ theme }) =>
    theme.type === 'dark' ? dark : light;

export default defineEcConfig({
  plugins: [pluginLineNumbers(), pluginCodeblocks()],
  defaultProps: { showLineNumbers: false },
  styleOverrides: {
    borderRadius: '6px',
    borderColor: pick('#31353d', '#dddede'),
    codeFontFamily: "'Roboto Mono', Menlo, monospace",
    frames: {
      frameBoxShadowCssValue: 'none',
      editorTabBarBackground: pick('#1e2127', '#f7f7f7'),
      editorTabBarBorderBottomColor: pick('#31353d', '#dddede'),
      editorActiveTabBackground: pick('#1e2127', '#f7f7f7'),
      editorActiveTabIndicatorTopColor: 'transparent',
      editorActiveTabIndicatorBottomColor: '#0dc09d',
      editorActiveTabIndicatorHeight: '2px',
      terminalTitlebarBackground: pick('#1e2127', '#f7f7f7'),
      terminalTitlebarBorderBottomColor: pick('#31353d', '#dddede'),
    },
  },
});
