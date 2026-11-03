import { pluginLineNumbers } from '@expressive-code/plugin-line-numbers';
import { defineEcConfig } from '@astrojs/starlight/expressive-code';
import { pluginCodeblocks } from 'starlight-codeblocks/expressive-code';

export default defineEcConfig({
  plugins: [pluginLineNumbers(), pluginCodeblocks()],
  defaultProps: { showLineNumbers: false },
});
