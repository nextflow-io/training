# nf-greeting plugin

## Building

To build the plugin:

```bash
make assemble
```

## Testing with Nextflow

The plugin can be tested without a local Nextflow installation:

1. Build and install the plugin to your local Nextflow installation: `make install`
2. Run a pipeline with the plugin: `nextflow run hello -plugins nf-greeting@0.1.0`

## Publishing

Plugins can be published to a central plugin registry to make them accessible to the Nextflow community.

Follow these steps to publish the plugin to the Nextflow Plugin Registry:

1. Create a file named `$HOME/.gradle/gradle.properties`, where $HOME is your home directory. Add the following properties:

   - `pluginRegistry.accessToken`: Your Nextflow Plugin Registry access token.

2. Use the following command to package and create a release for your plugin on GitHub: `make release`.

> [!NOTE]
> The Nextflow Pluging registry is currently avaialable as private beta technology. Contact info@nextflow.io to learn how to get access to it.
