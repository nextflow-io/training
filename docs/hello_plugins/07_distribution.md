# Part 7: Distribution

Once your plugin is working locally, you have several options for sharing it with others.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 6 to use as your starting point:

    ```bash
    cp -r solutions/6-configuration/* .
    ```

---

## 1. Distribution options

| Distribution method  | Use case                                           | Approval required          |
| -------------------- | -------------------------------------------------- | -------------------------- |
| **Public registry**  | Open source plugins for the community              | Yes (name must be claimed) |
| **Internal hosting** | Private/proprietary plugins within an organization | No                         |

---

## 2. Publishing to the public registry

The [Nextflow plugin registry](https://registry.nextflow.io/) is the official way to share plugins with the community.

!!! tip "Plugin registry"

    The Nextflow plugin registry is currently in public preview.
    See the [Nextflow documentation](https://www.nextflow.io/docs/latest/guides/gradle-plugin.html#publishing-a-plugin) for the latest details.

### 2.1. Claim your plugin name

Before publishing, claim your plugin name in the registry:

1. Go to the [Nextflow plugin registry](https://registry.nextflow.io/)
2. Sign in with your GitHub account
3. Claim your plugin name (e.g., `nf-greeting`)

You can claim a name before the plugin exists to reserve it.

### 2.2. Configure API credentials

Create a Gradle properties file to store your registry credentials:

```bash
touch ~/.gradle/gradle.properties
```

Add your API key (obtain this from the registry after signing in):

```properties title="~/.gradle/gradle.properties"
npr.apiKey=YOUR_API_KEY_HERE
```

!!! warning "Keep your API key secret"

    Don't commit this file to version control.
    The `~/.gradle/` directory is outside your project, so it won't be included in your repository.

### 2.3. Prepare for release

Before publishing, ensure your plugin is ready:

1. **Update the version** in `build.gradle` (use [semantic versioning](https://semver.org/))
2. **Run tests** to ensure everything works: `make test`
3. **Update documentation** in your README

```groovy title="build.gradle"
version = '1.0.0'  // Use semantic versioning: MAJOR.MINOR.PATCH
```

### 2.4. Publish to the registry

Run the release command from your plugin directory:

```bash
cd nf-greeting
make release
```

This builds the plugin and publishes it to the registry in one step.

??? info "What `make release` does"

    The `make release` command runs `./gradlew publishPlugin`, which:

    1. Compiles your plugin code
    2. Runs tests
    3. Packages the plugin as a JAR file
    4. Uploads to the Nextflow plugin registry
    5. Makes it available for users to install

### 2.5. Using published plugins

Once published, users can install your plugin without any local setup:

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting'        // Latest version
    id 'nf-greeting@1.0.0'  // Specific version (recommended)
}
```

Nextflow automatically downloads the plugin from the registry on first use.

---

## 3. Internal distribution

Organizations often need to distribute plugins internally without using the public registry.
This is useful for proprietary plugins, plugins under development, or plugins that shouldn't be publicly available.

!!! note "What internal distribution provides"

    Internal distribution uses a simple `plugins.json` file that tells Nextflow where to download plugin ZIP files.
    This is **not** a full self-hosted registry (no web UI, search, or automatic updates).
    A full self-hosted registry solution may be available in the future.

### 3.1. Build the plugin ZIP

The ZIP file will be at `build/distributions/`:

```bash
./gradlew assemble
ls build/distributions/
# nf-myplugin-0.1.0.zip
```

### 3.2. Host the files

Host the ZIP file(s) somewhere accessible to your users:

- Internal web server
- S3 bucket (with appropriate access)
- GitHub releases (for private repos)
- Shared network drive (using `file://` URLs)

### 3.3. Create the plugins.json index

Create a `plugins.json` file that describes available plugins:

```json
[
  {
    "id": "nf-myplugin",
    "releases": [
      {
        "version": "1.0.0",
        "url": "https://internal.example.com/plugins/nf-myplugin-1.0.0.zip",
        "date": "2025-01-09T10:00:00Z",
        "sha512sum": "5abe4cbc643ca0333cba545846494b17488d19d17...",
        "requires": ">=24.04.0"
      }
    ]
  }
]
```

Host this file alongside your plugin ZIPs (or anywhere accessible).

??? tip "Generating the checksum"

    ```bash
    sha512sum nf-myplugin-1.0.0.zip | awk '{print $1}'
    ```

??? info "plugins.json field reference"

    | Field | Description |
    |-------|-------------|
    | `id` | Plugin identifier (e.g., `nf-myplugin`) |
    | `version` | Semantic version string |
    | `url` | Direct download URL to the plugin ZIP |
    | `date` | ISO 8601 timestamp |
    | `sha512sum` | SHA-512 checksum of the ZIP file |
    | `requires` | Minimum Nextflow version (e.g., `>=24.04.0`) |

### 3.4. Configure Nextflow to use your index

Set the environment variable before running Nextflow:

```bash
export NXF_PLUGINS_TEST_REPOSITORY="https://internal.example.com/plugins/plugins.json"
```

Then use the plugin as normal:

```groovy title="nextflow.config"
plugins {
    id 'nf-myplugin@1.0.0'
}
```

!!! tip "Setting the variable permanently"

    Add the export to your shell profile (`~/.bashrc`, `~/.zshrc`) or set it in your CI/CD pipeline configuration.

---

## 4. Versioning best practices

Follow semantic versioning for your releases:

| Version change            | When to use                       | Example                                    |
| ------------------------- | --------------------------------- | ------------------------------------------ |
| **MAJOR** (1.0.0 → 2.0.0) | Breaking changes                  | Removing a function, changing return types |
| **MINOR** (1.0.0 → 1.1.0) | New features, backward compatible | Adding a new function                      |
| **PATCH** (1.0.0 → 1.0.1) | Bug fixes, backward compatible    | Fixing a bug in existing function          |

---

## 5. Advanced extension types

Some extension types require significant infrastructure or deep Nextflow knowledge to implement.
This section provides a brief conceptual overview.
For implementation details, see the [Nextflow plugin documentation](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html).

### 5.1. Executors

Executors define how tasks are submitted to compute resources:

- AWS Batch, Google Cloud Batch, Azure Batch
- Kubernetes, SLURM, PBS, LSF
- Creating a custom executor is complex and typically done by platform vendors

### 5.2. Filesystems

Filesystems define how files are accessed:

- S3, Google Cloud Storage, Azure Blob
- Custom storage systems
- Creating a custom filesystem requires implementing Java NIO interfaces

---

## Takeaway

You learned that:

- Use the public registry for open source plugins (requires claiming a name)
- For internal distribution, host plugin ZIPs and a `plugins.json` index, then set `NXF_PLUGINS_TEST_REPOSITORY`
- Use semantic versioning to communicate changes to users
- Executors and filesystems are advanced extension types typically created by platform vendors

---

## What's next?

You've completed the plugin development training!

[Continue to Summary :material-arrow-right:](summary.md){ .md-button .md-button--primary }
