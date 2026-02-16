# Part 7: Distribution

Once your plugin is working locally, you can share it with others through the Nextflow plugin registry.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 6 to use as your starting point:

    ```bash
    cp -r solutions/6-configuration/* .
    ```

---

## 1. Versioning

Follow [semantic versioning](https://semver.org/) for your releases:

| Version change             | When to use                       | Example                                   |
| -------------------------- | --------------------------------- | ----------------------------------------- |
| **MAJOR** (1.0.0 → 2.0.0) | Breaking changes                  | Removing a function, changing return types |
| **MINOR** (1.0.0 → 1.1.0) | New features, backward compatible | Adding a new function                     |
| **PATCH** (1.0.0 → 1.0.1) | Bug fixes, backward compatible    | Fixing a bug in existing function         |

Update the version in `build.gradle` before each release:

```groovy title="build.gradle"
version = '1.0.0'  // Use semantic versioning: MAJOR.MINOR.PATCH
```

---

## 2. Publishing to the registry

The [Nextflow plugin registry](https://registry.nextflow.io/) is the official way to share plugins with the community.

The publishing workflow:

1. Claim your plugin name on the [registry](https://registry.nextflow.io/) (sign in with your GitHub account)
2. Configure your API credentials in `~/.gradle/gradle.properties`
3. Run tests to verify everything works: `make test`
4. Publish with `make release`

For step-by-step instructions, see the [official publishing documentation](https://www.nextflow.io/docs/latest/guides/gradle-plugin.html#publishing-a-plugin).

Once published, users install your plugin without any local setup:

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting@1.0.0'
}
```

Nextflow automatically downloads the plugin from the registry on first use.

---

## Takeaway

You learned that:

- Use semantic versioning to communicate changes to users
- The public registry is the standard way to distribute plugins
- Publishing requires claiming a name and configuring API credentials

---

## What's next?

You've completed the plugin development training.

[Continue to Summary :material-arrow-right:](summary.md){ .md-button .md-button--primary }
