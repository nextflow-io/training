# Summary

You have completed the Hello Plugins training.
This page recaps what you built in each part, covers distribution, and provides guidance on where to go next.

---

## What you learned

### Part 1: Using plugins

You discovered how Nextflow plugins work from a user's perspective.
You installed nf-schema and nf-co2footprint, configured them via `nextflow.config`, and saw how plugins can validate inputs, add functions, and hook into pipeline lifecycle events.

### Part 2: Setting up

You set up a plugin development environment with Java 21+, created a new plugin project using the `nextflow plugin create` command, and learned the project structure that Nextflow expects: source files, build configuration, and the Makefile workflow.

### Part 3: Custom functions

You implemented your first extension point by creating `@Function`-annotated methods in a `PluginExtensionPoint` class.
You built `reverseGreeting` and `decorateGreeting`, then imported and called them from a pipeline script.

### Part 4: Testing

You wrote unit tests for your custom functions using the Groovy testing framework.
You learned how to run tests with `make test` and verify that your plugin behaves correctly before installing it.

### Part 5: Observers

You implemented the `TraceObserver` interface to hook into pipeline lifecycle events.
You built `GreetingObserver` (reacting to pipeline start and completion) and `TaskCounterObserver` (counting completed tasks), then registered them through a `TraceObserverFactory`.

### Part 6: Configuration

You made your plugin configurable through `nextflow.config` using `session.config.navigate()` to read values at runtime.
You added a `@ConfigScope` class to formally declare your plugin's options, eliminating the "Unrecognized config option" warnings and enabling IDE support.

---

## Distribution

Once your plugin is working locally, you can share it with others through the Nextflow plugin registry.

### Versioning

Follow [semantic versioning](https://semver.org/) for your releases:

| Version change             | When to use                       | Example                                    |
| -------------------------- | --------------------------------- | ------------------------------------------ |
| **MAJOR** (1.0.0 → 2.0.0) | Breaking changes                  | Removing a function, changing return types |
| **MINOR** (1.0.0 → 1.1.0) | New features, backward compatible | Adding a new function                      |
| **PATCH** (1.0.0 → 1.0.1) | Bug fixes, backward compatible    | Fixing a bug in existing function          |

Update the version in `build.gradle` before each release:

```groovy title="build.gradle"
version = '1.0.0'  // Use semantic versioning: MAJOR.MINOR.PATCH
```

### Publishing to the registry

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

## Plugin development checklist

- [ ] Java 21+ installed
- [ ] Create project with `nextflow plugin create <name> <org>`
- [ ] Implement extension class with `@Function` methods
- [ ] Build with `make assemble`
- [ ] Write unit tests and run with `make test`
- [ ] Install with `make install`
- [ ] Optionally add `TraceObserver` implementations for workflow events
- [ ] Optionally add `ConfigScope` for plugin configuration
- [ ] Enable in `nextflow.config` with `plugins { id 'plugin-id' }`
- [ ] Import functions with `include { fn } from 'plugin/plugin-id'`
- [ ] Version and publish to the registry

---

## Key code patterns

**Function definition:**

```groovy
@Function
String myFunction(String input, String optional = 'default') {
    return input.transform()
}
```

**Plugin configuration:**

```groovy
nextflowPlugin {
    provider = 'my-org'
    className = 'my.org.MyPlugin'
    extensionPoints = ['my.org.MyExtension']
}
```

**Using in workflows:**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { item -> myFunction(item) }
        .view()
}
```

---

## Extension point summary

| Type                | Class/Annotation | Purpose                                        |
| ------------------- | ---------------- | ---------------------------------------------- |
| Function            | `@Function`      | Callable from workflows                        |
| Trace Observer      | `TraceObserver`  | Hook into workflow lifecycle events            |
| Configuration Scope | `@ScopeName`     | Define plugin configuration in nextflow.config |

---

## What to do next

Here are some practical next steps for continuing your plugin development journey.

**Build something real.**
Pick a use case from your own work: a custom function that your team uses repeatedly, an observer that sends Slack notifications on pipeline completion, or a config scope that standardizes options across your organization's pipelines.
Starting from a real problem is the fastest way to deepen your understanding.

**Use nf-hello as a reference.**
The [nf-hello](https://github.com/nextflow-io/nf-hello) repository is the official minimal plugin example.
It is a good starting point for new projects and a useful reference when you need to check how something is structured.

**Read the official documentation.**
The Nextflow docs cover topics beyond this training, including channel factories, operator overloading, and advanced observer patterns.
The [developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html) guide is the most comprehensive reference.

**Study existing plugins.**
The [Nextflow plugins repository](https://github.com/nextflow-io/plugins) contains the source code for official plugins like nf-schema, nf-wave, and nf-tower.
Reading production plugin code is one of the best ways to learn patterns and conventions that go beyond introductory examples.

---

## Additional resources

**Official documentation:**

- [Using plugins](https://www.nextflow.io/docs/latest/plugins/plugins.html): comprehensive guide to installing and configuring plugins
- [Developing plugins](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html): detailed plugin development reference
- [Config scopes](https://nextflow.io/docs/latest/developer/config-scopes.html): creating configuration scopes for plugins

**Plugin discovery:**

- [Nextflow Plugin Registry](https://registry.nextflow.io/): browse and discover available plugins
- [Plugin registry docs](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html): registry documentation

**Examples and references:**

- [nf-hello](https://github.com/nextflow-io/nf-hello): simple example plugin (great starting point)
- [Nextflow plugins repository](https://github.com/nextflow-io/plugins): collection of official plugins for reference
