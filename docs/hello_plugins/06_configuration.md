# Part 6: Configuration & Distribution

Users should be able to control your plugin from `nextflow.config` without editing code.
In this section, you'll add configuration options to your plugin using two approaches, then learn how to share your plugin.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 5 to use as your starting point:

    ```bash
    cp -r solutions/5-observers/* .
    ```

!!! info "Official documentation"

    For comprehensive configuration details, see the [Nextflow config scopes documentation](https://nextflow.io/docs/latest/developer/config-scopes.html).

Nextflow provides two approaches for plugin configuration:

| Approach                    | Best for                           | Trade-offs                              |
| --------------------------- | ---------------------------------- | --------------------------------------- |
| `session.config.navigate()` | Quick prototyping, simple plugins  | No IDE support, manual type conversion  |
| `@ConfigScope` classes      | Production plugins, complex config | More code, but type-safe and documented |

We'll start with the simple approach, then upgrade to the formal approach.

---

## 1. Simple configuration with navigate()

The `session.config.navigate()` method reads nested configuration values:

```groovy
// Read 'greeting.enabled' from nextflow.config, defaulting to true
final enabled = session.config.navigate('greeting.enabled', true)
```

This lets users control plugin behavior:

```groovy title="nextflow.config"
greeting {
    enabled = false
}
```

---

## 2. Make the task counter configurable

Add configuration options to:

1. Enable/disable the entire greeting plugin
2. Control whether per-task counter messages are shown

### 2.1. Update TaskCounterObserver

Edit `TaskCounterObserver.groovy` to accept a configuration flag:

=== "After"

    ```groovy title="TaskCounterObserver.groovy" linenums="1" hl_lines="14 17-19 24-26"
    package training.plugin

    import groovy.transform.CompileStatic
    import nextflow.processor.TaskHandler
    import nextflow.trace.TraceObserver
    import nextflow.trace.TraceRecord

    /**
     * Observer that counts completed tasks
     */
    @CompileStatic
    class TaskCounterObserver implements TraceObserver {

        private final boolean verbose
        private int taskCount = 0

        TaskCounterObserver(boolean verbose) {
            this.verbose = verbose
        }

        @Override
        void onProcessComplete(TaskHandler handler, TraceRecord trace) {
            taskCount++
            if (verbose) {
                println "📊 Tasks completed so far: ${taskCount}"
            }
        }

        @Override
        void onFlowComplete() {
            println "📈 Final task count: ${taskCount}"
        }
    }
    ```

=== "Before"

    ```groovy title="TaskCounterObserver.groovy" linenums="1" hl_lines="19"
    package training.plugin

    import groovy.transform.CompileStatic
    import nextflow.processor.TaskHandler
    import nextflow.trace.TraceObserver
    import nextflow.trace.TraceRecord

    /**
     * Observer that counts completed tasks
     */
    @CompileStatic
    class TaskCounterObserver implements TraceObserver {

        private int taskCount = 0

        @Override
        void onProcessComplete(TaskHandler handler, TraceRecord trace) {
            taskCount++
            println "📊 Tasks completed so far: ${taskCount}"
        }

        @Override
        void onFlowComplete() {
            println "📈 Final task count: ${taskCount}"
        }
    }
    ```

The key changes:

- **Line 14**: A `verbose` flag controls whether per-task messages are printed
- **Lines 17-19**: Constructor that accepts the verbose setting
- **Lines 24-26**: Only print per-task messages if `verbose` is true

### 2.2. Update the Factory

Update `GreetingFactory.groovy` to read the configuration and pass it to the observer:

=== "After"

    ```groovy title="GreetingFactory.groovy" linenums="31" hl_lines="3-6 9"
    @Override
    Collection<TraceObserver> create(Session session) {
        final enabled = session.config.navigate('greeting.enabled', true)
        if (!enabled) return []

        final verbose = session.config.navigate('greeting.taskCounter.verbose', true) as boolean
        return [
            new GreetingObserver(),
            new TaskCounterObserver(verbose)
        ]
    }
    ```

=== "Before"

    ```groovy title="GreetingFactory.groovy" linenums="31"
    @Override
    Collection<TraceObserver> create(Session session) {
        return [
            new GreetingObserver(),
            new TaskCounterObserver()
        ]
    }
    ```

The factory now:

- **Lines 33-34**: Reads the `greeting.enabled` config and returns early if disabled
- **Line 36**: Reads the `greeting.taskCounter.verbose` config (defaulting to `true`)
- **Line 39**: Passes the verbose setting to the `TaskCounterObserver` constructor

### 2.3. Build and test

Rebuild and reinstall the plugin:

```bash
cd nf-greeting && make assemble && make install && cd ..
```

Update `nextflow.config` to disable the per-task messages:

=== "After"

    ```groovy title="nextflow.config" linenums="1" hl_lines="5-7"
    plugins {
        id 'nf-greeting@0.1.0'
    }

    greeting {
        // enabled = false        // Disable plugin entirely
        taskCounter.verbose = false  // Disable per-task messages
    }
    ```

=== "Before"

    ```groovy title="nextflow.config" linenums="1"
    plugins {
        id 'nf-greeting@0.1.0'
    }
    ```

Run the pipeline and observe that only the final count appears:

```bash
nextflow run main.nf -ansi-log false
```

??? example "Output"

    ```console
    N E X T F L O W  ~  version 25.10.2
    Launching `main.nf` [stoic_wegener] DSL2 - revision: 63f3119fbc
    WARN: Unrecognized config option 'greeting.taskCounter.verbose'
    Pipeline is starting! 🚀
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    [5e/9c1f21] Submitted process > SAY_HELLO (2)
    [20/8f6f91] Submitted process > SAY_HELLO (1)
    [6d/496bae] Submitted process > SAY_HELLO (4)
    [5c/a7fe10] Submitted process > SAY_HELLO (3)
    [48/18199f] Submitted process > SAY_HELLO (5)
    Decorated: *** Hello ***
    Decorated: *** Bonjour ***
    Decorated: *** Holà ***
    Decorated: *** Ciao ***
    Decorated: *** Hallo ***
    Pipeline complete! 👋
    📈 Final task count: 5
    ```

    Note the "Unrecognized config option" warning on line 3.

!!! note

    The "Unrecognized config option" warning appears because Nextflow doesn't know about the `greeting` scope yet.
    The config values still work via `session.config.navigate()`, but Nextflow flags them as unrecognized.
    This goes away in Section 4 when you register a formal config scope class.

??? exercise "Disable the plugin entirely"

    Try setting `greeting.enabled = false` in `nextflow.config` and run the pipeline again.
    What changes in the output?

    ??? solution

        ```groovy title="nextflow.config" hl_lines="6"
        plugins {
            id 'nf-greeting@0.1.0'
        }

        greeting {
            enabled = false
        }
        ```

        The "Pipeline is starting!", "Pipeline complete!", and task count messages all disappear because the factory returns an empty list when `enabled` is false.
        The pipeline itself still runs, but no observers are active.

        Remember to set `enabled` back to `true` (or remove the line) before continuing.

---

## 3. Make the decorator configurable

This exercise makes the `decorateGreeting` function use configurable prefix/suffix.

### 3.1. Add the configuration reading (this will fail!)

Edit `GreetingExtension.groovy` to read configuration in `init()` and use it in `decorateGreeting()`:

```groovy title="GreetingExtension.groovy" linenums="35" hl_lines="7-8 18"
@CompileStatic
class GreetingExtension extends PluginExtensionPoint {

    @Override
    protected void init(Session session) {
        // Read configuration with defaults
        prefix = session.config.navigate('greeting.prefix', '***') as String
        suffix = session.config.navigate('greeting.suffix', '***') as String
    }

    // ... other methods unchanged ...

    /**
    * Decorate a greeting with celebratory markers
    */
    @Function
    String decorateGreeting(String greeting) {
        return "${prefix} ${greeting} ${suffix}"
    }
```

Try to build:

```bash
cd nf-greeting && make assemble
```

### 3.2. Observe the error

The build fails:

```console
> Task :compileGroovy FAILED
GreetingExtension.groovy: 30: [Static type checking] - The variable [prefix] is undeclared.
 @ line 30, column 9.
           prefix = session.config.navigate('greeting.prefix', '***') as String
           ^

GreetingExtension.groovy: 31: [Static type checking] - The variable [suffix] is undeclared.
```

In Groovy (and Java), you must _declare_ a variable before using it.
We're trying to assign values to `prefix` and `suffix`, but we never told the class that these variables exist.

### 3.3. Fix by declaring instance variables

Add the variable declarations at the top of the class, right after the opening brace:

```groovy title="GreetingExtension.groovy" linenums="35" hl_lines="4-5"
@CompileStatic
class GreetingExtension extends PluginExtensionPoint {

    private String prefix = '***'
    private String suffix = '***'

    @Override
    protected void init(Session session) {
        // Read configuration with defaults
        prefix = session.config.navigate('greeting.prefix', '***') as String
        suffix = session.config.navigate('greeting.suffix', '***') as String
    }

    // ... rest of class unchanged ...
```

The `private String prefix = '***'` line does two things:

1. **Declares** a variable named `prefix` that can hold a String
2. **Initializes** it with a default value of `'***'`

Now the `init()` method can assign new values to these variables, and `decorateGreeting()` can read them.

### 3.4. Build again

```bash
make assemble
```

This time it should succeed with "BUILD SUCCESSFUL".

```bash
make install && cd ..
```

!!! tip "Learning from errors"

    This "declare before use" pattern is fundamental to Java/Groovy but unfamiliar if you come from Python or R where variables spring into existence when you first assign them.
    Experiencing this error once helps you recognize and fix it quickly in the future.

### 3.5. Test the configurable decorator

Update `nextflow.config` to customize the decoration:

=== "After"

    ```groovy title="nextflow.config" hl_lines="7-8"
    plugins {
        id 'nf-greeting@0.1.0'
    }

    greeting {
        taskCounter.verbose = false
        prefix = '>>>'
        suffix = '<<<'
    }
    ```

=== "Before"

    ```groovy title="nextflow.config"
    plugins {
        id 'nf-greeting@0.1.0'
    }

    greeting {
        taskCounter.verbose = false
    }
    ```

Run the pipeline:

```bash
nextflow run main.nf -ansi-log false
```

```console title="Output (partial)"
Decorated: >>> Hello <<<
Decorated: >>> Bonjour <<<
...
```

The decorator now uses your custom prefix and suffix.

---

## 4. Formal configuration with ConfigScope

The `session.config.navigate()` approach works, but has limitations:

- No IDE autocompletion for users writing `nextflow.config`
- Configuration options aren't self-documenting
- Manual type conversion with `as String`, `as boolean`, etc.

For production plugins, Nextflow provides a formal configuration system using annotations.
By creating a config scope class, you define a new top-level configuration block, equivalent to how built-in scopes like `process {}` or `docker {}` work.

### 4.1. Create the config class (minimal version)

Create a new file:

```bash
touch nf-greeting/src/main/groovy/training/plugin/GreetingConfig.groovy
```

Start with a minimal config class for the `enabled`, `prefix`, and `suffix` options:

```groovy title="GreetingConfig.groovy" linenums="1"
package training.plugin

import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description

/**
 * Configuration options for the nf-greeting plugin.
 *
 * Users configure these in nextflow.config:
 *
 *     greeting {
 *         enabled = true
 *         prefix = '>>>'
 *         suffix = '<<<'
 *     }
 */
@ScopeName('greeting')
class GreetingConfig implements ConfigScope {

    GreetingConfig() {}

    GreetingConfig(Map opts) {
        this.enabled = opts.enabled as Boolean ?: true
        this.prefix = opts.prefix as String ?: '***'
        this.suffix = opts.suffix as String ?: '***'
    }

    @ConfigOption
    @Description('Enable or disable the plugin entirely')
    boolean enabled = true

    @ConfigOption
    @Description('Prefix for decorated greetings')
    String prefix = '***'

    @ConfigOption
    @Description('Suffix for decorated greetings')
    String suffix = '***'
}
```

Key points:

- **`@ScopeName('greeting')`**: Maps to the `greeting { }` block in config
- **`implements ConfigScope`**: Required interface for config classes
- **`@ConfigOption`**: Each field becomes a configuration option
- **`@Description`**: Documents each option for language server support (imported from `nextflow.script.dsl`)
- **Constructors**: Both no-arg and Map constructors are needed

### 4.2. Add nested configuration

Now add the `taskCounter` nested scope for the verbose option:

```groovy title="GreetingConfig.groovy (final version)" linenums="1" hl_lines="29-31 46 48-57"
package training.plugin

import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description

/**
 * Configuration options for the nf-greeting plugin.
 *
 * Users configure these in nextflow.config:
 *
 *     greeting {
 *         enabled = true
 *         prefix = '>>>'
 *         suffix = '<<<'
 *         taskCounter.verbose = false
 *     }
 */
@ScopeName('greeting')
class GreetingConfig implements ConfigScope {

    GreetingConfig() {}

    GreetingConfig(Map opts) {
        this.enabled = opts.enabled as Boolean ?: true
        this.prefix = opts.prefix as String ?: '***'
        this.suffix = opts.suffix as String ?: '***'
        if (opts.taskCounter instanceof Map) {
            this.taskCounter = new TaskCounterConfig(opts.taskCounter as Map)
        }
    }

    @ConfigOption
    @Description('Enable or disable the plugin entirely')
    boolean enabled = true

    @ConfigOption
    @Description('Prefix for decorated greetings')
    String prefix = '***'

    @ConfigOption
    @Description('Suffix for decorated greetings')
    String suffix = '***'

    TaskCounterConfig taskCounter = new TaskCounterConfig()

    static class TaskCounterConfig implements ConfigScope {
        TaskCounterConfig() {}
        TaskCounterConfig(Map opts) {
            this.verbose = opts.verbose as Boolean ?: true
        }

        @ConfigOption
        @Description('Show per-task completion messages')
        boolean verbose = true
    }
}
```

For nested paths like `taskCounter.verbose`, use a nested class that also implements `ConfigScope`.

### 4.3. Register the config class

Every class that implements `ExtensionPoint` needs to be listed in `extensionPoints` in `build.gradle` for the Nextflow plugin system to find it.

Update `build.gradle`:

=== "After"

    ```groovy title="build.gradle" hl_lines="4"
    extensionPoints = [
        'training.plugin.GreetingExtension',
        'training.plugin.GreetingFactory',
        'training.plugin.GreetingConfig'
    ]
    ```

=== "Before"

    ```groovy title="build.gradle"
    extensionPoints = [
        'training.plugin.GreetingExtension',
        'training.plugin.GreetingFactory'
    ]
    ```

Note the difference between the factory and extension points registration:

- **`extensionPoints` in `build.gradle`**: Compile-time registration. Tells the Nextflow plugin system which classes implement extension points.
- **Factory `create()` method**: Runtime registration. The factory creates observer instances when a workflow actually starts.

### 4.4. Build and test

The `ConfigScope` class doesn't change how your code reads configuration at runtime; your factory and extension still use `session.config.navigate()`.
What it does is register the `greeting {}` scope with Nextflow, which eliminates the "Unrecognized config option" warning and enables tooling support.

```bash
cd nf-greeting && make assemble && make install && cd ..
nextflow run main.nf -ansi-log false
```

The pipeline behavior is the same, but the warning is gone and your configuration options are now formally documented within the plugin.

---

## 5. Distribution

Once your plugin is working locally, you can share it with others through the Nextflow plugin registry.

### 5.1. Versioning

Follow [semantic versioning](https://semver.org/) for your releases:

| Version change            | When to use                       | Example                                    |
| ------------------------- | --------------------------------- | ------------------------------------------ |
| **MAJOR** (1.0.0 → 2.0.0) | Breaking changes                  | Removing a function, changing return types |
| **MINOR** (1.0.0 → 1.1.0) | New features, backward compatible | Adding a new function                      |
| **PATCH** (1.0.0 → 1.0.1) | Bug fixes, backward compatible    | Fixing a bug in existing function          |

Update the version in `build.gradle` before each release:

```groovy title="build.gradle"
version = '1.0.0'  // Use semantic versioning: MAJOR.MINOR.PATCH
```

### 5.2. Publishing to the registry

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

- `session.config.navigate()` provides simple, quick configuration reading
- `@ScopeName` and `@ConfigOption` annotations create self-documenting configuration
- Configuration can be applied to both observers and extension functions
- Variables must be declared before use in Groovy/Java
- Use semantic versioning and the public registry to distribute plugins

| Use case                            | Recommended approach                      |
| ----------------------------------- | ----------------------------------------- |
| Quick prototype or simple plugin    | `session.config.navigate()` only          |
| Production plugin with many options | Add `ConfigScope` class for documentation |
| Plugin you'll share publicly        | Add `ConfigScope` class for documentation |

---

## What's next?

You've completed the plugin development training.

[Continue to Summary :material-arrow-right:](summary.md){ .md-button .md-button--primary }
