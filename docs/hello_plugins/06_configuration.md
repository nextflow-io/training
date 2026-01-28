# Part 6: Configuration

In this section, you'll make your plugin configurable by reading settings from `nextflow.config`.
Users will be able to customize plugin behavior without modifying code.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 5 to use as your starting point:

    ```bash
    cp -r solutions/5-observers/* .
    ```

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

This approach works well for quick prototyping and simple plugins.

---

## 2. Try it: Make the task counter configurable

This exercise adds configuration options to:

1. Enable/disable the entire greeting plugin
2. Control whether per-task counter messages are shown

### 2.1. Update TaskCounterObserver

First, edit `TaskCounterObserver.groovy` to accept a configuration flag:

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

- **Line 14**: Add a `verbose` flag to control whether per-task messages are printed
- **Lines 17-19**: Constructor that accepts the verbose setting
- **Lines 24-26**: Only print per-task messages if `verbose` is true

### 2.2. Update the Factory

Now update `GreetingFactory.groovy` to read the configuration and pass it to the observer:

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

Now update `nextflow.config` to disable the per-task messages:

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

```console title="Expected output"
N E X T F L O W  ~  version 25.10.2
Launching `main.nf` [stoic_wegener] DSL2 - revision: 63f3119fbc
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

---

## 3. Try it: Make the decorator configurable

This exercise makes the `decorateGreeting` function use configurable prefix/suffix.
We'll intentionally make a common mistake to understand how Groovy/Java handles variables.

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

Now try to build:

```bash
cd nf-greeting && make assemble
```

### 3.2. Observe the error

The build fails with an error like:

```console
> Task :compileGroovy FAILED
GreetingExtension.groovy: 30: [Static type checking] - The variable [prefix] is undeclared.
 @ line 30, column 9.
           prefix = session.config.navigate('greeting.prefix', '***') as String
           ^

GreetingExtension.groovy: 31: [Static type checking] - The variable [suffix] is undeclared.
```

**What went wrong?** In Groovy (and Java), you can't just use a variable.
You must _declare_ it first.
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

```console title="Expected output (partial)"
Decorated: >>> Hello <<<
Decorated: >>> Bonjour <<<
...
```

The decorator now uses our custom prefix and suffix.

---

## 4. Formal configuration with ConfigScope

The `session.config.navigate()` approach works, but has limitations:

- No IDE autocompletion for users writing `nextflow.config`
- Configuration options aren't self-documenting
- Manual type conversion with `as String`, `as boolean`, etc.

For production plugins, Nextflow provides a formal configuration system using annotations.
This creates a schema that documents available options.

### 4.1. Understanding the annotations

| Annotation              | Purpose                                               |
| ----------------------- | ----------------------------------------------------- |
| `@ScopeName('name')`    | Declares a configuration block (e.g., `greeting { }`) |
| `@ConfigOption`         | Marks a field as a configuration option               |
| `ConfigScope` interface | Must be implemented by config classes                 |

### 4.2. Create a configuration class

Create a new file `GreetingConfig.groovy`:

```bash
touch nf-greeting/src/main/groovy/training/plugin/GreetingConfig.groovy
```

Add the configuration class:

```groovy title="GreetingConfig.groovy" linenums="1"
package training.plugin

import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName

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

    /**
     * Enable or disable the plugin entirely.
     */
    @ConfigOption
    boolean enabled = true

    /**
     * Prefix for decorated greetings.
     */
    @ConfigOption
    String prefix = '***'

    /**
     * Suffix for decorated greetings.
     */
    @ConfigOption
    String suffix = '***'

    /**
     * Task counter configuration
     */
    TaskCounterConfig taskCounter = new TaskCounterConfig()

    static class TaskCounterConfig implements ConfigScope {
        TaskCounterConfig() {}
        TaskCounterConfig(Map opts) {
            this.verbose = opts.verbose as Boolean ?: true
        }

        @ConfigOption
        boolean verbose = true
    }
}
```

Key points:

- **`@ScopeName('greeting')`**: Maps to the `greeting { }` block in config
- **`implements ConfigScope`**: Required interface for config classes
- **`@ConfigOption`**: Each field becomes a configuration option
- **Nested class**: For nested paths like `taskCounter.verbose`, use a nested class
- **Constructors**: Both no-arg and Map constructors are needed
- **Default values**: Set directly on the fields

### 4.3. Register the config class

Update `build.gradle` to register the config class as an extension point:

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

### 4.4. Build and test

The config class provides documentation and schema validation.
Your code continues using `session.config.navigate()` for reading values:

```bash
cd nf-greeting && make assemble && make install && cd ..
nextflow run main.nf -ansi-log false
```

The behavior is identical, but now your configuration is:

- **Self-documenting**: The config class shows all available options
- **Structured**: Nested configuration is explicit

!!! note "Config class vs runtime access"

    The config class primarily serves as documentation and schema definition.
    Runtime value access still uses `session.config.navigate()`.
    This is the pattern used by most Nextflow plugins including nf-validation.

---

## 5. Summary

| Use case                            | Recommended approach                      |
| ----------------------------------- | ----------------------------------------- |
| Quick prototype or simple plugin    | `session.config.navigate()` only          |
| Production plugin with many options | Add `ConfigScope` class for documentation |
| Plugin you'll share publicly        | Add `ConfigScope` class for documentation |

For this training, the `navigate()` approach is sufficient.
Adding a config class helps users understand available options.

---

## Takeaway

You learned that:

- `session.config.navigate()` provides simple, quick configuration reading
- `@ScopeName` and `@ConfigOption` annotations create self-documenting configuration
- Configuration can be applied to both observers and extension functions
- Variables must be declared before use in Groovy/Java

---

## What's next?

The next section covers how to share your plugin with others.

[Continue to Part 7 :material-arrow-right:](07_distribution.md){ .md-button .md-button--primary }
