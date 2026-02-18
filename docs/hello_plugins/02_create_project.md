# Part 2: Create a Plugin Project

You've seen how plugins extend Nextflow with reusable functionality.
Now you'll create your own, starting with a project template that handles the build configuration for you.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 1 to use as your starting point:

    ```bash
    cp -r solutions/1-plugin-basics/* .
    ```

!!! info "Official documentation"

    This section and those that follow cover plugin development essentials.
    For comprehensive details, see the [official Nextflow plugin development documentation](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html).

---

## 1. Create the plugin project

The built-in `nextflow plugin create` command generates a complete plugin project:

```bash
nextflow plugin create nf-greeting training
```

The first argument is the plugin name, and the second is your organization name (used to organize the generated code into folders).

!!! tip "Manual creation"

    You can also create plugin projects manually or use the [nf-hello template](https://github.com/nextflow-io/nf-hello) on GitHub as a starting point.

---

## 2. Examine the project structure

Change into the plugin directory:

```bash
cd nf-greeting
```

List the contents:

```bash
tree
```

You should see:

```console
.
├── build.gradle
├── COPYING
├── gradle
│   └── wrapper
│       ├── gradle-wrapper.jar
│       └── gradle-wrapper.properties
├── gradlew
├── Makefile
├── README.md
├── settings.gradle
└── src
    ├── main
    │   └── groovy
    │       └── training
    │           └── plugin
    │               ├── GreetingExtension.groovy
    │               ├── GreetingFactory.groovy
    │               ├── GreetingObserver.groovy
    │               └── GreetingPlugin.groovy
    └── test
        └── groovy
            └── training
                └── plugin
                    └── GreetingObserverTest.groovy

11 directories, 13 files
```

---

## 3. Explore the build configuration

### 3.1. settings.gradle

This file identifies the project:

```bash
cat settings.gradle
```

```groovy title="settings.gradle"
rootProject.name = 'nf-greeting'
```

The name here must match what you'll put in `nextflow.config` when using the plugin.

### 3.2. build.gradle

The build file is where most configuration happens:

```bash
cat build.gradle
```

The file contains several sections.
The most important is the `nextflowPlugin` block:

```groovy title="build.gradle"
plugins {
    id 'io.nextflow.nextflow-plugin' version '1.0.0-beta.10'
}

version = '0.1.0'

nextflowPlugin {
    nextflowVersion = '24.10.0'

    provider = 'training'
    className = 'training.plugin.GreetingPlugin'
    extensionPoints = [
        'training.plugin.GreetingExtension',
        'training.plugin.GreetingFactory'
    ]

}
```

The `nextflowPlugin` block configures:

- `nextflowVersion`: Minimum Nextflow version required
- `provider`: Your name or organization
- `className`: The main plugin class (the entry point that Nextflow loads first, specified in `build.gradle`)
- `extensionPoints`: Classes that add features to Nextflow (your functions, monitoring, etc.)

### 3.3. Update nextflowVersion

The template generates a `nextflowVersion` value that may be outdated.
Update it to match your installed Nextflow version for full compatibility:

=== "After"

    ```groovy title="build.gradle" hl_lines="2"
    nextflowPlugin {
        nextflowVersion = '25.10.0'

        provider = 'training'
    ```

=== "Before"

    ```groovy title="build.gradle" hl_lines="2"
    nextflowPlugin {
        nextflowVersion = '24.10.0'

        provider = 'training'
    ```

---

## 4. Tour the source files

The plugin code lives in `src/main/groovy/training/plugin/`.
Examine each file to see what the template generated:

```bash
cat src/main/groovy/training/plugin/GreetingPlugin.groovy
```

`GreetingPlugin` is the plugin entry point that you specified in `build.gradle` via `className`.
When Nextflow loads your plugin, this is the first class it instantiates.
You won't need to modify this file; it's generated automatically.

```bash
cat src/main/groovy/training/plugin/GreetingExtension.groovy
```

The extension class holds functions marked with `@Function` that become callable from Nextflow workflows.
This is where you'll add most of your plugin's functionality (Part 3).

```bash
cat src/main/groovy/training/plugin/GreetingFactory.groovy
```

The factory is responsible for creating observer instances when a workflow starts.
You'll work with this in Part 5 when you add a task counter.

```bash
cat src/main/groovy/training/plugin/GreetingObserver.groovy
```

The observer runs code when things happen in the workflow, like when a task completes or the pipeline finishes.
The template includes messages that print "Pipeline is starting!" and "Pipeline complete!" (Part 5).

??? info "How the components relate"

    ```mermaid
    graph TD
        A[GreetingPlugin] -->|registers| B[GreetingExtension]
        A -->|registers| C[GreetingFactory]
        C -->|creates| D[GreetingObserver]

        B -->|provides| E["@Function methods<br/>(callable from workflows)"]
        D -->|hooks into| F["Lifecycle events<br/>(onFlowCreate, onProcessComplete, etc.)"]

        style A fill:#e1f5fe
        style B fill:#fff3e0
        style C fill:#fff3e0
        style D fill:#fff3e0
    ```

---

## 5. Build, install, and run

Compile the plugin and install it locally:

```bash
make assemble && make install
```

??? example "Build output"

    The first time you run this, Gradle will download itself (this may take a minute):

    ```console
    Downloading https://services.gradle.org/distributions/gradle-8.14-bin.zip
    ...10%...20%...30%...40%...50%...60%...70%...80%...90%...100%

    Welcome to Gradle 8.14!
    ...

    Deprecated Gradle features were used in this build...

    BUILD SUCCESSFUL in 23s
    4 actionable tasks: 4 executed
    ```

    **The warnings are expected.**

    - **"Downloading gradle..."**: This only happens the first time. Subsequent builds are much faster.
    - **"Deprecated Gradle features..."**: This warning comes from the plugin template, not your code. It's safe to ignore.
    - **"BUILD SUCCESSFUL"**: This is what matters. Your plugin compiled without errors.

??? info "What is `./gradlew`?"

    The `./gradlew` script is the **Gradle wrapper**, a small script included with the project that automatically downloads and runs the correct version of Gradle.

    This means you don't need Gradle installed on your system.
    The `make` commands in the Makefile are shortcuts that call `./gradlew` for you.

Go back to the pipeline directory:

```bash
cd ..
```

Add the nf-greeting plugin to `nextflow.config`:

=== "After"

    ```groovy title="nextflow.config" hl_lines="4"
    // Configuration for plugin development exercises
    plugins {
        id 'nf-schema@2.6.1'
        id 'nf-greeting@0.1.0'
    }
    ```

=== "Before"

    ```groovy title="nextflow.config"
    // Configuration for plugin development exercises
    plugins {
        id 'nf-schema@2.6.1'
    }
    ```

!!! note "Version required for local plugins"

    When using locally installed plugins, you must specify the version (e.g., `nf-greeting@0.1.0`).
    Published plugins in the registry can use just the name.

Run the pipeline:

```bash
nextflow run greet.nf
```

```console title="Output"
Pipeline is starting! 🚀
Output: Hello
Output: Bonjour
Output: Holà
Output: Ciao
Output: Hallo
Pipeline complete! 👋
```

(Your output order may differ.)

The "Pipeline is starting!" and "Pipeline complete!" messages come from the `GreetingObserver` in your plugin.
The pipeline itself is unchanged; the observer runs automatically because it's registered in the factory.

---

## Takeaway

You learned that:

- The `nextflow plugin create` command generates a complete starter project
- `build.gradle` configures plugin metadata, dependencies, and which classes provide features
- The plugin has four main components: Plugin (entry point), Extension (functions), Factory (creates monitors), and Observer (responds to workflow events)
- The development cycle is: edit code, `make assemble`, `make install`, run the pipeline

---

## What's next?

Now you'll implement custom functions in the Extension class and use them in the workflow.

[Continue to Part 3 :material-arrow-right:](03_custom_functions.md){ .md-button .md-button--primary }
