# Plugin development

Nextflow's plugin system allows you to extend the language with custom functions, operators, executors, and more.
In this side quest, you'll build a simple plugin from scratch, learning the fundamentals of plugin architecture along the way.

### Learning goals

In this side quest, you'll learn how to create, build, and use a custom Nextflow plugin.

By the end of this side quest, you'll be able to:

- Understand the Nextflow plugin architecture
- Create a new plugin project
- Implement custom functions using the `@Function` annotation
- Build and test your plugin locally
- Use your plugin in a Nextflow workflow
- Understand other extension points (operators, trace observers)

### Prerequisites

Before taking on this side quest, you should:

- Have completed the [Hello Nextflow](../../hello_nextflow/) tutorial or equivalent beginner's course
- Have Java 17 or later installed (check with `java -version`)
- Have basic familiarity with object-oriented programming concepts

!!! note "Development environment"

    This side quest requires Java and Gradle for building plugins.
    The training Codespace comes with Java pre-installed.

---

## 0. Get started

#### Open the training codespace

If you haven't yet done so, make sure to open the training environment as described in the [Environment Setup](../envsetup/index.md).

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/nextflow-io/training?quickstart=1&ref=master)

#### Verify Java installation

Check that Java is available:

```bash
java -version
```

You should see Java 17 or later.

#### Move into the project directory

```bash
cd side-quests/plugin_development
```

#### Review the materials

```console title="Directory contents"
.
├── greetings.csv
├── main.nf
└── nextflow.config
```

We have a simple greeting pipeline that we'll enhance with custom plugin functions.

#### What we'll build

We'll create a plugin called `nf-greeting` that provides functions to manipulate greeting strings:

- `reverseGreeting()` - Reverses a greeting
- `decorateGreeting()` - Adds decorative markers
- `friendlyGreeting()` - Creates a friendly greeting with a name

#### Readiness checklist

- [ ] My codespace is running
- [ ] Java is installed and working
- [ ] I understand we're building a plugin from scratch

---

## 1. Plugin architecture

### 1.1. How plugins extend Nextflow

Nextflow's plugin system is built on [PF4J](https://pf4j.org/), a lightweight plugin framework for Java.
Plugins can extend Nextflow in several ways:

| Extension Type  | Purpose                                  | Example                 |
| --------------- | ---------------------------------------- | ----------------------- |
| Functions       | Custom functions callable from workflows | `reverseString()`       |
| Operators       | Custom channel operators                 | `myFilter()`            |
| Factories       | Create new channel types                 | `mySource()`            |
| Executors       | Custom task execution backends           | AWS Batch, Kubernetes   |
| Filesystems     | Custom storage backends                  | S3, Azure Blob          |
| Trace Observers | Monitor workflow execution               | Custom logging, metrics |

### 1.2. Plugin project structure

A typical plugin project looks like this:

```
nf-greeting/
├── build.gradle          # Build configuration
├── settings.gradle       # Project settings
├── Makefile             # Convenience commands
└── src/
    ├── main/
    │   └── groovy/
    │       └── nextflow/greeting/
    │           ├── GreetingPlugin.groovy    # Main plugin class
    │           └── GreetingExtension.groovy # Extension with functions
    ├── test/
    │   └── groovy/
    │       └── nextflow/greeting/
    │           └── GreetingExtensionTest.groovy
    └── resources/
        └── META-INF/
            └── MANIFEST.MF   # Plugin metadata
```

### 1.3. Key components

**Plugin class**: The entry point that manages plugin lifecycle.

**Extension classes**: Contain the actual functionality (functions, operators, etc.).

**Build configuration**: Gradle scripts that compile and package the plugin.

### Takeaway

Plugins extend Nextflow through well-defined extension points.
The plugin system uses standard Java/Groovy tooling.

### What's next?

Let's create our plugin project.

---

## 2. Creating a plugin project

### 2.1. Using the Nextflow plugin create command

The easiest way to create a plugin is with the built-in command:

```bash
nextflow plugin create nf-greeting
```

This scaffolds a complete plugin project.

!!! tip "Manual creation"

    You can also create plugin projects manually or use the [nf-hello template](https://github.com/nextflow-io/nf-hello) on GitHub as a starting point.

### 2.2. Examine the generated project

Change into the plugin directory:

```bash
cd nf-greeting
```

List the contents:

```bash
ls -la
```

You should see:

```console
.
├── build.gradle
├── settings.gradle
├── Makefile
└── src/
    └── ...
```

### 2.3. Understand settings.gradle

```bash
cat settings.gradle
```

```groovy title="settings.gradle"
rootProject.name = 'nf-greeting'
```

This simply sets the project name.

### 2.4. Understand build.gradle

```bash
cat build.gradle
```

Key sections in the build file:

```groovy title="build.gradle (key parts)"
plugins {
    id 'io.nextflow.nextflow-plugin' version '1.0.0-beta.6'
}

version = '0.1.0'

nextflowPlugin {
    nextflowVersion = '25.04.0'
    provider = 'training'
    className = 'nextflow.greeting.GreetingPlugin'
    extensionPoints = [
        'nextflow.greeting.GreetingExtension'
    ]
}
```

The `nextflowPlugin` block configures:

- `nextflowVersion`: Minimum Nextflow version required
- `provider`: Your name or organization
- `className`: The main plugin class
- `extensionPoints`: Classes providing extensions (functions, operators)

### Takeaway

The `nextflow plugin create` command scaffolds a complete project.
The `build.gradle` file configures the plugin metadata and dependencies.

### What's next?

Let's implement our custom functions.

---

## 3. Implementing custom functions

### 3.1. The PluginExtensionPoint class

Functions are defined in classes that extend `PluginExtensionPoint`.
Open the extension file:

```bash
cat src/main/groovy/nextflow/greeting/GreetingExtension.groovy
```

The template includes a sample function. Let's replace it with our greeting functions.

### 3.2. Create our functions

Edit the file to contain:

```groovy title="src/main/groovy/nextflow/greeting/GreetingExtension.groovy" linenums="1"
package nextflow.greeting

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.PluginExtensionPoint

/**
 * Plugin extension providing custom functions for greeting manipulation
 */
@CompileStatic
class GreetingExtension extends PluginExtensionPoint {

    @Override
    void init(Session session) {
        // Initialization logic if needed
    }

    /**
     * Reverse a greeting string
     */
    @Function
    String reverseGreeting(String greeting) {
        return greeting.reverse()
    }

    /**
     * Decorate a greeting with celebratory markers
     */
    @Function
    String decorateGreeting(String greeting) {
        return "*** ${greeting} ***"
    }

    /**
     * Convert greeting to a friendly format with exclamation
     */
    @Function
    String friendlyGreeting(String greeting, String name = 'World') {
        return "${greeting}, ${name}!"
    }
}
```

### 3.3. Understanding the @Function annotation

The `@Function` annotation marks a method as callable from Nextflow workflows:

```groovy
@Function
String reverseGreeting(String greeting) {
    return greeting.reverse()
}
```

Key points:

- Methods must be public (default in Groovy)
- Can have parameters with default values
- Return type can be any serializable type
- Will be available via `include { reverseGreeting } from 'plugin/nf-greeting'`

### 3.4. The init() method

The `init()` method is called when the plugin loads:

```groovy
@Override
void init(Session session) {
    // Access session configuration
    // Initialize resources
    // Set up state
}
```

You can access configuration via `session.config`.

### Takeaway

Functions are defined with the `@Function` annotation in `PluginExtensionPoint` subclasses.
They become available to import in Nextflow workflows.

### What's next?

Let's build and test our plugin.

---

## 4. Building and testing

### 4.1. Build the plugin

The Makefile provides convenient commands:

```bash
make build
```

Or directly with Gradle:

```bash
./gradlew build
```

??? example "Build output"

    ```console
    > Task :compileJava NO-SOURCE
    > Task :compileGroovy
    > Task :processResources
    > Task :classes
    > Task :jar
    > Task :assemble
    > Task :compileTestJava NO-SOURCE
    > Task :compileTestGroovy
    > Task :processTestResources NO-SOURCE
    > Task :testClasses
    > Task :test
    > Task :check
    > Task :build

    BUILD SUCCESSFUL
    ```

### 4.2. Write unit tests

Good plugins have tests. Create or edit the test file:

```groovy title="src/test/groovy/nextflow/greeting/GreetingExtensionTest.groovy" linenums="1"
package nextflow.greeting

import org.junit.jupiter.api.Test
import static org.junit.jupiter.api.Assertions.*

class GreetingExtensionTest {

    @Test
    void testReverseGreeting() {
        def ext = new GreetingExtension()
        assertEquals('olleH', ext.reverseGreeting('Hello'))
        assertEquals('ruojnoB', ext.reverseGreeting('Bonjour'))
    }

    @Test
    void testDecorateGreeting() {
        def ext = new GreetingExtension()
        assertEquals('*** Hello ***', ext.decorateGreeting('Hello'))
    }

    @Test
    void testFriendlyGreetingDefault() {
        def ext = new GreetingExtension()
        assertEquals('Hello, World!', ext.friendlyGreeting('Hello'))
    }

    @Test
    void testFriendlyGreetingWithName() {
        def ext = new GreetingExtension()
        assertEquals('Hello, Alice!', ext.friendlyGreeting('Hello', 'Alice'))
    }
}
```

### 4.3. Run the tests

```bash
make test
```

Or:

```bash
./gradlew test
```

??? example "Test output"

    ```console
    > Task :test

    GreetingExtensionTest > testReverseGreeting() PASSED
    GreetingExtensionTest > testDecorateGreeting() PASSED
    GreetingExtensionTest > testFriendlyGreetingDefault() PASSED
    GreetingExtensionTest > testFriendlyGreetingWithName() PASSED

    BUILD SUCCESSFUL
    ```

### 4.4. Install locally

To use the plugin with Nextflow, install it to your local plugins directory:

```bash
make install
```

This copies the plugin to `~/.nextflow/plugins/`.

### Takeaway

Use `make build` to compile and `make test` to run tests.
Install with `make install` to use the plugin locally.

### What's next?

Let's use our plugin in a workflow.

---

## 5. Using your plugin

### 5.1. Configure the plugin

Go back to the pipeline directory:

```bash
cd ..
```

Edit `nextflow.config` to enable the plugin:

```groovy title="nextflow.config"
plugins {
    id 'nf-greeting@0.1.0'
}
```

!!! note "Version required for local plugins"

    When using locally installed plugins, you must specify the version (e.g., `nf-greeting@0.1.0`).
    Published plugins in the registry can use just the name.

### 5.2. Import and use functions

Edit `main.nf` to use our custom functions:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

// Import custom functions from our plugin
include { reverseGreeting; decorateGreeting } from 'plugin/nf-greeting'

params.input = 'greetings.csv'

process SAY_HELLO {

    input:
        val greeting

    output:
        path "${greeting}-output.txt"

    script:
    // Use our custom plugin function to decorate the greeting
    def decorated = decorateGreeting(greeting)
    """
    echo '$decorated' > '${greeting}-output.txt'
    """
}

workflow {

    greeting_ch = channel.fromPath(params.input)
                        .splitCsv(header: true)
                        .map { row -> row.greeting }

    // Demonstrate using reverseGreeting function
    greeting_ch
        .map { greeting -> reverseGreeting(greeting) }
        .view { "Reversed: $it" }

    SAY_HELLO(greeting_ch)

    SAY_HELLO.out.view()
}
```

### 5.3. Run the pipeline

```bash
nextflow run main.nf
```

??? example "Output"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `main.nf` [happy_euler] DSL2 - revision: abc123

    executor >  local (5)
    [12/abc123] SAY_HELLO (1) | 5 of 5 ✔

    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    /path/to/work/.../Hello-output.txt
    /path/to/work/.../Bonjour-output.txt
    ...
    ```

### 5.4. Check the output files

```bash
cat work/*/*/Hello-output.txt
```

You should see:

```
*** Hello ***
```

The `decorateGreeting()` function was applied in the process script.

### Takeaway

Import plugin functions with `include { function } from 'plugin/plugin-id'`.
Functions work in both workflow blocks and process scripts.

### What's next?

Let's explore other extension types.

---

## 6. Going further

### 6.1. Custom operators

Operators work with channels. Use `@Operator` annotation:

```groovy
import nextflow.plugin.extension.Operator
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel

@Operator
DataflowWriteChannel reverseAll(DataflowReadChannel source) {
    def target = CH.create()
    // Transform each item
    new SubscribeOp()
        .withSource(source)
        .withTarget(target)
        .withTransform { it.reverse() }
        .apply()
    return target
}
```

Usage in workflow:

```groovy
greeting_ch.reverseAll().view()
```

### 6.2. Trace observers

Monitor workflow events with `TraceObserverV2`:

```groovy
import nextflow.trace.TraceObserverV2

class MyObserver implements TraceObserverV2 {

    @Override
    void onFlowBegin() {
        log.info "Workflow started!"
    }

    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        log.info "Process ${trace.name} completed"
    }

    @Override
    void onFlowComplete() {
        log.info "Workflow finished!"
    }
}
```

### 6.3. Accessing configuration

Plugins can read custom configuration:

```groovy
// In nextflow.config
greeting {
    prefix = '>>>'
    suffix = '<<<'
}

// In your extension
@Override
void init(Session session) {
    def prefix = session.config.navigate('greeting.prefix', '***')
    def suffix = session.config.navigate('greeting.suffix', '***')
}
```

### 6.4. Publishing your plugin

To share your plugin:

1. Build a release: `make release`
2. Register with the Nextflow plugin registry
3. Users can then install via `plugins { id 'nf-greeting' }` without local installation

!!! tip "Plugin registry"

    The Nextflow plugin registry is in public preview.
    Contact Seqera for access to publish plugins.

### Takeaway

Plugins can provide operators, trace observers, and more.
Configuration is accessible via the Session object.
Published plugins can be shared with the community.

### What's next?

Let's summarize what we've learned.

---

## Summary

### Plugin development checklist

- [ ] Java 17+ installed
- [ ] Create project with `nextflow plugin create`
- [ ] Implement extension class with `@Function` methods
- [ ] Write unit tests
- [ ] Build with `make build`
- [ ] Install with `make install`
- [ ] Enable in `nextflow.config` with `plugins { id 'plugin-id' }`
- [ ] Import functions with `include { fn } from 'plugin/plugin-id'`

### Key code patterns

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
    id = 'my-plugin'
    className = 'my.package.MyPlugin'
    extensionClasses = ['my.package.MyExtension']
}
```

**Using in workflows:**

```groovy
include { myFunction } from 'plugin/my-plugin'

workflow {
    channel.of('a', 'b', 'c')
        .map { myFunction(it) }
        .view()
}
```

### Extension point summary

| Type     | Annotation  | Purpose                 |
| -------- | ----------- | ----------------------- |
| Function | `@Function` | Callable from workflows |
| Operator | `@Operator` | Transform channels      |
| Factory  | `@Factory`  | Create channels         |

### Additional resources

- [Nextflow plugin development docs](https://www.nextflow.io/docs/latest/plugins/developing-plugins.html)
- [Plugin registry](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html)
- [nf-hello example plugin](https://github.com/nextflow-io/nf-hello)
- [Existing plugins](https://github.com/nextflow-io/plugins) for reference

---

## What's next?

Congratulations on completing this side quest!
You've learned how to create Nextflow plugins and extend the language with custom functionality.

Plugin development opens up powerful possibilities for:

- Sharing reusable functions across pipelines
- Integrating with external services
- Custom monitoring and reporting
- Supporting new execution platforms

Return to the [Side Quests](./index.md) menu to continue your training journey.
