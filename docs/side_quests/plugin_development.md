# Plugin development

Nextflow's plugin system allows you to extend the language with custom functions, operators, executors, and more.
In this side quest, you'll build a simple plugin from scratch, learning the fundamentals of plugin architecture along the way.

!!! warning "Advanced topic"

    Plugin development is one of the more advanced Nextflow topics.
    It involves Java/Groovy programming, build tools, and software engineering concepts that may be unfamiliar if you come from a pure bioinformatics background.

    **Most Nextflow users will never need to develop plugins** - the existing plugin ecosystem and Nextflow's built-in features cover the vast majority of use cases.
    This side quest is for those who want to understand how plugins work or have specific needs that require custom extensions.

    If you find this material challenging, that's completely normal!
    Consider bookmarking it for later and focusing on the [Using existing plugins](#1-using-existing-plugins) section for now.

### Learning goals

In this side quest, you'll learn how to use existing Nextflow plugins and create your own custom plugin.

By the end of this side quest, you'll be able to:

- Install and use existing plugins in your workflows
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

??? info "What are Java, Groovy, and Gradle?"

    If these terms are unfamiliar, here's a quick primer:

    **Java** is a widely-used programming language. Nextflow itself is built with Java, and plugins must be compatible with the Java runtime.

    **Groovy** is a programming language that runs on Java and is designed to be more concise and flexible. Nextflow's DSL is based on Groovy, which is why Nextflow syntax looks the way it does. Plugin code is typically written in Groovy.

    **Gradle** is a build tool that compiles code, runs tests, and packages software. Think of it like `make` but for Java/Groovy projects. You don't need to understand Gradle deeply - we'll use simple commands like `./gradlew build`.

    The good news: you don't need to be an expert in any of these. We'll explain the relevant concepts as we go.

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

## 1. Using existing plugins

Before diving into plugin development, let's understand how to use plugins that others have created.
Nextflow has a growing ecosystem of plugins that extend its functionality.

### Why use plugins instead of local functions?

You can define custom functions directly in your Nextflow scripts, so why use plugins?

| Approach | Best for | Limitations |
| -------- | -------- | ----------- |
| **Local functions** | Project-specific logic | Copy-paste between pipelines, no versioning |
| **Plugins** | Reusable utilities | Requires Java/Groovy knowledge to create |

Plugins are ideal when you have functions that:

- Are useful across multiple pipelines
- Need to be shared with the community
- Require versioning and dependency management
- Need access to Nextflow internals (channels, sessions, etc.)

!!! tip "This is the most important section for most users"

    Even if you never develop your own plugin, knowing how to use existing plugins is valuable.
    Many powerful features - like input validation with nf-schema - come from plugins.
    If plugin development seems daunting, focus on mastering this section first.

### 1.1. Installing plugins

Plugins are declared in your `nextflow.config` file using the `plugins {}` block:

```groovy title="nextflow.config"
plugins {
    id 'nf-schema@2.1.1'
}
```

Key points:

- Use the `id` keyword followed by the plugin name
- Specify a version with `@version` (recommended for reproducibility)
- Nextflow automatically downloads plugins from the plugin registry

### 1.2. Importing plugin functions

Once a plugin is installed, you can import its functions using the familiar `include` syntax with a special `plugin/` prefix:

```groovy title="main.nf"
include { samplesheetToList } from 'plugin/nf-schema'
```

This imports the `samplesheetToList` function from the nf-schema plugin, making it available in your workflow.

### 1.3. Example: Using nf-schema for validation

The nf-schema plugin is widely used in nf-core pipelines for input validation.
Here's how it works in practice:

```groovy title="main.nf" linenums="1"
#!/usr/bin/env nextflow

include { samplesheetToList } from 'plugin/nf-schema'

params.input = 'samplesheet.csv'

workflow {
    // Validate and parse input samplesheet
    ch_samples = Channel.fromList(
        samplesheetToList(params.input, "assets/schema_input.json")
    )

    ch_samples.view { "Sample: $it" }
}
```

The `samplesheetToList` function:

1. Reads the input CSV file
2. Validates it against a JSON schema
3. Returns a list of validated entries
4. Throws helpful errors if validation fails

This pattern is used extensively in nf-core pipelines to ensure input data is valid before processing begins.

### 1.4. Popular community plugins

Here are some useful plugins available in the Nextflow ecosystem:

| Plugin        | Purpose                                       |
| ------------- | --------------------------------------------- |
| nf-schema     | Input validation and samplesheet parsing      |
| nf-prov       | Provenance reporting (RO-Crate format)        |
| nf-wave       | Container provisioning with Wave              |
| nf-amazon     | AWS integration (S3, Batch)                   |
| nf-google     | Google Cloud integration (GCS, Batch)         |
| nf-azure      | Azure integration (Blob Storage, Batch)       |
| nf-cloudcache | Cloud-based caching for distributed execution |

!!! tip "Finding plugins"

    Browse available plugins in the [Nextflow plugin registry](https://www.nextflow.io/docs/latest/plugins/plugin-registry.html) or search GitHub for repositories with the `nf-` prefix.

### 1.5. Plugin configuration

Some plugins accept configuration options in `nextflow.config`:

```groovy title="nextflow.config"
plugins {
    id 'nf-schema@2.1.1'
}

// Plugin-specific configuration
validation {
    monochromeLogs = true
    ignoreParams = ['custom_param']
}
```

Each plugin documents its configuration options.
Check the plugin's documentation for available settings.

### 1.6. Try it: From local function to plugin

Let's see the difference between a local function and a plugin function in practice.

Create a new directory for this exercise:

```bash
mkdir -p /tmp/plugin-test && cd /tmp/plugin-test
```

#### Using a local function

First, create a workflow with a locally defined `sayHello` function:

```groovy title="main.nf"
#!/usr/bin/env nextflow

// Local function - defined in this file
def sayHello(name) {
    return "Hello, ${name}!"
}

workflow {
    Channel.of('Alice', 'Bob', 'Carol')
        | map { name -> sayHello(name) }
        | view
}
```

Run it:

```bash
nextflow run main.nf
```

```console title="Output"
Hello, Alice!
Hello, Bob!
Hello, Carol!
```

This works fine, but if you wanted to use `sayHello` in another pipeline, you'd have to copy the function definition.

#### Using a plugin instead

Now let's replace our local function with the [nf-hello](https://github.com/nextflow-io/nf-hello) plugin, which provides the same functionality.

Create a `nextflow.config` to enable the plugin:

```groovy title="nextflow.config"
plugins {
    id 'nf-hello@0.5.0'
}
```

Update `main.nf` to import the function from the plugin instead of defining it locally:

=== "After (plugin)"

    ```groovy title="main.nf" hl_lines="3-4"
    #!/usr/bin/env nextflow

    // Import function from plugin - no local definition needed
    include { sayHello } from 'plugin/nf-hello'

    workflow {
        Channel.of('Alice', 'Bob', 'Carol')
            | map { name -> sayHello(name) }
            | view
    }
    ```

=== "Before (local)"

    ```groovy title="main.nf" hl_lines="3-6"
    #!/usr/bin/env nextflow

    // Local function - defined in this file
    def sayHello(name) {
        return "Hello, ${name}!"
    }

    workflow {
        Channel.of('Alice', 'Bob', 'Carol')
            | map { name -> sayHello(name) }
            | view
    }
    ```

Run it again:

```bash
nextflow run main.nf
```

The first run will download the plugin automatically. The output is the same:

```console title="Output"
Hello, Alice!
Hello, Bob!
Hello, Carol!
```

The key difference: now the `sayHello` function comes from a versioned, shareable plugin rather than copy-pasted code. Any pipeline can use `nf-hello@0.5.0` and get the exact same function.

### Takeaway

Using plugins is straightforward: declare them in `nextflow.config`, import their functions, and use them in your workflows.
The plugin ecosystem extends Nextflow with powerful features like validation, cloud integration, and provenance tracking.

### What's next?

Now that you understand how to use plugins, let's explore how they work under the hood.

---

## 2. Plugin architecture

### 2.1. How plugins extend Nextflow

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

### 2.2. Plugin project structure

A typical plugin project looks like this:

```
nf-greeting/
├── build.gradle          # Build configuration
├── settings.gradle       # Project settings
├── gradlew               # Gradle wrapper script
├── Makefile              # Convenience commands
└── src/
    ├── main/
    │   └── groovy/
    │       └── training/plugin/
    │           ├── NfGreetingPlugin.groovy    # Main plugin class
    │           ├── NfGreetingExtension.groovy # Extension with functions
    │           ├── NfGreetingFactory.groovy   # Channel factory (optional)
    │           └── NfGreetingObserver.groovy  # Trace observer (optional)
    └── test/
        └── groovy/
            └── training/plugin/
                └── NfGreetingObserverTest.groovy
```

The package name (`training/plugin`) comes from the organization name you provide when creating the plugin.

### 2.3. Key components

**Plugin class**: The entry point that manages plugin lifecycle.

**Extension classes**: Contain the actual functionality (functions, operators, etc.).

**Build configuration**: Gradle scripts that compile and package the plugin.

### Takeaway

Plugins extend Nextflow through well-defined extension points.
The plugin system uses standard Java/Groovy tooling.

### What's next?

Let's create our plugin project.

---

## 3. Creating a plugin project

### 3.1. Using the Nextflow plugin create command

The easiest way to create a plugin is with the built-in command:

```bash
nextflow plugin create nf-greeting training
```

This scaffolds a complete plugin project.
The first argument is the plugin name, and the second is your organization name (used for the package namespace).

!!! tip "Manual creation"

    You can also create plugin projects manually or use the [nf-hello template](https://github.com/nextflow-io/nf-hello) on GitHub as a starting point.

### 3.2. Examine the generated project

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
    │               ├── NfGreetingExtension.groovy
    │               ├── NfGreetingFactory.groovy
    │               ├── NfGreetingObserver.groovy
    │               └── NfGreetingPlugin.groovy
    └── test
        └── groovy
            └── training
                └── plugin
                    └── NfGreetingObserverTest.groovy

11 directories, 13 files
```

### 3.3. Understand settings.gradle

```bash
cat settings.gradle
```

```groovy title="settings.gradle"
rootProject.name = 'nf-greeting'
```

This simply sets the project name.

### 3.4. Understand build.gradle

```bash
cat build.gradle
```

Key sections in the build file:

```groovy title="build.gradle"
plugins {
    id 'io.nextflow.nextflow-plugin' version '0.0.1-alpha4'
}

version = '0.1.0'

nextflowPlugin {
    nextflowVersion = '24.10.0'

    provider = 'training'
    className = 'training.plugin.NfGreetingPlugin'
    extensionPoints = [
        'training.plugin.NfGreetingExtension',
        'training.plugin.NfGreetingFactory'
    ]

    publishing {
        registry {
            url = 'https://nf-plugins-registry.dev-tower.net/api'
            authToken = project.findProperty('pluginRegistry.accessToken')
        }
    }
}
```

The `nextflowPlugin` block configures:

- `nextflowVersion`: Minimum Nextflow version required
- `provider`: Your name or organization
- `className`: The main plugin class (uses your package name)
- `extensionPoints`: Classes providing extensions (functions, factories, etc.)
- `publishing`: Configuration for publishing to the plugin registry (optional)

### Takeaway

The `nextflow plugin create` command scaffolds a complete project.
The `build.gradle` file configures the plugin metadata and dependencies.

### What's next?

Let's implement our custom functions.

---

## 4. Implementing custom functions

### 4.1. The PluginExtensionPoint class

Functions are defined in classes that extend `PluginExtensionPoint`.
Open the extension file:

```bash
cat src/main/groovy/training/plugin/NfGreetingExtension.groovy
```

The template includes sample functions. Let's replace them with our greeting functions.

### 4.2. Create our functions

The template includes a simple `sayHello` function.
Let's replace it with our greeting manipulation functions.

Edit the file to replace the `sayHello` function with our three new functions:

=== "After"

    ```groovy title="src/main/groovy/training/plugin/NfGreetingExtension.groovy" hl_lines="28-50" linenums="1"
    /*
     * Copyright 2025, Seqera Labs
     *
     * Licensed under the Apache License, Version 2.0 (the "License");
     * you may not use this file except in compliance with the License.
     * You may obtain a copy of the License at
     *
     *     http://www.apache.org/licenses/LICENSE-2.0
     *
     * Unless required by applicable law or agreed to in writing, software
     * distributed under the License is distributed on an "AS IS" BASIS,
     * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
     * See the License for the specific language governing permissions and
     * limitations under the License.
     */

    package training.plugin

    import groovy.transform.CompileStatic
    import nextflow.Session
    import nextflow.plugin.extension.Function
    import nextflow.plugin.extension.PluginExtensionPoint

    @CompileStatic
    class NfGreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
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
         * Convert greeting to a friendly format with a name
         */
        @Function
        String friendlyGreeting(String greeting, String name = 'World') {
            return "${greeting}, ${name}!"
        }

    }
    ```

=== "Before"

    ```groovy title="src/main/groovy/training/plugin/NfGreetingExtension.groovy" hl_lines="28-42" linenums="1"
    /*
     * Copyright 2025, Seqera Labs
     *
     * Licensed under the Apache License, Version 2.0 (the "License");
     * you may not use this file except in compliance with the License.
     * You may obtain a copy of the License at
     *
     *     http://www.apache.org/licenses/LICENSE-2.0
     *
     * Unless required by applicable law or agreed to in writing, software
     * distributed under the License is distributed on an "AS IS" BASIS,
     * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
     * See the License for the specific language governing permissions and
     * limitations under the License.
     */

    package training.plugin

    import groovy.transform.CompileStatic
    import nextflow.Session
    import nextflow.plugin.extension.Function
    import nextflow.plugin.extension.PluginExtensionPoint

    /**
     * Implements a custom function which can be imported by
     * Nextflow scripts.
     */
    @CompileStatic
    class NfGreetingExtension extends PluginExtensionPoint {

        @Override
        protected void init(Session session) {
        }

        /**
         * Say hello to the given target.
         *
         * @param target
         */
        @Function
        void sayHello(String target) {
            println "Hello, ${target}!"
        }

    }
    ```

??? info "Understanding the Groovy syntax"

    If this code looks unfamiliar, here's a breakdown of the key elements:

    **`package training.plugin`** - Declares which package (folder structure) this code belongs to. This must match the directory structure.

    **`import ...`** - Brings in code from other packages, similar to Python's `import` or R's `library()`.

    **`@CompileStatic`** - An annotation (marked with `@`) that tells Groovy to check types at compile time. This catches errors earlier.

    **`class NfGreetingExtension extends PluginExtensionPoint`** - Defines a class that inherits from `PluginExtensionPoint`. The `extends` keyword means "this class is a type of that class."

    **`@Override`** - Indicates we're replacing a method from the parent class.

    **`@Function`** - The key annotation that makes a method available as a Nextflow function.

    **`String reverseGreeting(String greeting)`** - A method that takes a String parameter and returns a String. In Groovy, you can often omit `return` - the last expression is returned automatically.

    **`String name = 'World'`** - A parameter with a default value, just like in Python.

### 4.3. Understanding the @Function annotation

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

### 4.4. The init() method

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

## 5. Building and testing

### 5.1. Build the plugin

The Makefile provides convenient commands:

```bash
make assemble
```

Or directly with the Gradle wrapper:

```bash
./gradlew assemble
```

??? info "What is `./gradlew`?"

    The `./gradlew` script is the **Gradle wrapper** - a small script included with the project that automatically downloads and runs the correct version of Gradle.

    This means you don't need Gradle installed on your system.
    The first time you run `./gradlew`, it will download Gradle (which may take a moment), then run your command.

    The `make` commands in the Makefile are just shortcuts that call `./gradlew` for you.

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

    **Don't worry about the warnings!**

    - **"Downloading gradle..."**: This only happens the first time. Subsequent builds are much faster.
    - **"Deprecated Gradle features..."**: This warning comes from the plugin template, not your code. It's safe to ignore.
    - **"BUILD SUCCESSFUL"**: This is what matters! Your plugin compiled without errors.

### 5.2. Write unit tests

Good plugins have tests.
Tests verify that your code works correctly and help catch bugs when you make changes later.

??? info "What are unit tests?"

    **Unit tests** are small pieces of code that automatically check if your functions work correctly.
    Each test calls a function with known inputs and checks that the output matches what you expect.

    For example, if you have a function that reverses strings, a test might check that `reverse("Hello")` returns `"olleH"`.

    Tests are valuable because:

    - They catch bugs before users do
    - They give you confidence to make changes without breaking things
    - They serve as documentation showing how functions should be used

    You don't need to write tests to use a plugin, but they're good practice for any code you plan to share or maintain.

The generated project includes a test for the Observer class, but we need to create a new test file for our extension functions.

Create the new test file:

```bash
touch src/test/groovy/training/plugin/NfGreetingExtensionTest.groovy
```

Open it in your editor and add the following content:

```groovy title="src/test/groovy/training/plugin/NfGreetingExtensionTest.groovy" linenums="1"
package training.plugin

import spock.lang.Specification

/**
 * Tests for the greeting extension functions
 */
class NfGreetingExtensionTest extends Specification {

    def 'should reverse a greeting'() {
        given:
        def ext = new NfGreetingExtension()

        expect:
        ext.reverseGreeting('Hello') == 'olleH'
        ext.reverseGreeting('Bonjour') == 'ruojnoB'
    }

    def 'should decorate a greeting'() {
        given:
        def ext = new NfGreetingExtension()

        expect:
        ext.decorateGreeting('Hello') == '*** Hello ***'
    }

    def 'should create friendly greeting with default name'() {
        given:
        def ext = new NfGreetingExtension()

        expect:
        ext.friendlyGreeting('Hello') == 'Hello, World!'
    }

    def 'should create friendly greeting with custom name'() {
        given:
        def ext = new NfGreetingExtension()

        expect:
        ext.friendlyGreeting('Hello', 'Alice') == 'Hello, Alice!'
    }
}
```

!!! note "Spock testing framework"

    The generated project uses [Spock](https://spockframework.org/), a testing framework for Groovy.
    Spock tests use descriptive method names in quotes (like `'should reverse a greeting'`) and a `given`/`expect` structure that reads almost like plain English.

### 5.3. Run the tests

```bash
make test
```

Or:

```bash
./gradlew test
```

??? example "Test output"

    ```console
    BUILD SUCCESSFUL in 5s
    6 actionable tasks: 2 executed, 4 up-to-date
    ```

    **Where are the test results?** Gradle hides detailed output when all tests pass - "BUILD SUCCESSFUL" means everything worked! If any test fails, you'll see detailed error messages.

### 5.4. View the test report

To see detailed results for each test, you can view the HTML test report that Gradle generates.

Start a simple web server in the test report directory:

```bash
cd build/reports/tests/test
python -m http.server
```

VS Code will prompt you to open the application in your browser.
Click through to your test class to see individual test results:

![Test report showing all tests passed](./img/test_report.png)

The report shows each test method, its duration, and whether it passed or failed.
This confirms that all four of our greeting functions are being tested correctly.

Press ++ctrl+c++ in the terminal to stop the server when you're done, then return to the plugin directory:

```bash
cd ../../../..
```

!!! tip "If the build fails"

    Build errors can be intimidating, but they usually point to a specific problem.
    Common issues include:

    - **Syntax errors**: A missing bracket, quote, or semicolon. The error message usually includes a line number.
    - **Import errors**: A class name is misspelled or the import statement is missing.
    - **Type errors**: You're passing the wrong type of data to a function.

    Read the error message carefully - it often tells you exactly what's wrong and where.
    If you're stuck, compare your code character-by-character with the examples.

### 5.4. Install locally

To use the plugin with Nextflow, install it to your local plugins directory:

```bash
make install
```

This copies the plugin to `~/.nextflow/plugins/`.

### Takeaway

Use `make assemble` to compile and `make test` to run tests.
Install with `make install` to use the plugin locally.

### What's next?

Let's use our plugin in a workflow.

---

## 6. Using your plugin

### 6.1. Configure the plugin

Go back to the pipeline directory:

```bash
cd ..
```

Edit `nextflow.config` to enable the plugin by uncommenting the plugins block and adding the version:

=== "After"

    ```groovy title="nextflow.config" hl_lines="3"
    // Plugin development example config
    plugins {
        id 'nf-greeting@0.1.0'
    }
    ```

=== "Before"

    ```groovy title="nextflow.config"
    // Plugin development example config
    // Uncomment the following when your plugin is ready:
    // plugins {
    //     id 'nf-greeting'
    // }
    ```

!!! note "Version required for local plugins"

    When using locally installed plugins, you must specify the version (e.g., `nf-greeting@0.1.0`).
    Published plugins in the registry can use just the name.

### 6.2. Import and use functions

Edit `main.nf` to import and use our custom functions:

=== "After"

    ```groovy title="main.nf" hl_lines="4-5 18-19 30-33" linenums="1"
    #!/usr/bin/env nextflow

    // Import custom functions from our plugin
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'

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

=== "Before"

    ```groovy title="main.nf" linenums="1"
    #!/usr/bin/env nextflow

    /*
     * Simple greeting pipeline
     * Will be enhanced to use custom plugin functions
     */

    params.input = 'greetings.csv'

    process SAY_HELLO {

        input:
            val greeting

        output:
            path "${greeting}-output.txt"

        script:
        """
        echo '$greeting' > '${greeting}-output.txt'
        """
    }

    workflow {

        greeting_ch = channel.fromPath(params.input)
                            .splitCsv(header: true)
                            .map { row -> row.greeting }

        SAY_HELLO(greeting_ch)

        SAY_HELLO.out.view()
    }
    ```

The key changes:

- **Lines 4-5**: Import our plugin functions using `include { function } from 'plugin/plugin-name'`
- **Lines 18-19**: Use `decorateGreeting()` in the process script to wrap the greeting
- **Lines 30-33**: Use `reverseGreeting()` in the workflow to demonstrate channel operations with plugin functions

### 6.3. Run the pipeline

```bash
nextflow run main.nf
```

??? example "Output"

    ```console
    N E X T F L O W   ~  version 25.04.3

    Launching `main.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    /workspaces/.../work/.../Holà-output.txt
    /workspaces/.../work/.../Hello-output.txt
    ...
    Pipeline complete! 🎉
    ```

    The "Pipeline is starting!" and "Pipeline complete!" messages come from the `NfGreetingObserver` trace observer that was included in the generated plugin template.

### 6.4. Check the output files

```bash
cat work/*/*/Hello-output.txt
```

```console title="Output"
*** Hello ***
```

The `decorateGreeting()` function was applied in the process script, wrapping the greeting with decorative markers.

### Takeaway

Import plugin functions with `include { function } from 'plugin/plugin-id'`.
Functions work in both workflow blocks and process scripts.

### What's next?

Let's explore other extension types.

---

## 7. Going further

The `@Function` annotation covers most common use cases, but plugins can do much more.
Let's explore some advanced capabilities using the code already in your generated plugin.

### 7.1. Trace observers: How the startup messages work

Remember the "Pipeline is starting! 🚀" message when you ran the pipeline?
That came from the `NfGreetingObserver` class in your plugin.

Look at the observer code:

```bash
cat src/main/groovy/training/plugin/NfGreetingObserver.groovy
```

You'll see something like:

```groovy
class NfGreetingObserver implements TraceObserver {

    @Override
    void onFlowCreate(Session session) {
        println "Pipeline is starting! 🚀"
    }

    @Override
    void onFlowComplete() {
        println "Pipeline complete! 🎉"
    }
}
```

Trace observers let you hook into workflow lifecycle events - useful for custom logging, metrics, or notifications.

#### Try it: Customize the messages

Let's change the messages to something more descriptive.

Edit `src/main/groovy/training/plugin/NfGreetingObserver.groovy` and change the `println` statements:

=== "After"

    ```groovy hl_lines="4 9"
    @Override
    void onFlowCreate(Session session) {
        // Custom startup message
        println "🔬 Starting analysis pipeline..."
    }

    @Override
    void onFlowComplete() {
        // Custom completion message
        println "✅ Analysis complete! Check your results."
    }
    ```

=== "Before"

    ```groovy hl_lines="3 8"
    @Override
    void onFlowCreate(Session session) {
        println "Pipeline is starting! 🚀"
    }

    @Override
    void onFlowComplete() {
        println "Pipeline complete! 🎉"
    }
    ```

Rebuild and reinstall the plugin:

```bash
make assemble && make install
```

Run the pipeline again from the parent directory:

```bash
cd .. && nextflow run main.nf
```

You should see your custom messages in the output:

```console
🔬 Starting analysis pipeline...
executor >  local (5)
[ab/123456] process > SAY_HELLO (5) [100%] 5 of 5 ✔
...
✅ Analysis complete! Check your results.
```

### 7.2. Available observer hooks

Trace observers can respond to many events:

| Method | When it's called |
| ------ | ---------------- |
| `onFlowCreate` | Workflow starts |
| `onFlowComplete` | Workflow finishes |
| `onProcessCreate` | A process is defined |
| `onProcessStart` | A task begins execution |
| `onProcessComplete` | A task finishes |
| `onProcessCached` | A cached task is reused |
| `onFilePublish` | A file is published |

This enables powerful use cases like custom reports, Slack notifications, or metrics collection.

### 7.3. Custom operators (reference)

Operators work with channels directly. Use the `@Operator` annotation:

```groovy
import nextflow.plugin.extension.Operator
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel

@Operator
DataflowWriteChannel reverseAll(DataflowReadChannel source) {
    def target = CH.create()
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

!!! note "Operators are more complex"

    Creating operators requires understanding Nextflow's internal channel APIs.
    For most use cases, `@Function` is simpler and sufficient.

### 7.4. Accessing configuration

Plugins can read custom configuration from `nextflow.config`:

```groovy title="nextflow.config"
greeting {
    prefix = '>>>'
    suffix = '<<<'
}
```

```groovy title="In your extension's init() method"
@Override
void init(Session session) {
    def prefix = session.config.navigate('greeting.prefix', '***')
    def suffix = session.config.navigate('greeting.suffix', '***')
}
```

This lets users customize plugin behavior without modifying plugin code.

### 7.5. Publishing your plugin

To share your plugin with others:

1. Build a release: `make release`
2. Push to GitHub with a release tag
3. Register with the Nextflow plugin registry

Once published, users can install without local builds:

```groovy
plugins {
    id 'nf-greeting'  // No version needed for registry plugins
}
```

!!! tip "Plugin registry"

    The Nextflow plugin registry is in public preview.
    See the [Nextflow documentation](https://www.nextflow.io/docs/latest/plugins/publishing-plugins.html) for publishing details.

### Takeaway

Beyond functions, plugins can provide trace observers for lifecycle hooks, custom operators for channel transformations, and configuration-driven behavior.
The generated plugin template includes working examples of each.

### What's next?

Let's summarize what we've learned.

---

## Summary

### Plugin development checklist

- [ ] Java 17+ installed
- [ ] Create project with `nextflow plugin create <name> <org>`
- [ ] Implement extension class with `@Function` methods
- [ ] Write unit tests
- [ ] Build with `make assemble`
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
