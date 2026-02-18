# Part 3: Custom Functions

By the end of this section, you'll have custom functions in your plugin, built and installed locally, running in a real workflow.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 2 to use as your starting point:

    ```bash
    cp -r solutions/2-create-project/* .
    ```

---

## 1. See what the template generated

Change into the plugin directory:

```bash
cd nf-greeting
```

Functions are defined in classes that extend `PluginExtensionPoint`.
Open the extension file to see what the template created:

```bash
cat src/main/groovy/training/plugin/GreetingExtension.groovy
```

The template includes a sample `sayHello` function.
The `@Function` annotation is what makes a method callable from Nextflow workflows.
Without it, the method exists only inside the plugin code.

??? info "Understanding the Groovy syntax"

    If the code looks unfamiliar, here's a breakdown of the key elements:

    **`package training.plugin`**: Declares which package (folder structure) this code belongs to.
    This must match the directory structure.

    **`import ...`**: Brings in code from other packages, similar to Python's `import` or R's `library()`.

    **`@CompileStatic`**: An annotation (marked with `@`) that tells Groovy to check types at compile time.
    This catches errors earlier.

    **`class GreetingExtension extends PluginExtensionPoint`**: Defines a class that inherits from `PluginExtensionPoint`.
    The `extends` keyword means "this class is a type of that class."

    **`void`** means the method doesn't return a value.
    Our functions use `String` as their return type because they return text.

---

## 2. Replace sayHello with reverseGreeting

Replace the template's `sayHello` function with a function that reverses a greeting string.

Edit `src/main/groovy/training/plugin/GreetingExtension.groovy` to replace the `sayHello` method:

=== "After"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="8-14"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

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

    }
    ```

=== "Before"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="12-20"
    /**
     * Implements a custom function which can be imported by
     * Nextflow scripts.
     */
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

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

Key parts of this function:

- **`@Function`**: Makes the method callable from Nextflow workflows
- **`String reverseGreeting(String greeting)`**: Takes a String, returns a String
- **`greeting.reverse()`**: Groovy's built-in string reversal method

!!! tip "Public and private methods"

    Methods without `@Function` are not exposed to Nextflow workflows.
    You can add helper methods to your class without worrying about them leaking into the workflow namespace.

---

## 3. Build and install your plugin

The plugin code must be compiled before Nextflow can use it.

??? info "Why do we need to build?"

    If you're used to scripting languages like Python, R, or even Nextflow's DSL, you might wonder why we need a "build" step.
    In those languages, you write code and run it directly.

    Nextflow plugins are written in Groovy, which runs on the Java Virtual Machine (JVM).
    JVM languages need to be **compiled** before they can run.
    The build tools handle all this automatically.
    Run `make assemble` and let Gradle do the work.

Build the plugin:

```bash
make assemble
```

Install the plugin:

```bash
make install
```

!!! tip "If the build fails"

    Build errors usually point to a specific problem.
    Common issues include:

    - **Syntax errors**: A missing bracket, quote, or semicolon. The error message usually includes a line number.
    - **Import errors**: A class name is misspelled or the import statement is missing.
    - **Type errors**: You're passing the wrong type of data to a function.

    Read the error message carefully.
    It often tells you exactly what's wrong and where.
    If you're stuck, compare your code character-by-character with the examples.

??? warning "Common runtime issues"

    Even if the build succeeds, you might encounter issues when running:

    - **"Plugin not found"**: Did you run `make install`? The plugin must be installed locally before Nextflow can use it.
    - **"Unknown function"**: Check that you've imported the function with `include { functionName } from 'plugin/nf-greeting'`.
    - **Wrong directory**: Make sure you're in the right directory. Use `pwd` to check, and `cd ..` or `cd nf-greeting` as needed.
    - **IDE showing errors**: The VS Code Nextflow extension may show warnings for plugin imports. If the build succeeds and Nextflow runs correctly, you can ignore these.

---

## 4. Use your function in a workflow

Go back to the pipeline directory:

```bash
cd ..
```

Edit `greet.nf` to import and use `reverseGreeting`:

=== "After"

    ```groovy title="greet.nf" hl_lines="4 23-25" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    include { reverseGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
                            .map { row -> row[0] }

        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

=== "Before"

    ```groovy title="greet.nf" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
                            .map { row -> row[0] }
        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

Run the pipeline:

```bash
nextflow run greet.nf
```

??? example "Output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `greet.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    Output: Hello
    Output: Bonjour
    Output: Holà
    Output: Ciao
    Output: Hallo
    Pipeline complete! 👋
    ```

Your first custom plugin function is working in a real workflow.
The same `include { ... } from 'plugin/...'` pattern you used with nf-hello and nf-schema in Part 1 works with your own plugin.

??? exercise "Try a different transformation"

    Change the `map` closure to pass a different string to `reverseGreeting`, like your own name.
    No rebuild is needed since the plugin function is already installed; just edit `greet.nf` and rerun.

    ??? solution

        ```groovy
        greeting_ch
            .map { greeting -> reverseGreeting('Nextflow') }
            .view { reversed -> "Reversed: $reversed" }
        ```

        ```console title="Output (partial)"
        Reversed: wolftxeN
        ```

---

## 5. Add decorateGreeting

Add a second function that wraps a greeting with decorative markers.

Edit `GreetingExtension.groovy` to add `decorateGreeting` after `reverseGreeting`, before the closing brace of the class:

=== "After"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="16-22"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

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

    }
    ```

=== "Before"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @CompileStatic
    class GreetingExtension extends PluginExtensionPoint {

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

    }
    ```

This function uses Groovy string interpolation (`"*** ${greeting} ***"`) to embed the greeting variable inside a string.

Build, install, and update the workflow:

```bash
cd nf-greeting && make assemble && make install && cd ..
```

Update `greet.nf` to also import and use `decorateGreeting`:

=== "After"

    ```groovy title="greet.nf" hl_lines="4-6 16-17 19 33" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    // Import custom functions from our plugin
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        // Use our custom plugin function to decorate the greeting
        def decorated = decorateGreeting(greeting)
        """
        echo '$decorated'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
                            .map { row -> row[0] }

        // Demonstrate using reverseGreeting function
        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Decorated: ${result.trim()}" }
    }
    ```

=== "Before"

    ```groovy title="greet.nf" linenums="1" hl_lines="4 15 28"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    include { reverseGreeting } from 'plugin/nf-greeting'

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            stdout
        script:
        """
        echo '$greeting'
        """
    }

    workflow {
        greeting_ch = Channel.fromList(samplesheetToList(params.input, 'greetings_schema.json'))
                            .map { row -> row[0] }

        greeting_ch
            .map { greeting -> reverseGreeting(greeting) }
            .view { reversed -> "Reversed: $reversed" }

        SAY_HELLO(greeting_ch)
        SAY_HELLO.out.view { result -> "Output: ${result.trim()}" }
    }
    ```

```bash
nextflow run greet.nf
```

??? example "Output"

    ```console
    N E X T F L O W   ~  version 25.10.2

    Launching `greet.nf` [elated_marconi] DSL2 - revision: cd8d52c97c

    Pipeline is starting! 🚀
    executor >  local (5)
    [fe/109754] process > SAY_HELLO (5) [100%] 5 of 5 ✔
    Reversed: olleH
    Reversed: ruojnoB
    Reversed: àloH
    Reversed: oaiC
    Reversed: ollaH
    Decorated: *** Hello ***
    Decorated: *** Bonjour ***
    Decorated: *** Holà ***
    Decorated: *** Ciao ***
    Decorated: *** Hallo ***
    Pipeline complete! 👋
    ```

Plugin functions work in both process scripts (like `decorateGreeting` inside `SAY_HELLO`) and workflow operations (like `reverseGreeting` in a `map`).

??? exercise "Add a function with default parameters"

    Add a third function, `friendlyGreeting`, that combines a greeting with a name.
    Use a default parameter value so the name is optional.

    ??? solution

        Add this method after `decorateGreeting` in `GreetingExtension.groovy`:

        ```groovy
        /**
         * Convert greeting to a friendly format with a name
         */
        @Function
        String friendlyGreeting(String greeting, String name = 'World') {
            return "${greeting}, ${name}!"
        }
        ```

        The `String name = 'World'` syntax provides a default value, just like in Python.
        Users can call `friendlyGreeting('Hello')` or `friendlyGreeting('Hello', 'Alice')`.

        Build and install, then try it in a `map`:

        ```bash
        cd nf-greeting && make assemble && make install && cd ..
        ```

---

## Takeaway

You learned that:

- Functions are defined with the `@Function` annotation in `PluginExtensionPoint` subclasses
- Plugin functions imported with `include` work identically whether from your own plugin or an existing one
- Plugin functions work in both process scripts and workflow operations
- Methods can accept parameters with default values

---

## What's next?

Now you'll write unit tests to verify your functions work correctly.

[Continue to Part 4 :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
