# Part 3: Custom Functions

By the end of this section, you'll have custom functions in your plugin, built and installed locally, running in a real workflow.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 2 to use as your starting point:

    ```bash
    cp -r solutions/2-create-project/* .
    ```

---

## 1. See what the template generated

Before writing your own functions, look at the example function the template created to understand the pattern.

Change into the plugin directory:

```bash
cd nf-greeting
```

The template created a file called `GreetingExtension.groovy` where plugin functions are defined.
Open it to see the starting point:

```bash
cat src/main/groovy/training/plugin/GreetingExtension.groovy
```

```groovy title="Output" hl_lines="29 40-43"
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
class GreetingExtension extends PluginExtensionPoint { // (1)!

    @Override
    protected void init(Session session) {             // (2)!
    }

    /**
     * Say hello to the given target.
     *
     * @param target
     */
    @Function                                          // (3)!
    void sayHello(String target) {
        println "Hello, ${target}!"
    }

}
```

1. The class your extension builds on. Nextflow requires this to recognize your functions.
2. Called when the plugin loads; use for initialization
3. Makes this method callable from workflows via `include`

The template includes a sample `sayHello` function.
The `@Function` annotation is what makes a method callable from Nextflow workflows.
Without it, the method exists only inside the plugin code.

In Groovy (and Java), methods declare what type they return and what types their parameters are.
For example, `String reverseGreeting(String greeting)` declares a method that takes a `String` parameter and returns a `String`.
The keyword `void` means the method returns nothing, as with `sayHello` above.
This is different from Python or R, where types do not need to be declared explicitly.

---

## 2. Replace sayHello with reverseGreeting

The template's `sayHello` function is a placeholder.
Replace it with your own function to see the full cycle of writing, building, and using a plugin function.

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
        @Function                                // (1)!
        String reverseGreeting(String greeting) { // (2)!
            return greeting.reverse()             // (3)!
        }

    }
    ```

    1. Makes the method callable from Nextflow workflows
    2. Takes a String, returns a String
    3. Groovy's built-in string reversal method

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

Build and install the plugin:

```bash
make install
```

!!! tip "If the build fails"

    Read the error message carefully; it usually includes a line number and describes the problem.
    Common causes are syntax errors (missing bracket or quote), misspelled class names, and type mismatches.
    If you are stuck, compare your code character-by-character with the examples.

---

## 4. Use your function in a workflow

The plugin is built and installed.
The next step is to use `reverseGreeting` in a workflow to verify it works end-to-end.

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

---

## 5. Add decorateGreeting

A plugin can provide multiple functions.
Add a second one that wraps a greeting with decorative markers; you'll make it configurable in Part 6.

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
            return "*** ${greeting} ***"             // (1)!
        }

    }
    ```

    1. Groovy string interpolation: `#!groovy ${...}` inserts the variable's value into the string

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
cd nf-greeting && make install && cd ..
```

Update `greet.nf` to also import and use `decorateGreeting`:

=== "After"

    ```groovy title="greet.nf" hl_lines="4-6 14 16-17 19 33" linenums="1"
    #!/usr/bin/env nextflow

    include { samplesheetToList } from 'plugin/nf-schema'
    // Import custom functions from our plugin
    include { reverseGreeting } from 'plugin/nf-greeting'
    include { decorateGreeting } from 'plugin/nf-greeting'  // (1)!

    params.input = 'greetings.csv'

    process SAY_HELLO {
        input:
            val greeting
        output:
            path 'greeting.txt'
        script:
        // Use our custom plugin function to decorate the greeting
        def decorated = decorateGreeting(greeting)  // (2)!
        """
        echo '$decorated' > greeting.txt
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
        SAY_HELLO.out.view { file -> "Decorated: ${file.text.trim()}" }
    }
    ```

    1. Multiple functions from the same plugin need separate `include` statements
    2. Plugin functions work inside process `script:` blocks too

=== "Before"

    ```groovy title="greet.nf" linenums="1" hl_lines="4 12 15 28"
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

---

## Takeaway

You learned that:

- Functions are defined with the `@Function` annotation in `PluginExtensionPoint` subclasses
- Plugin functions imported with `include` work identically whether from your own plugin or an existing one
- Plugin functions work in both process scripts and workflow operations

---

## What's next?

Your functions work, but so far you've only verified that by running the full pipeline and checking the output by eye.
That approach doesn't scale: as you add more functions, you need a faster way to check that each one behaves correctly, especially after making changes.
The next section introduces unit tests, which let you verify individual functions automatically without running a pipeline.

[Continue to Part 4 :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
