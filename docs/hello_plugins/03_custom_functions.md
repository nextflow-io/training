# Part 3: Custom Functions

In this section, you'll implement custom functions that can be called from Nextflow workflows.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 2 to use as your starting point:

    ```bash
    cp -r solutions/2-create-project/* .
    ```

    Then change into the plugin directory:

    ```bash
    cd nf-greeting
    ```

---

## 1. The PluginExtensionPoint class

Functions are defined in classes that extend `PluginExtensionPoint`.
Open the extension file:

```bash
cat src/main/groovy/training/plugin/GreetingExtension.groovy
```

The template includes a sample `sayHello` function.
We'll replace it with our own functions.
The goal is to create a small library of string manipulation functions: one that reverses text, one that decorates text with markers, and one that formats a friendly greeting.

??? info "Understanding the Groovy syntax"

    If the code looks unfamiliar, here's a breakdown of the key elements:

    **`package training.plugin`**: Declares which package (folder structure) this code belongs to.
    This must match the directory structure.

    **`import ...`**: Brings in code from other packages, similar to Python's `import` or R's `library()`.

    **`@CompileStatic`**: An annotation (marked with `@`) that tells Groovy to check types at compile time.
    This catches errors earlier.

    **`class GreetingExtension extends PluginExtensionPoint`**: Defines a class that inherits from `PluginExtensionPoint`.
    The `extends` keyword means "this class is a type of that class."

    **`@Override`**: Indicates we're replacing a method from the parent class.

    **`@Function`**: The key annotation that makes a method available as a Nextflow function.

    **`String reverseGreeting(String greeting)`**: A method that takes a String parameter and returns a String.
    In Groovy, you can often omit `return`; the last expression is returned automatically.

---

## 2. Add the first function: reverseGreeting

Start by replacing the template's `sayHello` function with something more interesting: a function that reverses a greeting string.
This demonstrates the basic pattern of defining a plugin function.

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

The key parts of this function:

- **`@Function`**: This annotation makes the method callable from Nextflow workflows
- **`String reverseGreeting(String greeting)`**: Takes a String, returns a String
- **`greeting.reverse()`**: Groovy's built-in string reversal method

---

## 3. Add the second function: decorateGreeting

With the basic pattern established, add a second function.
This one wraps a greeting with decorative markers, demonstrating string interpolation.

Add this method after `reverseGreeting`, before the closing brace of the class:

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

---

## 4. Add the third function: friendlyGreeting

The final function demonstrates default parameter values, a feature that makes functions more flexible without requiring callers to provide every argument.

Add this method after `decorateGreeting`:

=== "After"

    ```groovy title="GreetingExtension.groovy" linenums="24" hl_lines="24-30"
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

        /**
         * Decorate a greeting with celebratory markers
         */
        @Function
        String decorateGreeting(String greeting) {
            return "*** ${greeting} ***"
        }

    }
    ```

The `String name = 'World'` syntax provides a default value, just like in Python.
Users can call `friendlyGreeting('Hello')` or `friendlyGreeting('Hello', 'Alice')`.

---

## 5. Understanding the @Function annotation

All three functions share a common pattern: the `@Function` annotation.
This annotation is what makes a method callable from Nextflow workflows.

Key requirements:

- **Methods must be public**: In Groovy, methods are public by default
- **Return type**: Can be any serializable type (`String`, `List`, `Map`, etc.)
- **Parameters**: Can have any number of parameters, including default values

Once defined, functions are available via the `include` statement:

```groovy
include { reverseGreeting; decorateGreeting } from 'plugin/nf-greeting'
```

---

## 6. The init() method

You may have noticed the `init()` method in the extension class.
This method is called when the plugin loads:

```groovy
@Override
void init(Session session) {
    // Access session configuration
    // Initialize resources
    // Set up state
}
```

You can access configuration via `session.config`.
We'll use this in Part 6 to make our plugin configurable.

---

## Takeaway

You learned that:

- Functions are defined with the `@Function` annotation in `PluginExtensionPoint` subclasses
- Methods can have any return type and accept parameters with default values
- Once defined, functions become available to import in Nextflow workflows

---

## What's next?

Now we build and test our plugin.

[Continue to Part 4 :material-arrow-right:](04_build_and_test.md){ .md-button .md-button--primary }
