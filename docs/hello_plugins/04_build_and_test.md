# Part 4: Testing

Unit tests verify your functions work correctly before you use them in pipelines.
In this section, you'll write and run tests using the Spock testing framework.

!!! tip "Starting from here?"

    If you're joining at this part, copy the solution from Part 3 to use as your starting point:

    ```bash
    cp -r solutions/3-custom-functions/* .
    ```

    Then change into the plugin directory:

    ```bash
    cd nf-greeting
    ```

---

## 1. Why test?

A successful build means the code compiles, but not that it works correctly.
Tests verify your functions behave as expected and help catch bugs when you make changes later.

??? info "What are unit tests?"

    **Unit tests** are small pieces of code that automatically check if your functions work correctly.
    Each test calls a function with known inputs and checks that the output matches what you expect.

    For example, if you have a function that reverses strings, a test might check that `reverse("Hello")` returns `"olleH"`.

    Tests are valuable because:

    - They catch bugs before users do
    - They give you confidence to make changes without breaking things
    - They serve as documentation showing how functions should be used

---

## 2. Understanding Spock tests

The plugin template uses [Spock](https://spockframework.org/), a testing framework for Groovy that reads almost like plain English.
Here's the basic structure:

```groovy
def 'should reverse a greeting'() {   // (1)!
    given:                             // (2)!
    def ext = new GreetingExtension()

    expect:                            // (3)!
    ext.reverseGreeting('Hello') == 'olleH'
}
```

1. **Test name in quotes**: Describes what the test checks. Use plain English.
2. **`given:` block**: Set up what you need for the test (create objects, prepare data)
3. **`expect:` block**: The actual checks. Each line should be `true` for the test to pass

This structure makes tests readable: "Given an extension object, expect that `reverseGreeting('Hello')` equals `'olleH'`."

---

## 3. Create the test file

```bash
touch src/test/groovy/training/plugin/GreetingExtensionTest.groovy
```

Open it in your editor and add the following content:

```groovy title="src/test/groovy/training/plugin/GreetingExtensionTest.groovy" linenums="1"
package training.plugin

import spock.lang.Specification

/**
 * Tests for the greeting extension functions
 */
class GreetingExtensionTest extends Specification {

    def 'should reverse a greeting'() {
        given:
        def ext = new GreetingExtension()

        expect:
        ext.reverseGreeting('Hello') == 'olleH'
        ext.reverseGreeting('Bonjour') == 'ruojnoB'
    }

    def 'should decorate a greeting'() {
        given:
        def ext = new GreetingExtension()

        expect:
        ext.decorateGreeting('Hello') == '*** Hello ***'
    }

    def 'should create friendly greeting with default name'() {
        given:
        def ext = new GreetingExtension()

        expect:
        ext.friendlyGreeting('Hello') == 'Hello, World!'
    }

    def 'should create friendly greeting with custom name'() {
        given:
        def ext = new GreetingExtension()

        expect:
        ext.friendlyGreeting('Hello', 'Alice') == 'Hello, Alice!'
    }
}
```

---

## 4. Run the tests

```bash
make test
```

??? example "Test output"

    ```console
    BUILD SUCCESSFUL in 5s
    6 actionable tasks: 6 executed
    ```

    **Where are the test results?** Gradle hides detailed output when all tests pass.
    "BUILD SUCCESSFUL" means everything worked.
    If any test fails, you'll see detailed error messages.

??? exercise "Add an edge case test"

    Add a test that checks `reverseGreeting` handles an empty string.
    What should `reverseGreeting('')` return?
    Add the test, run `make test`, and verify it passes.

    ??? solution

        Add this test method to `GreetingExtensionTest.groovy`:

        ```groovy
        def 'should handle empty string'() {
            given:
            def ext = new GreetingExtension()

            expect:
            ext.reverseGreeting('') == ''
        }
        ```

        An empty string reversed is still an empty string.

---

## 5. View the test report

Gradle generates an HTML test report with detailed results for each test.

Start a simple web server in the test report directory:

```bash
pushd build/reports/tests/test
python -m http.server
```

VS Code will prompt you to open the application in your browser.
Click through to your test class to see individual test results:

![Test report showing all tests passed](./img/test_report.png)

The report shows each test method, its duration, and whether it passed or failed.

Press ++ctrl+c++ in the terminal to stop the server when you're done, then return to the plugin directory:

```bash
popd
```

---

## Takeaway

You learned that:

- Spock tests use a readable `given:`/`expect:` structure
- Use `make test` to run tests and `build/reports/tests/test/` for the HTML report
- Tests verify behavior and serve as documentation for how functions should be used

---

## What's next?

The next section explores trace observers for hooking into workflow lifecycle events.

[Continue to Part 5 :material-arrow-right:](05_observers.md){ .md-button .md-button--primary }
