<!-- fixture: lang=en rel=side_quests/plugin_development/03_custom_functions.md -->
=== "After"

    ```groovy {3-4} title="GreetingExtension.groovy" showLineNumbers startLineNumber=24
    class GreetingExtension extends PluginExtensionPoint {

        @Function                                // [!annotate] Makes the method callable from Nextflow workflows
        String reverseGreeting(String greeting) { // [!annotate] Takes a String, returns a String
            return greeting.reverse()
        }
    }
    ```

=== "Before"

    ```groovy title="GreetingExtension.groovy" showLineNumbers startLineNumber=24
    class GreetingExtension extends PluginExtensionPoint {
    }
    ```
