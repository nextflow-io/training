```groovy title="build.gradle"
nextflowPlugin {
    nextflowVersion = '25.10.0'       // (1)!
    provider = 'training'             // (2)!
}
```

1. **`nextflowVersion`**: Minimum Nextflow version required
2. **`provider`**: Your name or organization

Text after the list.

=== "After"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @Function                                // (1)!
    String decorateGreeting(String greeting) {
        return "*** ${greeting} ***"         // (2)!
    }
    ```

    1. Makes the method callable from Nextflow workflows
    2. Groovy string interpolation: `#!groovy ${...}` inserts the variable's value
       into the string

=== "Before"

    ```groovy title="GreetingExtension.groovy" linenums="24"
    @Function
    ```

```groovy title="greet.nf"
include { decorateGreeting } from 'plugin/nf-greeting'  // (1)!
def decorated = decorateGreeting(greeting)  // (2)!
```

1. Multiple functions from the same plugin need separate `include` statements
