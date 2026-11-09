```groovy
channel.of('Hello', 'Bonjour', 'Hola')
    .view()
```

```groovy
def parseSample(filename) {
    return filename.tokenize('_')
}
```

```groovy
include { myFunction } from 'plugin/my-plugin'
```

```groovy
def 'should reverse a greeting'() {
    given:
    def ext = new GreetingExtension()
}
```

```groovy
// Read 'greeting.prefix' from nextflow.config
final prefix = session.config.navigate('greeting.prefix', '***') as String
```

```groovy
@Function
String myFunction(String input) {
    return input
}
```

```groovy
nextflowPlugin {
    provider = 'my-org'
}
```

```groovy title="Syntax example"
input:
    tuple val(meta), path(file)
```
