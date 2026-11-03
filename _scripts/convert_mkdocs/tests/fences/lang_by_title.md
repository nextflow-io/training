```groovy title="main.nf"
workflow {
}
```

```groovy title="modules/local/collectGreetings.nf (excerpt)"
process collectGreetings {
}
```

```groovy title="main.nf - Combined handlers"
workflow.onError = {
}
```

```groovy title="nextflow.config"
docker.enabled = true
```

```groovy title="molkart/nextflow.config (excerpt)"
profiles {
}
```

```groovy title="GreetingExtension.groovy"
class GreetingExtension extends PluginExtensionPoint {
}
```

```groovy title="build.gradle"
plugins {
    id 'io.nextflow.nextflow-plugin' version '1.0.0'
}
```

```groovy title="tests/main.sayhello.nf.test"
nextflow_process {
}
```

```groovy title="Output"
package training.plugin
```

```groovy title="Output"
#!/usr/bin/env nextflow
params.input = 'greetings.csv'
```
