```groovy title="Built-in configuration"
process {
    executor = 'local'
}
```

```groovy
docker.enabled = true
process {
    resourceLimits = [cpus: 2]
}
```

```groovy
outputDir = params.outdir
```

```groovy
plugins {
    id 'nf-schema@2.7.2'
}
```

```groovy title="Code to add"
withName: 'COWPY' {
    ext.args = { "-c ${params.character}" }
}
```

```groovy title="Code to add"
ext.prefix = { "cowpy-${meta.id}" }
```

Process directives and old-style params are script code:

```groovy
container 'ghcr.io/schapirolabor/molkart-local:v0.0.4'
```

```groovy
params.input = 'greetings.csv'
```
