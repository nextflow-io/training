# More complex file formats

## JSON

We can also easily parse the JSON file format using the following groovy schema:

```groovy linenums="1"
import groovy.json.JsonSlurper

def f = file('data/meta/regions.json')
def records = new JsonSlurper().parse(f)


for( def entry : records ) {
  log.info "$entry.patient_id -- $entry.feature"
}
```

!!! warning

    When using an older JSON version, you may need to replace `parse(f)` with `parseText(f.text)`

## YAML

This can also be used as a way to parse YAML files:

```groovy linenums="1"
import org.yaml.snakeyaml.Yaml

def f = file('data/meta/regions.yml')
def records = new Yaml().load(f)


for( def entry : records ) {
  log.info "$entry.patient_id -- $entry.feature"
}
```

## Storage of parsers into modules

The best way to store parser scripts is to keep them in a Nextflow module file.

See the following Nextflow script:

```groovy linenums="1"
include{ parseJsonFile } from './modules/parsers.nf'

process foo {
  input:
  tuple val(meta), path(data_file)

  """
  echo your_command $meta.region_id $data_file
  """
}

workflow {
  Channel.fromPath('data/meta/regions*.json') \
    | flatMap { parseJsonFile(it) } \
    | map { entry -> tuple(entry,"/some/data/${entry.patient_id}.txt") } \
    | foo
}
```

For this script to work, a module file called `parsers.nf` needs to be created and stored in a modules folder in the current directory.

The `parsers.nf` file should contain the `parseJsonFile` function.

Nextflow will use this as a custom function within the workflow scope.

!!! tip

    You will learn more about module files later in section 8.1 of this tutorial.
