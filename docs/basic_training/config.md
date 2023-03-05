---
title: Configuration
description: Basic Nextflow Training Workshop
---

# Nextflow configuration

A key Nextflow feature is the ability to decouple the workflow implementation by the configuration setting required by the underlying execution platform.

This enables portable deployment without the need to modify the application code.

## Configuration file

When a pipeline script is launched, Nextflow looks for a file named `nextflow.config` in the current directory and in the script base directory (if it is not the same as the current directory). Finally, it checks for the file: `$HOME/.nextflow/config`.

When more than one of the above files exists, they are merged, so that the settings in the first override the same settings that may appear in the second, and so on.

The default config file search mechanism can be extended by providing an extra configuration file by using the command line option: `-c <config file>`.

### Config syntax

A Nextflow configuration file is a simple text file containing a set of properties defined using the syntax:

```groovy linenums="1"
name = value
```

!!! info

    Please note that string values need to be wrapped in quotation characters while numbers and boolean values (`true`, `false`) do not. Also, note that values are typed, meaning for example that, `1` is different from `'1'`, since the first is interpreted as the number one, while the latter is interpreted as a string value.

### Config variables

Configuration properties can be used as variables in the configuration file itself, by using the usual `$propertyName` or `${expression}` syntax.

```groovy linenums="1"
propertyOne = 'world'
anotherProp = "Hello $propertyOne"
customPath = "$PATH:/my/app/folder"
```

!!! tip

    In the configuration file it’s possible to access any variable defined in the host environment such as `$PATH`, `$HOME`, `$PWD`, etc.

### Config comments

Configuration files use the same conventions for comments used in the Nextflow script:

```groovy linenums="1"
// comment a single line

/*
   a comment spanning
   multiple lines
 */
```

### Config scopes

Configuration settings can be organized in different scopes by dot prefixing the property names with a scope identifier or grouping the properties in the same scope using the curly brackets notation. This is shown in the following example:

```groovy linenums="1"
alpha.x  = 1
alpha.y  = 'string value..'

beta {
    p = 2
    q = 'another string ..'
}
```

### Config params

The scope `params` allows the definition of workflow parameters that override the values defined in the main workflow script.

This is useful to consolidate one or more execution parameters in a separate file.

```groovy linenums="1" title="Config file"
params.foo = 'Bonjour'
params.bar = 'le monde!'
```

```groovy linenums="1" title="Workflow script"
params.foo = 'Hello'
params.bar = 'world!'

// print both params
println "$params.foo $params.bar"
```

!!! exercise

    Save the first snippet as `nextflow.config` and the second one as `params.nf`. Then run:

    ```bash
    nextflow run params.nf
    ```

    ??? result

        ```console
        Bonjour le monde!
        ```

    Execute is again specifying the `foo` parameter on the command line:

    ```bash
    nextflow run params.nf --foo Hola
    ```

    ??? result

        ```console
        Hola le monde!
        ```

    Compare the result of the two executions.

### Config env

The `env` scope allows the definition of one or more variables that will be exported into the environment where the workflow tasks will be executed.

```groovy linenums="1"
env.ALPHA = 'some value'
env.BETA = "$HOME/user/some/path"
```

Save the above snippet as a file named `my-env.config`. Then save the snippet below in a file named `foo.nf`:

```groovy linenums="1"
process foo {
  echo true
  '''
  env | egrep 'ALPHA|BETA'
  '''
}
```

Finally, execute the following command:

```bash
nextflow run foo.nf -c my-env.config
```

??? result

    ```console
    BETA=/home/some/path
    ALPHA=some value
    ```

### Config process

Process [directives](https://www.nextflow.io/docs/latest/process.html#directives) allow the specification of settings for the task execution such as `cpus`, `memory`, `container`, and other resources in the pipeline script.

This is useful when prototyping a small workflow script.

However, it’s always a good practice to decouple the workflow execution logic from the process configuration settings, i.e. it’s strongly suggested to define the process settings in the workflow configuration file instead of the workflow script.

The `process` configuration scope allows the setting of any `process` [directives](https://www.nextflow.io/docs/latest/process.html#directives) in the Nextflow configuration file. For example:

```groovy linenums="1"
process {
    cpus = 10
    memory = 8.GB
    container = 'biocontainers/bamtools:v2.4.0_cv3'
}
```

The above config snippet defines the `cpus`, `memory` and `container` directives for all processes in your workflow script.

The [process selector](https://www.nextflow.io/docs/latest/config.html#process-selectors) can be used to apply the configuration to a specific process or group of processes (discussed later).

!!! info

    Memory and time duration units can be specified either using a string-based notation in which the digit(s) and the unit **can** be separated by a blank or by using the numeric notation in which the digit(s) and the unit are separated by a dot character and are not enclosed by quote characters.

| String syntax     | Numeric syntax | Value                 |
| ----------------- | -------------- | --------------------- |
| `'10 KB'`         | `10.KB`        | 10240 bytes           |
| `'500 MB'`        | `500.MB`       | 524288000 bytes       |
| `'1 min'`         | 1.min          | 60 seconds            |
| `'1 hour 25 sec'` | \-             | 1 hour and 25 seconds |

The syntax for setting `process` directives in the configuration file requires `=` (i.e. assignment operator), whereas it should not be used when setting the process directives within the workflow script.

??? example

    ```groovy linenums="1"
    process foo {
        cpus 4
        memory 2.GB
        time 1.hour
        maxRetries 3

        script:
        """
            your_command --cpus $task.cpus --mem $task.memory
        """
    }
    ```

This is especially important when you want to define a config setting using a dynamic expression using a closure. For example:

```groovy linenums="1"
process foo {
    memory = { 4.GB * task.cpus }
}
```

Directives that require more than one value, e.g. [pod](https://www.nextflow.io/docs/latest/process.html#pod), in the configuration file need to be expressed as a map object.

```groovy linenums="1"
process {
    pod = [env: 'FOO', value: '123']
}
```

Finally, directives that are to be repeated in the process definition, in the configuration files need to be defined as a list object. For example:

```groovy linenums="1"
process {
    pod = [ [env: 'FOO', value: '123'],
            [env: 'BAR', value: '456'] ]
}
```

### Config Docker execution

The container image to be used for the process execution can be specified in the `nextflow.config` file:

```groovy linenums="1"
process.container = 'nextflow/rnaseq-nf'
docker.enabled = true
```

The use of unique "SHA256" docker image IDs guarantees that the image content does not change over time, for example:

```groovy linenums="1"
process.container = 'nextflow/rnaseq-nf@sha256:aeacbd7ea1154f263cda972a96920fb228b2033544c2641476350b9317dab266'
docker.enabled = true
```

### Config Singularity execution

To run a workflow execution with Singularity, a container image file path is required in the Nextflow config file using the container directive:

```groovy linenums="1"
process.container = '/some/singularity/image.sif'
singularity.enabled = true
```

!!! info

    The container image file must be an absolute path: it must start with a `/`.

The following protocols are supported:

-   `library://` download the container image from the [Singularity Library service](https://cloud.sylabs.io/library).
-   `shub://` download the container image from the [Singularity Hub](https://singularity-hub.org/).
-   `docker://` download the container image from the [Docker Hub](https://hub.docker.com/) and convert it to the Singularity format.
-   `docker-daemon://` pull the container image from a local Docker installation and convert it to a Singularity image file.

!!! warning

    Singularity hub `shub://` is no longer available as a builder service. Though existing images from before 19th April 2021 will still work.

!!! tip

    By specifying a plain Docker container image name, Nextflow implicitly downloads and converts it to a Singularity image when the Singularity execution is enabled.

    ```groovy linenums="1"
    process.container = 'nextflow/rnaseq-nf'
    singularity.enabled = true
    ```

    The above configuration instructs Nextflow to use the Singularity engine to run your script processes. The container is pulled from the Docker registry and cached in the current directory to be used for further runs.

    Alternatively, if you have a Singularity image file, its absolute path location can be specified as the container name either using the `-with-singularity` option or the `process.container` setting in the config file.

!!! exercise

    Try to run the script as shown below, changing the `nextflow.config` file to the one above using `singularity`:

    ```bash
    nextflow run script7.nf
    ```

    !!! note

        Nextflow will pull the container image automatically, it will require a few seconds depending on the network connection speed.

### Config Conda execution

The use of a Conda environment can also be provided in the configuration file by adding the following setting in the `nextflow.config` file:

```groovy linenums="1"
process.conda = "/home/ubuntu/miniconda2/envs/nf-tutorial"
```

You can specify the path of an existing Conda environment as either **directory** or the path of Conda environment YAML file.
