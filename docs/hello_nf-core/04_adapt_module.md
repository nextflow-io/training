# Part 4: Adapt local modules to nf-core conventions

In this fourth part of the Hello nf-core training course, we show you how to adapt your local modules to follow nf-core conventions.

Now that we've successfully integrated the nf-core `CAT_CAT` module in [Part 3](./03_use_module.md), let's adapt our local `cowpy` module to follow the same nf-core patterns. We'll do this incrementally, introducing one pattern at a time:

1. First, we'll update `cowpy` to accept and propagate metadata tuples
2. Then, we'll simplify its interface using `ext.args`
3. Finally, we'll add configurable output naming with `ext.prefix`

!!! note

    This section assumes you have completed [Part 3: Use an nf-core module](./03_use_module.md) and have integrated the `CAT_CAT` module into your pipeline.

    If you didn't complete Part 3 or want to start fresh for this section, you can use the `core-hello-part3` solution as your starting point:

    ```bash
    cp -r hello-nf-core/solutions/core-hello-part3 core-hello
    cd core-hello
    ```

    This gives you a pipeline with the `CAT_CAT` module already integrated.

---

## 1. Adapt local modules to nf-core conventions

### 1.1. Update cowpy to use metadata tuples

Currently, we're extracting the file from `CAT_CAT`'s output tuple to pass to `cowpy`. It would be better to have `cowpy` accept metadata tuples directly, allowing metadata to flow through the entire workflow.

Open [core-hello/modules/local/cowpy.nf](core-hello/modules/local/cowpy.nf) and modify it to accept metadata tuples:

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="12 16"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

=== "Before"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="12 16"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            path input_file
            val character

        output:
            path "cowpy-${input_file}"

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

Key changes:

1. **Input**: Changed from `path input_file` to `tuple val(meta), path(input_file)` to accept metadata
2. **Output**: Changed to emit a tuple with metadata: `tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output`
3. **Named emit**: Added `emit: cowpy_output` to give the output channel a descriptive name

Now update the workflow to pass the tuple directly instead of extracting the file. Open [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf):

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        cowpy(CAT_CAT.out.file_out, params.character)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2-4"
        // generate ASCII art of the greetings with cowpy
        // Extract the file from the tuple since cowpy doesn't use metadata yet
        ch_for_cowpy = CAT_CAT.out.file_out.map{ meta, file -> file }
        cowpy(ch_for_cowpy, params.character)
    ```

Also update the emit block to use the named emit:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="58" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out.cowpy_output
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="58" hl_lines="2"
        emit:
        cowpy_hellos   = cowpy.out
        versions       = ch_versions                 // channel: [ path(versions.yml) ]
    ```

Test the workflow to ensure metadata flows through correctly:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

The pipeline should run successfully with metadata now flowing from `CAT_CAT` through `cowpy`.

### 1.2. Simplify the interface with ext.args

Now let's address another nf-core pattern: simplifying module interfaces by using `ext.args` for optional command-line arguments.

Currently, our `cowpy` module requires the `character` parameter to be passed as a separate input. While this works, nf-core modules follow a convention of keeping interfaces minimal - only essential inputs (metadata and files) should be declared. Optional tool arguments are instead passed via configuration.

#### Understanding ext.args

The `task.ext.args` pattern is an nf-core convention for passing optional command-line arguments to tools. Instead of adding multiple input parameters for every possible tool option, nf-core modules accept optional arguments through the `ext.args` configuration directive.

Benefits of this approach:

- **Minimal interface**: The module only requires essential inputs (metadata and files)
- **Flexibility**: Users can specify any tool arguments via configuration
- **Consistency**: All nf-core modules follow this pattern
- **Portability**: Modules can be reused in other pipelines without expecting specific parameter names
- **No workflow changes**: Adding new tool options doesn't require updating workflow code

#### Update the module

Let's update the cowpy module to use `ext.args` instead of the `character` input parameter. We'll also remove the local `publishDir` directive to rely on the centralized configuration in `modules.config`.

!!! note "Why remove the local publishDir?"

    nf-core modules should not contain hardcoded `publishDir` directives. Instead, publishing is configured centrally in `conf/modules.config`. This provides several benefits:

    - **Single source of truth**: All output paths are configured in one place
    - **Flexibility**: Users can easily customize where outputs are published
    - **Consistency**: All modules follow the same publishing pattern
    - **No conflicts**: Avoids having two separate publishing locations (local and centralized)

Open [core-hello/modules/local/cowpy.nf](core-hello/modules/local/cowpy.nf):

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="16 18"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

=== "Before"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="6 13"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        publishDir 'results', mode: 'copy'

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)
            val character

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        """
        cat $input_file | cowpy -c "$character" > cowpy-${input_file}
        """
    }
    ```

Key changes:

1. **Removed character input**: The module no longer requires `character` as a separate input
2. **Removed local publishDir**: Deleted the `publishDir 'results', mode: 'copy'` directive to rely on centralized configuration
3. **Added ext.args**: The line `def args = task.ext.args ?: ''` uses the Elvis operator (`?:`) to provide an empty string as default if `task.ext.args` is not set
4. **Updated command**: Changed from hardcoded `-c "$character"` to using the configurable `$args`

The module interface is now simpler - it only accepts the essential metadata and file inputs. By removing the local `publishDir`, we follow the nf-core convention of centralizing all publishing configuration in `modules.config`.

#### Configure ext.args

Now we need to configure the `ext.args` to pass the character option. This allows us to keep the module interface simple while still providing the character option at the pipeline level.

Open [core-hello/conf/modules.config](core-hello/conf/modules.config) and add the cowpy configuration:

=== "After"

    ```groovy title="core-hello/conf/modules.config" linenums="13" hl_lines="6-8"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        ]

        withName: 'cowpy' {
            ext.args = { "-c ${params.character}" }
        }
    }
    ```

=== "Before"

    ```groovy title="core-hello/conf/modules.config" linenums="13"
    process {
        publishDir = [
            path: { "${params.outdir}/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}" },
        ]
    }
    ```

This configuration passes the `params.character` value to cowpy's `-c` flag. Note that we use a closure (`{ "-c ${params.character}" }`) to allow the parameter to be evaluated at runtime.

Key points:

- The **module interface stays simple** - it only accepts the essential metadata and file inputs
- The **pipeline still exposes `params.character`** - users can configure it as before
- The **module is now portable** - it can be reused in other pipelines without expecting a specific parameter name
- Configuration is **centralized** in `modules.config`, keeping workflow logic clean

!!! note

    The `modules.config` file is where nf-core pipelines centralize per-module configuration. This separation of concerns makes modules more reusable across different pipelines.

#### Update the workflow

Since the cowpy module no longer requires the `character` parameter as an input, we need to update the workflow call.

Open [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf) and update the cowpy call:

=== "After"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        cowpy(CAT_CAT.out.file_out)
    ```

=== "Before"

    ```groovy title="core-hello/workflows/hello.nf" linenums="39" hl_lines="2"
        // generate ASCII art of the greetings with cowpy
        cowpy(CAT_CAT.out.file_out, params.character)
    ```

The workflow code is now cleaner - we don't need to pass `params.character` directly to the process. The module interface is kept minimal, making it more portable, while the pipeline still provides the explicit option through configuration.

#### Test

Test that the workflow still works with the ext.args configuration. Let's specify a different character to verify the configuration is working:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false --character cow
```

The pipeline should run successfully. In the output, look for the cowpy process execution line which will show something like:

```console title="Output (excerpt)"
[f3/abc123] process > CORE_HELLO:HELLO:cowpy [100%] 1 of 1 ✔
```

Now let's verify that the `ext.args` configuration actually passed the character argument to the cowpy command. Use the task hash (the `f3/abc123` part) to inspect the `.command.sh` file in the work directory:

```bash
cat work/f3/abc123*/command.sh
```

You should see the cowpy command with the `-c cow` argument:

```console title="Output"
#!/usr/bin/env bash
...
cat test.txt | cowpy -c cow > cowpy-test.txt
```

This confirms that `task.ext.args` successfully passed the character parameter through the configuration rather than requiring it as a process input.

### 1.3. Add configurable output naming with ext.prefix

There's one more nf-core pattern we can apply: using `ext.prefix` for configurable output file naming.

#### Understanding ext.prefix

The `task.ext.prefix` pattern is another nf-core convention for standardizing output file naming across modules while keeping it configurable.

Benefits:

- **Standardized naming**: Output files are typically named using sample IDs from metadata
- **Configurable**: Users can override the default naming if needed
- **Consistent**: All nf-core modules follow this pattern
- **Predictable**: Easy to know what output files will be called

#### Update the module

Let's update the cowpy module to use `ext.prefix` for output file naming.

Open [core-hello/modules/local/cowpy.nf](core-hello/modules/local/cowpy.nf):

=== "After"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 17 19"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("${prefix}.txt"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        prefix = task.ext.prefix ?: "${meta.id}"
        """
        cat $input_file | cowpy $args > ${prefix}.txt
        """
    }
    ```

=== "Before"

    ```groovy title="core-hello/modules/local/cowpy.nf" linenums="1" hl_lines="13 18"
    #!/usr/bin/env nextflow

    // Generate ASCII art with cowpy (https://github.com/jeffbuttars/cowpy)
    process cowpy {

        container 'community.wave.seqera.io/library/cowpy:1.1.5--3db457ae1977a273'
        conda 'conda-forge::cowpy==1.1.5'

        input:
            tuple val(meta), path(input_file)

        output:
            tuple val(meta), path("cowpy-${input_file}"), emit: cowpy_output

        script:
        def args = task.ext.args ?: ''
        """
        cat $input_file | cowpy $args > cowpy-${input_file}
        """
    }
    ```

Key changes:

1. **Added ext.prefix**: `prefix = task.ext.prefix ?: "${meta.id}"` provides a configurable prefix with a sensible default (the sample ID)
2. **Updated output**: Changed from hardcoded `cowpy-${input_file}` to `${prefix}.txt`
3. **Updated command**: Uses the configured prefix for the output filename

Note that the local `publishDir` has already been removed in the previous step, so we're continuing with the centralized configuration approach.

#### Configure ext.prefix

To maintain the same output file naming as before (`cowpy-<id>.txt`), we can configure `ext.prefix` in modules.config.

Update [core-hello/conf/modules.config](core-hello/conf/modules.config):

=== "After"

    ```groovy title="core-hello/conf/modules.config" linenums="21" hl_lines="3"
        withName: 'cowpy' {
            ext.args = { "-c ${params.character}" }
            ext.prefix = { "cowpy-${meta.id}" }
        }
    ```

=== "Before"

    ```groovy title="core-hello/conf/modules.config" linenums="21"
        withName: 'cowpy' {
            ext.args = { "-c ${params.character}" }
        }
    ```

Note that we use a closure (`{ "cowpy-${meta.id}" }`) which has access to `meta` because it's evaluated in the context of the process execution.

!!! note

    The `ext.prefix` closure has access to `meta` because the configuration is evaluated in the context of the process execution, where metadata is available.

#### Test and verify

Test the workflow once more:

```bash
nextflow run . --outdir core-hello-results -profile test,docker --validate_params false
```

Check the outputs:

```bash
ls results/
```

You should see the cowpy output files with the same naming as before: `cowpy-test.txt` (based on the batch name). This demonstrates how `ext.prefix` allows you to maintain your preferred naming convention while keeping the module interface flexible.

If you wanted to change the naming (for example, to just `test.txt`), you would only need to modify the `ext.prefix` configuration - no changes to the module or workflow code would be required.

### Takeaway

You now know how to adapt local modules to follow nf-core conventions:

- Update modules to accept and propagate metadata tuples
- Use `ext.args` to keep module interfaces minimal and portable
- Use `ext.prefix` for configurable, standardized output file naming
- Configure process-specific parameters through `modules.config`

### What's next?

Clean up by optionally removing the now-unused local module.

---

### 1.4. Optional: Clean up unused local modules

Now that we're using the nf-core `cat/cat` module, the local `collectGreetings` module is no longer needed.

Remove or comment out the import line for `collectGreetings` in [core-hello/workflows/hello.nf](core-hello/workflows/hello.nf):

```groovy title="core-hello/workflows/hello.nf" linenums="10"
include { sayHello               } from '../modules/local/sayHello.nf'
include { convertToUpper         } from '../modules/local/convertToUpper.nf'
// include { collectGreetings    } from '../modules/local/collectGreetings.nf'  // No longer needed
include { cowpy                  } from '../modules/local/cowpy.nf'
include { CAT_CAT                } from '../modules/nf-core/cat/cat/main'
```

You can optionally delete the `collectGreetings.nf` file:

```bash
rm modules/local/collectGreetings.nf
```

However, you might want to keep it as a reference for understanding the differences between local and nf-core modules.

---

## 2. Contributing modules back to nf-core

Now that you've learned how to create modules following nf-core conventions, you might develop modules that could benefit the wider community. The [nf-core/modules](https://github.com/nf-core/modules) repository welcomes contributions of well-tested, standardized modules.

### Why contribute?

Contributing your modules to nf-core:

- Makes your tools available to the entire nf-core community through the modules catalog at [nf-co.re/modules](https://nf-co.re/modules)
- Ensures ongoing community maintenance and improvements
- Provides quality assurance through code review and automated testing
- Gives your work visibility and recognition

### Getting started

Before contributing a new module:

1. Check if it already exists at [nf-co.re/modules](https://nf-co.re/modules) or in [open PRs](https://github.com/nf-core/modules/pulls)
2. Create an issue on [nf-core/modules](https://github.com/nf-core/modules/issues) to notify the community of your plans
3. Use `nf-core modules create <tool>/<subtool>` in the nf-core/modules repository to generate the module structure
4. Follow the patterns you've learned: metadata tuples, `ext.args`, `ext.prefix`, and comprehensive testing
5. Submit a pull request and request review from `@nf-core/modules-team`

### Resources

- **Comprehensive guide**: [nf-core components tutorial](https://nf-co.re/docs/tutorials/nf-core_components/components)
- **Module specifications**: [Module guidelines](https://nf-co.re/docs/guidelines/components/modules)
- **Community support**: Join the `#modules` channel on [nf-core Slack](https://nf-co.re/join)

Contributing to nf-core is a rewarding way to give back to the community while ensuring your tools follow best practices and reach researchers worldwide.

---

## Takeaway

You now know how to leverage pre-built nf-core modules in your pipeline and adapt your local modules to follow nf-core conventions. You learned how to use metadata tuples to track sample information through the workflow, simplify module interfaces with `ext.args` for configurable arguments, and use `ext.prefix` for standardized output file naming. These patterns keep modules portable and reusable while centralizing configuration in `modules.config`, making your pipeline more maintainable and consistent with nf-core best practices.

Finally, you learned how to contribute your modules back to the nf-core community, making them available to researchers worldwide and ensuring they benefit from ongoing community maintenance.

## What's next?

Continue to [Part 5: Input validation](./05_input_validation.md) to learn how to add schema-based input validation to your pipeline, or explore other nf-core modules you might add to enhance your pipeline further.
