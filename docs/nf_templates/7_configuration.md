# Configuration files

The nf-core pipeline template comes with...

Each pipeline comes with a set of “sensible defaults” for the resource requirements, found in `config/base.config`. These are always loaded and can be overwritten by subsequent layers of configuration.

In addition to this base config, pipelines have configuration “profiles” that can be enabled with the command line flag `-profile`. Multiple profiles can be specified in a comma-separated list (e.g., `-profile test,docker`).

!!! note

    The order of arguments is important! They are loaded in sequence, so later profiles can overwrite earlier profiles.

## Software dependencies

The nf-core template pipeline comes with profiles for the management of software dependencies (e.g., `docker`, `singularity`, & `conda`). For pipelines with modules that have been shipped with containers/images/recipes, these profiles can be used to change the way dependencies are handled when you execute your pipeline.

!!! note

    Conda should only be used as a last resort, i.e., when it’s not possible to run the pipeline with Docker or Singularity.

If `-profile` for managing software dependencies is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. **This is not recommended.**

## Test profiles

The nf-core pipeline template also comes with config files called `test` and `test_full` for running the pipeline with a minimal and full test dataset, respectively. The test profile is located in the `conf/test.config` file and is included as a profile in the `nextflow.config` file in your repository.

The test profile is versatile and can be used to test your pipeline during development, test your pipeline is working on different infrastructures, and is utilized by CI testing.

It is recommended that you maintain a working test profile in your pipeline.

Using the `test` profile with one of the profiles for managing your software dependencies the nf-core template pipeline will execute two processes, `FASTQC` and `MULTIQC`.

!!! question "Exercise"

    From within the `nf-core-mypipeline` folder, execute your pipeline template using the `test` and `docker` profiles using the following command.

    ```bash
    nextflow run main.nf -profile test,docker --outdir results
    ```

## Adding new test profiles

Additional test profiles can be added to the template to test different implementations of your pipeline.

To add a new test profile you can create a new test configuration file in the `conf/` folder and include it in the `nextflow.config` file.

For example, you could create a new test configuration file named `test2.conf` in the `conf/` folder.

```console title="conf/test2.config"
params {
    config_profile_name        = 'Test profile 2'
    config_profile_description = 'A second minimal test dataset to check pipeline function'

    // Limit resources so that this can run on GitHub Actions
    max_cpus   = 2
    max_memory = '6.GB'
    max_time   = '6.h'

    // Input data
    input  = 'https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv'

    // Genome references
    genome = 'R64-1-1'
}
```

!!! note

    The example above is a modified version of the `test.conf` file created by the `nf-core create` command.

The `nextflow.config` would need to be updated to include it as an optional profile.

```console title="nextflow.config"
test      { includeConfig 'conf/test.config'      }
test2     { includeConfig 'conf/test2.config'     }
test_full { includeConfig 'conf/test_full.config' }
```

Finally, you can customize CI pipeline run tests as required.

For example, adding new test profiles runs, with or without different parameters.

```console title=".github/workflows/ci.yml"
- name: Run pipeline with a second test profile
        run: |
          nextflow run ${GITHUB_WORKSPACE} -profile test2,docker --outdir ./results
```

!!! question "Exercise"

    Create a new test profile and add it to your CI workflows
