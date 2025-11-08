# Hello Nextflow Quiz

Great stuff, you've completed all of the training material in Hello Nextflow!
But, did it stick?

Test yourself with the quiz below to see how many Nextflow fundamentals you can remember :think:

## Part 1: Hello World

<?quiz?>

What is the minimum requirement for a Nextflow process?

- [ ] Input and output blocks
- [x] Output and script blocks
- [ ] Only a script block
- [ ] Input block
<?/quiz?>

<?quiz?>

What does the output block in a process do?

- [ ] Creates the output files automatically
- [x] Declares what output files to expect from the process
- [ ] Determines the format of the output
- [ ] Publishes outputs to a directory
<?/quiz?>

<?quiz?>

What command is used to run a Nextflow workflow?

- [ ] `nextflow execute script.nf`
- [x] `nextflow run script.nf`
- [ ] `nextflow start script.nf`
- [ ] `nextflow launch script.nf`
<?/quiz?>

<?quiz?>

What is stored in the `work/` directory?

- [ ] Only the final output files
- [ ] Configuration files
- [x] Temporary files for each process execution
- [x] Execution metadata like `.command.sh` and `.exitcode`
- [x] Output files from processes
- [ ] Container images
<?/quiz?>

<?quiz?>

What does the `-resume` flag do?

- [ ] Restarts the entire workflow from scratch
- [x] Skips processes that already ran with the same code, settings, and inputs
- [ ] Resumes only failed processes
- [ ] Saves the workflow state for later
<?/quiz?>

<?quiz?>

What is the default mode for the `publishDir` directive?

- [ ] copy
- [x] symlink
- [ ] move
- [ ] hardlink
<?/quiz?>

<?quiz?>

How do you pass a parameter to a Nextflow workflow from the command line?

- [ ] `-parameter value`
- [x] `--parameter value`
- [ ] `params.parameter = value`
- [ ] `-p parameter=value`
<?/quiz?>

<?quiz?>

Which quote type is required for variable interpolation in Nextflow?

- [x] Double quotes `"`
- [ ] Single quotes `'`
- [ ] Backticks `` ` ``
- [ ] Triple quotes `"""`
<?/quiz?>

## Part 2: Hello Channels

<?quiz?>

What is a channel in Nextflow?

- [ ] A file directory
- [x] A queue-based mechanism for passing data between processes
- [ ] A configuration parameter
- [ ] A container registry
<?/quiz?>

<?quiz?>

Which channel factory creates a simple queue channel with values?

- [ ] `channel.fromPath()`
- [x] `channel.of()`
- [ ] `channel.from()`
- [ ] `channel.create()`
<?/quiz?>

<?quiz?>

What happens when you pass a channel with multiple values to a process?

- [ ] Only the first value is processed
- [x] The process is called once for each value in the channel
- [ ] All values are processed together in a single call
- [ ] The workflow fails with an error
<?/quiz?>

<?quiz?>

What does the `flatten()` operator do?

- [x] Unpacks array contents into individual items
- [ ] Removes duplicate items from a channel
- [ ] Sorts channel contents alphabetically
- [ ] Converts items to lowercase
<?/quiz?>

<?quiz?>

What is the purpose of the `view()` operator?

- [ ] To visualize the workflow structure
- [x] To inspect and print channel contents for debugging
- [ ] To create graphical output
- [ ] To validate channel data types
<?/quiz?>

<?quiz?>

What does the `splitCsv()` operator do?

- [ ] Creates a CSV file from channel data
- [x] Parses CSV-formatted content, reading each line into an array
- [ ] Splits a channel into multiple channels
- [ ] Validates CSV formatting
<?/quiz?>

<?quiz?>

Which operator transforms channel contents element-by-element?

- [ ] `transform()`
- [ ] `apply()`
- [x] `map()`
- [ ] `modify()`
<?/quiz?>

<?quiz?>

Why do you need dynamic filenames when processing multiple inputs?

- [ ] To save disk space
- [ ] To improve performance
- [x] To prevent output files from overwriting each other
- [ ] To make resume work correctly
<?/quiz?>

## Part 3: Hello Workflow

<?quiz?>

How do you access the output of a process in Nextflow?

- [ ] `processName.output`
- [x] `processName.out`
- [ ] `processName.result`
- [ ] `output.processName`
<?/quiz?>

<?quiz?>

What determines the execution order of processes in a workflow?

- [ ] The order they are defined in the file
- [ ] Alphabetical order by process name
- [x] The dependencies between processes (output → input)
- [ ] The order they are called in the workflow block
<?/quiz?>

<?quiz?>

What does the `collect()` operator do?

- [ ] Unpacks arrays into individual items
- [x] Packages many channel items into a single array element
- [ ] Removes null values from a channel
- [ ] Collects error messages
<?/quiz?>

<?quiz?>

When should you use the `collect()` operator?

- [ ] Before processing items individually in parallel
- [x] Before a process that needs all results at once
- [ ] To speed up channel operations
- [ ] When publishing outputs
<?/quiz?>

<?quiz?>

How do you access named outputs from a process?

- [ ] `processName.out[outputname]`
- [x] `processName.out.outputname`
- [ ] `processName.outputname`
- [ ] `output.processName.outputname`
<?/quiz?>

<?quiz?>

What is the syntax for naming a process output?

- [x] `output: path "file.txt", emit: myoutput`
- [ ] `output: path "file.txt", name: myoutput`
- [ ] `output: path "file.txt" as myoutput`
- [ ] `output: path "file.txt" -> myoutput`
<?/quiz?>

<?quiz?>

When passing multiple inputs to a process, what must match?

- [x] The order of inputs in the workflow call must match the input block order
- [ ] The data types must all be the same
- [ ] The variable names must match exactly
- [ ] The number of inputs must equal the number of outputs
<?/quiz?>

## Part 4: Hello Modules

<?quiz?>

What is a Nextflow module?

- [ ] A configuration file
- [ ] A workflow definition
- [x] A standalone file containing a single process definition
- [ ] A container image
<?/quiz?>

<?quiz?>

What is the recommended filename convention for modules?

- [ ] `module_processName.nf`
- [x] `processName.nf`
- [ ] `processName.module`
- [ ] `processName_module.nf`
<?/quiz?>

<?quiz?>

Where should module files typically be stored?

- [ ] In the root directory
- [ ] In the `lib/` directory
- [x] In the `modules/` directory
- [ ] In the `processes/` directory
<?/quiz?>

<?quiz?>

What is the correct syntax to import a module?

- [ ] `import { processName } from './modules/processName.nf'`
- [x] `include { processName } from './modules/processName.nf'`
- [ ] `require { processName } from './modules/processName.nf'`
- [ ] `use { processName } from './modules/processName.nf'`
<?/quiz?>

<?quiz?>

What happens to resume functionality when you convert processes to modules?

- [ ] Resume stops working
- [ ] Resume only works for the main workflow
- [x] Resume continues to work normally
- [ ] Resume works but cache is reset
<?/quiz?>

<?quiz?>

What are key benefits of using modules?

- [ ] Faster execution time
- [ ] Reduced memory usage
- [x] Reuse processes in multiple workflows without duplication
- [x] Easier maintenance - improve module once, all workflows benefit
- [x] Better code organization
- [ ] Automatic parallelization
<?/quiz?>

## Part 5: Hello Containers

<?quiz?>

What is a container?

- [ ] A directory structure for organizing code
- [x] A lightweight, standalone executable unit with application and dependencies
- [ ] A Nextflow configuration file
- [ ] A type of channel
<?/quiz?>

<?quiz?>

What is the difference between a container image and instance?

- [ ] They are the same thing
- [x] Image is a template; instance is a running container created from the image
- [ ] Image is local; instance is remote
- [ ] Instance is smaller than image
<?/quiz?>

<?quiz?>

What does the `-v` flag do in Docker?

- [ ] Enables verbose output
- [ ] Sets the Docker version
- [x] Mounts a volume to share directories between host and container
- [ ] Validates the container image
<?/quiz?>

<?quiz?>

Why do you need to mount volumes when using containers?

- [ ] To improve performance
- [ ] To reduce disk space usage
- [x] Because containers are isolated from the host filesystem by default
- [ ] To enable networking
<?/quiz?>

<?quiz?>

How do you specify a container for a Nextflow process?

- [ ] `docker 'image_uri'` directive
- [x] `container 'image_uri'` directive
- [ ] `image 'image_uri'` directive
- [ ] `containerImage 'image_uri'` directive
<?/quiz?>

<?quiz?>

How do you enable Docker in a Nextflow workflow?

- [ ] Use the `-with-docker` flag only
- [x] Set `docker.enabled = true` in the config file
- [ ] Docker is enabled by default
- [ ] Use the `--docker` parameter
<?/quiz?>

<?quiz?>

What does Nextflow automatically do when using containers?

- [ ] Builds the container from a Dockerfile
- [x] Pulls the image if not already local
- [x] Mounts the work directory into the container
- [x] Runs the process script inside the container
- [x] Handles cleanup after execution
- [ ] Installs Docker on the system
- [ ] Creates a new container image
<?/quiz?>

## Part 6: Hello Config

<?quiz?>

What file does Nextflow automatically load for configuration?

- [ ] `config.nf`
- [ ] `workflow.config`
- [x] `nextflow.config`
- [ ] `settings.config`
<?/quiz?>

<?quiz?>

What is the configuration precedence in Nextflow?

- [ ] Config file overrides CLI arguments
- [x] CLI arguments override config file values
- [ ] Default values override everything
- [ ] Profile settings always take precedence
<?/quiz?>

<?quiz?>

Can you enable both Docker and Conda in the same configuration?

- [ ] No, only one can be enabled at a time
- [x] Yes, both can be enabled simultaneously
- [ ] Only if using different profiles
- [ ] Only on certain operating systems
<?/quiz?>

<?quiz?>

What happens when both Docker and Conda are enabled for a process that has both directives?

- [ ] Conda is prioritized
- [x] Docker is prioritized
- [ ] The workflow fails with an error
- [ ] The first one defined is used
<?/quiz?>

<?quiz?>

What is the default memory allocation for a Nextflow process?

- [ ] 1.GB
- [x] 2.GB
- [ ] 4.GB
- [ ] It depends on the system
<?/quiz?>

<?quiz?>

How do you set resources for a specific process in the config?

- [x] `process { withName: 'processName' { cpus = 4 } }`
- [ ] `process.processName.cpus = 4`
- [ ] `processName { cpus = 4 }`
- [ ] `resources.processName.cpus = 4`
<?/quiz?>

<?quiz?>

What command generates a resource profiling report?

- [ ] `nextflow run script.nf -profile`
- [ ] `nextflow run script.nf -report`
- [x] `nextflow run script.nf -with-report filename.html`
- [ ] `nextflow run script.nf --resource-report`
<?/quiz?>

<?quiz?>

What is the purpose of `resourceLimits` in the configuration?

- [ ] To set minimum resource requirements
- [x] To cap the maximum resources that can be requested
- [ ] To automatically allocate optimal resources
- [ ] To track resource usage over time
<?/quiz?>

<?quiz?>

What is the default executor in Nextflow?

- [ ] `slurm`
- [x] `local`
- [ ] `cloud`
- [ ] `batch`
<?/quiz?>

<?quiz?>

How do you use a parameter file with Nextflow?

- [ ] `nextflow run script.nf -config file.json`
- [x] `nextflow run script.nf -params-file file.json`
- [ ] `nextflow run script.nf --parameters file.json`
- [ ] Parameter files are loaded automatically
<?/quiz?>

<?quiz?>

What can profiles be used for in Nextflow?

- [ ] User account settings
- [x] Grouping infrastructure settings (executor, container technology)
- [x] Setting different resource limits for different environments
- [x] Providing test parameter values
- [ ] Performance benchmarking
- [ ] Container registry authentication
<?/quiz?>

<?quiz?>

How do you use multiple profiles simultaneously?

- [ ] You can only use one profile at a time
- [ ] `nextflow run script.nf -profile profile1 -profile profile2`
- [x] `nextflow run script.nf -profile profile1,profile2`
- [ ] Profiles are automatically merged
<?/quiz?>

## General Concepts

<?quiz?>

Which of these is part of the Nextflow language?

- [x] Process
- [ ] Rule
- [x] Workflow
- [x] Module
- [ ] Switch
- [x] Channel
- [x] Include
<?/quiz?>

<?quiz?>

What does DSL2 refer to?

- [ ] Docker Specification Language version 2
- [x] Domain-Specific Language version 2
- [ ] Data Structure Library version 2
- [ ] Distributed System Layer version 2
<?/quiz?>

<?quiz?>

What is the shebang required for Nextflow scripts?

- [ ] `#!/bin/bash`
- [x] `#!/usr/bin/env nextflow`
- [ ] `#!/usr/bin/nextflow`
- [ ] No shebang is required
<?/quiz?>

<?quiz?>

What is the naming convention for process names in DSL2?

- [ ] lowercase
- [ ] camelCase
- [x] UPPERCASE
- [ ] snake_case
<?/quiz?>

<?quiz?>

Which tool can generate container images from Conda packages?

- [ ] DockerHub
- [x] Seqera Containers
- [ ] BioContainers
- [ ] Conda Cloud
<?/quiz?>

<?quiz?>

What type of quotes should be used for multi-line script blocks?

- [ ] Single quotes `'...'`
- [ ] Double quotes `"..."`
- [x] Triple quotes `"""..."""`
- [ ] Backticks `` `...` ``
<?/quiz?>

<?quiz?>

What is the file extension for Nextflow scripts?

- [ ] `.nextflow`
- [x] `.nf`
- [ ] `.nx`
- [ ] `.workflow`
<?/quiz?>

<?quiz?>

Which files in the work directory are useful for debugging? (Select all that apply)

- [x] `.command.sh` - shows the actual command executed
- [x] `.command.err` - shows stderr output
- [x] `.command.out` - shows stdout output
- [x] `.command.log` - shows complete log
- [x] `.exitcode` - shows the exit status
- [ ] `.command.config` - shows process configuration
<?/quiz?>

<?quiz?>

What happens if you use single quotes for a string that needs variable interpolation?

- [ ] Nextflow automatically converts them to double quotes
- [x] The variable name appears literally; interpolation doesn't work
- [ ] The workflow fails with a syntax error
- [ ] The variable is interpolated anyway
<?/quiz?>

<?quiz?>

What does the `nextflow clean` command do?

- [ ] Cleans up container images
- [x] Removes work subdirectories from past runs
- [ ] Deletes configuration files
- [ ] Clears the Nextflow cache
<?/quiz?>
