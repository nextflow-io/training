# Development Best Practices

Writing Nextflow has never been easier, thanks to improvements in the available tooling and in the Nextflow language itself. Modern development practices, IDE integration, and automated tooling can dramatically improve your productivity, code quality, and debugging experience.

Whether you're a seasoned developer or just starting with Nextflow, adopting these best practices will help you write more maintainable, robust, and efficient workflows. In this side quest, we'll explore the tools and techniques that the Seqera Scientific Development team uses daily to build high-quality Nextflow pipelines.

The benefits of following development best practices include:

- **Faster development**: IDE features like syntax highlighting, auto-completion, and inline documentation help you write code faster
- **Fewer bugs**: Linting catches common errors before runtime, and proper debugging techniques help you identify issues quickly
- **Better collaboration**: Consistent formatting and clear code structure make it easier for team members to understand and contribute
- **Improved maintainability**: Well-structured, properly formatted code is easier to modify and extend over time
- **Professional workflow**: Using modern tooling gives you confidence and makes development more enjoyable

---

## 0. Warmup

Let's move into the project directory and explore the materials we'll be working with:

```bash
cd side-quests/development_best_practices
```

The `development_best_practices` directory contains example workflows and configuration files that we'll use throughout this side quest:

```console title="Directory contents"
development_best_practices/
├── examples/
│   ├── buggy_workflow.nf
│   ├── messy_workflow.nf
│   └── basic_workflow.nf
├── data/
│   └── sample_data.csv
└── nextflow.config
```

These files represent common scenarios you'll encounter when developing Nextflow workflows, including code that needs debugging, formatting, and optimization.

---

## 1. IDE Setup and Configuration

A properly configured Integrated Development Environment (IDE) can transform your Nextflow development experience. We'll focus on Visual Studio Code (VS Code) as it has excellent Nextflow support, but many of these principles apply to other editors as well.

### 1.1. Installing the Nextflow Extension

First, let's ensure you have the Nextflow extension installed in VS Code. If you're working in GitHub Codespaces, this may already be configured.

The Nextflow extension provides:
- **Syntax highlighting**: Makes your code easier to read and understand
- **Auto-completion**: Suggests process names, channel operators, and built-in functions
- **Error detection**: Highlights syntax errors as you type
- **Documentation**: Inline help for Nextflow functions and operators
- **Formatting**: Automatic code formatting and indentation

To install the extension manually:

1. Open VS Code
2. Go to the Extensions view (`Ctrl+Shift+X` or `Cmd+Shift+X`)
3. Search for "Nextflow"
4. Install the official Nextflow extension

TODO: Screenshot/ vid

### 1.2. Using VS Code with Nextflow

Once you have the Nextflow extension installed, you're ready to go! The extension automatically:
- Associates `.nf` files with Nextflow syntax highlighting
- Provides error detection and syntax checking
- Offers auto-completion for Nextflow functions and operators

### 1.3. Using AI Assistance Effectively

Modern AI coding assistants like GitHub Copilot, Codeium, or Cursor can significantly speed up Nextflow development when used properly. Here are some best practices:

**Good AI prompts for Nextflow:**
- "Create a process that runs FastQC on paired-end reads"
- "Write a channel operator to group files by sample name"
- "Generate a config file with resource requirements for different processes"

**Tips for effective AI assistance:**
- Be specific about input/output types
- Mention resource requirements (CPU, memory)
- Specify container or conda environments
- Ask for error handling and edge cases

Let's try this with the basic workflow. Open `examples/basic_workflow.nf` and examine the structure:

```bash
cat examples/basic_workflow.nf
```

**Exercise 1.1**: Use your AI assistant to help expand this basic workflow by adding a quality control process. Ask it to:
1. Add a FastQC process after the input channel
2. Include proper resource requirements
3. Add a container or conda directive

!!! tip

    When using AI assistance, always review the generated code for:
    - Proper channel handling
    - Resource requirements
    - Error handling
    - Best practices compliance

---

## 2. Debugging Techniques

Debugging is an essential skill for any Nextflow developer. Nextflow provides several built-in tools and techniques to help you identify and resolve issues in your workflows.

### 2.1. Common Debugging Scenarios

Let's examine a workflow with intentional bugs:

```bash
cat examples/buggy_workflow.nf
```

This workflow contains several common issues:
- Channel handling errors
- Process input/output mismatches
- Configuration problems
- Resource requirement issues

### 2.2. Using Nextflow's Built-in Debugging Features

**Trace and Timeline Reports**

Run the workflow with tracing enabled:

```bash
nextflow run examples/buggy_workflow.nf -with-trace -with-timeline -with-report
```

This generates:
- `trace.txt`: Detailed execution trace
- `timeline.html`: Visual timeline of process execution
- `report.html`: Execution summary and resource usage

**Dry Run Mode**

Before running a complex workflow, use dry run mode to check the logic:

```bash
nextflow run examples/buggy_workflow.nf -preview
```

This shows you what processes would be executed without actually running them.

### 2.3. Channel Debugging

One of the most common sources of bugs in Nextflow is incorrect channel handling. Use the `.view()` operator to inspect channel contents:

```groovy
input_ch = Channel.fromPath(params.input)
    .view { "Input file: $it" }
    .splitCsv(header: true)
    .view { "CSV row: $it" }
```

### 2.4. Process Debugging

For debugging individual processes, you can:

1. **Add debug statements in the script block:**
```groovy
script:
"""
echo "Processing: ${sample_id}"
echo "Input files: ${input_files}"
# Your actual commands here
"""
```

2. **Use the `echo` directive:**
```groovy
process myProcess {
    echo true

    script:
    """
    echo "Debug info from process"
    """
}
```

3. **Inspect work directories:**
```bash
ls -la work/
```

**Exercise 2.1**: Debug the `buggy_workflow.nf` file by:
1. Running it with trace enabled
2. Adding `.view()` operators to inspect channels
3. Identifying and fixing the bugs
4. Running the corrected version

You can find the corrected version in `solutions/development_best_practices/debugged_workflow.nf`.

### 2.5. Advanced Debugging Techniques

**Resume and Caching**

When debugging, use Nextflow's resume feature to avoid re-running successful processes:

```bash
nextflow run examples/buggy_workflow.nf -resume
```

**Stub Running**

For testing workflow logic without running actual processes:

```groovy
process myProcess {
    stub:
    """
    touch output.txt
    """

    script:
    """
    # Actual process logic
    """
}
```

---

## 3. Workflow Linting and Formatting

Consistent code style and adherence to best practices are crucial for maintainable workflows. Nextflow provides built-in linting capabilities to help you identify and fix common issues.

### 3.1. Using Nextflow Lint

The `nextflow config` command can validate your configuration:

```bash
nextflow config -profile test
```

For workflow linting, use:

```bash
nextflow run examples/messy_workflow.nf -profile test -preview
```

### 3.2. Common Linting Issues

Let's examine a poorly formatted workflow:

```bash
cat examples/messy_workflow.nf
```

Common issues include:
- Inconsistent indentation
- Missing documentation
- Poor variable naming
- Inefficient channel operations
- Missing error handling

### 3.3. Formatting Best Practices

**Consistent indentation** (use 4 spaces):
```groovy
process myProcess {
    input:
        path input_file

    output:
        path "output.txt"

    script:
    """
    command --input ${input_file} --output output.txt
    """
}
```

**Clear variable naming:**
```groovy
// Good
samples_ch = Channel.fromPath(params.samples)
fastqc_results_ch = fastqc(samples_ch)

// Bad
ch1 = Channel.fromPath(params.samples)
ch2 = fastqc(ch1)
```

**Proper documentation:**
```groovy
/*
 * Quality control process using FastQC
 * Input: Raw sequencing reads
 * Output: Quality control reports
 */
process fastqc {
    // Process definition
}
```

### 3.4. Automated Formatting

While Nextflow doesn't have a built-in formatter like `black` for Python or `prettier` for JavaScript, you can use VS Code's formatting features:

1. Install the Nextflow extension
2. Use `Shift+Alt+F` (or `Shift+Option+F` on Mac) to format
3. Configure auto-formatting on save

**Exercise 3.1**: Format the `messy_workflow.nf` file by:
1. Fixing indentation issues
2. Adding proper documentation
3. Improving variable names
4. Organizing imports and configurations

### 3.5. Configuration Best Practices

**Organize your `nextflow.config`:**
```groovy
// Global parameters
params {
    input = null
    output = 'results'
    publish_dir_mode = 'copy'
}

// Process-specific configurations
process {
    withName: 'fastqc' {
        cpus = 2
        memory = '4.GB'
        container = 'biocontainers/fastqc:v0.11.9_cv8'
    }
}

// Profile configurations
profiles {
    test {
        params.input = 'data/test_samples.csv'
    }

    cluster {
        process.executor = 'slurm'
        process.queue = 'compute'
    }
}
```

---

## 4. Advanced Development Practices

### 4.1. Modular Development

Break your workflows into reusable modules:

```groovy
// modules/fastqc.nf
process fastqc {
    publishDir "${params.output}/fastqc", mode: params.publish_dir_mode

    input:
        tuple val(sample_id), path(reads)

    output:
        tuple val(sample_id), path("${sample_id}_fastqc_report.html")

    script:
    """
    fastqc ${reads} --outdir . --noextract
    mv *_fastqc.html ${sample_id}_fastqc_report.html
    """
}
```

### 4.2. Error Handling

Implement proper error handling in your processes:

```groovy
process robustProcess {
    errorStrategy 'retry'
    maxRetries 3

    script:
    """
    # Your command here
    if [ \$? -ne 0 ]; then
        echo "Process failed, will retry"
        exit 1
    fi
    """
}
```

### 4.3. Testing and Validation

Always test your workflows with different inputs:

```bash
# Test with minimal dataset
nextflow run main.nf -profile test

# Test with different parameters
nextflow run main.nf --input data/large_dataset.csv --output results_large

# Test resume functionality
nextflow run main.nf -resume
```

**Exercise 4.1**: Create a well-structured workflow that includes:
1. Modular processes
2. Proper error handling
3. Comprehensive documentation
4. Multiple test profiles

---

## 5. Best Practices Summary

### 5.1. Development Workflow

1. **Start with a plan**: Define your inputs, outputs, and process flow
2. **Use version control**: Track changes with git
3. **Write tests**: Create test profiles and sample data
4. **Document as you go**: Add comments and documentation
5. **Lint and format**: Use consistent style and check for errors
6. **Debug systematically**: Use tracing and channel inspection

### 5.2. Code Quality Checklist

- [ ] Consistent indentation (4 spaces)
- [ ] Clear variable and process names
- [ ] Proper input/output declarations
- [ ] Resource requirements specified
- [ ] Error handling implemented
- [ ] Documentation included
- [ ] Test data provided
- [ ] Configuration organized

### 5.3. Performance Considerations

- Use appropriate resource requirements
- Optimize channel operations
- Consider process parallelization
- Monitor resource usage with reports
- Use caching and resume effectively

---

## Takeaway

You now have a comprehensive toolkit for developing high-quality Nextflow workflows. By incorporating these practices into your development routine, you'll write more maintainable, efficient, and robust pipelines.

Key takeaways:
- **IDE tooling** dramatically improves development speed and reduces errors
- **Systematic debugging** helps you identify and fix issues quickly
- **Consistent formatting** and linting improve code quality and collaboration
- **Modern development practices** make Nextflow development more enjoyable and professional

### What's next?

Continue practicing these techniques in your own workflows. Consider exploring:
- Advanced testing with nf-test
- Continuous integration for Nextflow workflows
- Performance optimization techniques
- Advanced debugging and profiling tools

Remember: good development practices are habits that improve over time. Start with a few techniques and gradually incorporate more as they become natural.

Happy coding!
