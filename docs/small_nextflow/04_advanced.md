# Part 4: Advanced Topics

In this final part, we'll explore version control integration, cloud execution, and extension exercises to deepen your Nextflow skills.

---

## Version control

One of Nextflow's most powerful features is its deep integration with version control systems.
This allows you to share workflows, track changes, and ensure reproducibility by pinning to specific versions.

### Create a GitHub repository

First, let's create a new repository on GitHub to store your workflow.

1. Go to [github.com](https://github.com) and log in
2. Click the "+" icon in the top right and select "New repository"
3. Name it something like `cat-classifier` (or any name you prefer)
4. Make it **public** (so Nextflow can access it easily)
5. **Don't** initialize with a README, .gitignore, or license
6. Click "Create repository"

GitHub will show you some commands to push an existing repository.
Keep this page open - we'll use those commands in a moment.

### Initialize and push your workflow

Now let's version control your workflow.
From your workshop directory:

```bash
# Initialize a git repository
git init

# Add your workflow files
git add main.nf
git add bin/

# If you have a nextflow.config, add that too
git add nextflow.config

# Create your first commit
git commit -m "Initial commit of cat classifier workflow"

# Connect to your GitHub repository (replace with your username and repo name)
git remote add origin https://github.com/YOUR-USERNAME/cat-classifier.git

# Push to GitHub
git branch -M main
git push -u origin main
```

Your workflow is now on GitHub!
Visit your repository URL to see your code online.

### Running remote workflows

Here's where it gets interesting: **you don't need a local copy of a workflow to run it**.

Nextflow can pull workflows directly from GitHub and run them.
For example, to run the nf-core RNA-seq pipeline:

```bash
nextflow run nf-core/rnaseq --help
```

This command pulls the workflow from `github.com/nf-core/rnaseq` (the `nf-core` organization, `rnaseq` repository), downloads it to `$HOME/.nextflow/assets/`, and runs it.

Let's try a simpler example:

```bash
nextflow run hello
```

This pulls and runs `github.com/nextflow-io/hello`.
Notice we didn't specify the full path - Nextflow uses sensible defaults:

- If no provider is specified, it defaults to `github.com`
- If no organization is specified, it defaults to `nextflow-io`

!!! tip "Other Git providers"

    Nextflow also supports:

    - **GitLab**: `nextflow run gitlab.com/user/repo`
    - **Bitbucket**: `nextflow run bitbucket.org/user/repo`
    - **Gitea**: With custom configuration
    - **Azure Repos**: With custom configuration
    - **AWS CodeCommit**: With custom configuration

### Running specific versions with revisions

Now you can run your own workflow from anywhere:

```bash
# Run from GitHub (replace with your username and repo name)
nextflow run YOUR-USERNAME/cat-classifier
```

But what about version control?
What if you want to continue developing while also maintaining a stable version?

Nextflow allows you to specify a **revision** - a specific branch, tag, or commit:

```bash
# Run a specific branch
nextflow run YOUR-USERNAME/cat-classifier -revision dev-branch

# Or use the short form
nextflow run YOUR-USERNAME/cat-classifier -r dev-branch
```

### Using Git tags for stable versions

**Git tags** are named references to specific commits, typically used to mark release versions.
They're like bookmarks in your repository's history - they don't change, making them perfect for reproducible pipelines.

Let's create a `1.0` tag for your workflow:

```bash
# Create an annotated tag
git tag -a 1.0 -m "First stable release of cat classifier"

# Push the tag to GitHub
git push origin 1.0
```

Now you can run this exact version forever:

```bash
nextflow run YOUR-USERNAME/cat-classifier -r 1.0
```

This will always run the code as it existed when you created the tag, even if you continue developing on the `main` branch.

### Testing with different revisions

Let's see this in action.
Create a new branch and make a change:

```bash
# Create and switch to a new branch
git checkout -b experimental

# Make a small change (e.g., modify a parameter default in main.nf)
# Then commit it
git add main.nf
git commit -m "Experimental feature"
git push origin experimental
```

Now you can run different versions:

```bash
# Run the stable 1.0 release
nextflow run YOUR-USERNAME/cat-classifier -r 1.0

# Run the main branch
nextflow run YOUR-USERNAME/cat-classifier -r main

# Run the experimental branch
nextflow run YOUR-USERNAME/cat-classifier -r experimental

# Run a specific commit (use any commit hash from git log)
nextflow run YOUR-USERNAME/cat-classifier -r abc123def
```

This is incredibly powerful: you can have a stable, reproducible pipeline (using a tag) while actively developing new features (on branches), all from the same repository.

### Takeaway

Nextflow's git integration allows you to version your workflows, share them easily, and run specific commits, branches, or tags from anywhere - no local copy required.

### What's next?

Let's explore running workflows on cloud infrastructure.

---

## Cloud executors

Nextflow supports running workflows on cloud infrastructure through various executors including AWS Batch, Google Cloud Batch, and Azure Batch.

### Benefits of cloud execution

Running workflows in the cloud offers several advantages:

- **Scalability**: Automatically scale to hundreds or thousands of parallel tasks
- **Cost efficiency**: Pay only for the compute resources you use
- **No local infrastructure**: No need to maintain local HPC clusters
- **Global accessibility**: Run workflows from anywhere with internet access

### Cloud executor configuration

Each cloud provider has its own executor configuration.
Here's a basic example for AWS Batch:

```groovy title="nextflow.config for AWS Batch" linenums="1"
process {
    executor = 'awsbatch'
    queue = 'my-batch-queue'
    container = 'my-default-container'
}

aws {
    region = 'us-east-1'
    batch {
        cliPath = '/home/ec2-user/miniconda/bin/aws'
    }
}
```

### Data staging with cloud storage

When running in the cloud, you'll typically stage your data in cloud storage:

```bash
# Upload model to S3
aws s3 cp data/models/b32_400m.pt s3://my-bucket/models/

# Run workflow with cloud-based data
nextflow run main.nf --model s3://my-bucket/models/b32_400m.pt
```

!!! tip

    For detailed cloud setup instructions, see the Nextflow documentation for [AWS](https://www.nextflow.io/docs/latest/aws.html), [Google Cloud](https://www.nextflow.io/docs/latest/google.html), or [Azure](https://www.nextflow.io/docs/latest/azure.html).

### Takeaway

Nextflow's cloud executors enable scalable, cost-effective workflow execution without managing local infrastructure.

### What's next?

Try the extension exercises to practice your new skills!

---

## Extension Exercise 1: Find extreme scores

Our team is interested in which cat is the cutest cat and which cat is the ugliest cat.
Can you extend the workflow to identify (for each label) which picture scores the highest?

### Hints

**Hint 1:** You can use the `--json` flag on the `classify.py` script to output structured data instead of plain text.

**Hint 2:** You can parse a JSON file in a closure by using the JsonSlurper class, part of the standard library:

```groovy
| map { meta, jsonFile -> new groovy.json.JsonSlurper().parseText(jsonFile.text) }
```

**Hint 3:** You can use the `min` and `max` operators to return a channel containing the minimum or maximum item, and you can pass a closure to those operators to describe how the elements in the channel should be compared ([docs](https://www.nextflow.io/docs/latest/reference/operator.html#min)).

### Exercise goals

- Modify the Classify process to output JSON
- Parse the JSON to extract scores
- Use operators to find the highest and lowest scoring images
- Publish these special images separately

### Takeaway

This exercise demonstrates how to work with structured data and use comparison operators to filter channel contents.

### What's next?

Try making the classification labels configurable!

---

## Extension Exercise 2: Configurable labels

We've decided that "bad" and "good" are too cruel a classification system for the cats.
Can you modify the workflow to add a `--labels` parameter?
The parameter should take a comma-separated list of labels and use those labels in preference to the default "good cat" and "bad cat".

### Example usage

```bash
nextflow run main.nf --labels 'red cat','orange cat','black cat'
```

### Hints

**Hint 1:** You'll need to modify how the labels are passed to the `classify.py` script.

**Hint 2:** The classify.py script accepts multiple labels via the `--labels` argument:

```bash
classify.py --labels "label1" "label2" "label3" image.jpg
```

**Hint 3:** You can split a comma-separated string in Groovy:

```groovy
params.labels = 'good cat,bad cat'
labels = params.labels.split(',')
```

### Exercise goals

- Add a `--labels` parameter to the workflow
- Parse the comma-separated labels
- Pass them to the classification script
- Ensure the grouping and collage steps still work with custom labels

### Takeaway

This exercise shows how to make workflows more flexible by parameterizing key values.

### What's next?

Congratulations! You've completed the Small Nextflow workshop.

---

## Final notes

### Additional topics to explore

TODO: explain difference between `path(img)` and `path("inputs/*.png")`

TODO: Add in resources directive memory, cpus, etc. (note: this is partially covered in section 8)

### What you've learned

Congratulations on completing the Small Nextflow workshop!
You've built a complete image classification workflow from scratch and learned:

- **Fundamentals**: Channels, processes, parameters, and metadata
- **Data transformation**: Operators, resource management, and grouping
- **Publishing**: Organized outputs with indexes and custom paths
- **Portability**: Filesystem independence and containerization
- **Advanced patterns**: Version control and cloud execution

### Where to go from here

Now that you understand the basics, you can:

- Explore the [Hello Nextflow](../../hello_nextflow/) course for more detailed explanations
- Learn about [nf-core](../../hello_nf-core/) for production-ready pipeline templates
- Try the [Side Quests](../../side_quests/) for advanced topics
- Build your own scientific workflows!

### Takeaway

You now have the foundational knowledge to build, run, and share reproducible scientific workflows with Nextflow.
