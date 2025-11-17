# Part 4: Advanced Topics

In this final part, we'll explore version control integration, cloud execution, and extension exercises to deepen your Nextflow skills.

---

## 1. Version control

TODO: Create git repository at project root

TODO: Commit current state, create git tag, and create branch, and then change and re-commit.

TODO: Change directories and then run using revision argument, pointing to branch, tag, and then specific commit.

### Takeaway

Nextflow's git integration allows you to version your workflows and run specific commits, branches, or tags from anywhere.

### What's next?

Let's explore running workflows on cloud infrastructure.

---

## 2. Cloud executors

Nextflow supports running workflows on cloud infrastructure through various executors including AWS Batch, Google Cloud Batch, and Azure Batch.

### 15.1. Benefits of cloud execution

Running workflows in the cloud offers several advantages:

- **Scalability**: Automatically scale to hundreds or thousands of parallel tasks
- **Cost efficiency**: Pay only for the compute resources you use
- **No local infrastructure**: No need to maintain local HPC clusters
- **Global accessibility**: Run workflows from anywhere with internet access

### 15.2. Cloud executor configuration

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

### 15.3. Data staging with cloud storage

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

## 3. Extension Exercise 1: Find extreme scores

Our team is interested in which cat is the cutest cat and which cat is the ugliest cat.
Can you extend the workflow to identify (for each label) which picture scores the highest?

### 16.1. Hints

**Hint 1:** You can use the `--json` flag on the `classify.py` script to output structured data instead of plain text.

**Hint 2:** You can parse a JSON file in a closure by using the JsonSlurper class, part of the standard library:

```groovy
| map { meta, jsonFile -> new groovy.json.JsonSlurper().parseText(jsonFile.text) }
```

**Hint 3:** You can use the `min` and `max` operators to return a channel containing the minimum or maximum item, and you can pass a closure to those operators to describe how the elements in the channel should be compared ([docs](https://www.nextflow.io/docs/latest/reference/operator.html#min)).

### 16.2. Exercise goals

- Modify the Classify process to output JSON
- Parse the JSON to extract scores
- Use operators to find the highest and lowest scoring images
- Publish these special images separately

### Takeaway

This exercise demonstrates how to work with structured data and use comparison operators to filter channel contents.

### What's next?

Try making the classification labels configurable!

---

## 4. Extension Exercise 2: Configurable labels

We've decided that "bad" and "good" are too cruel a classification system for the cats.
Can you modify the workflow to add a `--labels` parameter?
The parameter should take a comma-separated list of labels and use those labels in preference to the default "good cat" and "bad cat".

### 17.1. Example usage

```bash
nextflow run main.nf --labels 'red cat','orange cat','black cat'
```

### 17.2. Hints

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

### 17.3. Exercise goals

- Add a `--labels` parameter to the workflow
- Parse the comma-separated labels
- Pass them to the classification script
- Ensure the grouping and collage steps still work with custom labels

### Takeaway

This exercise shows how to make workflows more flexible by parameterizing key values.

### What's next?

Congratulations! You've completed the Small Nextflow workshop.

---

## 5. Final notes

### 18.1. Additional topics to explore

TODO: explain difference between `path(img)` and `path("inputs/*.png")`

TODO: Add in resources directive memory, cpus, etc. (note: this is partially covered in section 8)

### 18.2. What you've learned

Congratulations on completing the Small Nextflow workshop!
You've built a complete image classification workflow from scratch and learned:

- **Fundamentals**: Channels, processes, parameters, and metadata
- **Data transformation**: Operators, resource management, and grouping
- **Publishing**: Organized outputs with indexes and custom paths
- **Portability**: Filesystem independence and containerization
- **Advanced patterns**: Version control and cloud execution

### 18.3. Where to go from here

Now that you understand the basics, you can:

- Explore the [Hello Nextflow](../../hello_nextflow/) course for more detailed explanations
- Learn about [nf-core](../../hello_nf-core/) for production-ready pipeline templates
- Try the [Side Quests](../../side_quests/) for advanced topics
- Build your own scientific workflows!

### Takeaway

You now have the foundational knowledge to build, run, and share reproducible scientific workflows with Nextflow.
