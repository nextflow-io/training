# Part 2: Hello Containers

NOTE: THIS IS A PLACEHOLDER FOR MATERIAL THAT IS COMING SOON

In Part 1, you learned how to use the basic building blocks of Nextflow to assemble a simple pipeline capable of processing some text and parallelizing execution if there were multiple inputs. 

However, you were limited to utilizing only basic UNIX tools that were already installed in your computing environment. Real-world work typically requires you to use all sorts of tools and packages that don't come standard. Usually you'd have to figure out how to install all of the necessary tools and their software dependencies, as well as manage any conflicts that may arise between dependencies that aren't compatible with each other.

That is all very tedious and annoying, so we're going to show you how to use **containers** to solve this problem much more conveniently.

TODO [shortest possible summary of what are containers]

---

## 1. Use a container directly [basics]

TODO

### 1.1. Pull the container

TODO

```bash
docker pull <container>
```

### 1.2. Spin up the container interactively

TODO

```bash
docker run -it -v ./data:/data <container>
```

### 1.3. Run the command

TODO

```bash
<example command>>
```

This should complete immediately, and you should now see a file called `<file>` in your working directory.

#### 1.4. Exit the container

```bash
exit
```

### Takeaway

You know how to pull a container and run it interactively, and you know how to use that to try out commands without having to install any software on your system.

### What's next?

Learn how to use containers from within a workflow.

---

## 2. Use the container in a workflow

TODO

### 2.1. [step 1]

TODO

### 2.2. [step 2]

TODO

### 2.3. [step 3]

TODO

### Takeaway

You know how to use containers from within a workflow.

### What's next?

Celebrate, take a stretch break and drink some water!

When you are ready, move on to Part 3 of this training series to learn how to apply what you've learned so far to a more realistic data analysis use case.
