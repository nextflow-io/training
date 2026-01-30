# Part 5: Hello Containers - Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/5PyOWjKnNmg?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important notes"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../05_hello_containers.md).

    The section numbers shown in the transcript are provided for indicative purposes only and may not include all section numbers in the materials.

## Welcome

Hi, welcome to Part Five of the Hello Nextflow training course.

This chapter's called Hello Containers. We're going to talk about how Nextflow integrates with tools such as Docker and Singularity to use software containers to provision software to the users of your pipeline.

This means that when people run your pipeline, they don't have to go and install all the different tools themselves. Nextflow will do it for them.

Containers are an extremely powerful technology and crucial in reproducibility and ease of use. We're going to start off by doing a brief introduction to containers themselves, running some docker commands manually, and then we'll take those same containers and put them into our Nextflow pipeline.

Okay. Let's get started.

So just as before, let's start off by loading up the training material. Go to training.nextflow.io. Hello Nextflow, Chapter Five, Hello Containers.

I'm going to jump into my Codespaces environment and on the left here we see hello containers dot nf.

Just as before, this is the same script that we finished previous chapter four on, so it should look familiar.

We've got our command line parameters to specify the input file and the batch name. We are including our three modules, and we have our workflow where we run the three processes.

## 0. Warmup: Run hello-containers.nf

Feel free to run this workflow again and double check that it is producing the outputs that you expect. For now, I'm actually going to close it and dive into the terminal.

## 1. Use a container 'manually'

To start off this chapter, we're going to do a bit of a recap over container technology. If you're very used to docker or singularity or other container techs, then treat this as a refresher , or feel free to skip over it completely.

Nextflow supports many different types of container technologies. That includes Docker, Singularity, Podman, Shifter, Charliecloud, and more.

In this training, we're going to be focusing on Docker. That comes pre-installed in the code spaces and is one of the most popular container technologies, especially if you're developing on your own computer or your own laptop.

If you're working in an academic environment on a shared HPC, you might find that Singularity is available and not Docker. That's fine. All of the concepts are exactly the same. A few of the manual commands are different, but if you understand Docker, you'll also understand singularity.

In fact, Singularity is also installed in the Code Spaces environment. So if you like, you can try and do the same tasks using Singularity instead of Docker.

Okay, so what is container technology? The idea behind Docker is that it can fetch an image from a remote source. Pull it to your local machine and then create a container based on that image.

This running container is a little bit like a virtual machine running on your computer. It's isolated from your environment, and it comes prepackaged with an operating system and a set of available software.

## 1.1. Pull the container image

The syntax that we need to fetch a preexisting image is "docker pull". So I'm going to type that into my terminal, but now we need an image to play with.

You can build images yourself. You can find them on public registries like Docker Hub or quay.io. But a really good way to get images quickly is using Seqera Containers.

This is a free to use community service that we've built in 2024, which you can use without login or anything.

If you go to seqera.io/containers or click containers at the top here, you're presented with a search interface and you can type in the name of any tool available in Conda or on the Python package Index.

By default, it searches the Bioconda and Conda Forge channels, but you can prefix any Conda channel. I'm here if you want to.

For a bit of fun, let's use cowpy. I'm going to type in cowpy. It gives me results from Python Package Index and Conda Forge. I'm going to click that to add it to my container. I could add multiple packages here if I wanted to. Select Docker, select linux/amd64, and click Get Container.

This builds the image for me on demand if it hasn't already been created, and gives me a URL that I can copy.

If you're interested, you can click view Build Details, and that takes you to a page, which shows the conda environment file that was used and the complete build log for the build, along with the security scan results.

If I go back to my code spaces, I can now paste this container name and hit enter.

Docker now downloads all the different layers within this container image, and now tells us that this image is available to use.

## Pulling a Singularity image

If you are using singularity, the process is basically the same. We select our image packages, select cowpy. Now we choose Singularity instead of Docker and click Get Container. That gives us an image URL using oras://. Or if you prefer, you can use https:// by checking that box. Copy that URL. Now go to Code Spaces. We actually have Apptainer installed in this space, which is the same as Singularity, but they're aliased to one another. So I'm going to do apptainer pull and then I'm going to call it cowpy sif, but you can call it whatever you want. Paste URL. And that's going to download that image for me.

I could do ls -lh and see cowpy.sif

Singularity is different to Docker, that singularity stores all images in flat files, whereas Docker has a registry where it keeps all the layers separately on your host machine, and it has a running demon to keep track of all of that.

## 1.2. Use the container to run cowpy as a one-off command

Okay, let's go back to Docker. We can now try running this image that we created by doing docker run.

I'm going to do dash dash rm, which just does a one-off execution of the image. And I'm going to paste the image URL. And then finally, you finish this with a command that you want to run.

The image we generated had cowpy installed, so let's try cowpy.

There you go. It ran our command. I don't have cowpy installed locally. You can see if I try and run it, it doesn't exist. However, in this command, I ran it using Docker and it correctly generated this output.

## 1.3. Use the container to run cowpy interactively

We can go further than this if we like and spin up a container interactively and look around inside. Again, I do "docker run dash dash rm". Now I'm going to do dash it, which tells Docker that we want an interactive terminal. I do the image URL again, and this time, instead of doing cowpy, I'm going to do bin bash because the command that we want to run is bash.

This takes us into this running container and you can see that the prompt has changed now.

If I do LS slash you can see the directories here are different.

If I open a second terminal here on the right hand side, which is just running in GitHub Code Spaces and do LS slash, you see we have directories like workspaces and temp, whereas over here in Docker it's different.

So this environment is completely separate within Docker and isolated from my host environment. That's a good thing, because that isolates the execution of this command into the Docker image and keeps it reproducible between different people on different host systems.

If you want to use data from your host system within the Docker image, you have to explicitly mount that into the container.

We're going to do that in a second.

## 1.3.2. Run the desired tool command(s)

First though, let's see if we can run cowpy. There again, the command is available now directly on the command line, and we can start to do more complex things and pass arguments. Hello containers and instead of the cow, let's do the tux penguin. Let's see what else we have.

Let's do cheese. Wonderful. How about Dragon and Cow? Pretty good.

## 1.3.3. Exit the container

Okay. I can't do much more because I don't have any data in this container. So let's exit out this running image and see if we can mount some data into the container. I can do that by doing control D or typing exit. Okay, I am now back in my regular GitHub code space.

## 1.3.4. Mount data into the container

In order to mount some data into the Docker container, I need to use dash V. So I'm going to take my previous docker command, go back to the start do dash v. I'm going to do "." for the current local working directory, and then a colon to say where that should be mounted in the host directory and do slash data. So that's mounting this particular directory into the container at slash data.

Now if I do LS slash we can see we have a new directory called data, and if I do LS data, you can see all of the files that we have in the sidebar here. Fantastic.

## 1.3.5. Use the mounted data

Now we can start to use some of the files which are on the host system within the Docker image. So I can say cat data greetings csv. If you remember, this is our CSV file with our different greetings from before, and I can pipe that to cowpy. Fantastic. Now we're getting somewhere.

Okay. That's enough for running Docker interactively. Hopefully you now have a feel for roughly what Docker is and how to use it both to run a command in a one-off way, and also to use an image interactively. If you are using singularity. The commands are all very similar except you do things like apptainer exec or apptainer run, or singularity exec or singularity run.

## 2. Use containers in Nextflow

Next we're going to go back to our Nextflow workflow and see how to use this technology within the Nextflow pipeline.

Let's close the terminal and open up Hello Containers again.

## 2.1. Write a cowpy module

To stick with our cowpy example, let's create a new process in our workflow, which uses cowpy. Let's go up to modules, create a new file and call it cowpy nf. I'm now going to cheat a little bit and copy in the code from the training material and hit save. And let's have a look.

So this is a simple process. Hopefully now you understand what the building blocks of a process look like. We've got our publishDir again, going to results. We have two inputs, an input file and a string called character. We have an output cowpy input file, and we have a script which looks exactly the same as what we ran manually inside our docker image a second ago: cat to print a file, piping it onto cowpy, saying which type of cowpy character we want to use, and outputting that to the output file, which we pass as the output here.

## 2.2. Add cowpy to the workflow

Okay, let's go back to our workflow, import this new process. So cowpy from modules cowpy nf. Let's create a new parameter so that we can specify which character we wanted. Let's say Turkey by default. And then let's call this new process at the end of the workflow,

cowpy. And let's use the output here from Collect Greetings. So collect greetings out, out file here. And then we need a second argument, which is the new params that we've just made. params dot character.

## 2.2.4. Run the workflow to verify that it works

Okay, let's see if our new process works. Nextflow run hello containers. This should run those first three processes and then try and run cowpy at the end.

We've got an error. What it's saying here, cowpy had an error and it had an exit status 127 and sure enough, command sh cowpy command not found.

We didn't tell Nextflow that we have a Docker image available for cowpy, so it tried to run it on our host system and we don't have cowpy installed on our host system, so it triggered an error.

## 2.3. Use a container to run it

So what we need to do is we need to tell Nextflow that we have a container available. Let's go to our cowpy process and we're going to add a new directive at the top of the process called container.

We then find our image, copy the URL, and put that in a string.

This isn't enough by itself because an X Flow pipeline can have several ways to specify software. I could also do conda conda-forge cowpy, for example. And Nextflow needs to know which of these technologies you want to use.

## 2.3.2. Enable use of Docker via the nextflow.config file

So in order to run with Docker enabled, we're going to get ahead of ourselves slightly and use the Nextflow config file, which is something we're going to cover in more detail in the next chapter. You can see in this directory that we have a file called Nextflow Config, and in here you already have docker.enabled False.

We're going to change that to True to enable Docker, and then we can try and run the workflow again.

## 2.3.3. Run the workflow with Docker enabled

Nextflow run hello containers nf and this time cowpy ran successfully. Let's look in Results. cowpy collected test and there's our Turkey. Wonderful.

So in the background there, Nextflow knew that it had a container available for that process.

It fetched the image and it ran the commands for us.

## 2.3.4. Inspect how Nextflow launched the containerized task

If you're curious, we can actually see exactly what it did by looking in the work directory. If I do code work, and then the hash and then command run, which if you remember is the actual file that's executed for that task, we can go in and we can look for a function called NXF launch. And here you can see the exact docker command that Nextflow used, which looks a lot like what we were doing manually in the terminal earlier. Docker run. Binding this host directory into the container, and then specifying the container URL.

So there's no magic here. It's just that Nextflow is automatically doing the heavy lifting for you in a way that means you can easily specify containers in your pipeline, which are then readily available for anyone else to use who runs your workflow. And those people no longer have to think about managing software to run your analysis pipeline.

Very, very simple, very convenient, and also really reproducible. Good all around.

Okay, great work. That's the end of Chapter Five. Join us in the next video for part six, which is the final part of this Hello Nextflow training, where we'll talk about Nextflow configuration in more detail.

See you in the next video.

[Next video transcript :octicons-arrow-right-24:](06_hello_config.md)
