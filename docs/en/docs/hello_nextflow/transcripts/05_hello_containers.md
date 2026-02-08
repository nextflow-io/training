# Part 5: Hello Containers - Video Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xqr--bKEN9U?si=QinuAnFwFj-Z8CrO&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important notes"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../05_hello_containers.md).

    The section numbers shown in the transcript are provided for indicative purposes only and may not include all section numbers in the materials.

## Welcome and background

Hi, and welcome back to Hello Nextflow. This is part five called Hello Containers. And in this part of the course, we're going to talk all about how to encapsulate the software requirements for a pipeline so that people running the pipeline don't have to think about installing the software.

If you've been working in bioinformatics as long as I have, you might remember what I often call the bad old days, where when you wanted to run someone else's pipeline or replicate their work, you would spend hours or days trying to install all the different software tools that they used, at the same versions, trying to compile them on your machine, and it was a nightmare. It was really difficult.

If you're running on an HPC, you might've used environment modules where the sysadmins tried to install software for you, which was, okay, but still imperfect.

But now we have better ways to do this. Nextflow has support built in for different software container technologies. Docker is the most common one. That's the one we're gonna use today. It works well in Codespaces. It works well on your local computer and it works well in the cloud.

But also Singularity or Apptainer, which are very common on HPC systems and effectively work in exactly the same way. Or Poman, Shifter, there's a bunch of other ones which are all very similar.

The one extra one, which is kind of similar but not quite which Nextflow supports is Conda. And Nextflow can manage Conda environments for you on a per process basis, which is much better than doing your own Conda environments. And again, can ship with a pipeline.

We are gonna start off this chapter by talking a little bit about container technologies and Docker and how they work. And we're gonna do the first half just manually in Docker so that you understand what's going on under the hood and how this works. 'Cause that's really important to understand what Nextflow is doing and how to understand what your workflow is doing when it's being executed.

So. Let's jump over into our Codespaces. Now I've cleaned everything up again, but if we jump into Hello Containers, you should see that all of our scripts and everything are there the same as the end of the modules chapter. So we've got our different modules here, which I created in the modules directory.

They are still there. They need to be there so that it can run. and the workflow and the output are all the same except we've changed the output publishing path to Hello Containers, so that your files end up in that directory.

We can run this now to check it works if you like, or we can crack on with the terminal.

## 1. Use a container 'manually'

We're going to be using Docker to manage our containers, and I can check that it's installed on my Codespaces by doing "docker -v", which shows me the, the version that's installed and everything, and that it's working properly.

Now containers and Docker have two concepts which are really important. One is called image, and one's called container. The image is the snapshot, if you like, of the whole file system that you'll be using, and the container is the running environment. So you create a container using an image.

Once you're in that container, it typically works like a whole operating system. It is cut off from the outside world. It is separated from everything else, and that's a good thing. That's how we get such good reproducibility with Nextflow.

Because for tasks run inside a container, they're not tainted by any config files on your local system. Any other external influences, they run in their own little sandbox. The files are then produced in a very, very reproducible way because you're using the same underlying libraries, all the same dependencies, exactly the same software for every person running on every different computing environment. Which frankly I think is fantastic and amazing that it works. And even, even to this day still kind of blows my mind that this is possible.

## 1.1. Pull the container image

So we are going to try out using some Docker images and Docker, when you run it on your system, has a docker registry on your computer, or in this case, in the code space, which keeps track of all the different images which have been downloaded and used in the past, and the different layers of which they're built up on.

We can see what images we have locally with Docker by doing "docker image ls". And in this case you can see there's a bunch of Docker images here, which all to do with setting up this Codespaces. All to do with dev containers and things. So you don't need to worry about them too much, but as we add in more images and download them, as this course goes on, you can check that list and you'll see that the local registry is keeping track of all these things that we've pulled.

But we're gonna grab a new one by doing "docker pull". And that tells Docker to fetch a new image from the web.

We then put in the URI for that container. Now this could be a docker image that you have built locally and then pushed to the internet. it could be an image that someone else has made. There's many, many, many different ways to make Docker images, but arguably one of the simplest ways is to outsource that, and get someone else to do it for you.

And what we're going to use in this, tutorial is a service from Seqera called Seqera Containers.

Now, Seqera Containers is totally free, and it uses a piece of open source software that we develop called Wave, which was built to manage containers in a complimentary way to Nextflow. And it handles many of the common use cases that we find ourselves dealing with, with Nextflow.

It is very common that the software we need is packaged in Conda, in Bioconda, or conda-forge channels or other more domain specific channels. And Wave and Seqera Containers is really good at building images from that.

So I can go to this web UI and we're going to mess around with the package called "cowpy". So I type in the name of the package I want. It searches, it's found it on the Python package index, so I can use that. Or if I wait a little bit longer, it's searching bioconda and conda-forge. And you can see, I can specify any con channel here. So if you want to find an Nvidia channel or anything else, that should work as well.

And then I can specify whether I want it to build a docker image for me or a singularity image and also what CPU architecture I want to use. So amd64 or arm64.

And once bioconda results are listed, I can now see all the different versions that are available as well. I'm gonna pop that in. And now I could keep searching and get more packages from Conda if I want to and compose this container however I want, but I just want that one. So I'm gonna click Get Container.

Now, someone else has already requested the same container before and it's returned from a registry, so we just get it immediately. But if no one else had ever asked for this software package or this combination of software packages, Wave and Seqera Containers would build it on the fly for us.

We can copy this URL and we can also see the view build details. And this shows us what the service did on the backend. It created a conda environment file. A docker file, and then this is it, running the docker build process. It also ran a scan, a security scan, so you can see any CVEs. And it tells you when this was created.

Wave and Seqera Containers can do much more than this, but this is kind of a simple use case, which is most common. And I should say that these images are hosted for at least five years. So you can build these URLs into your pipelines and know that they're not gonna go away anytime soon.

So I've got my URL for my docker image for cowpy.

I can now do "docker pull" that URL, and it will fetch all the different layers and download this image so it's available for me locally.

## 1.2. Use the container to run cowpy as a one-off command

Okay, now let's try and actually use it. So now I'm now gonna use a "docker run" command instead of "docker pull", and I'm gonna use the flag "--rm", which just tells Docker to shut down this container once it's finished doing what I've asked it to. And then I put in identifier for the container, which is just a URI.

And then at the end, I specify the command that I want Docker to run inside the container generated from this image. I'm just gonna say cowpy, which is the name of the tool that's installed from Conta Forge, which is available inside the image.

I'm gonna hit enter and there you go. We've run cowpy on a system. We have a little cow giving us some information.

Now note cowpy is not installed on my local system. So if I run it just without all the Docker stuff, it says, command not found. So this is pulled an image. It has created a container using Docker, and then it's gone into that container and run this command for us and given us the output back to our terminal. Very, very cool.

## 1.3. Use the container to run cowpy interactively

Okay, we're gonna go one step further now and run this container interactively and have a poke around, so we can see what's happening inside the container.

So if I go back up and I take my run command and I'm gonna get rid of cowpy at the end there, 'cause I actually don't wanna run cowpy. I want to run a Bash terminal.

And then I'm gonna go back here and I'm gonna do "-it", which stands for Interactive and Terminal or TTY, and I'm gonna press enter.

And now you can see the prompt, the bit before I type, has changed. This was the Codespaces prompt where it said the directory, and now it says base and roots and tmp. So I'm now inside the container, and if I do "ls", you'll see that the files I see in this directory are different to the files I have in my workspace.

And in fact, I can't see any of the files from my local codespaces workspace or my local drive inside the Docker container. The docker container runtime, is completely isolated and it can't write or read any files from a host file system outside.

I can, however, see the software which is installed inside the container and run it. So I can run cowpy and we can see a little bit more about how to use cowpy. Here I can do "cowpy 'Hello World'" and that tells, tells it to actually put my quote inside a little speech bubble. And you can also run different types of cows, so it doesn't have to be a cow. You can do a "-c". And I'm in Sweden, so I'm gonna choose a moose. Very nice. Given him some antlers.

And there's a whole bunch of different ones you can play around with, which you can see described in the training docs.

## 1.3.4. Mount data into the container

Okay. It would be nice if we could run cowpy on the files in our file system.

Of course, it's not super useful just to have the container and no access to anything at all. It might be safe and reproducible, but it's not very useful.

So how do we do that? I'm gonna get out of this Docker container by typing exit, and you can see the prompt tells us that we're now back in our regular Codespaces again.

And I'm gonna run the same command again. But this time I'm gonna add some additional flags back here. And the important one is "-v", which stands for mounting a volume, which is like basically a part, part of a disk space.

The "-v" takes two parts: there's like a string and then a colon and a string. And the first part is the local file system, which should be mounted into the container. And then the second part is where that should end up inside the container.

Now I just wanna load my whole local file system here. So "." is the current working directory. So I'm just gonna do "." and then ":", and then we're gonna put that into a new directory inside the container called "my_project". This could really be called anything.

And then I'm gonna run again.

In the working directory where I'm dumped, which is /tmp, the files are not there. But if I do "ls my_project", there we have it: all of the same files that we had locally on our Codespaces are now available inside the container at that path.

This is read and write access so I can create new files in this directory and they will show up on my host file system. This particular directory, then behaves exactly as if I was outside the container so I can now read and write and do things.

## 1.3.5. Use the mounted data

Okay, let's just prove that we can do this. I do "cat /my_project/data/greetings.csv". If you remember this file contents look like this. I can now pipe that to cowpy and the cow will print out the different outputs of that file in its little speech bubble, which is kind of fun.

So you can see, we can now use the software in the container to interact with the files on our host system.

Okay, let's drop back out and we will get on with the rest of the training material.

## 2. Use containers in Nextflow

So that's really cool using containers. Hopefully that makes sense. And you can see the value of these containers and why that's useful for running analysis software.

But how do we do this whole same process inside Nextflow? We don't wanna be running loads of Docker commands ourselves. We want to just let Nextflow handle all of that for us.

So let's work through this. We're gonna add a new process to our pipeline, to run cowpy. Okay, so let's create a new module for our new process. So go into modules, let's call it cowPy.nf, and then I'm gonna copy the code from the training material here.

But you can see the process is very simple. It looks much like the ones we've done so far, we have an input block with a path, which is our input file, and also a value here so that this will be a character, so we could use a moose again if we want to.

And then an output, which is a single file here, a path and then a script. And we're doing the same thing we did interactively inside the container: we're doing "cat" to read the, the input file. We're piping that contents onto cowpy. We're choosing a specific character based on that input, we're writing to an output file called cowpy input file, which is then echoed to the output.

Great. Let's include that. So include \{ cowpy \} from "./modules/cowpy.nf", did I call it cowpy? Yep.

And then let's call our new process down here in the main block of the workflow. So let's run cowpy. And we'll take our new cowpy process and we're gonna say collectGreetings.out.

And then if you remember, there were two outputs for this module. One called outfile and one called report. The VS Code extension is auto- suggesting these for us and we want the .outfile.

You can always hop into this process here. You either hover on it and it should show you the quickly what the outputs were. And we can also command click in it and it will open the module file if you wanna see in more detail.

So here we go. That's the outfile there, and that's the path. So that will now be the input file for our cowpy process. Fantastic.

Now if you remember, a cowpy process has two inputs. We also had the the value channel for the character. So we can add in "params.character" here. I could have hard coded this if I wanted to, but let's make it a CLI option so we can do dash, dash character.

Right. I now need to define the input parameter that we've just called and give it a default. So character, String. And I like the moose, so I'm gonna set it to moose by default.

Right, let's try and run it. So if I do Nextflow run hello containers, we'll see what happens.

I could have used dash resume if I had the old work directories kicking around. And again, these first processes would've been, cached and it would've been a bit faster, but it should be basically the same.

Now we can see straight away it has thrown an error when it got to our new process, it's telling us here that there was an error executing the, the cowpy process and it exited with an exit status 127. This is the command it tried to run. It, it looks right, it looks how we expected. It's taking that output file name, which looks about right, it's running it with a moose character and trying to save it.

But you can see the command error here is saying cowpy command's not found. And that makes sense because we haven't actually told Nextflow to use a container yet. We've just given it the cowpy command. And like I said before, cowpy is not installed on our local system. So when it tried to run it, it failed.

## 2.3.1. Specify a container for cowpy

We need to tell Nextflow there is a container available and it can use it. So how do we do that?

If we pop into our module here, we're gonna add a new declaration up at the top called "container". And we're gonna then set that to a string.

Now, if you remember, over in Seqera Containers, I can copy that URL and I just drop that into quotes here.

Now go back and try and run it again.

Lemme see if it works this time.

Unfortunately, it fails in exactly the same way, even though now we've defined a container for the process to run. So in order to use our docker image, we need to tell Nextflow to enable Docker usage when we run the workflow.

And we're gonna do that by creating a new config file. So I'm going to say touch nextflow.config.

This is a special file name where if it's in the working directory while I launch the pipeline, it will be loaded automatically. So if I go into this Nextflow dot config file, you can see it actually already exists, which I had forgotten. And we have docker.enabled in here already, but it's set to false, which is the default.

So if I change that to equals True instead, docker.enabled. And there's reference docs for all of these config scopes in the Nextflow docs. And also you can see that when I hover over with a VS Code extension, it pulls in the docs specific to this and tells me what it means and how to set it.

So now we've set it to true, and if I run Nextflow again, Nextflow will now know to fetch that docker image for us if we don't already have it locally, and then execute that process with that container environment.

And so we can see that it has run successfully and we have a little tick next to a cowpy. Fantastic. If I go up and look in the results directory, the file isn't there yet. And that's because we still need to, publish this output file just as the same as the, all the others.

So we go to the published block within the workflow, say mycowpy equals cowpy.out.

And then down here in the output block, mycowpy, squiggly brackets path. Oops. Hello containers. Mode, copy.

If I run again now, it should run in exactly the same way. I could have used dash resume and I forget every time. And then I go up and now we have a new file created called cowpy-COLLECTED, and there is my moose saying BONJOUR, HELLO, HOLà Fantastic.

Now of course I could also pass now "--character". What are the different options? I think there's a Turkey? So I can use character Turkey. It's gonna run exactly the same way. I missed another opportunity to use dash resume, and now if we load our file up and now we have a Turkey. Fantastic.

## 2.3.4. Inspect how Nextflow launched the containerized task

Okay. Final little thing. Let's just quickly run this command again, resume this time, and have a quick look in the work directory to see what it is that Nextflow is doing under the hood to make all of this work for us.

This time it's super fast, let's go into this work directory, cd work/. Now if you remember we have a bunch of dot files here and the one that we're interested in this case is the one that I said we almost never need to look at, called .command.run.

If I do code dot command run, it's gonna open it up in the editor. And I can search in this file and if I scroll down I should see Docker run. And you can see Nextflow is doing the docker run command for us, when Docker is enabled in a config. It's got a whole bunch of different, flags and things here, but you can see the "-v" flag that we used ourselves when we were running. And you can see it's mounting the local, workspace directory into the container, so that the container can access our input files and save the outputs. And then at the end, it is also running .command.sh, which is the generated script, which has got the cowpy command in.

And so you can see that Nextflow is taking the workflow logic, which is the stuff that we actually care about, which is specific to our analysis, and it's doing all the clever behind the scenes stuff to make Docker work on our system.

And it's doing that in a really portable way so that an end user of the pipeline can switch out the technology they're using: Docker, Singularity, Apptainer, Conda. That doesn't really matter to the pipeline logic, but Nextflow will handle all the underlying infrastructure needs, so that it runs anywhere.

And that is really the superpower of Nextflow. Is reproducibility and portability. And with Nextflow you can actually share your workflow and other people can run it on their systems and it will just work.

That's a really, really difficult thing to do, and now you know how to do it as well with your workflows.

Okay, that's it for this chapter. if you go down to the end of a course, you'll find a, a quiz again about some containers. Hopefully that all made sense. It's a really cool way to work with analysis. And if you're new to containers, I hope I've convinced you that it's the way to go, and you'll never look back.

But with that, have a little break maybe, and you join me in a couple of minutes to go through the final part six of the Hello Nextflow, which is all about configuration.

Thanks very much.
