# Part 1: Hello World - Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/8X2hHI-9vms?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important notes"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../01_hello_world.md).

    The section numbers shown in the transcript are provided for indicative purposes only and may not include all section numbers in the materials.

## Welcome

Hi, welcome to Chapter One of Hello Nextflow.

In this first part of a six part course, we're going to go into the very basics of Nextflow. We're going to start off by running some commands in a terminal, and then we'll take those Bash commands and see how to build them into a Nextflow script.

We'll try running that first Nextflow pipeline, see what Nextflow does, where it runs, what files it creates, and what the purpose of where those files is.

All right, let's get started.

## training.nextflow.io

First things first, go to training.nextflow.io. Just as before, all of the material is written here, and I'll be working through it step by step. I'll be showing my screen as I do the steps of a training, but. Everything I'm saying is in the training material so you can follow it at your own speed, and you can find it all written there.

This video also has video subtitles enabled, so feel free to put those up and track exactly what I'm saying as I say it.

Okay, let's go to Hello Nextflow. That's the course we're going to be doing today, and we've done the orientation already in the first video, so we're going to go straight into part one. Hello World.

Okay, I'm going to leave this training material now and hop into my Code Spaces environment. This is what we set up in the first video. Hopefully you have something that looks very similar to this in your own system. I'm using VS Code and I'm looking at the training material and I've changed directories into the hello Nextflow directory.

## 0. Warmup: Run Hello World directly

Okay. Let's start off with a couple of basics, which hopefully will feel familiar to everybody. I'm going to start off just by writing very basic command in the terminal. Down here going to say 'echo Hello World!"' press enter and, no surprises, the terminal does what I ask it and returns that string. Hello world.

Okay, then I'm going to press up to get that command and edit it a bit more. Let's this time redirect that output to a file. I'm going to write it instead to output.txt and press enter nothing on the terminal this time because the output didn't come to the terminal. It went into that file.

I can then read that file by doing 'cat output.txt' hit tab there to auto expand the file name and there you go. The file's there.

I can also see that file in the sidebar over in the file explorer in VS code. I can double click it and open it here. If you want to open it in VS Code without clicking anything, you can also do "code" and then "output.txt" and it does the same thing.

Great. That's the first step. Very simple.

## 1. Examine the Hello World workflow starter script

Okay. We're now going to do exactly the same thing, but in Nextflow, instead of directly in the terminal.

We're going to use the first example script to start with, this file is called Hello World. I can do "ls" to view it in a terminal, and I'm on Mac, so I can do command click to open that file, or I could have just double clicked in the sidebar over here.

There are a few things we can see in this file. Right at the top, there's a hash statement saying that this is a Nextflow file and that's how it could be executed. There are some comments here, just regular code comments in light gray, which don't affect the execution, and just help us read the script.

And then there are two main structures. There's a process here and a workflow.

Processes in Nextflow are the steps of the pipeline. They're the parts which actually do the logic and do the processing.

The workflow then at the bottom stitches these processes together and governs the logic of the workflow, how everything connects to one another.

We're going to start off looking at a process. We'll come back to the workflow in a moment.

## 1.2 The process definition

So every process starts with a key word process. Has a name and then has some curly brackets and everything within those curly brackets is that single process.

A process must have a script section, and contained here is a bash snippet in a multi-line string, which is the part of the code which is actually executed in the compute environment.

We also have an output statement here, which tells Nextflow, which files are expected to be created by the script. Note that the output here has a keyword path, which tells Nextflow that this is a file, not a value, or a string.

Within the script block, this is just a regular bash statement, and it's exactly the same as what we wrote in the terminal. We're echoing hello world to a file called output.txt. This output.txt is then picked up by the output definition. The output definition isn't actually doing anything. It's just telling Nextflow what to expect, and if this file wasn't created, Nextflow would throw an error.

Note that this example is not a great one because we've hardcoded the file name here, output.txt and output.txt. If either of these were changed, that would cause an error in our workflow.

There is better way to do this with variables, which we'll cover in a minute.

## 1.3 The workflow definition

Okay. Moving down to the workflow, we can see that we have a comment and then we run the process called sayHello. This is the same keyword that's up here. This is about as simple as a workflow can get. We're just calling a single process with no variable input, so we're not connecting it to anything else. In the later part of this course, we'll talk about how to make this more powerful by using variable inputs and connecting things with channels.

## 2. Run the workflow

Okay, this is all we need. Let's see if we can run it and see what happens. I am going to just clear the terminal and then I'm going to do "nextflow run", and I'm going to call the file name, which is hello-world.nf. That's all we need to run a Nextflow pipeline. This pipeline doesn't take any input, so we don't need any other arguments.

Let's hit enter and see what happens.

Okay. Hopefully you should have some output, which looks like this. We have a few bits of information telling us that Nextflow ran and what version it was using. Tells us which script was launched and it gives us a randomly generated name for this particular workflow execution. In this case, mine was called "gloomy_crick".

The most important part of this though, is it tells us which steps ran in the pipeline. You can see that our process called sayHello ran, and it ran once and it was a hundred percent complete.

This part here is the hash for that particular workflow task. Each process runs one or more times, and each one of those executions is called a task.

## 2.2. Find the output and logs in the work directory

Every task gets its own isolated directory where it runs, so it's separate from the rest of the execution of the workflow. This hash corresponds to the file structure Within the work directory. If I do "tree work", we can see a0, and then a longer version of a short hash, and then our output.txt file. You can also see it in a sidebar.

You can see in the sidebar there are some additional files here. The reason that these didn't show up in a terminal is because they are hidden files, they start with a dot. And indeed, if I do "tree -a" for all, and "work", we can see them here.

These dot files are present in every single work directory. That Nextflow creates, and each one has a slightly different task. Firstly .command.begin just includes some instructions for Nextflow that sets up the task before it runs. .command.run are the actual instructions executed by Nextflow itself. Then .command.sh is probably the most interesting one. This is the script that was resolved from our process block script.

If I open it, you can see we've got our "echo Hello World" to output.txt file. This is exactly the same as our process in this case, but if we have any variables within our Nextflow code, every task will have a different .command.sh, and you can see how those variables were resolved.

The other files are to do with how the task executed. So .command.err, .log and .out are the standard error, standard output and the two combined. And .exitcode tells Nextflow how this task executed with what exit code, whether it's successful or not.

Finally, we have our output.txt file and sure enough, "Hello World" this is what we're expecting and this is what was created.

Okay, great. That was your first ever Nextflow run. Congratulations. It really is that simple.

Next, we're going to go onto how to do this a little bit more conveniently so that we don't have to edit the code every time we want to make a change to how the pipeline runs.

## 3. Manage workflow executions

This directory structure is great for keeping all the tasks separated and everything organized, but of course, it's not very convenient to find your output files. You don't want to be digging through loads of nested directories trying to find the results of your pipeline.

## 3.1. Publish outputs

The good news is you're not meant to. The work directories are really just for Nextflow to use itself. So what we're going to do is we're going to use a function for Nextflow called "publishDir".

We go back to our workflow, go to the process. We can add a new statement here called a directive. This is what Nextflow calls these things at the top of processes which augment how the functionality works, and the one we're going to use is called publishDir.

You can see I've started typing here and the Nextflow extension for VS Code has suggested to directive for me, so I can just hit enter.

Okay, I am going to follow this with a directory called "results" and we're going to tell it to copy the output files there . So I'm going to say mode copy. Great. going to hit save and let's run the workflow again.

nextflow run hello-world.nf

It runs exactly the same. Though note we have a slightly different hash this time. Nextflow, will use a different hash every time you run the workflow. And we have a different set of work directories as a result. Areas, ones called EB instead, but you can see all the files are the same . However, what's new this time is that we also have a directory called "results".

Within "results" here we have our output file. That's what we told Nextflow to do. We said, save the results files in a directory called "results" and copy them there. And so this is now much easier to find. It's just there alongside where we launched a workflow and all the different files can be organized there however, we wish, irrespective of where or how Nextflow ran the actual execution.

Note that publishDir can handle symlinks, which is good if you're working on a shared file system and you want to save on space. And also you don't have to define all the files which are created by a process as an output.

Nextflow will only copy the things which are defined in this output block. So if you have intermediate files created by the step, which are not needed downstream of this process, you just don't define them in output and they won't turn up in the publishDir. So this is a way of keeping your output files from a pipeline clean and easily deleting intermediate files once the workplace finished.

A quick note here. There's some new Nextflow syntax coming called workflow output definitions, which will eventually replace publishDir. This gives us a way to define all the outputs from a workflow at pipeline level down in the workflow block. This is described in the Nextflow docs if you wanna give it a try. But for now, publishDir will be around for a while, so still have that in a training for 2025.

## 3.2. Re-launch a workflow with -resume

Okay. I mentioned that the work directory here now has two sets of results with a different hash from each time we run the workflow. That's good. However, sometimes we don't want to recompute steps every time if we don't need to.

Maybe you are iteratively building your workflow and you're adding steps in and you want the first steps just to reuse the cached versions. Or maybe something went wrong on your compute system halfway through your workflow and you want it to carry on from where it left off, but skip the steps it had already completed.

Nextflow has built-in functionality for this called resume. Let's try it out. So first off, I'm going to just have a look at the work directory so we can remember what was there.

And then I'm going to do "nextflow run hello-world.nf" and I'm going to add a single command here, "-resume".

Note, single dash, that's really important. I'm going to run it and the output's going to look basically exactly the same, with a couple of small differences.

Note here it says "cached" in gray. That means that Nextflow didn't run the task. This time it found something that matched what were requirements and it reused those outputs directly rather than rerunning the step.

And sure enough, if you look at the hash here, you can see this corresponds to the existing hash that we had from a previous run.

## 3.3. Delete older work directories

Okay. But if you are developing iteratively, you're going to build up a lot of these workflow files. That can be a problem if you might be short on space.

Nextflow can help us clean up these work directories with a couple of helper commands. If I do "nextflow log". That will give me a list of all the different workflow runs that I've done in this directory, and they have the run names here. You can see the gloomy quick one, which was the first one we ran, and then these two new ones.

We can now take that name and use those with the "nextflow clean" command. I can specify a single run name. Or even better, I can tell Nextflow to delete everything from before a single workflow name with "-before", and I'm going to put in "stupefied_shaw". That was my most recent run, "-n".

The "-n" command told Nextflow to do it as a dry run without actually deleting anything for real, and it tells us which of the hash directories it would've been removed. Sure enough, it's just that one from the first execution. Both of the second executions use the same hash directory.

I am going to run it again, but now instead of "-n" for dry run, I'm going to do "-f" for force and it has removed that hash directory. Now if I do "tree work", we can see, we just have this output file left.

Great. So we've managed to clean up a whole bunch of disc space there.

A couple of things to note when deleting work directories, if you symlink stuff to your results directory, those symlink sources will now be deleted and your results will be gone forever. So that's why using copy mode is a safer thing to do, and generally what we recommend.

Secondly, Nextflow's resume functionality relies on these work directories. So if you delete them and you run Nextflow again, the resume functionality will no longer work. So it's up to you to keep track of which things you may need or may not need, and only delete things when you're sure that it's safe to do so.

The other thing we can do is we can just delete the entire work directory if we've finished our workflow run and we're sure we don't need it anymore.

So I can do "rm -r work". I know there was nothing important in there. I've got my results that I care about in the results directory where we copied them. And so it was safe to delete the work directory. It's up to you which of these approaches you use.

## 4. Use a variable input passed on the command line

Okay, what's next? I mentioned that we had hardcoded some of the values in our workflow script here, the output.txt file, and that there might be a better way to do that.

Let's make a start on this. What we're going to do is three things. We're going to add a new input to the process. We're going to tell the process script how to use that input, and then we're going to wire it up in the workflow so that we can use it dynamically with a command line flag when running Nextflow.

So first things first. Let's add an input block here. Just the same as output. This is a new section for the process, and I'm going to say, "val greeting".

Note here, I'm saying "val", which says that this is a variable, not a path.

I can then go down into the script and then I can take out this hardcoded text here and do $greeting. This works just like any other programming language. We're defining a variable here and we're referencing it within this script block. When Nextflow runs this process, the variable will be interpolated. And when we go and look at that .command.sh file, we'll see the actual hard coded string here instead.

## 4.1.3. Set up a CLI parameter and provide it as input to the process call

Okay, but where do we provide the variable? Next we go down to the workflow section, and you can see that the extension here is saying, we now expect an input, and it's given me a warning.

Now, the simplest thing we could do is just hard code it . I could write "Hello World" and provide that string input to the process. But again, that wouldn't really solve any problems. We'd still have to go back and edit the pipeline code every time we wanted to change something, which is no good.

The good news is that Nextflow has a built in system to handle command line arguments called parameters. So instead, I can use one of these special variables called paras and I can call it whatever I want, but I'm going to say greeting so that it matches the workflow logic.

Hit save and let's see what we can do with this.

So if I go back to the terminal. So we do "nextflow run hello-world.nf". Just as before, but the key difference is we do --greeting

Note, there are two dashes here because this is a parameter. When we resumed the workflow before, that was a single dash . That's because resume is a core Nextflow option, and this is a parameter which is specific to our pipeline.

Don't mix the two up. It's easy to do that. If you did --resume instead of just one dash, then that would be "params.resume", which wouldn't do anything. Likewise, if you did a single dash here, Nextflow wouldn't recognize it as a key argument.

So it's --greeting, which corresponds to parameters greeting.

I can now follow that with whatever text I want. So I'm in Sweden at the moment, so I'm going to say, "Hej världen".

So let's run it, see what happens, moment of truth.

Okay, so you can see that the process ran again, just as before, sayHello with a single execution.

This will have overwritten the file that was in the publishDir "results" directory. And so be careful when you're rerunning the files because things in the published air will be overwritten.

I can now do "code results/output.txt", and sure enough, our output's been updated and now says "Hej världen".

## 4.2. Use default values for command line parameters

Okay, that's great. But the problem now is our workflow relies on us always defining this parameter, and it's nice to have sensible defaults so that things will run in a sensible way for your workflow unless you override the defaults.

So the way we do that is by setting a default value for the parameter in our workflow script.

So if I go back to my hello-world.nf file, I can go into the script just above workflow, type "prams.greeting" and define it like any other variable. So let's put a string here and let's say "Holà mundo!"

Now this parameter has got a default defined, which will be used here , or we can still override it on the command line with --greeting, just as we did before.

So let's check it works. "nextflow run hello-world.nf"

No command-line arguments this time, and check whether it did the right thing.

"code results/output.txt". And there it is. We got our default.

Okay, let's try again, just check I'm not telling you any lies. Let's run it again, but do --greeting, and use the example from a training material, let's say "Konnichiwa!"

Reruns, the workflow, and sure enough, our output file up at the top is just updated with the new value which we provided on the command line.

Great. This is a real central aspect to writing any Nextflow workflow. Defining sensible defaults in your pipeline code, but making it very easy to configure for the end user by having command line arguments on the terminal.

Note that the end user can overwrite the config in multiple different places. You can have a config file in your home directory, which is applied to every single Nextflow run that you do. You can have a config file in a launch directory. You can have a config file in a pipeline directory. All of these different config locations are loaded in a specific order, which is described in the Nextflow docs.

Okay, that's the end of section one. We've had our first ever workflow script in Nextflow with a process and a workflow. We've looked at inputs, outputs, scripts, and publishing, and how to wire up parameters and an input channel into our process.

Congratulations, your first step towards writing Next low code is complete.

Have a little break and I'll see you back in a few minutes for chapter two.

[Next video transcript :octicons-arrow-right-24:](02_hello_channels.md)
