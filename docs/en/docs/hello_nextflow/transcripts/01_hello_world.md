# Part 1: Hello World - Video Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/tOukLxWCHiA?si=F0t9LFYLjAWoyRXj&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important notes"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../01_hello_world.md).

    The section numbers shown in the transcript are provided for indicative purposes only and may not include all section numbers in the materials.

## Welcome

Hi, and welcome back.

You are now in Part One of the "Hello Nextflow" course called "Hello World". In this chapter, we're gonna start to build up some understanding of the very basics of Nextflow.

So hopefully you are now set up in Codespaces or somewhere equivalent with VS Code running, and you have your Hello Nextflow folder in the workspace in the Explorer with all these different files here.

We're gonna start off by doing just some very basic things in the terminal using Bash, and then we'll see if we can do the same things within Nextflow so you get a feel for what the syntax looks like.

## 0. Warmup

So let's start really, simple. Let's just start with "echo", to print something to a terminal. "Hello World". I press enter and that goes to a terminal. Hello World. Hopefully that's not a surprise to anyone watching this course.

Okay, let's do something with this. Rather than just printing it to the terminal, let's write it to a file. I'm going to press up my keyboard up cursor, which cycles through the Bash history, so it gives me my last command, and I'm gonna tack onto the end of it there, little greater than symbol, which redirects the output from this command to a file, and I'm gonna call it output.txt.

Enter again, to run that command, nothing in the terminal this time, but we can see on the left hand side, the new file has popped up here, called output.txt.

We can view that in a terminal with something like cat. So cat output.txt and sure enough it says "Hello World". We can also double click it and it opens it up in the code editor in VS Code.

## 1.1. Examine the code

All right. I told you it was simple. What's next? Let's try and take this process and do it again, but this time, let's do it inside Nextflow.

As I said, all of the different chapters in this course start with a script and this one is called, Hello World. So I'm gonna find Hello World. It previews it if I single click it, I'm gonna double click it to open it in the editor here. And I'm gonna just quickly get rid of the terminal.

Now this is a very, simple script, so about as simple as it gets. It's only 22 lines long, and it does basically the same thing. In fact. Some of this should look familiar. It is what we just typed. We can see our bash command redirecting to a file there.

Okay. What else? Also, in this file, we can start to see some of the core concepts of Nextflow. We have a process in red here and a workflow. These the special keywords and special terminology in Nextflow.

## 1.1.1. The process definition

Different processes within a workflow wrap different logical units of your workflow. Each process does one thing.

When we run it, it generates a task or multiple tasks, which are an actual doing steps of a pipeline. All of the processes are then orchestrated within a workflow block, which we see at the bottom, and in this case, just runs that one process.

The process name follows this keyword here, and this can be basically anything. And then the contents of the process are within these curly brackets.

There's only really one requirement for process, which is that it includes some kind of a script or an exec block. This is in the triple quotes here, and this is the bash script which gets written to the working directory when we run the pipeline and there's a thing which actually runs on your computer or server.

This is bash typically, but you can also put in a different hash bang here at the top, and it could be a Python script or a R script. It doesn't matter. Whatever is in this script will be executed.

There's one other thing we've added into this process here, which is the output declaration. This tells Nextflow that this process is expecting an output file called output.txt. It says it's a path, so it should be handled like a file, not say, if this was val, it would say it's like a variable or value.

Note that this is not creating this file. It's not actually generating it. That's done by the script down here. It's just telling Nextflow to expect an output file with this filename.

## 1.1.2. The workflow definition

Okay. And then at the bottom we have a workflow here, and again, we have a declaration. This one's called Main. This is the workflow equivalent of a script block, if you like. It's the part of the workflow that does something. And in this case, we are saying, call the process called sayHello.

Normally, of course, your pipeline will look a lot more complex than this. You'll have probably more than one process, and you'll use channels to orchestrate the data flow between them. We're gonna come onto that in the next parts of this course, but for now, this is enough. This is a valid pipeline, which should work.

I can even click preview DAG here in VS Code. The DAG or DAG is a representation of a data flow structure in the pipeline, and we can see it rendered on the side as a mermaid diagram. In this case it's very, simple. There's one box, which is the workflow and one process, which is called sayHello, but that might look more interesting as we go on.

## 1.2. Run the workflow

Okay, let's try and run this workflow and see what happens.

I'm gonna bring up the terminal again at the bottom, clear the output, and I'm gonna type Nextflow Run. And then I'm just gonna type in the script name, which is hello-world.nf. And I'm gonna press enter.

Okay, it's got some standard stuff at the top, which tells us that Nextflow ran and which version was running and what the script name and everything was.

And really the important thing that we're looking for here is _here_, which is a summary of the different tasks which were executed.

If yours looks like this with a little green tick, then well done. You've just run your first pipeline. Fantastic.

It tells us here the name of the process, which Ran, which was called Say Hello, and it told us that it ran once and that it was successful. This updates as you go along, so when you're running a bigger pipeline, you'll see the progress represented here. But because this is so tiny, it runs basically immediately.

## 1.2.2. Find the output and logs in the work directory

Now when you run a Nextflow pipeline, each one of those processes is stitched together, and each process, like I said before, can generate tasks one or multiple. So in this case, we had a single task from this process. It just ran once and that was done under this task _hash_.

Nextflow doesn't deal with the files in your working directory directly, it creates a special folder called work. And if I do "ls", we'll see it has appeared here: _work_, and within here there are sub directories for every single task which runs. And that matches this hash. So you can see if I go to "ls work/c4", and then it's truncated, but it starts 203, and that's the working directory, which was created by this process when we ran the pipeline. And you can see it on the side as well.

When I list those files, you can see the output.txt file was generated. You can see it here as well. And there's a bunch of hidden files , which are not showing up with my regular "ls".

If I click on output.txt, sure enough, we have our output. Fantastic. So the pipeline worked.

It might seem like quite a lot of boilerplate for running what was essentially a one line bash script, but it will make more sense as our processes get more complicated. And this work directory with Nextflow and these files, which are created is really the backbone of what makes Nextflow so powerful.

Each task, each element of a pipeline is isolated from every other task. It's reproducible. They don't conflict with one another, and everything can run in parallel. It's actually a really nice way when you get used to it because of this isolation that you can go in and see exactly what happened for a single task and debug.

Let's have a quick look at these other files in the work directory. From top to bottom, we have a file called _.command.begin_. This is empty. It is just at what's called a sentinel file, created by Nextflow saying, okay, I'm starting the task. Nothing interesting there.

Then there's _.command.error_, _.command.log_ and _.command.out_. These are all outputs from the bash command or this script which ran. This is standard error. This is standard out, and this is the two of them combined as they came out. So you get the logical order.

Okay, those were all empty for this as well, so not very interesting, but things get more interesting when you get to _.command.run_.

This is typically a very long script. And this is what Nextflow actually executes. If you go in here, you'll start to see all the inner logic of Nextflow and see what it's doing and how it's executing your process. This will depend on where you are running, whether we're running locally or submitting it as a job to SLURM, in which case we'll have SLURM headers at the top. All these different setups.

Generally, you don't really need to ever look in this file, though. It is autogenerated by Nextflow and there's nothing really particularly unique to your pipeline, which is in it. But that's really the core of what's running.

The next one is much more interesting. _.command.sh_ is the generated script, which came from your process, and here you can see that Nextflow added the Bash header, and then it executed our command, which was in our script block.

And that's all the _.command.run_ file does is it just runs this _.command.sh_ file.

This is a really useful one, which is the one you usually end up looking at the most when you're trying to debug something and check that the logic of your Nextflow pipeline is doing what you expect it to do.

Finally, we have a file called _.exitcode_, and this just captures the exit code from a task , which in this case was successful. So the exit code was zero.

If something goes wrong, you run outta memory or something else and it fails, then this very useful to understand what went wrong.

## 1.3. Run the workflow again

One more thing to understand about work directories is that if I keep running this pipeline repeatedly, so if I _"nextflow run hello-world.nf"_, it's gonna do exactly the same thing, but this time it will have a new task id. You can see that this hash here is different, and now if I look in work, there's two hash directories. And these are, again, separate from one another.

So every time you run a Nextflow workflow, unless you use the resume, which uses the cache, we'll touch on later, it's gonna rerun those processes in new work directories, which are separate from one another. You won't get any file name collisions, you won't have any problems like that. Everything is isolated and clean.

And if we go in this directory, you can see all the same files and the same _output.txt_, which has been recreated from scratch.

## 2. Publish outputs

Okay, that's great for Nextflow for itself, whilst it's running your pipeline so that all the things are separate from one another and clean and can be managed.

But it's not super useful if you're a person trying to explore your results. You don't really want to be digging around through thousands and thousands of different work directories trying to find your result files. And you're not really meant to. The work directories are not meant to be the final state of where your files are created.

We do this by publishing our files.

## 2.1.1. Declare the output of the sayHello process

So if I go back to our script, we're gonna work in our workflow block here. We're going to tell it what files to expect, which files we care about, and then we're gonna create a new block underneath called the output block.

This is the new syntax, which came with a the syntax parser and be default in version 26.04 of Nextflow. So if you have used Nextflow a little bit before, this is one of the things which is new.

So we've got the main block, and next I'm gonna say publish and I'm gonna tell Nextflow what to expect from the publishing. We're gonna call it _first_output_, and we're gonna call it, _sayHello.out_.

I accidentally made a typo there, but this is a good opportunity to also point out some of the features of the Nextflow VS Code extension. You can see that right away it gave me a little wiggly red line underneath this saying something's wrong. And if I hover over it, it's gonna tell me this variable is not defined. I dunno what it is.

It's pretty obvious in this case, I made a typo. I meant to type, sayHello, and then the squiggly line goes away.

Now it's purple. The Nextflow syntax parser knows that this is a process and when I hover over it, it gives me a reduced representation of what this process looks like. So I can see very quickly at a glance it doesn't take any inputs and it gives us this output. So working in VS Code with this extension gives you lots of contextual information as you're writing code.

Note that we can refer to the output from this process with the _.out_ syntax. And at the moment we can call this whatever we like, it's just an arbitrary variable name.

## 2.1.2. Add an output: block to the script

Where it becomes important is when we do our new block here, and this is below the workflow block now, we're not no longer inside workflow. Squiggly brackets again. And this is where we just tell Nextflow where to put all of the files, which are created by the workflow.

Now I'm gonna take this variable name, which I created here, and I'm gonna put that there and put some squiggly brackets for this. And I'm gonna tell Nextflow to use a path. Oops. Path, in quote marks. And I'm gonna use dot. That just tells Nextflow to put the file in the root of the results directory. So not any sub directories or anything.

Let's try running our workflow again. If I do _"nextflow run hello-world.nf"_, then hopefully it should look basically exactly the same. Nothing has really changed with Nextflow here. It's running the same things. It's just doing them in work directories again.

But now if I do _"ls results/"_, you'll see there's a new directory here that's been created called results, which is the default base directory for workflow publishing. And in there's a file called _output.txt_.

If I do _"ls -l results"_, you'll see this is actually soft linked to the work directory. So this is not a real file, it's linked to the work directory and it's collected all of the files there for us.

## 2.2. Set a custom location

"Results" is the default name for this path. If I run the workflow again, and this time I do _dash_ single hyphen, this is, 'cause it's a core Nextflow option. _" -Output-dir **my** results"._ Could also just do _"-o"_ for short. Then it's gonna set a different base directory for where the files are stored and once again, up here in _myresults/_, now we have an _output.txt_.

That's great, but we probably don't want all the files just in the root. We want some organization, so we can also create a subdirectory here called whatever we want. Let's say _"path 'hello_world'"_, and I just run this again. _"nextflow run hello-world.nf"_. It should go into the results directory into a subdirectory and sure enough, now under results here at the top we have _hello_world/_ and we have _output.txt_.

Important thing to notice, the old _output.txt_ file is still there. The results directory is not wiped when you do this. Just new files are copied in there. They'll overwrite files which are already there if they have the same file name, but they won't clear out old ones. So you need to be a little bit careful about when you rerun pipelines. If you don't want them to be on top of the files that are already there. Make sure you use a blank empty directory.

## 2.3. Set the publish mode to copy

Okay, I mentioned that these files are soft links, so if I do _"ls -l results/hello_world/"_, you can see it's soft linking to the work directory. That's generally a good thing if you're working on something like HPC, and these are really huge files and you don't want to duplicate them, because it means that the files only stored once on the file system.

However, it does mean that if you delete the work directory: if I do _"rm -r work"_ and clear out all those intermediate files which were created. Now, if I try and read this file _"results/hello_world/"_. It's gonna be pointing as a soft link to a file that no longer exists and the data is gone forever and is irretrievable, which is maybe not great.

So generally we, I say it's good practice to copy the files instead of soft linking if you can, because it's safer. Just be aware that it will use twice as much disk space unless you delete those work directories.

To do that with the output block, I'm gonna go to the first output here. I set the path before and now I'm gonna set the mode and you can see as I type, the VS code extension is, suggesting stuff it knows it's an output directive here. And I'm gonna say copy. I hit save.

Let's rerun the workflow. It is gonna create the files again, new work directory.

Now, if I go to _"ls -l results/hello_world/"_ you can see this is a real file and it is not a soft link anymore, and Nextflow copied that. Good to know. So path and mode are things you'll find yourself writing quite a lot.

Now, of course, this is very simple. We will make this more complex and powerful as we go along, and you'll see how to make these things dynamic and not too verbose.

## 2.4. Note on process-level publishDir directives

Now, I said as we started on this, that this is a fairly new form of syntax. It's only available in the latest versions of Nextflow as I record this, and it's called Workflow Outputs.

If you use this, it's great. It unlocks lots of other cool features within Nextflow, such as, Nextflow Lineage to help track the heritage of these files as they're created, and soon will be the default in 26.04. And at a later date in the future, this will be the only way to write your workflows.

However, as we're in this transition phase right now, you might well see pipelines in the wild, which you use something called publishDir, which is the old way to do it, and this is defined not at the workflow and output level, but this is defined at process level.

And this declaration says basically the same thing. It says, publish the results files into a directory called results, and use a copy mode. So you can see the syntax is very similar. But when you're writing new pipelines now, try not to use this publishDir directive, even if you see it, in AI results or in documentation or other pipelines, because that's the old way to do it.

In 2026 we should all be using workflow outputs.

This is all documented, if you are doing this and you've used Nextflow before, you can go to the Nextflow docs here, nextflow.io/docs/. And if I scroll down to tutorials, there's a tutorial called _Migrating to Workflow Outputs_.

It is really good. It goes through all the syntax, how it's equivalent to the old syntax, why we changed it, and, have a timeline and everything. And it goes through all the different scenarios with loads and lots of examples. So you can easily convert existing Nextflow code over to the new syntax.

## 3.1. Change the sayHello process to expect a variable input

Okay, so we've got our simple script, which is running a process, creating a file, telling Nextflow it's an output, and then we're telling Nextflow where to save that file. That's a good start.

But it'd be more interesting if it wasn't all hardcoded. So next, let's think about how to tell Nextflow that this process can take a variable input, which is something we can control at runtime when we launch a workflow.

We need to do a few different things to make this happen.

Firstly, we need to tell this process that it can accept an input variable and we type _input_ here as a new declaration block. And we're gonna call this _"val greeting"_.

The val bit is the equivalent of a path down here. It tells Nextflow this is a variable, like a string in this case. And if you hover over it again, it tells you from the extension of what this means.

Next we're gonna tell Nextflow what to do with this. It's not, enough just to say that there is a variable. You have to say in the script how to use that variable. And so I'm gonna get rid of this hardcoded string here, and I'm gonna put in a variable.

I'm gonna quickly do it without squiggly brackets just to show you that this is, allowed, and this is the old, style way of doing it. But now with the new syntax, we really recommend putting it inside squiggly brackets like this, and it makes it really clear that this is being interpolated by Nextflow here.

Great. So _"input greeting"_ goes into _$\{greeting\}._ Last thing is we need to tell Nextflow at the workflow level that this process now takes an input. And to do that, we're basically going to give it a variable.

## 3.2. Set up a command-line parameter to capture user input

We could hard code it again, like Hello World, and that would work fine, but obviously it doesn't really give us any advantage. We wanted to be able to configure this at run time, so we want to be able to do it on the CLI, when you launch Nextflow.

And the way we do that is a special Nextflow concept called _params_. We're gonna call that _params.input_.

What this does is it exposes this input variable on the CLI and that's where we use a double dash when we launch Nextflow.

I can call this whatever I like, I can call it _hello, greeting_. Doesn't matter. Whatever I do there will be exposed as a CLI option when we launch a pipeline. And this is a real magic trick by Nextflow 'cause it means you can build your workflow script very quickly with these parameters, and you're essentially building out a custom CLI for your pipeline, making it really easy to customize different options on the fly when you launch.

So. Let's try it out. Go back to our terminal. We have our _"nextflow run"_ command here. And now I'm gonna do _"--input"_, which matches the _"params.input"_ we saw before. I think in the docs it's in French. Geraldine likes to speak French. I'm gonna do it in Swedish 'cause I live in Sweden. so I'm gonna say, "_Hej Världen_" and hit enter.

Can use single quotes or double quotes, it just affects how Bash interprets it.

It runs the Nextflow pipeline exactly the same way. You can see the working directory and everything is the same. But now if I go up to _"results/hello_world/output"_. We can see our nice Swedish here instead.

So we have dynamically passed an input from a CLI to a parameter. We've passed that as an input to the process and the process interpreted that and put it into a script block, which has then dynamically changed the output of that script result. Pretty cool.

Quite complex logic with very, little syntax here. And you can hopefully see how this now starts to scale. And this is how we really build the logic and the customizability of our pipelines into the Nextflow script.

## 3.4. Use default values for command line parameters

Okay, that's great. The problem though now is, every single time I run this pipeline, I need to do dash, input for it to run.

If I try and run without this parameter, now Nextflow is gonna throw an error saying it needed this parameter and it wasn't set. and so it didn't know what to do.

This is a cool new thing, by the way. In the past, Nextflow would've just run with an empty string, and you'd have had all kinds of weird errors, which would've been difficult to understand. But in the new Nextflow syntax parser, it is a bit more careful and it tells you right away.

So we don't always want to specify every single option. It's good practice to specify sensible defaults. So how do we do that in our script?

You'll notice that when we wrote this, we just put _params.input_ straight into where we're using it. So the obvious solution is we define a default, and we do that up at the top of the script here in a special params block in the workflow. This is in the workflow script here.

Again, some new syntax here, so pay attention. This is really cool stuff. We've got the name of the parameter, which will be expected here.

And then after this colon character, we're defining a type of the variable. You don't have to do this, you can just leave it blank, but it's really nice. It tells Nextflow that we're expecting a string and treat it as such.

If we want a number instead, for example, we could write float, and that would say we want a floating point number. And if we try and run with that, then it will throw an error. If we give it a string, which is not a float. And it'll also pass it as such. As if we do string, then it knows it's a string. And even if it has leading zeros and is all numeric, it will still pass it as an actual string.

So that type safety is a very new feature of Nextflow, but really powerful to make your code safer to write and to run.

Then after that we've got an equal symbol and then the default value here. Nextflow is written in Barcelona originally, so it seems appropriate that we've got some, Spanish here, _"Holà mundo!"_ as a default.

Right i'm gonna save that script, go back, run the script again without _--input_. And this time it should run and it will create our new file up in _results_. And in this file now it says _"Holà mundo!"_.

This is just a default though, so it doesn't mean that we can't still do the same thing as before. If I go back and find my old script here, _"Hej Världen"_, because I do _--input_ on the command line, that will overwrite that default and use that again in the output.txt file.

So this in the script is only the default value that I'm setting.

As we build up our workflow to be more complex and include more parameters, this params block at the top of the script will start to collect all of them in one place.

And you end up with this quite nice symmetry in your script, where you effectively have all your workflow inputs here and your workflow outputs down at the bottom. And it's very clear what the interface of your workflow is to the outside world. So you can pick up a new pipeline very quickly with the new syntax and understand how to use it.

One last cool thing. We don't have to set a default value with this. If we do params input but don't set a default value, then it tells Nextflow that this parameter is required, and again, the pipeline will fail to run without it, but it'll give you a more useful error message rather than something about it being null.

So it says we are expecting its input is required, but it wasn't specified on the command line. Very nice.

Okay, so hopefully now it's clear about how to set up your Nextflow pipeline with variable inputs and parameters, how to set the default, set, the types, it could be a Boolean true false flag or an integer or different types here. How to pass them into your workflow, where it goes through, and then interpolates into your process. And then you also know how to customize those on the command line when you launch Nextflow. This is starting to look more interesting than our simple bash command.

## 4. Manage workflow executions

Okay. What's next? For the final part of this chapter, we're gonna talk a little bit about how to manage all different workflow executions. If you look in my sidebar here and the Explorer underneath work, you'll see I've run a bunch of different pipelines and these work directories are getting quite long, there's a lot of them.

And the other thing is, as I said before, every time I rerun this pipeline, it's creating a new set of work directories, and it's rerunning all the processes from scratch, which is a good thing. That's intended behavior. It's reproducible and it's regenerating everything fresh. But it obviously, if you're running very long running processes, it's annoying to always have to start your pipeline from the beginning if it crashed halfway through, or if you change something at the end of the pipeline.

## 4.1. Re-launch a workflow with -resume

Luckily, Nextflow is really, good at knowing what has previously been run and what's available, and to reuse those old results is very, simple. We just add a, new flag at the end of the command _"-resume"_.

Now, note there are two hyphens on input 'cause that's the parameter. There's only one hyphen on resume because that's a core Nextflow option.

It trips people up all the time, even if you've been using Nextflow for a long time. So always remember one or two hyphens. Depends if it's a core Nextflow option.

Okay, so now I do _-resume_ and I run exactly the same workflow again. And this time it should look pretty much exactly the same with one key difference.

In the output here, you can see that the results were cached. And in fact, this task hash here is exactly the same as the previous run, and it's just reused that work directory in its entirety. The inputs and the outputs and the script were all unmodified. And so it just takes that file from that and if there are downstream steps in the process, it would pass them onto the next step in the pipeline.

So it's still running the whole pipeline from start to end, but it's using cached results for each one of those tasks, where it can.

Now, when you do _-resume_, it just resumes the last pipeline run in your working directory, whatever that was. But you can actually resume from any previous run that you've done there. And we've done quite a lot now.

## 4.2. Inspect the log of past executions

To look at all of them, we can do _"nextflow log"_ instead of _"nextflow run"_, and that will give us a nice output showing all of these different.. I need to make my screen a bit smaller so we can see it, all of these different runs when we did them, the session id, the command and everything.

And we can look in here and we can take the run name of any of these and then resume one of those specific ones. So I can go back and I can resume that one called _hungry_ekeblad_. And I just put that after the _resume_.

If you're curious, by the way, all of these adjectives and scientists names are in the Nextflow source code. It's a really good way to get your first ever pull request to Nextflow by going and finding it and adding your favorite scientists.

And anyway, so I did that and it went back and it looked at the cached results from this workflow run, realized it could still reuse them, and it did. So I got the cached results again.

## 4.3. Delete older work directories

That's great. What about if I want to clean up these work directories? There are loads of them here. There's loads of files. Maybe I know for a fact that I want to resume from the last couple of pipeline runs, but I don't care about all the ones before that.

Then I can pick one here and I can use another Nextflow command, which is _"nextflow clean"_, and I can do _"nextflow clean"_, I'm gonna do _"-before"_, and the particular run name, which in this case was _reverent_pike_ and I'm gonna do _"-n"_, which tells Nextflow just to do a dry run. So it just tells me what it will delete. Without actually doing anything, so it would remove these work directories.

That looks sensible. So I'm gonna do the same command again, but instead of _"-n"_ I'll do _"-f"_ to actually do the cleanup. And this time it's actually removed all these directories. And if I go in and look at the work directories, it's now looking a lot lighter. Fantastic.

So that's how to clean up all of your local work directories in a pretty safe way without completely destroying the cache. So you can still resume if you want to.

If ever you forget what these flags are for every Nextflow command you can do _"nextflow help"_, and then the name of the command. So if I do _"nextflow help clean"_, you can see all the different options: _-after, -before, -but_, all different ways to configure this cleanup behavior. Pretty cool.

## Takeaway

Okay, that's the end of part one of Hello Nextflow. It's quite an intense start to the course, but hopefully now you have a pretty good understanding of what a Nextflow script looks like; with different key parts, the processes, the workflows, the outputs, and the parameters. You know how to configure them with basic overrides from the command line, how to make a dynamic input block with a dynamic script and you know how to manage all of your workload executions: seeing what you've already run, resuming, cleaning up. There's a lot of stuff. You've come a long way. So if you want to take a break and have a quick walk around and a cup of tea, now is probably a good time. You've earned it.

From here on, we are basically building on this foundation. How can we make this more complex, more powerful? How can we make it more flexible? Do the things we want to do our analysis at scale.

## Quiz

Now if you scroll down to the part one, hello world, on the webpage you'll see a little quiz and this is something new that we've done for this version of the Nextflow training. And you can go through and quiz yourself to check that you've understood all the material that we've done in this chapter.

This isn't sent to us or anything, it's just stored in your browser. So we don't know what your answers are, but it's just a little self check to make sure that you haven't missed anything or misunderstood anything. And you can try it as many times as you like.

If you're like me, maybe you want to stay in the terminal in your VS Code instance, in which case you can type the _quiz_ command and then just tell it which, chapter you're on. So we do _"Hello World"_, and then you can do exactly the same, quiz questions, which are in the web browser, but just in your terminal.

Cool. Okay. Hope you enjoy that. Have a bit of fun and, we'll see you in the next chapter in just a minute to talk all about Nextflow channels.

​
