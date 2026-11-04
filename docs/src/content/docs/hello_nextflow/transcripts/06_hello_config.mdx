# Part 6: Hello Config - Video Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/FcZTiE25TeA?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important notes"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../06_hello_config.md).

    The section numbers shown in the transcript are provided for indicative purposes only and may not include all section numbers in the materials.

## Welcome

Hi, and welcome back to Part Six of Hello Nextflow. This section's all about configs, and it's the last part of this course.

Nextflow is particularly good at two things, reproducibility and portability. Configs is where we see the second of these really shine. The ability to configure a Nextflow pipeline to run in different ways and work on different systems, without having to edit the underlying pipeline code.

This superpower allows Nextflow pipelines to be reused by other people in different places, or across different infrastructures that you might have access to yourself.

It means you can develop pipeline code on your laptop, push it to the cloud, run it on your HPC, and it's the same pipeline code and it runs everywhere.

In this section, we're gonna go through a few topics. We will start with how Nextflow handles config files, where it loads them from, and how you write them and how you structure them, and that separation between the pipeline itself and what should go in a config file.

Then we'll go onto some common use cases such as changing where output files are stored, and also how to get the pipeline to work on different infrastructures, both using different types of software packaging or submitting jobs to different infrastructures.

## Config file hierarchies

Okay, let's get started. When it comes to loading config files, Nextflow can pull from many different places, which is a good thing and also can be a slightly risky thing because sometimes it can be a little bit difficult to know where it's getting a config file from and what order it loads things in.

So I really recommend that you do click on this link here, which takes us to the Nextflow docs. And on this configuration page, it lists the key places that config is loaded from, and importantly, the order in which these things are loaded.

So you can see, you can put a config file into your Nextflow home directory, which is typically ".nextflow" in your home dir. And that file will always be loaded by every Nextflow run on your system.

The next place to look is a file in the root, of your pipeline, repository or directory called "nextflow.config".

Then after that, another file called "nextflow.config", but this time in the directory that you're launching Nextflow from: the launch directory.

Finally, you can provide config file paths on the command line with a "-c" argument, and you can do that multiple times. And they're applied in the order that you specify them.

You can, provide config files in all of these locations if you want to, and they'll be loaded iteratively, each overwriting the previous in just the config scopes where they clash.

This is a really, powerful system because it means you can set sensible defaults and then get gradually more and more specific as you narrow down on that config.

## 0. Warmup: Run hello-config.nf

Okay, let's close this and jump into our Codespaces and get started. As before I've cleaned up here, I've removed my previous results directories, my Nextflow, and my work directories and so on. Don't worry if you still have those files kicking around. It's just 'cause I'm very zoomed in and so things get messy very quickly otherwise.

We're gonna be working with hello-config.nf, the last file in our directory, and this should follow on from where we left off in the previous section.

So we've got our four different processes, which are being included from module files. We have our pipeline parameters, our workflow block where we're calling the different processes and stitching the channels together, publishing the output channels, and then the output block at the bottom where we define where those files should be stored and how they should be copied.

We also already have a "nextflow.config" file from the last chapter, where we enable Docker, and we're gonna be building into this file today.

As before, we've changed the output path in this main script to hello config, just so it doesn't clash with previous results that you've generated.

Okay, let's just quickly check everything is still working as we expect. Bring up a terminal and we do nextflow run hello-config.nf. Nextflow loads up. Should run our four different processes. Generate some nice asci artwork using cowpy and then save our results to our results files in that directory.

I can have a quick look in here just to make sure that these files look as we expect, and sure enough, there is our giant Turkey. Great.

## 1.1. Move default values to nextflow.config

Now the first thing we're going to do is start to move some things from our script into our config file.

And what we care about is mostly the parameters at this stage. We want to take the default values, into the config file, so it's clearer what the defaults are and it's easier for people to overwrite them.

I am gonna just take this params block here from the script and put it into the config file. And we need to be a bit careful here, 'cause right now the syntax is slightly different between config and scripts. The config file cannot take type declarations 'cause we're not really defining these params, we're just referencing them. So I'm gonna get rid of those.

But otherwise it's very much the same. We have a params block and then we have our different input parameters, batch parameter, character parameter.

I can now go back to my script and I don't need to define these defaults anymore because these values are now in my Nextflow config file.

However, I do leave the parameter names and their types, so that Nextflow knows that information and can still do all of the type safety and everything.

Okay. We save those files and quickly check that everything still works the same as it did before. There shouldn't be any changes here. We have kept the values the same. We've just moved where they've been defined.

Great.

## 1.2. Use a run-specific configuration file

Now, so far we've been launching Nextflow from the same directory where we have our pipeline script. So our launch dir and our pipeline dir are kind of the same thing.

To show how we can have different config files with different launch directories, we're gonna create a new subdirectory now.

So I'm gonna say mkdir, and we're gonna call it tux-run.

And then I'm going to cd, change directory into tux-run. And note that we are now in our working directory is now no longer in the same directory as the pipeline scripts.

Okay, let's create a new "nextflow.config" file. So touch Nextflow config, and let's just open it up in VS Code. You can see also in the sidebar over here that we're now in this subdirectory.

Now we can take the same params block that we had in our top level nextflow.config, copy this over and now we can change these values.

Firstly, the data is now a different relative path because we're in a subdirectory, so we need to update that. And then we're gonna change batch to experiment, and we're gonna change the character from Turkey to tux.

Now click save there, and let's try it out. Just as with data, I need to now say ../ to get to the script. So it's Hello config. And I press enter.

The pipeline code hasn't changed at all, but now we're gonna have two sets of config loading, and the launch dir config file should overwrite the defaults, which were set in the pipeline nextflow.config, and we should get different sets of results.

Sure enough, within our directory here, within tux-run, you can see we have a dot Nextflow directory and a work directory and that's because these are created always in your launch directory. So these are different to the work and results ones that we had from earlier runs.

Now, if I look in results, we can see our collected and there is our little tux character. So you can see those parameters were properly interpreted.

## 1.3. Use a parameter file

Okay. Before when I was talking about the different config files that could be loaded, I missed out one other place we can get config.

You can get it from a command line as we've seen with dash dash parameter names, but we can also supply a YAML or a JSON file, just of params.

The config file can have all different types of scopes, but these files are just parameters, and it's a nice user-friendly way to supply many parameters at once, and perhaps a little bit more reproducible way because you write them to file, so it's easy to get them at a later stage.

So let's go back to our terminal and just before we forget, make sure we move back up a directory, so I'm no longer in the subdirectory, and i'm going to look at the YAML file that we have here called test-params.yaml.

So if I just do code test-params.yaml, you can see this is just a regular YAML file. Nothing special about it. With the keys being our parameter names, with the YAML formatting so a colon here, and then a value.

Note that this is not Nextflow code, so we can't put things like variables in here. These are just static values.

Also because JSON actually parses as YAML, we all can also have a test-params.json file, which looks very similar. It is just, different data format.

So we've got two different test files here and we've got slightly different variables.

Okay, so how do we give these to Nextflow? It is very simple. We do Nextflow run hello config, as before. And instead of "-c" for config file, or loading those default file names, we do -params-file. Single hyphen because it's a core Nextflow option.

And then pass the path for that file. So I'm gonna do "-params-file test-params.yaml", and we'll see if those are properly loaded.

Okay. It ran. Let's just remind ourselves what was in this YAML file. So the batch was set to YAML, so that's what it should be called, and it should have a stegosaurus. So let's go up and look in results. And we have COLLECTED-yaml. So let's see if we have a Stegosaurus. Fantastic, a Stegosaurus wearing a hat. That's what we like.

So that's worked really well, and is exactly the same with the JSON file. We just switch out the file extension here and Nextflow knows how to, read that.

And in this case, we should have a batch called JSON and we should have a turtle. Let's have a look. Wonderful. One of my favorite CLI tools.

## 2.1. Customize the output directory with -output-dir

Okay, so that's been mostly thinking about inputs to the pipeline and changing parameters. what about the outputs?

Now, although we've been changing the sub directories using params, you might have noticed that all of our files are still going to results.

We can change that base directory that all files are published to with a command line flag called -output-dir. So if I do Nextflow run hello config, and then I do -output-dir, and we're gonna call it "custom-outdir-cli". Can't type. Just so we remember where these files came from.

This is a core Nextflow option and it's a very new one. This was only recently added, and this is one of the things that we can do with the new language parser and everything.

It's a bit of a mouthful to type. You can also just call it "-o" if you want to. So if I just go back. Can just shorten that to "-o", which is a bit simpler.

Okay. We run that. We haven't changed anything in our pipeline or even in our config at this point, and it should hopefully save all of our results into a different top level directory. And you can imagine you can set this to basically any path that you want.

It is just arrived up at the top. We've got a custom-outdir-cli, and all the files are organized in there in exactly the same way, with their same sub directories and file names. So this is a really easy way just to change where the pipeline publishes its results to, without thinking too much about how those results are organized.

## 2.1.2. Remove hardcoded paths from the output block

If I look into this directory, we can see that we still have a subdirectory called Hello Config, which feels a little redundant now.

So let's just load up our script again and we can now remove that subdirectory from the output block at the bottom. 'cause we don't really need it anymore. So we can just do that now, delete that from here. And then if it's just this, you can either delete that completely or leave it as an empty string. I'm gonna leave it as an empty string for now, because we're gonna come back and put some different things in its place in the future. But if you don't care about sub directories, it's cleanest just to completely remove the path declaration there.

Okay, let's hit save. Just quickly try it again. I'm actually gonna remove my "custom-outdir-cli" directory so we're not confused by any existing files there. 'cause remember, when you publish things, it doesn't remove the files that were there already. It just adds new ones. Let's run that command again, custom-outdir-cli.

And now if you do "ls custom-outdir-cli", there's no more directory there called Hello Config.

## 2.2.1. Set outputDir in the configuration file

Okay, the command line flag here, "-o" or "-output-dir" is good. But how about setting defaults for this in config? How do we do that?

I open up the "nextflow.config" file, close everything else and get rid of that. We can add in a new, config option here, which is I've just copied from the training material website, and it's called outputDir.

It's not under any scope. It's not under params or anything. It's top level, and we can set this to a string. Now a simple thing to do is just to change it to anything other than results as a hard-coded string. But because this is in a Nextflow config file, we can be a little bit clever here and also include variables.

And you can see here that we've included a params variable, params.batch, which is part of this string. This means that we can reuse variables which are coming in from other places. And in this case, if we do --batch, when we run Nextflow Pipeline, we're gonna get a subdirectory in our custom path based on what the batch name was.

Okay, so let's try this out and just have a quick look to see how it, how the results look. So if I do Nextflow run hello config and --batch my_run. Remind ourselves what the config looked like. So it's custom-outdir-config.

Tree custom-outdir-config. And you can see the batch was called my_run. And then we have that subdirectory called my_run. So that dynamic file path worked.

And not only that, it didn't go into a default results directory anymore, and I didn't have to specify anything on the command line to change the base directory. So we've successfully reset the default value for the default outputDir.

## 2.2.2. Subdirectories with batch and process names

Okay, let's take that a little bit further. That's a dynamic variable within the config file. How about the, the script itself? Now, so far we've had these paths here and these can also be dynamic. So instead of just hard coding something, we can put in some squiggly brackets and put something dynamic.

So for example, we have our processes called sayHello. We could do sayHello.name, which is an attribute of the process, which is kind of boring 'cause it's just "sayHello" in this case. But it's variable.

So this gives you an idea. So we can put this in here and say convertToUpper.name, collectGreetings.name, collectGreetings.name again, and cowpy.

Now when we run, the base directory is still gonna be custom-outdir-config. And it's gonna be in a subdirectory called params.batch, but the sub directories under that should be organized by process name.

Let's just try that out and see if it works. So I'm gonna remove the previous directory so we don't get confused, and just use exactly the same Nextflow Run command.

It should run in the same way. I could be using dash resume on all of these to make it a little bit faster and use the previously calculated results. Now, if I do tree custom-outdir-config, you can see it's not in results, it's in our base directory with, the batch name. And you can see all of the results are now organized within sub directories named after the process. So we've got two different places where we're defining dynamic output paths here.

Okay. Final thing, let's add back those intermediate folders, which we had before 'cause they were kind of nice. Intermediates.

And we can also think a little bit about this params.batch, maybe as a pipeline developer I really liked having that in the subdirectory, but if end users of the pipeline are setting "-o" or -output-dir on the CLI, it's completely overwriting this entire statement, and we lose that subdirectory.

So what we can do is we can take that dynamic path out of the outputDir config, which would be clobbered, and put it into the output path, which is not clobbered.

So we can do params.batch slash intermediates slash sayHello.name, and do all of this in a double quoted string, so it's interpolated by Nextflow.

Can now copy, whoops. Copy these down to the other processes. Remember to put them all in quotes. And remove intermediates from these particular outputs.

Okay? It's looking slightly more complex now, but you can see we are really starting to build up nice organized output directory structure in our code.

And what's really nice is that this extra complexity in the code that doesn't pass through to the CLI. So we can run our command with -output-dir and whatever batch variables, just thinking about the how to run the pipeline and not really thinking too much about what's in the code. And our output files are gonna be constructed really nicely in a very well organized way, which is nice for people using the pipeline basically.

Great. As I write this, I realize that I made a mistake. See if anyone caught me out here.. We have collectGreetings.name, so something's gone a bit wrong. And yeah, sure enough, I accidentally forgot to put these in squiggly brackets.

So remember, be careful when you're writing your code and make sure that you tell Nextflow what is a variable and what is just a string. 'Cause it will do exactly what you tell it to do. And nothing more. Like all good computers. Okay, that should fix it.

## 2.3. Set the publish mode at the workflow level

There's one bit of this script, which I don't love still, which is the fact that we're writing mode copy again and again, and if there's one thing we don't like, it's repeating ourselves.

So we can clean this up a bit by taking this and moving it into the config. And in fact, we can set it for the whole pipeline in one go. So we don't have to say it multiple times.

We go over to our config file and we've got a new scope here called workflow. And we can either do squiggly brackets or we can do dot notation. Doesn't make any difference. The training material website uses dot notation. I can say output and we can mix and match, so mode equals copy. Great.

And now we can go back over here and delete these. Now we could leave them in place. The config is basically overwriting what's written here, but as we have it in the pipeline level config, and these two files ship together, there's no reason to really do it twice.

Okay. Just a sanity check ourselves, because apparently we do make mistakes. Let's, run that again and just check that we are correctly using the copy mode for publishing files. So we're gonna run the script again and this time we've put the results into a directory called config-output-mode, see what the files look like in there.

And then if I do "ls -l" to look at batch, and we can look at cowpy, for example. And we should see, yeah, that this is a proper file here, which is not a soft link, so that config attribute has been applied properly.

## 3. Select a software packaging technology

Okay. So far we've been focusing on the inputs and the outputs, the files that the workflow is running with. But how about infrastructure? I said at the start that Nextflow allows you to run the same pipeline on different computing setups. So how does that look?

To show this, we're going to switch from using Docker to run cowpy, and instead we'll use Conda to do the same thing.

I can do this very simply. If I go to code, " nextflow.config". If you remember at the top, we defined docker.enabled earlier on, and the last chapter so that we could use the container with cowpy in.

I'm gonna tell Nextflow not to use Docker. Set that to false. And I'm gonna say Conda enabled equals true. So tell Nextflow, please use Conda.

Now just enabling Conda is not enough by itself. Just as we did with Docker, we have to tell Nextflow where it can get the software it needs.

So if we hop into the modules here. And open up the cowpy script. We can see we have a container declaration up at the top. And the container is used by Docker, but also Singularity, Apptainer, and many of the other software tools.

But it can't be used for Conda, so we have a separate declaration called "conda", and we could just write "cowpy". And that will leave it up to the conda package resolution to figure out the best way to solve, that according to your local conda environment.

Or it's good practice to do what the, the training material website says to do, which is to define a specific conda channel with its, double colon notation, and definitely define a specific version of the software so that every person who runs the pipeline will get the same version.

Note that containers are a bit superior in this respect, because when you install something with Conda, it's still gonna figure out all of the dependencies for that package, and they can change over time. Called dependency drift.

So containers, however, lock the entire stack of a whole software dependencies all the way down, so you can be a bit more confident that A, it's gonna work, and B, it'll be reproducible.

So if you're able to use Docker or Singularity or Apptainer, I would definitely recommend that.

Now what's nice about this is the module file, which is written by the pipeline developer, now has both Container and Conda, and so we're telling the person who's running this pipeline, we don't mind what software packaging solution you use. It will work with both Docker and with Conda, and this is where to get the software in both cases.

We can pull up the terminal and let's give this a try. So Nextflow run hello config --batch conda. And the first time this runs with conda, it's gonna be a little bit slow when it gets to that particular process, because it has to run "conda install".

And it's creating a special conda environment just for this one process. So it's not using my global conda environment, which I have on my terminal. It's creating one just for that one process. This is good because it avoids things like dependency clashes between different processes in your workflow. If your processes have tools which need different versions of Python or things like that, that's okay because they're using different conda environments.

Nextflow caches these conda environments locally, you can see it tells you where that path is, it's in the work directory here. And so the next time I run this script with Conda, it'll be much faster 'cause it will find that existing conda environment and just reuse it. But the first time we do it, it has to go and fetch it, resolve it, download all the dependencies, and set everything up.

Okay, great, it ran. We can just remind ourselves what the pipeline is currently configured to use. If we look in the config file, it was "custom-outdir-config" right now for me. See if I go up to that base directory. And I did --batch conda. There's our conda subdirectory. So it worked and there's our cowpy output.

So it fetched cowpy, installed it on my local system using conda, and ran the process. And what's great is, as that end user, I didn't have to think at all about any of the software management there. Nextflow just sorted it for me. I said, I need to use conda on this system. The pipeline developer said which packages I needed. And Nextflow did the rest. Very powerful.

Note that you can actually use a mixture of different technologies. So I can enable Docker for specific processes, and conda for other processes, or say that some processes should just use whatever local software I had installed. This is pretty unusual, but it is possible, and in some cases, for example, if you're using certain software that might be difficult to package in Docker, you do have an escape

## 4. Select an execution platform

So that is software packaging. The other part of portability to other systems is where the jobs actually run. At the moment, I'm running on basically my_laptop or in this Codespaces, which is a single computer. There's nothing fancy. Nextflow is being a little bit clever about parallelizing the jobs as best it can, but it's all on one system.

Now, if you are running on an HPC, you probably have some kind of job scheduler such as SLURM or PBS or something, and you will submit jobs to that scheduler and it will farm all the jobs out to different compute nodes.

Another way of running is on the cloud. So maybe you're using AWS Batch, or Azure Cloud, or Google. And these all work in a similar system where you have a scheduler and you submit jobs and they are submitted to different places to be computed.

Now in the long distant past when I started doing bioinformatics, everyone's software for running analysis was very tied to their computational infrastructure, which made it almost impossible to replicate.

But with this config separation in Nextflow, and with Nextflow's ability to interact with very many different compute infrastructure backends, it's very, simple to take our pipeline without modifying the pipeline code at all and just switch that out.

## 4.1. Targeting a different backend

So if we go to our "nextflow.config" file, and we can now put in some process level config. So if I put up at the top process scope and I can set the executor, and here it's set to local, which is the default.

Notice because this is process level, we can target things to different processes. And so you can actually set up executors to be process specific and have a hybrid execution, where some jobs might run locally, wherever the Nextflow job is being executed. Some are submitted to, different HPC and some might be submitted to the cloud. You can be as clever as you like.

Now, it's very difficult to demo this in a training environment like this because I don't have an HPC to submit to. But I can do is if I type in slurm, we can cheat a little bit and you can get a feel for this.

And this is really most interesting for people who are used to running on SLURM and know what the SLURM headers look like. But if I do Nextflow run, hello config. It is gonna fail because it's going to try and submit jobs to a cluster which doesn't exist. So we'll get some kind of error about sbatch not being available.

Yeah, written. That's the tool. That's the CLI tool that you use to submit jobs to a slurm cluster. But what we can do is we can go and look in our work directory here by command, click, open that directory and look at the .command.run. And you can see at the top of the .command.run file, we have our sbatch headers, telling a theoretical SLURM cluster how to handle this job submission.

And so you can see that Nextflow is being clever, it's doing all the right things. It's just that we didn't have a cluster to submit to.

## 5. Control compute resource allocations

What else is different between different computing infrastructures? Another thing is how much available resources you have, and in fact, in many compute environments, it's a requirement that you have to specify how many CPUs and how much memory a job needs.

Again, Nextflow abstracts this for us, so that it's no longer specific to a single compute environment type, and we can type in the process level scope here. CPUs equals one, memory equals two gigabytes. Our pipeline's not very demanding, so that should be fine.

Now, I've just guessed these numbers here, but how do you know what is a sensible amount of resources to use? It's quite a difficult job to go off and dig through all these different processes of a big pipeline of many samples and understand what the resource utilization was.

So a good approach for this is to set these values to high numbers to start off with, just so that your pipeline runs without any errors, and then ask Nextflow to generate a usage report for you.

This is super easy to do, so I'm gonna go back to a terminal. Oh, I'm need to remember to set that back to local so that my pipeline actually runs. And I'm gonna say Nextflow run, and I'm gonna use a command line flag -with-report.

And I can leave that blank and it will give a default file name, but I'm gonna give it a specific file name so that, that's saved to a specific place.

Hit Enter, and the pipeline runs exactly as normal, but when it finishes, it's gonna generate a nice HTML report for me.

So in the sidebar here, I've got this HTML file. If I was running this locally, I'd just open it. I'm, 'cause I'm in Codespaces, I'm gonna right click on that and click download, which is gonna download it to my local computer. And I can just easily open it in the web browser.

Nextflow can generate a report like this for any pipeline and it's got some really nice information. So it's good practice to always save these things. It tells us when we ran, where we ran, whether it's successful or not, what parameters were used, what the CLI command was, things like this.

And there are also these plots about resource usage. So it tells us what percentage of, CPU calls were used for each process as a box plot here, because there are many tasks for each process, so we can see the distribution.

You can see our processes here, cowpy and collectGreetings only had a single task, so it's just a single line. And we have both CPU and memory and job duration, and they were very fast.

If you're using Seqera Platform, by the way, you get the same plots built into the Platform interface without having to do anything. So you always get this information at your fingertips.

Okay, so we can use this report and on a real run, and get a feel for how many CPUs and how much memory is being used by our pipeline and come back and, and put those values back into our config file, so that next time maybe we don't request quite so much. And we can be a little bit more lean.

Now you can get really clever about configuring pipeline config files. And again, if you're using Seqera Platform, look out for a little button that looks like a light bulb. 'cause if you click that, it will generate a highly optimized config file, which is tailored specifically for your data, your run and your pipeline. To run it in the most efficient way possible.

But for now, I'm gonna say that actually the default number of CPUs that Nextflow was giving was fine and but only need one gigabyte of memory.

## 5.3. Set resource allocations for a specific process

Now, in real life, it's quite unusual that all of the processes in your pipeline are going to need the same requirements. You might have something like MultiQC as a reporting tool, which needs very little in terms of resources and runs quite quickly.

And then maybe you have something which is indexing a reference genome or doing some alignment or doing some other job. It doesn't matter what it is, which takes a lot of resources. And so for these different job submissions to a scheduler, you want to give different amounts of resources.

Under this process scope, we can define a config, which targets specific processes in different ways.

Here we're using withName, we can also use labels, and these can use a pattern to target one or multiple processes. Here we are just saying any processes that have got a name cowpy set to two gigabytes memory and two CPUs, and because this is more specific selector than the top level process one, this is overwritten in these cases, so you can build up a nice config file here, which really tailors all of your different processes in your pipeline to make them really efficient.

## 5.5. Add resource limits

Now as a pipeline developer, I probably know the tools quite well, and I want everything to run as fast and as well as possible. So it might be that I put in pretty high numbers for some of these because I know it'll run much faster if I give cowpy 20 CPUs instead of two.

That's fine until you go to run on your laptop or on GitHub Actions Continuous Integration test, or some other system, which maybe doesn't have 20 CPUs available.

Now when you try and run the pipeline, it will crash because Nextflow will say, I can't submit this job anywhere. I don't have the available resources.

Now to avoid that hard crash, we can add a bit more config, which is specific to our system now, called resource limits. And that looks like this. It's under the process scope again.

And resource limits, you can specify basically the ceiling of what you have available. It's a map here, and you can, within this map, you can set the memory, the CPUs, and the time.

Now what happens is when Nextflow submits a task from a process, it looks at what's requested and it basically just does a minimum between that and that. So if we requested 20 CPUs, but only four are available, it will request four. The pipeline doesn't crash and it uses as close to what it was designed by the pipeline developer as possible.

## 6. Use profiles to switch between preset configurations

Okay. I said that the resource limits here might be system specific, and maybe I have an Nextflow config file in my pipeline, and I know that people are going to be using this in a range of different places. Now, instead of forcing everybody to create their own Nextflow config file every single time, what I can do is I can group different presets of configuration together into config profiles.

I am gonna scroll down a little bit here and then just past params, because the order of the config file here is important, the config file is loaded sequentially, so I'm gonna put these profiles after everything else so that it overrides the previously defined params. And I'm gonna paste in these profiles from the training material.

So there's new top, top level scope called profiles. We can have arbitrary names here. So we have my_laptop and univ_hpc. And here we can see we're setting the other same config parameters that we were before. Now within just a profile. So we've got a local executor for running on my_laptop and I'm submitting to a SLURM cluster on the HPC.

I'm using Docker locally, conda on the HPC, and the HPC system has much higher resource limits.

Now I can run the pipeline with the -profile CLI option, say which profile I want to use. So I'm gonna use my_laptop, and Nextflow will apply all of the config within that profile scope. So I can try that now. It's the same command as before. Nextflow run hello config, and I do dash profile, single dash 'cause it's the core Nextflow option, dash profile my_laptop.

It's now gonna batch apply all that config option. Oh, and you can see, I said before this might happen that the process requirement, it asked for four CPUs and I've only got two on this Codespaces instance.

So this is a good opportunity just to try out the process resource limits, and say that I only have two CPUs on my_laptop, or in this Codespaces. Now if we run it again, it should cap that, requirement to two and hopefully the pipeline will run. Great.

## 6.2. Create a profile of test parameters

Note that these profiles don't have to only have configuration about their infrastructure. You can have groupings of any config here, including parameters.

So another thing you'll see very often in people's pipelines is a test profile, which includes parameters, which you would normally submit on a per user basis. But here we have, basically different sensible defaults for when I want to run test cases.

And this is great because I don't have to necessarily go and specify all these things, which might be required parameters. Otherwise I can just say dash profile test and it will just run out of the box.

Now something to note is that profiles can also be combined more than one. So I can do profile my_laptop here, and then also add on test. I don't do profile twice. I just do a comma separated list here without spaces. And it's going to apply these profiles in order. So it will take the config from my_laptop profile, and then it will apply the test config on top.

Really convenient and you can see how you can set up lots of sensible default groups here to make it easy to run your pipeline.

## 6.3. Use nextflow config to see the resolved configuration

Hopefully, I've convinced you that Nextflow config resolution is powerful, but I wouldn't blame you if you're going a little bit cross-eyed at this point after I've said about 20 different ways to provide config and give all these different layers like an onion skin.

So if ever you are feeling unsure about what the final resolved config is for Nextflow, know that there's a command called "nextflow config", and we can run that and it will tell us what the resolved configuration is at our current location.

So when I run it here, it finds the "nextflow.config" file in the current working directory, and it processes all the different config, and it gives me the resolved output.

Note that the Nextflow config file can also take the profile CLI option. So if I tell it to resolve in my_laptop and test profiles, and you can see it's also applied, the resource limits here from my_laptop config option and also set the params, which were in the test.

So this is a nice way just to explore how the config resolution is working, if you're at all unsure.

## Wrap up

Okay, that's it. That is Nextflow config in a nutshell. You can do a lot of stuff with config. It's really powerful. But these are most of the common use cases that you'll find yourselves doing, and these concepts apply to all the different options.

Give yourself a pat on the back 'cause this is the end of the Hello Nextflow training course. You're hopefully now confident at both writing your own Nextflow pipeline from scratch, configuring it and running it, and you know all the ins and outs and the things to look out for.

There's one more quiz you can try out on the config training page. So go down and try that out and make sure you've understood all these parts about the config.

And, join us in the last video just for a quick wrap up about some of the next steps that might be good to do after this training course.

Thanks for sticking with us. Well done and I'll see you in the next video.
