# Part 6: Hello Config - Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/IuDO2HeKvXk?si=tnXTi6mRkITY0zW_&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important notes"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../06_hello_config.md).

    The section numbers shown in the transcript are provided for indicative purposes only and may not include all section numbers in the materials.

## Welcome

Hi, welcome to part six of the Hello Nextflow training course.

This chapter is called Hello Config, and it's the final part of our training course.

In this chapter, we're going to be talking about Nextflow configuration. Nextflow configuration is really powerful. It allows us to run the same pipeline on multiple different compute infrastructures with different software provisioning. and different options in the pipeline itself.

This means that you can take Nextflow pipelines built by other people and run them on your system, even though they may have been built for an entirely different infrastructure. This ability to configure Nextflow makes workflows truly portable and shareable.

In this chapter, we'll be using the workflow we've built in previous parts, but we're not going to edit the workflow code at all. We're just going to look at our Nextflow config file and see how changing the config alters the way that Nextflow runs.

Okay, let's get started.

Just as before, let's start off by going to training.nextflow.io. Go over to the left on Hello Nextflow and chapter six. Hello config. I'm now going to go into my GitHub code Spaces environment and double check the script that we'll be using.

## 0. Warmup: Check that Docker is enabled and run the Hello Config workflow

This one's called Hello Config, and it's starting off from where we were before. So looking exactly the same with our three parameters. Greetings for the CSV file, batch for the output batch name and the character for the cowpy name. We have our four imports of the different processes, and then we have a workflow where we chain them together.

I am actually going to close this file now because we're not going to touch the Nextflow file at all in this chapter. We're going to work purely within configuration file. If I look into the Nextflow dot config file that we briefly looked at in the previous chapter five, we can see that we have a single statement here: Docker enabled equals true, which which is telling Nextflow to use Docker when it executes this workflow.

I am using Nextflow dot config in the pipeline root here, which is loaded automatically when I run Nextflow. But remember, Nextflow can load config files from multiple places.

If I check with Nextflow docs go to Configuration, you can see a list of these places and a priority which they load in.

Okay. Let's check that our workflow is executing as we expect it to. Bring up a terminal. Do Nextflow. Run. Hello, config. And hit enter. We should have those four processes running, ending up with a cowpy command. Sure enough, this worked properly. I had Docker enabled it, pulled Docker and ran cowpy for me, just as it did at the end of chapter five.

## 1. Determine what software packaging technology to use

Okay. Let's say I'm running on an HPC and I don't have Docker installed. The best thing to do in this scenario would be to use Singularity or Apptainer. If I was going to do that, I'd go into the module cowpy and change this container to use the singularity image like I showed in the previous chapter, with an oras://, which you can also get from Seqera Containers.

I'd then go to Nextflow dot config set Docker enabled to false and do singularity enabled equals true. Or, if using Apptainer, apptainer enabled equals true and that would work.

Nextflow, it does support other technologies as well apart from containers though, something you might be familiar with is conda. Here we can do conda enabled equals true and set Docker to false. conda doesn't use the same container directive. Instead, we can add a new one here called conda. We then specify the conda package we want to use. It's good practice to be as specific as possible to try and make the pipeline as reproducible as possible. So I'm going to specify the conda channel, conda- forge, and then cowpy, and the exact version, which was 1.1.5.

I could also just write cowpy if I wanted, but that might resolve to a different version of cowpy on different executions of the pipeline.

The nice thing about this is that I haven't touched the docker directive at all. This Docker image is still there. I'm just providing two alternatives now, and these can be switched on or off by using a config file alone.

## 1.3. Run the workflow to verify that it can use Conda

Conda's now enabled, so let's try it out.

Great. It's running and you can see there's a message from Nextflow here saying that Nextflow is creating a conda environment for me, and it's using this cache location.

In the background, Nextflow is running "conda create" commands for me to create a new isolated conda environment with just the packages I want, and then installing and fetching those conda packages so that it can run the process.

You can see it took a little bit of time there because it was creating the environment and installing the software for the first time. However, it's cached this environment, so if I run the same Nextflow command again, it should be a lot quicker because it will reuse the same conda environment.

One of the cool things about this is that these directives can be specified at process level, not just the entire workflow. So if you want to, you can mix and match what technology is used for different processes.

## 2. Allocate compute resources with process directives

The Nextflow configuration file can do a lot more than just software packaging. We can also tell Nextflow how to actually run the steps in the pipeline. One example is telling a host system what resources should be made available to each executing task.

By default, Nextflow doesn't give very much. It gives a single CPU and only two gigabytes of memory to each process.

This is probably something we'd want to change, so that processes which take a long time to run can have more resources and run more quickly, but it can be difficult to know what to allocate to a process. Nextflow has some nice tricks up its sleeve to help you with this.

## 2.1. Run the workflow to generate a resource utilization report

Let's run the workflow again. This time, I'm going to add an additional argument, which is dash with reports. It's a core Nextflow option, so it's a single hyphen. And then whatever file name I like. In this case, I'm going to call it report config one html.

I'm going to run the workflow again. It's going to run exactly as before, but it's going to give me an additional helper report, which you can see has now popped up here in the sidebar.

I am going to right click on this file, click download, which downloads it from GitHub code Spaces to my local system, so that I can easily then view it in the web browser up here.

This report can be generated for any Nextflow run, and it has a lot of information in. It starts off at the top with some metadata about what command was used, when the workflow ran, how long it took, but as you scroll down, we get more detailed information about the resources, which were used by every step in the pipeline.

Because each process runs multiple times for different tasks. We have a box plot showing the variation of the resources that we used for each process.

If I scroll down a bit further, I see similar information about memory used and job duration. Also disk read write.

You can imagine for a large pipeline with long running tasks, this can be very informative about how to fine tune the configuration of the resources which you're requesting so that you don't over request, but also so that you can provide enough that it runs quickly.

If I keep scrolling down the report, we also see a task table, which shows us detailed information about every single task that was run in the workflow. This includes information such as the resolved script, which was run.

Okay, let's go back to our config file. We saw that we really didn't need much for our workflow, so let's tell Nextflow that we only need one gigabyte of memory for every process in the workflow.

Now when we define it like this at process level, this is applied to every single process in the pipeline.

## 2.3. Set resource allocations for an individual process

For the sake of argument, let's pretend that cowpy is really doing a lot of heavy lifting and it needs more resources than the other tasks. We can define an extra block of config here, which applies to just that process by using, with name cowpy.

This is called a config selector, and we can define different patterns here to match different processes. For example, I could do cow star. I then follow that with some curly brackets and let's give it two gigabytes of memory instead of one and let's say two CPUs.

Now Nextflow will be giving every process in the workflow one gigabyte apart from this request, which is more specific. So it overrides it. And just for any processes which are called cowpy, will get two gigs of memory and two CPUs.

Note that Nextflow is clever about resource utilization. So if you start putting these numbers to higher values, you'll see that Nextflow starts to queue job submissions one after another, rather than running all of them in parallel, so that it doesn't over request the resources which are available.

## 2.4. Run the workflow with the modified configuration

Let's try running a workflow again and let's save a new report this time.

Okay, we can download this file and take a look.

Yeah, unsurprisingly, it looks basically exactly the same because this is a dummy workflow, which is not doing anything real. But you can imagine how this iterative approach of defining limits and doing real life workflows with this kind of reporting allows you to do an evidence-based approach to setting appropriate configuration and really making the most of the computational resources that you have available to you.

You can start to be really clever about this. Nextflow has a built-in ability to retry failures, and you can take advantage in your config file by using a closure like this and dynamically setting the resources which are made available. So here I've told Nextflow to multiply that two gigabyte by the retry attempt. So the second retry will get four gigs, the third retry will get six gigs and so on. This is a bit beyond the scope of this training course, but if you're interested, check out the Nextflow docs, which has a nice section about dynamic retry logic.

## 2.5. Add resource limits

Now, one thing you might notice about this is this kind of thing may make it quite easy to accidentally go beyond the resources available on your system. If you request more resources than are available Nextflow will throw an error about your configuration and halt the run. To avoid that, you can use something called resource limits.

Under process scope, in our workflow, we can define resource limits like this, which takes an array, and we can specify the maximum memory CPUs and time which are available on this system.

Setting high values here doesn't increase the amount of resources which are requested. We're still going to be using one gigabyte in our requests, but it means that if any of these requests get to 750, they'll hit that ceiling and nothing more than that will be requested, which means that Nextflow will continue to run and won't crash because of unavailable resources.

So this is a nice safeguard to use, especially if you're using dynamic logic with your resource allocation.

The other situation where this is really useful is if you're using pipelines which are public and not controlled by you. They might come with configuration defaults, and Nextflow will automatically take the right approach of thresholding any resource requests to run on your system.

Okay, great. We've talked about software. We've talked about resource allocation, and we've described different scopes of config, both for all processes and specific processes.

## 3. Use a parameter file to store workflow parameters

Okay, next we're going to turn our attention to parameters. We can define parameters in the config file just as we did before in the Nextflow script. So params dot greeting equals hello or or use params scope and set foo equals bar.

And that's great for setting defaults for your workflow. However, when you're running pipelines, it can be nice to specify parameters in a JSON or a YAML file.

Using a file like this is much better than specifying command line options with dash dash. As when you run a workflow, you might have to specify many parameters and it can be tedious to write them all on a single CLI and error prone. Also, it is unlikely that you'll remember all the parameters that you used, so if you code that into a file, it's easier to launch the workflow again, using the same parameters in the future.

We've got an example file here called test params, and you can see this specifies the three parameters we've got in our workflow with three different values. Personally, I find YAML easier to write than JSON. So just to demonstrate that it works, I'm going to create a new file called Test yaml and copy these in, get rid of the quotes. And hit save.

These JSON and YAML files can be easier to write as they're more familiar syntax. But note that these are only for parameters and they only take key value syntax like this.

## 3.1. Run the workflow using a parameter file

Let's try it out. Do same command as before. Get rid of the report and I'm going to do dash params file test params yaml.

No, this is a core Nextflow option, so it's a single hyphen.

Okay. It ran the workflow and it used the parameters in that YAML file instead of me specifying them all on the command line. Might seem like overkill just for this simple example, but you can imagine if you have 10 or 20 different parameters, it can be a pain to type in manually, and this is just much easier to edit in a code editor and keep hold of for reproducibility sake.

## 3. Determine what executor(s) should be used to do the work

Okay. We've talked about software packaging with Docker and conda. We've talked about process resource requirements with CPUs and memory. And we talked a little bit about how to specify parameters when running workflows.

The final parts of the configuration really is the execution, the underlying compute infrastructure itself, and this is the real jewel in the crown of Nextflow: that we can run these same workflow across multiple different compute infrastructures.

I'm actually going to switch over to the written training material for a second. Under this part of the training, we can see a few different examples of how different executors, in this case, HPC schedulers, define the resource requirements needed to submit a job.

So for Slurm, you have these SBATCH headers, which define dash dash mem and the CPU number. If you're using PBS, you have different headers, and if you use Grid Engine, you have different headers again.

You can imagine it's even more different if you want to run on the cloud, be it AWS batch, Google Cloud, Azure, or more.

Each of these underlying compute infrastructures is called an executor and Nextflow knows how to talk to all of these different executors in order to submit jobs with the correct syntax.

The good news is you don't have to know about this. All you have to do is tell Nextflow, which executor to use.

## 3.1. Targeting a different backend

We go back to our config file and the process we do executor, and I'm going to type local.

Local's actually the default, if you don't specify any other executor, local is what will be used, and that just means your host system, wherever you launched Nextflow,

I could specify instead, Slurm. And that would submit Slurm jobs, or I could say AWS batch, and that would submit jobs to AWS batch.

You need some additional configuration in some cases, for example, running on cloud will need certain credentials, but really this is the core of it, and it can be as simple as one or two lines of config to run your workflow in a completely different compute environment.

Even though we are running on a simple system within code spaces, I can still play around with this a bit and pretend that we're running on Slurm. If I then launch the workflow again, Nextflow run, hello config. It will fail because it won't be able to submit jobs to Slurm. But we can still go into the work directories and see what Nextflow did. So if we go to this work directory and look at Command Run. You can see at the top of this file, we now have these sbatch header lines, which tried to specify the resources needed for the Slurm job.

## 4. Use profiles to select preset configurations

Okay, we're nearly there. Final part of this chapter is talking about configuration profiles. If you're running your pipeline on several different systems, it could be annoying to have all these different Nextflow config files, which you need to specify every time.

Instead, you can encode groupings of configuration within your Nextflow config file, and switch those groups on and off by using a profile flag. Let's see how that looks.

## 4.1. Create profiles for switching between local development and execution on HPC

We're going to create two profiles in our example here, one for my laptop and one for a heavier HPC system. I'm going to cheat a little bit and just copy the code from the training material and put it in here.

We have a new scope called profiles, and then we have a name for each profile, which can be anything. And within that we have configuration, which looks exactly the same as the top level config that we already wrote. So again, we have process scope. Docker scope.

On the profile called my laptop. I'm saying to run using the local executor, so on my host system and to use Docker.

On the university HPC profile here I'm saying to use Slurm to submit jobs, to use conda instead of Docker, and I'm specifying different resource limits, which may match for system size of a nodes on the HPC I'm using.

By default, none of this configuration will be used when I run Nextflow, I have to specify that I want to use one of these profiles.

## 4.2. Run the workflow with a profile

Let's do nextflow run hello config. And I'm going to do dash profile, single hyphen because it's a core Nextflow option. And then the name I gave it, which is my laptop. Nextflow should now use the block of config that was specified within that configuration profile, and apply it when it runs Nextflow. If I wanted to use the other config block, I just have to switch that profile name. Much easier to remember. Much easier to use.

## 4.3. Create a test profile

Note, the profiles can have any kind of configuration, so it doesn't have to be related to your execution environment. For example, let's create a new profile here, which has a set of parameters. We can change this to tux and change to my profile, and now when we do profile test, it's going to specify these parameters, which will overwrite the parameters which are specified at the top level of the workflow.

When you run Nextflow, you can chain multiple profiles and they'll be applied in sequence.

## 4.4. Run the workflow locally with the test profile

So I can take the previous command and do comma test. That will apply the, my laptop config first, and then it will apply the test config. If there's any overlap, then the profile on the right will overwrite any configuration in previous profiles. If I hit enter, let's see what happens.

Okay, we've got a new results file here. You can see the My Profile, which I specified as one of the options. And we can also see cowpy, my profile, and sure enough, there's tux. So that's worked.

## Wrap up

Okay! Amazing. That's it. You've made it to the end of the course. You get a little bit of celebration confetti. Well done for finishing this chapter.

[Next video transcript :octicons-arrow-right-24:](07_next_steps.md)
