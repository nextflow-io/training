# Part 3: Hello Workflow - Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/zJP7cUYPEbA?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important notes"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../03_hello_workflow.md).

    The section numbers shown in the transcript are provided for indicative purposes only and may not include all section numbers in the materials.

## Welcome

Hi, welcome to part three of the "Hello Nextflow" training course.

This chapter is called "Hello Workflow".

In chapter two, we built a simple workflow of one process, but in reality, pipelines are useful because they can chain multiple steps of analysis together.

In this chapter, we're going to take that initial example and extend it to be a little bit more realistic.

We're going to add some additional steps and we're going to look at how we use channels to connect those steps.

We are going to look at multiple tasks, which can collapse into a single process and we're going to look at processes which can have multiple inputs and multiple outputs.

Okay, let's get started.

So let's start off. Same as before. Let's go to training.nextflow.io. Hello Nextflow, chapter three. Hello Workflow. And let's open up our workspace. I've cleaned up all my work files from my previous chapters and I'm going to open Hello Workflow.

Now this is the same file we've been working on until now so this should look familiar. We've got our say hello process. We've got our params.greeting with its greetings CSV file, and we've got our workflow at the bottom, which loads that CSV file, creates the channel and passes it to our process.

## 0. Warmup: Run hello-workflow.nf

If you like, we can try this out and double check it is working as we expect. Load up a terminal to nextflow run hello workflow nf and click enter.

Okay, great. Our three processes run. We have our results directory with our three outputs. Bonjour. Hello. Holà. So let's close those files, close the terminal, go back to the script.

## 1. Add a second step to the workflow

Okay. For our example, we're staying basic and we're trying to stay domain agnostic. So our second process is just going to manipulate these strings, these words, in a simple way. We're going to use the translate Unix command to take these files and make them all uppercase. We do that with the "tr" command.

## 1.1. Define the uppercasing command and test it in the terminal

We can try this just in the bash terminal, and see if it works. So you do echo, Hello World, and then pass that on with the pipe character to tr, and we give it a recognition pattern, a to z and what it should translate to. A to Z in uppercase.

This is very simple because it's literally doing the A to Z characters. So it won't work on anything that's accented or anything like that. But for the purposes of the example, you should get the picture.

Going to hit enter and it prints to a terminal, HELLO WORLD in capitals. And just as before, we could redirect this to a file if we wanted to. Outfile.

Okay. Let's clean this up.

## 1.1. Write the uppercasing step as a Nextflow process

Let's go back to our script and write a new process to handle this bash command. I'm going to copy the previous process, paste it underneath, and call it convert to upper. For uppercase. I am going to use the same publishDir results, but I'm going to make a few changes here. Instead of taking a val, I'm going to take a path input file, and I'm going to have a prefix here upper, so that our output files don't clobber the output. And I'm going to use the variable name from the input.And then I'm going to change a script down here, and instead I'm going to use cat on the input file and just like we did in Bash TR, a-z, upper input file .txt. Okay, let's click save.

## 1.2. Add a call to the new process in the workflow block

Now if I scroll down, we need to actually call this process. Just adding the process into a script isn't enough. We have to tell Nextflow that we need to run this process and where to do that.

So I'm going to right here, convert to upper and

okay, we're getting an error here saying it expects an argument. Sure enough, we need to pass something to this process so that it actually has something to do.

## 1.3. Pass the output of the first process to the second process

What we're going to do is we're going to take the output from this process. So I take the name, say hello, and when I do dot out.

For a simple example like this, where we have a process that has just one output and we're passing that to a new process, so it has one input that should be all we need. So I'm going to click save, bring up the terminal, and let's try and run this again.

## 1.4. Run the workflow again

Now, I haven't cleared up my work directory from the last time I run this workflow. I'm going to run it again and I'm going to use this as an opportunity to show how partial caching works. So if I do single dash resume. Hopefully it should reuse the outputs from that first process, which were exactly the same as last time I ran. But now we have a new process here which hasn't run before, that runs from scratch. And sure enough, you can see the first process used, the cache outputs, and the second output ran three of three. You can also see that we have both of our processes here now, our first process, say hello, run three times, and our second process convert to upper run three times.

If I run this again, as a reminder, with -ansi-log false, we should see that six different process tasks run three for each of them. So this is doing exactly what we hoped it would. The first process is running three times, passing those outputs onto a second process, which is then running three times.

So let's have a look inside the work directory and see how Nextflow is handling these file inputs. If I take this hash directory here from the second process we can use a tree command again with -a just to look at these files. You can see in here that we have our input file, which is the Bonjour-output.txt file, and that's actually a symlink. That's what this arrow is showing us, and it's pointing to the file in the previous work directory.

This makes sense. Nextflow handles the execution of each task in its own encapsulated directory, so it's completely self enclosed. However, it needs to provide the files from a previous steps as an input. Rather than reaching outside of the work directory to get those files, Nextflow stages them into the work directory.

If we have a shared file system like here, it does that using a symlink so that it doesn't use any additional file space. If we use cloud storage with buckets in different locations, it would fetch those files and actually copy them into the work directory.

Let's have a look at the command sh file. If I do code work, command sh, you can see, sure enough, it's accessing that file from the local directory. So everything is very self-contained and clean.

We can also check the results directory and make sure that these files were outputted properly. And sure enough, in results, we can see all the output files from the first process and all the output files from the second. And they're all in uppercase as we hoped.

This is where the power of Nextflow starts to shine. With some very minimal code and Nextflow handled execution in parallel of these tasks with clean encapsulation within separate work directories and staging input and output files and file publishing all automatically for us just out of the box. So you can see how, as we scale this complexity of our analysis workflows, this functionality is really, really valuable.

## 2. Add a third step to collect all the greetings

Okay. These steps were one-to-one. We had one output from the first process going to one input for the second process. Next, we're going to talk about how to collect these different outputs into a single process task, which is again, a very common thing to do. So let's quickly bring up the terminal and do a dry run of this.

## 2.1. Define the collection command and test it in the terminal

I'm going to cheat and copy the example bash code from a training material and just hit enter.

What we can see here is we ran this echo command three times to three different output files, which I can see here. And then used the cat command to print the output of each of these three different files, and redirect that to a single collected file.

And if I do "cat COLLECTED-output", you can see it's got the contents of those three different files, now in a single file.

## 2.2. Create a new process to do the collection step

So let's see if we can replicate the same thing within our Nextflow pipeline.

Let's scroll up and create a third process. I'm going to copy this previous one, and this time I'm going to call it Collect Greetings.

In the bash terminal, we called it collected output txt. So I'm going to say the same path output here. And I'm going to do the redirection here, so it's saved in the same way.

Okay. We need to change what happens at the start of that command, and we need to think about what the input file is here. In fact, this process is going to take multiple input files. I'm going to keep path and I'm going to change this to a new variable called input files, plural.

I'm then going to again, cat them like we did in our bash script. And I'm going to use the variable here.

Now, you might think this wouldn't work. We've seen previously failures where an array of strings or an array of paths has been passed to a process and that caused an error. But in fact, here Nextflow is going to handle this automatically for us in the right way. It's going to take several different input files, and it's just going to print the different file paths here.

Of course it helps that the cat command can take a series of file names like this. If I was using a different command that required an argument before each file path or something, we'd have to have a bit more code here and logic to be able to handle the iteration of these file paths. But in this case, it should just work.

## 2.3. Add the collection step to the workflow

Okay, let's go down to the workflow and add in our new process. Collect greetings. And again, let's take the output from convert to upper out. Let's save this.

Give it a try. nextflow run hello workflow.

Okay, the workflow ran, but something's a bit strange here. We've got three executions of the first step, which we expect. Three tasks for the second, but we also have three tasks at the end when we expected to only have a single task here merging all the outputs.

If we go into our results directory. We also see the collected output only has a single value rather than all three. This is because that output file was overwritten three times with three different values.

This makes sense because we passed one output to one input here in the same way as we did in the previous step.

## 2.4. Use an operator to collect the greetings into a single input

So we need an operator here to take this channel with three elements and collapse them to a single element, so that that final process only runs once.

To do that, we're going to use the collect operator. I can do this directly within the workflow. I can do .out and chain onto an operator here at the end .collect.

Hit save. And then for the purposes of this training, I'm also going to do some view operators like we did before, so we can take a look at this channel before and after we use the collect operator, so we can understand what's happening.

I am going to take this channel, get rid of the collect and dot view greetings, and then I'm going to duplicate this line, add in the collect operator. And change that to after.

This is separate to where we're calling this, but that's fine because we're using the same operator calls on the same output channel.

Okay, let's hit save and let's try it out in the terminal. Going to lgoing nextflow run. Hello, workflow. Rerun our script.

Okay. This is looking better. As before we can see the first two processes run three times and now our final process only ran once.

If we look at what was printed by the view operator, down here, we said before collect, which is this output here, and that's printed three times. And you can see there's a single path for each one of those. And then after collect, you can see that we have this array of three paths. So that's as we expect.

Okay, let's check the results file and see if it's what we expect this time. Sure enough, there are now three lines in the file - that's successfully concatenated these three outputs into a single output file. Fantastic.

Okay, I am going to clean up and let's go on to the next step. And I'm going to delete these view statements just to keep things clean.

## 3. Pass more than one input to a process in order to name the final output file uniquely

Okay. So far, all of our processes have only taken a single input. We're now going to do an exercise where we add more than one input to a process to see how this works. To do this, we're going to use this collect greetings example.

Each time I ran the workflow, it overwrote that file in the results directory, which may not be what we want.

## 3.1. Modify the collector process to accept a user-defined name for the output file

So for this example, we're going to pass an additional parameter so that we can customize the output file name.

Adding a second input to a process is very simple. I just add a second line in the input block. This time it's going to be a value, rather than a path, because we want to pass a string and I'm going to call it batch underscore name.

I can now use this variable in the script block, and I'm going to say collected dash dollar batch name.

I'm using squiggly brackets here around the variable name. That's just to keep it separate from the rest of a string, and it probably isn't needed in this case, but I think it makes it easier to read.

Okay. Finally, remember to update the output path because now the file name has changed, so I'm going to do the same thing and put the batch name into the output of path as expected.

## 3.2. Add a batch command-line parameter

We now need to pass in a batch name from somewhere, and I'm going to create a second parameter to do this so that we can do it on the command line when we run the workflow.

So I'm going to do params batch name, and by default, let's call this test batch. Now I can use this special parameters variable down, where we call the process.

And sure enough VS Code is telling us that there's not enough arguments to this process now, and that it expects a second input.

Simply do comma and pass our new variable and the error goes away.

Note that the order of inputs here is really important. The first process input was the path, and the second input is the name. If I change the order here, I must also change the order when I call the process. Otherwise. Next, we will pass the wrong channel to the wrong input.

## 3.3. Run the workflow

Okay, let's try it and see if it works. Let's do "nextflow run hello- workflow. Okay, it ran as before. Let's have a look in the results directory.

Sure enough, our file name here is now called " collected test batch output txt". Fantastic.

And now let's see if we can overwrite that by running again. This time I'm going to do --batch_name to match that special parameter variable name here. And I'm going to call it demo output.

Run the workflow again and we'll see if something happens.

Okay, we now we have a collected demo output .txt. And because this file name is different to that one, it hasn't overwritten it. Both are now present in the results directory.

## 4. Add an output to the collector step

Okay, so there we showed giving multiple inputs to a process, but how about multiple outputs? For this example, we're going to calculate the number of greetings which are processed and output that as a secondary output for this collect greeting step.

## 4.1. Modify the process to count and output the number of greetings

We're going to do a bit of a trick here. Nextflow processes have this script block with a multi-line string, and that is passed as bash output to the dot command dot sh. But we can actually write any custom code above that, and that will be executed as part of a task but not included within the bash script.

One of the built-in functions in the Nextflow syntax is called size. So I'm going to take the path input, and I'm going to say count underscore greetings, just to define a variable name. I'm going to take the input files and I'm going to call "size" on it.

This function will count the size of this input channel and assign it to a variable.

We can now return that variable as part of the output block. So we say, val, because it is value, not a file. And count greetings.

Now this is enough by itself, and we could now access these different outputs from this process. However, we'd have to access them in a positional manner. So using an index key such as zero and one.

To make it a little bit easier to get at the outputs, we can name them and we do that by using an emit statement.

So we do comma emit out file or whatever I want to call this. And I do here emit count. This is basically just a decorator, which just helps us write slightly cleaner code so that we can easily reference the specific outputs later in the workflow block.

## 4.2. Report the output at the end of the workflow

Okay. If I scroll down to the workflow block, I can now take the outputs of collect greetings, do collect greetings, dot out, and we can see our two named outputs are suggested here by the VS Code extension. Very handy.

So I'm going to do dot count to get the count value that we've just created, and I'm going to do view, so that it prints it in the command line. So we can see it when we run the workflow.

Let's write something in the closure here just to make it a bit nicer. num greetings, there were greetings greetings.

And we don't actually care about the other output because we're not using that as an input for any other processes. But you can see how we could easily pass this as an input to another process if we wanted to, downstream.

## 4.3. Run the workflow

We're going to click save. Let's have a look at the terminal and try it out.

Okay, fantastic. Here we go. There are three greetings. That's exactly right.

Okay, great stuff. That's the end of this chapter. We're all done for making it this far. You're now starting to build up quite a realistic workflow, where we're able to handle inputs and outputs and logic within our workflow.

As these workflow files get longer, they start to become a little bit unwieldy. So in the next chapter, we'll look into how we can modularize Nextflow code into separate files so that it's easier to find and maintain the code within the workflow.

Join us in the next video for chapter four. Hello Modules.

[Next video transcript :octicons-arrow-right-24:](04_hello_modules.md)
