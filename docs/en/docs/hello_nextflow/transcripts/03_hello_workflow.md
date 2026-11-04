# Part 3: Hello Workflow - Video Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/_aO56V3iXGI?si=Irl9nAQniDyICp2b&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important notes"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../03_hello_workflow.md).

    The section numbers shown in the transcript are provided for indicative purposes only and may not include all section numbers in the materials.

## Welcome and recap

Hi, and welcome back to part three of Hello Nextflow. This part is called Hello Workflow, and it's in this part of the course where we really start to justify the name pipeline or workflow.

We're going to take our simple, pipeline script so far with its one process, and we're gonna start adding in additional processes and see how Nextflow handles this orchestration and the data flow through the pipeline.

Let's go back to our code spaces. You'll see I've deleted all my .nextflow\* directories and the work directories and everything to try and keep it clean. Don't worry if you still have those files kicking around from previous parts of the course.

We are gonna work from a file called hello-workflow.nf. As before, this basically represents the script that we have built up to this point, and gives us a clean starting place. And again, down in the output we can see that the path is now hello_workflow. So the published files should be going into a different subdirectory in your results folder.

To recap where we are so far, we have a single process in here, with one input greeting, one output greeting file. And then the simple Bash script, which just does an echo command to a file.

We have a single workflow input, the params block here, where we say that it's expecting a path, and the default is data/greetings.csv, which is this file up here.

Then in the workflow itself, we have a main block. We're creating a channel. We are parsing the CSV into rows and then taking the first element of each array, and we're passing that channel into that process, which is then generating three tasks, and we are publishing from the workflow, the outputs from that process.

And then finally, in the output block, we're telling Nextflow to publish these files from this channel to the directory called hello_workflow. And to copy those files rather than soft linking them.

## 1. Add a second step to the workflow

Okay, in this part we're going to add a second process to our workflow. We'll take the outputs of the sayHello process, and process them in a second step, which is gonna convert all the letters within those files convertToUppercase.

This is just a silly example, it's just some simple string processing again, but it shows you how we can take the logic, within the workflow.

We are going to do a bash command called "tr" for this, which is short for translate. It's a Unix command that's been around forever. If you're not familiar with it, I don't blame you. I don't think I ever used it before training, but you can try it out very quickly on the terminal. If I do "echo 'hello world'" and then pipe to 'tr' and then in quotes you say character range, so A to Z, lowercase, and then you want to do A to Z uppercase. And it just says, translate these letters into these letters.

And when I hit enter, you can see it's now capitalized everything. Very nice if you like shouting at people.

So that's a very simple style of bash command that we're going to use in our second process.

## 1.2. Write the uppercasing step as a Nextflow process

So if I go back to my script, I'm going to cheat a little bit and just copy the code from the, from the docs for the training. But you can see exactly what's going on.

We have a new process here. This one we've called convertToUpper, but we could call it whatever we want.

We have a single input path, as we did before. It's not a value channel, it's a path channel. And then a single output.

In the script block we do "cat" on the input file. And we can put this in squiggly brackets if we want to. and which takes that variable. And we run that same bash command in the pipe and we write the results to a file with this file name, and that's picked up by the output path.

We now need to do something with this new process. So we're gonna go down to the workflow where we build up the different logic of a workflow, and after that first process, we're going to run our second process. So convertToUpper is the name of the process here.

It takes an input so we can't just call it by itself. We want to process the output of the first process. So just like we did with this, sayHello out where we're publishing those results. We want to use those same results here as the input, so we can copy those and put them in there.

We want the sayHello process ".out", and Nextflow knows that this means a simple single output record here, which is this file. So that will then be passed as an input to a second process.

## 1.5. Set up the workflow output publishing

Okay. And finally, so that we actually save the results of this second process, we also need to publish them from the workflow, and then define them in the output block, same syntax as before. So we can copy this and say second outputs, or whatever you want to call it.

Take the process name that we're interested in, convertToUpper out, and then down here in the output block. Add this and we could do the same attributes here. So we also want these files in the Hello Workflow subdirectory, and we also want to copy them.

Great. Let's try and run it. So if I bring up the terminal and I do "nextflow run hello-workflow.nf", and we'll see what it does. See if it looks any different to the previous parts.

So it launches Nextflow. In the docs, it says to do this with "-resume", but I deleted all my work directory, so it wouldn't have made any difference here. But if you did, then that will work as well.

And it looks almost exactly the same. But you can see now there's a second line of output here, where you can see the name of the second process we've just added. And sure enough, you can see it ran three times successfully.

Brilliant. If I did have my previous work directories around and I'd done this with "-resume", these would've been, cached just the first step in the pipeline. 'Cause those outputs were exactly the same, so Nextflow would've known to reuse those again.

And so you can see how you can use -resume to iteratively build up your workflow, step by step, if you need to.

Okay, let's have a look in the results directory up here and see if it's worked. We can see we've got some more files up here. We've got our original files like we did before from the first process. And sure enough, we have our upper files and the letters are all uppercase, so it's worked. It is really nice to see.

It's also interesting just to check inside these work directories. As before the, hash here corresponds to the, work directories. So if I look into "ls work", and then expand that, we'll see the different files here.

We see the output file from the first process, which has been pulled in here as the input. And we can see the new output file that was generated.

Now if I do this with "-la" to list and show all files, we'll see a few more things. Firstly, you'll see that this file is actually a soft link to the first process. This is basically always a soft link if it can be, to save file space. We're not publishing the files here and it just references that file from a first task into a second task so that everything is encapsulated within that one working directory, and safe and isolated from everything else.

And that needs to be there because if we look at the .command.sh file, so if I do "cat work/b8/56\*", you can see that the, file parts here are relative, so it's catting that input file, which has been soft linked into the same working directory.

So that's how every work directory will look. When you look at it in Nextflow, you'll have all of the input files there staged into that work directory. And then you'll also have any output files that were created. So that's great. That looks how we expect.

## 2.1. Define the collection command and test it in the terminal

Okay, let's go back to our workflow. What's the next step we want to do?

We have two processes now and they are taking this one CSV file, parsing it and splitting it up. And then we have three tasks for each of these processes and Nextflow handles the parallelization of all of that, so it all runs side by side where possible.

That kind of way of splitting up work to run things in parallel is very common. And the inverse of that is then to gather everything back. So that's what we're gonna do with our final process in the workflow is we'll have a third one here, which takes these three different outputs and combines them all into a single file.

We can do this quite simply in a terminal, just to get a feel for what this will look like.

If I go to the results folder. So, "cd results/hello_workflow/", and we have all the UPPER files here. I can just use "cat", which we use to print the contents of that file, and you can give multiple files to "cat" and it'll read one after another.

So I can say "UPPER-\*", which gives me the same list of three file names with Bash expansion. And I can say combined.txt. I think in the docs, it lists the exact file names, but it's doing the same thing.

Now, if I use "cat combined.txt", we can see that we have the file contents of all three of those files.

So that's basically all that this process is going to do is we're going to try and give it all of the different output files from a previous process in a single process task, and then we're going to "cat" them together and save the output file.

## 2.2. Create a new process to do the collection step

Okay, so let's add our new process. I'm gonna paste this from the training materials, and you can see it's left us a bit of an exercise for the reader here with these question marks. But you can see the general outline of the process is basically what we just did in the terminal, where we're doing "cat" of a bunch of input files and writing it to an output file here called collected, and then the output expects that single path again.

So we need some kind of input here and they're gonna be a set of paths. So again, we define an input path channel and let's call it input_files. Now, this previously has given us a single path here, but a path can also have multiple files here, even though it's still a single declaration.

I'm gonna copy that down here because we want to "cat" these files. And you might think that we have some issues here with printing an array or things like that, but Nextflow is generally pretty kind of sensible when it comes to this. And if it's given a channel with multiple files in it like this, it will, put them all together with space separators. So this will give us the correct syntax.

That's great. So now let's wire up our new process. I go down to the workflow. I'm gonna say combine the outputs, the new process name, and just the same as before. I'm gonna take this previous process, convertToUpper and do ".out".

Great. Let's try it out and see if it works in the terminal. If I just go back up a couple of the directories and then rerun the Nextflow command, and we'll see what happens.

So the workflow has launched and now you can see that we have three different process names, which is great. The first two both look the same as before, and the third new one runs, which is good.

However, there's something a bit odd here. We wanted to combine those output files into a single file, and yet this process we can see has run three times, not once.

Sure enough, if we go into one of these work directories. And do "cat work/" "collected", then we'll see. There's only a single word in here, not three.

And so what's happened is that Nextflow has continued that parallelization just as it did in the previous steps. And this process gave us a channel with three elements, and those three channel elements were passed to that our downstream process, which generated three process tasks.

It basically tried to collect three separate times and each time it just had a single file, so it just did cat single file to an output, and in fact, we can see that in the .command.sh file as well.

If I do .command.sh, we can see it just has a single file name here and only a single file was staged into that working directory.

## 2.3. Add the collection step to the workflow

So somehow we need to tell Nextflow to bring all of those outputs together from a previous process and give them to this downstream process as a single channel element, rather than three.

We do that with a channel operator called _collect_.

This is a super useful operator, which you'll see in Nextflow pipelines all the time. This is a channel here, this output channel, just the same as the one we created up at the top. And so we can append channel operators to it just like we did before. We can just do dot, and then in this case, collect, brackets.

And that's all we need. That's gonna then manipulate this channel before it's passed into this process.

If you wanna see what's happening to it, we can also view it here. So here, this is not related to running this process at all, so I could put it at any point after running that process. But we take the same, output channel, and we're looking at it with .view, and then we're looking at it again with .collect.view.

And when we run this, it will show us the two different structures of that channel, before and after collect. So let's try that now. Okay, I've just zoomed out a little bit because some of the outputs are quite long, but if I run the pipeline, we'll see if it works.

I'm hoping a third process will run just once, because it's collecting the outputs and sure enough, you can see collectGreetings as one of one. So that's run just one task.

And then if we look at the view statements, we have three view statements for the three elements of before, with one file path in each one.

And then after that collect statement, that's just triggered once because there's a single element in that channel. And now we have this, list of three different file paths.

That's exactly what we hoped for. And you can see hopefully, this is basically the inverse of that "map" operator that we did to go from the CSV arrays into separate channel elements. Now we're taking separate channel elements and putting back into a single array.

Great, we can clear up these view statements. We don't need these anymore. We can go onto the next step.

Before I go any further, and before I forget, I'm gonna add a new publish statement here. Third output. You can call this something more semantic and descriptive in your workflow. And then I'm gonna add that down to the output block again and say path 'hello_workflow' mode 'copy'. Just so that the output file generated by this process is saved to our results folder up here.

Just to quickly double check that works. Should be a bit cleaner now 'cause we don't have those view statements. And, we'll see if we get our new output file up here. One of, one task ran, got a new file called collected, and now we have all three of those words. Fantastic. What's next?

## 3. Pass additional parameters to a process

Okay. Next we're going to look at handling multiple inputs into a single process. So far you can see that all of our processes are just taking one thing as an input. They all have a single line under their input.

We're going to demonstrate this by allowing Nextflow specify a different batch identifier so that maybe you run this, workflow multiple times and you can give it a different batch ID each time.

I am gonna simply add a second line in the input here for collectGreetings. And I'm gonna call it "val" ,'cause this is a string. Now it's a value, not a path, and I'm gonna call it "batch_name".

Then I'm gonna edit the script down here to use this variable, and I'm gonna try and put it in the same place as the training material. So I put it in the middle of this file path COLLECTED-$\{batch_name\}-output.

Not quite done yet. Remember that we have to tell Nextflow what the output file names are going to be. So we have to also do the same thing up here: COLLECTED-$\{batch_name\}-output.txt".

Fantastic. Nextflow now is getting a second variable input and it's interpolating that into the script and the output.

One last thing, we now have to find where this is being called, and we have to pass the second input to the process. This is just like any other input into a function in any other language.

Just like we did earlier in the training, I'm going to use the special "params" here, and we're gonna call it "params.batch" so that we can have a -- batch CLI option. And now you can see that our process here has two separate inputs just comma separated, which are being passed in.

It's really important to get the order right, so the order of arguments here for channel and then the param must match. The, channel and the batch name there. This is just positional matching.

Okay. I can run this pipeline now straight away with --batch, but let's first do the right thing and define it in the input here in Params. So I'm gonna add it to batch and then we're gonna say it is a string and let's give it a default. So let's just call it batch. Okay? Now let's try running the workflow.

--batch Trio. I think it says in the training material, but we could use any string we want there. And hopefully we will see that results output file come up here.

And sure enough, COLLECTED-trio-output - that has worked properly. It's renamed our file. And you can imagine now this is useful 'cause if I run that again with a different batch name, like replicate_two, then it's gonna give us a different batch name up here.

And and it won't then clobber the output files in this case. So that's nice.

## 4. Add an output to the collector step

Okay, so we've now got multiple inputs to our process here. But what happens if we want to create multiple outputs? Our example here then is we're going to create a report for this process, just saying this is how many files were collected.

And we'll do that with an echo command here. So we can say echo. There were, I'm gonna copy this from a training material, so you don't have to watch me type it out.

There were $\{count_greetings\} greetings in this batch, and save that to a new file now called $\{batch_name\}, so same variable, we can reuse that as many times as we want, report.txt.

## 4.1.1. Count the number of greetings collected

We need to actually calculate that somehow. We could do that logic in the Bash script if we wanted to, using Bash logic. However, we can also just do scripting directly within the Nextflow code, as long as it's within the script block in the process and above the quoted section.

Anything here won't be included in the final rendered script, and it will just be executed by Nextflow when it renders a task.

So here we're just doing some logic. We're creating a new variable called count_greetings. We take the input files channel here, and we're calling .size() on it.

Okay, that function is going to give me a number here into this variable, and now our warning's gone because this variable is being defined.

Okay, so we're creating that second file in the work directory, but we need to tell Nextflow to expect it as a published output of this process. So we do that by exactly the same syntax as we did for the first file.

We say path 'cause it's, again, we could be publishing a variable here if we wanted to with "val", but we're gonna say "path". And then the expected, file name. Notice it's not highlighted here. That's 'cause I used single quotes. I have to use double quotes.

## 4.1.2. Emit the report file and name outputs

Okay, that's great. And we could now start to out access these outputs down here just as I did here. But now it's an array of different objects, so I could do collectGreetings.out[0] to get the first, or one to get the second, which is our new report.

But I don't really like doing that very much 'cause it's quite easy to mess up the index counting. And you sit there counting lines a lot and you add in a new output and suddenly everything breaks. So

it's much nicer to reference everything by name instead. And we can do that with a special key here called "emit".

So we can call this whatever we want. Let's call it emit outfile, and emit reports. If you define these and you can do it on one or many, it's up to you. Now I can go down here and instead I can go dot out dot reports and just call it by name, which is much easier to understand your code when you read it, and it's safer to changes in the code.

I've, added the .out.report here, but actually I need to have two different outputs being published. So I'm gonna rename as something more interesting like collected and report and is that what I called it? I called it out file, sorry. So that emit name here outfile and report. 'cause we're publishing two different output channels and so we need to reference both of them in the publish block.

Then we also need to define these in the output block. So I renamed that collected, and again, for reports, little bit verbose here, but it's really useful when you come in to read a new workflow, to see all of the different outputs here, all the different channels listed side by side, and there are ways to make this less verbose, which we'll touch on later.

Okay, let's try it out and run our workflow and see what happens.

Hopefully now it should run basically the same as it did before. And we're gonna get a new output file up here called replicate_two, report. And there we go. It's opened and it says there are three greetings in the batch, which is what we expected, so it's perfect.

If I go into the work directory here just to prove to use was executed in the Nextflow, code rather than the bash script, I can go to cat work/ command.sh, and you'll see here that it's just echoing this string directly. There were three greetings in this batch, and so that variable was interpolated by Nextflow. It was calculated in the script block before it wrote the .command.sh file. So the resulting variable calculation is basically hard coded into this before it's executed on your compute environment in this case.

And so you can see that separation between the script. Block here and anything above it. I \_hope that makes sense.

## Takeaway and quiz

Okay, that's the end of this part of Hello Nextflow. So as before, go and check out the quiz. Do it on the webpage or in the CLI, go through some of the questions and just check you've understood some of the material we've covered. See if there's anything there which highlights anything you might not have understood. Not too many questions. Nice and easy to do. Or you can do it on the webpage down here as well.

And have a little break, little walk around and come back and join us in part four of Hello, Nextflow, where we'll talk about modules. Thanks very much.
