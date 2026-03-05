# Part 2: Hello Channels - Video Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/yDR66fzAMOg?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important notes"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../02_hello_channels.md).

    The section numbers shown in the transcript are provided for indicative purposes only and may not include all section numbers in the materials.

## Welcome

​Hello and welcome back to Part 2 of Hello Nextflow. This chapter is called Hello Channels.

Channels are like the glue in your Nextflow pipeline. They're the bits which hold all the different processes together, which Nextflow uses to pass all the information around and orchestrate your workflow.

There's another part to channels which are operators. These are basically functions that we can use on channels to modify the contents. Let's dive into VS code and see where we are.

I'm very zoomed in on this VS code, so to keep things clean and tidy, I have removed all the _.nextflow\*_ files and the* work/* directory and the results/ and everything from Chapter One. And I'm just starting fresh here. But don't worry too much about that. If you don't want to, you can leave those files around. They won't cause any problems.

We're gonna start off working on _hello-channels.nf_ for this chapter, and if I open this up, it should look very similar to the file we were working on previously. It might be that different parts are in different parts of the script, but everything should be basically the same.

One thing that is different is that the path in output block here is now _hello_channels_ for this part, which means that the result files will be stored in a different subdirectory in your results if you still have that there. So it should be nice and clean a place to start without being confused about outputs.

Okay, so let's quickly remember what this script does when we run this workflow. We do _"nextflow run hello-channels.nf"_. We can do _"--input myinput"_, and when we run this, it's going to use this parameter, params.input, which was passed as the variable for the sayHello process up here, which goes into greeting and gets saved to output.txt. And we can see that in the results file. Great.

## 1. Provide variable inputs via a channel explicitly

That's nice. But it is, it's quite simplistic. We have one variable in this parameter, which goes into a process which runs once, and doesn't really scale. And we can't give it lots of different files to create here. We can't give it lots of different greetings. We have just one.

In reality, Nextflow is all about scaling up your analysis. So you probably want it to do more than one thing. And we do that with _channels_.

Channels are a bit of a unique concept to many people picking up Nextflow. It comes from this kind of concepts of functional programming, and it can take a little bit of time to get your head around, but once you click, they really unlock the power of Nextflow and it's key to how you write your workflows.

## 1.1. Create an input channel

Let's start by taking this script and making it use a _channel_ instead of just a _param_.

We go down to the workflow, which is where all of our workflow logic is about stringing things together. And I'm gonna go in here and I'm gonna create a new channel.

Create a new channel.

And I'm gonna call it "_greeting_ch"_. This is convention to do "_\_ch"_ like this, just so that you can remember that this variable is a channel. But you can call it whatever you want.

And then I'm gonna say equals, and I'm gonna do _"channel.of"._

Channel is like the name space for everything to do with channels. Lower case "c" if you've been using Nextflow before. And the _".of"_ is something called a Channel factory, which is basically a way to create a channel.

There's lots of different channel factories. If I do just "." here, you can see that VS Code is suggesting loads of them, but _".of"_ is the simplest and just takes, an input here.

So I can do some brackets and I'm gonna say _"Hello Channels!"_.

Great. I have a channel. Fantastic. I can hit save, I could run it again, but nothing interesting is going to happen. VS Code has given me a orange warning line under here and told me that this is set up: you've created this, but you've never actually used it for anything. This channel's not being consumed.

Okay, so how do we use it? Very simple. I'm gonna take this, copy it, and I'm gonna delete _params.input_ and I'm gonna put in _"greeting_ch"_ here instead. So we're gonna pass this channel as the input to sayHello.

Note that I've hard coded this string for now. This is a bit of a backward step after our nice pram that we used at the end of the last chapter, but it just keeps things simple to start off with so you can see the logic.

Okay, I'm gonna go into my terminal and I'm gonna run the workflow again. Without any _"--input"_ this time, and it's gonna run and it's gonna use that channel that we've created and hopefully we should have a file up here in _results/hello_channels/_ and it now says "Hello* Channels!"*. Fantastic. So that is what we're hoping from, our channel here. Great.

## 1.4. Use view() to inspect the channel contents

One more thing to add on here, just a quick introduction to another function we can use on channels called "_.view"_.

This is analogous to the _print_ command in Python or other languages that you might be used to , and it just dumps the contents of this channel to the terminal when we run it.

So do "_.view"_, and then if I rerun the workflow again, it should print to the terminal what the contents of that channel is, at the time that we created it.

Sure enough, you can see it's printed to the terminal here. _"Hello Channels!"_.

Note that you can break these things across lines if you want to, and in fact, the Nextflow, automatic formatter will try and do that for you. White space is not really important here, so you can chain these things one after another.

## 2. Modify the workflow to run on multiple input values

Okay, so our channel has one thing in which is nice, but it's basically the same as it was before. So let's make it a bit more complicated. Let's add a few more things into our channel.

The "_.of()"_ channel factory can take multiple items, so let's write a few more. We'll do _Hello, Bonjour, Hej_. And then we can run this workflow again and we'll see what happens.

Should run again. And we've printed now. _"Hello", "Bonjour"_ and _"Hej"_ to the terminal with our view statement. Fantastic.

## 2.1.2. Run the command and look at the log output

You might think that we're done at this point. But actually there's a bit of a gotcha here, which is gonna trip us up. If we look at our output file here. You can see it's got _"Hello"_ in, but it doesn't have any of the other outputs. In fact, it's just this one.

If we run this workflow multiple times, we might even see that sometimes it has _"Bonjour"_, sometimes it has _"Hej"_. It's a bit random.

If we look at the terminal, we can see it ran three times and we can see the different view outputs. But if I go to the work directory, I can do _"cat work"_. Put this hash and expand that and _output.txt_. You can see that this file in the work directory is different than the results directory, and this one is, _"Hej"._ So there's something not quite working right here.

And the key is that, we have three tasks that ran. The Nextflow output tries to summarize that as the processing goes on, so that it doesn't completely take over your entire terminal, and that ANSI Logging uses ANSI escape codes, has basically overwritten the other tasks. So it just shows you the last one which happened to be updated.

## 2.1.3. Run the command again with the -ansi-log false option

There's a few things we can do to actually understand this a bit better. We can look into the work directory itself and you can see all the different work dirs there, but that's a bit confusing 'cause it will be mixed up with different Nextflow execution runs.

Or we can tell Nextflow not to use the ANSI escape codes.

So if I run the command again, but this time I say _"-ansi-log false"_ to turn it off, I could also use the environment variables _$NO_COLOR_ or _"$NXF_ANSI_LOG=false"_. Then it uses the kind of more old fashioned style of Nextflow logging without any of these escape codes. It just prints directly to a terminal with no clever updates.

And now we can see all three of these processes that ran. And each one of them its own task hash. And if we go into these work directories, we'll see the three different greetings that we specified.

So that makes a bit more sense now. Hopefully you understand that Nextflow was doing this, it was just being a bit clever with what it showed you in the terminal with those work directories.

However, this is fixed for one problem with the work directories, but it hasn't fixed a problem with the output file. We still just have one output file which says _"Hello"_.

## 2.2. Ensure the output file names will be unique

Now to understand this, we need to go back to our workflow script. We're generating our channel here, we're passing it to our process, and if we look at the process, we are writing the greeting to a file called _"output.txt"_ and passing that output file back out to the output block down here, publishing it.

However, each three times this process runs these three different tasks. They all generate a file called _"output.txt"_, that all of those output files are published into the results directory, and they all overwrite one another. So whatever result file you get there is just the last one that was generated, but clobbered all the others. That's not really what we want.

## 2.2.1. Construct a dynamic output file name

There are different ways to handle this, but the simplest for now is just to create different unique file names. So each time the task runs with a different greeting, it will generate a different output file, which will no longer clash when published. And then we'll get three unique output files.

We do this in exactly the same way. We can use this variable anywhere within the script block and we can use it multiple times.

So I can paste it here, _"$\{greeting\}\_output.txt"_, and then I also need to paste it up here because we're no longer creating a file called _output.txt_. So if I don't update this, Nextflow will crash with an error saying it expected a file, which was never generated.

So I need to do the same there and I need to use double quotes, not single quotes, so that this variable is understood.

Okay, let's try it out and see if it worked. We are gonna run the workflow again. Hopefully it'll show us the three different tasks within the three different work directories. And sure enough, you can see up in the results folder up here on the left. We now have three different files with three different file names and each with the different contents that we expect. So the files are no longer clobbered one another, and everything is there as we expect.

This is a bit of a kind of trivial setup that we've gone through here, but it underscores some of the key concepts you need to understand about how file publishing works, and some of the things that you might fall into as traps. So hopefully you can avoid that in your own workflows.

It is worth noting also that what we've done here is a bit impractical in real life situations. We've taken some input data and we're using that data, but we're also naming the file after that data, which you can't usually do.

So in real more mature Nextflow pipelines, you will often pass around a meta object with all the metadata associated with a given sample. You can then create dynamic file names based on that, which is a lot more practical.

If you're interested in how to do this with best practices, there's a side quest on _training.nextflow.io_, which is all about specifically metadata and meta maps, so you can dig in there for more detail.

## 3. Provide multiple inputs via an array

Okay. Next we're gonna explore a little bit about how channels are structured and how they differ to other kinds of data structures in the coding language. And I'm gonna think a little bit about how I could potentially use an array, which might be a familiar concept if you've come from other languages.

Can I use an array in a channel? Let's try it. I'm going to create an array, and I've copied this from the docs, _"greetings_array"_ and _"Hello", "Bonjour" \_and_ "Holà"_. And then I'm gonna put that here instead of my hardcoded strings. So I'm gonna say "channel.of" _"greetings_array",\_ passing this array into a channel. Let's try it.

Bring up the terminal, and run the pipeline.

Okay. You can see that the view statement here did print our array as expected, but then all of this red text, or it won't be red if you have still have _"-ansi-log"_ off, but all of this red text is telling us that something went wrong.

We don't have a nice green tick here anymore. We have a red cross, and if I just make this a little bit wider so it's easier to read, Nextflow is telling us what went wrong.

So let's break this down section by section. It says the error was caused by, and then the reason for the error, which is missing output files. So basically that output block said that this file should be created and it wasn't. Next it says this is the command which was executed. So this is basically the contents of that _.command.sh_ file. This is what it looked like after all those variables had been put in.

And you can see here our echo command is actually only been run once and it's used the entire array, but in a string representation, which is not really what we wanted.

And then the command exited like that, and that was the work directory where we can go and see the files to understand a bit more.

Okay. So what happened then was. Nextflow just passed this entire array as a single channel element to the process, which meant that the process only ran once. It had one task and it didn't use the data in a structure we expected.

## 3.2. Use an operator to transform channel contents

So we need to do something to this channel first, before it can be used. And this is setting a stage for using operators, which is special functions we can use on channels to manipulate channel contents.

In this case, we are going to use something called _flatten_. Which we pass on the end of the channel here. So we create the channel and then we run _flatten_. And again, if we hover over it, it shows us the documentation for this command straightaway in VS Code, which is very helpful. You can also find all of these docs on the Nextflow website, the documentation.

I could just run this code now and see if it works, but it's also a nice opportunity to introduce how to do dynamic code within operators and within Nextflow code, which are called closures.

So I'm gonna add back in a view command here before we run _flatten_. And here this one has got this squiggly brackets, which is the, dynamic closure. And there's just some arbitrary code within here which will be executed, within the context of a view operator.

Here, this is saying take the greeting, which is the inputs of the view operator, and that's here. I could call this whatever I wanted to, I could call this _"foo"_ and I just need to refer to it as _"foo"_ later. And then I say with this, return this.

And then set returning a string which says before the flatten for a variable. very simple.

I'm now gonna add another one of these exactly the same, but I'm gonna say after _flatten_.

So what this does, because this runs in sequence, you're gonna see what the channel looks like before we run _flatten_, and then again after we run _flatten_.

And then this greeting channel is still created, so it's still gonna be passed into the process. And hopefully now the workflow will run. Let's try it out.

Great. So first things first is that the pipeline didn't crash this time. We had three processes that ran properly and we've got a little tick mark. And then we can see our view statements did work.

We have before _flatten_, which is that array that we saw before from the failure, and then we have three times the after _flatten_ was called where we have _"Hello", "Bonjour",_ and all that other three separate elements in the array, which are now as we hoped, three separate elements in the channel.

And you can see that the _view_ operator was run three times. And that's because this channel after _flatten_ now has three elements. And so the operator gets called three times.

Very quickly, I would just mention that when I was creating channel factories before, I did _"."_, and then we saw there were lots of different ways to create channels, and one of them is called "_fromList"_. And that's actually specifically designed to do this same operation. So we could have just done from list greetings away, and that will work. It's slightly clean and nicer syntax. But for the purposes of this demonstration, we wanted to make it a bit more step-by-step so you could see how the channel is being manipulated and how different operators can change what's in the content of a channel.

## 4. Read input values from a CSV file

Okay, how can we make this a bit more realistic? You're probably not going to want to be creating lots of code in your Nextflow pipeline with hard coded arrays. You're probably going to want to take the data from outside when you launch, and that data is almost certainly going to be in files.

So the next thing we're gonna do is we're gonna replicate this, but instead of taking the data from a single CLI parameter or from a hardcoded string or array, we're gonna take it from a file.

So let's get rid of our greetings away. And now we're gonna change this channel factory again. I just said there were a bunch to choose from and there's one called _".fromPath"_. And I'm going to tell it to, in this case, take _params.input_, which is going back to our input that we were using earlier.

Now that parameter isn't really ready to be used yet. We're still saying that it's a string and it's hard coded here with a default, but we could overwrite that string. We now want this to be a file instead. So the type is different. It's no longer a _String_. It's a _Path_.

And then we can set the default if we want to, again, to a Path. And if I look in explore on the left, you can see in this repository, in this working directory, I have a directory called data. I have a file in there called _"greetings.csv"._

So I can just set the default here to _"data/greetings.csv"_. Now, when I run this pipeline again without any command line options, it will use this default value. It knows it's a path, so it knows it should handle that as a path and not a string.

And then it's gonna pass that into a channel factory from this _params.input_ and create our channel, which is then gonna be used in this process called _sayHello_. Let's try it out.

Okay. Failed. Don't worry. This was expected. And if you're following the training material, you'll see it was expected there as well. Let's see what's happening here.

It's tried to run the pipeline. It's tried to execute the process, and it's got quite a similar error to the one we saw before.

Here it says: we tried to run \_echo, \_but instead of echoing the contents of this CSV file, it just echoed the path. And you can see it's the full absolute path here to this CSV file.

And then sure enough, because it tried to write that to this really complicated path, it didn't really know what to do. And it was outside the scope of the process work directory.

I mentioned in the start that Nextflow encapsulates every executed task within a special work directory. And if you try and write to data, which is outside of that work directory, Nextflow will stop you as a safety precaution. And that's what's happened here. We try to write to an absolute path and Nextflow failed and prevented us.

## 4.2. Use the splitCsv() operator to parse the file

Okay, let's have a look at this, channel and see what it looks like. We can do _".view",_ and I've copied this from the website. So _.view_, and we have a, dynamic closure here and we say a variable name "_csv"_ as the input. So that's the channel contents, and we say before splitCsv, and this is what it looks like.

If I run it again, it will still fail, but it will show us what's inside this channel. It is not particularly exciting. It is that _path_ variable. So you can see it's just a string here because it's being printed to a terminal, but it is a _path_ object, which contains the information and metadata about this file.

We don't wanna pass the metadata of the file to the input. We want to pass the contents of that file. If we look at the _greetings.csv_ file, you can see here that it's got these different variables here. _Hello, Bonjour, Holà_ again. And these are the things really we want to be passing to our process, not just the file itself as a single object.

So we need to parse this CSV file. We need to unpack it, get at the contents of the CSV file, and then pass the contents within the channel to the process.

As you can probably tell from the log message, we want to use the _splitCsv_, which is another operator, another channel operator. So if I do "_dot" "s"_, and then you can see it's auto suggested. Oops, _splitCsv_ and some brackets.

And then after _splitCsv_, I'm gonna put another _view_ statement just so we can see how it looks afterwards. Let's run the pipeline and see what we've got.

Okay. It still failed, but in a new and exciting way, which is progress.

This time again, we have some problem with our script, which has been rendered. Now. We haven't got the final path anymore, but we've got an array of variables, which looks very much like the error we had earlier on when we were passing an array as a fixed input.

With our logging from the view operator, we can see before _splitCsv_ was the path. And sure enough, after _splitCsv_, we have three different outputs and each of those outputs looks an awful lot like each of the rows from the _greetings.csv_ file, which makes sense.

So what's happened here is that Nextflow has parsed, this CSV file given us three objects, one array for each line of the CSV file. So then three times we've passed an array of variables to the channel instead of a single string value.

Okay, so last time we had this problem, we used _flatten_. Let's just very quickly. Try flatten and see what happens.

I can call these variables, whatever. So I'm gonna call it _myarray_ 'cause it's no longer really a CSV. Let's try and run it again and see what happens with _flatten_.

So this time we're gonna run, we parsed the CSV into three array objects, and then we flattened it. And this time it, it passed. And the Nextflow pipeline ran. However you can see that _flatten_ really goes to town and flattens everything. And so we get three independent array, entries for each row. And so it ran the process three times every row of a CSV. And now we have a whole bunch of results files, and 123, 456, and all kinds of things, not just that first column of the CSV, which is what we really wanted.

## 4.3. Use the map() operator to extract the greetings

So how do we get at just the first column? If flatten is too simplistic here, we need a more complex operator where we can actually customize and tell it what we want from the CSV.

To do that, we're going to use _map_. Basically _map_ just says, run some code, some function over every element that I get given and do some kind of transformation on it. And because it's so flexible, you'll see it come up in Nextflow code all the time.

By itself, it doesn't do anything. So we don't want regular brackets, we want a closure here and we need to tell it what to do. So I'm gonna say _"row"_, 'cause that's being given rows from the CSV, so it's a logical variable name. Is the input. And I want to return just the first element of that array.

Arrays in Nextflow are zero based, so we're gonna say just the first element, which is row zero. If we wanted to the second column, I could be one or the third column be two, and so on. We can return whatever we want here, but I'm gonna return just the first value.

And now, we can run the pipeline again and see if it does what we expect.

Sure enough, after _splitCsv_ we've got our arrays, and then after the _map,_ we have our nice clean strings, just _"Hello", "Bonjour"_ and _"Holà"_. And the pipeline is now doing what we want it to. Fantastic.

So we can get rid of all these view commands now. We don't need them anymore.

## Recap

We finished our kind of debugging and this is the code we end up with. Taking our CLI parameter called _input_, which is classed as a _Path_. Nextflow finds the path, loads it, and understands the CSV file. Returns all the different rows. And then we map just the first element of that row into the channel that kind of gives us the channel contents, which is passed to the process.

And the process runs over each element in the channel, which is three. And it runs the process three times, giving it three tasks. And those results are then published from the workflow, picked up by the process output. Published from a workflow and saved in the output block to a subdirectory called _"hello_channels"_.

Pretty cool. We are getting now to something that more closely resembles a real life Nextflow pipeline that you might run for some real analysis.

## Takeaway

Okay. Hopefully you're now getting a feel for what Nextflow channels and operators are and how operators work on channels and how you can create them.

Channels, like I said at the start of this video, are the glue of Nextflow. And you can see here that we can take different inputs and manipulate them and take that data and then pass 'em into downstream workflow logic.

And this workflow block here is really where you build up all that parallelization and all the clever logic, and explain to Nextflow how to build your workflow DAG, and how to orchestrate your pipeline.

Channels are not the easiest concept to get your head around. So take a break, have a little think about this, maybe read over the material again, and really make sure that you've got these concepts down because this is key to your understanding of Nextflow and the better you understand channels and the different channel operators and the different channel factories. The more fun you'll have writing Nextflow and the more powerful your pipelines will be.

This is not the same as regular programming in Python or other languages. We're not using _if_ statements here, this is functional flow programming using channels and operators. So it is a bit different, but it's also super powerful.

That's the end of this chapter. Go and have a quick break and I'll see you in the next video for part three where we're gonna go through Hello Workflow, and talk a bit more about the workflows.

Just like the previous chapter, there's a few quiz questions at the bottom of the webpage here, so you can have a quick run through these and make sure you understand all the different parts of the material we've just done. And aside from that, I will see you in the next video. Thanks very much.

Okay.

​
