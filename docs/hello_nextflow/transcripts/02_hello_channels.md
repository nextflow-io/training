# Part 2: Hello Channels - Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/lJ41WMMm44M?si=xCItHLiOQWqoqBB9&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important notes"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../02_hello_channels.md).

    The section numbers shown in the transcript are provided for indicative purposes only and may not include all section numbers in the materials.

## Welcome

Hi, welcome to part two of Hello Nextflow.

This chapter is called Hello Channels. We're gonna be talking all about this fundamental part of Nextflow.

Channels are the things which connect with different steps in your pipeline, the way that your data and logic flows through your workflow.

Okay, let's dig in.

Let's start off by going to training.nextflow.io

Hello Nextflow in the sidebar and clicking part two. Hello Channels.

All of the material is written down here so you can follow at your own pace and catch anything that you might have missed.

Once you've got the website open, you can load up Codespaces and we'll continue from where we were at the end of the last chapter.

## 0. Warmup: Run hello-channels.nf

For this chapter, we're going to be editing a different file. This one's called Hello Channels, so you can find that in a side bar, double click it to open.

Now if you've just come from chapter one, this file will look very familiar. The starting point here is basically where we finish chapter one, with our process called sayHello, our input, output, our publishDir and our params.greeting, and our simple workflow.

We're starting with a new file, so it's a level playing ground for everybody, but you can continue with your previous file if you prefer.

Note, I've also deleted all of the .nextflow\* files and the work directories here, just so it's a clean starting point. It doesn't matter if you do that or not, it's up to you.

Okay. Let's start off by checking that this pipeline still works as we expect. I'm going to bring up the terminal here.

Do "nextflow run hello-channels.nf" and hit enter.

It's going to run that little workflow, runs our sayHello step, generates a work directory with that hash, and here's our results folder and there's our output file, just as we expected from our default params.greeting.

So that's great. Exactly the same as chapter one, working as we expect.

## 1. Provide variable inputs via a channel explicitly

In chapter one, you were actually already using channels, you just didn't realize it. When we specified a string here, Nextflow automatically created a channel around that string for us, just because it knew that we were calling a process, so we needed an input channel.

The first thing we're going to do is make it explicit by actually typing out the channel itself.

## 1.1. Create an input channel

So I'm going to go to the workflow here at the bottom of the script, and I'm going to say greeting_ch. This is a convention we often use in Nextflow code to have a underscore ch at the end of a variable name when it's a channel, just so it's easy to identify that it is a channel, but you don't have to do that. Equals channel of Hello Channels.

What we've just used is something called a "Channel Factory" in Nextflow language. This is this thing here, we're setting this variable to a new channel, and this channel factory here is creating a channel for us in a particular way.

There are a handful of different channel factories that Nextflow has, to create channels from different types of inputs. Dot of is the simplest one, and just takes any strings that we give it.

Notice that when I hover over these words in VS Code, the Nextflow extension is giving me a popup explaining what this syntax does, and there's also a read more text at the bottom of that popup window.

If I click that, it'll open the Nextflow docs. In a new tab and take me straight to the documentation for this specific thing. In this case for channel.of.

## 1.2. Add the channel as input to the process call

Note that the extension is also giving us a warning, saying that we've created a new channel here, but it's not being used by anything.

So, let's fix that. I'm going to take the new channel name and I'm going to replace this params.greeting with our new channel.

Note that we're no longer using the command line flag --greeting now, params.greeting isn't being used, we're going back to hard coding this string. That's okay. I'm just trying to keep things simple. We'll come back later and use the params again.

## 1.3. Run the workflow command again

Okay, let's just double check this works. Bring up the terminal and note again. Nextflow run hello channels. Check output.txt, and there it is.

Great bit of a boring example, doing exactly the same thing as we did before, but now at least the logic is a bit clearer. We're being explicit about writing a new channel.

We've effectively just written more code to do the same thing. But this will start to make more sense as we become a bit more complicated with how we create our channels.

## 2. Modify the workflow to run on multiple input values

Okay, let's make this a bit more interesting. It's very rare that you want to run a Nextflow pipeline on a single input, so let's give it several inputs.

## 2.1. Load multiple greetings into the input channel

From the docs here. I'm going to copy in these different strings, three of them. Hello, Bonjour, Olà. Oh, get Hope. Copilot is suggesting a couple of others. So let's tab enter those.

The Nextflow docs here tells us that we can give multiple values to this operator, so it should work, but let's try it out and see what happens.

## 2.1.2. Run the command and look at the log output

Well. Yes and no. Let's see. It says that five of five tasks have run here, but it only shows us one hash, which is a bit odd. That's okay. Everything is it's expected here. By default. Nextflow uses a special type of output to a terminal called ANSI control codes, which means it overwrites the certain lines to give a nice compressed view of all the different processes which are being run.

This makes much more sense when you have larger workflows and are running hundreds or thousands of different samples. You can just generate so much output on the terminal, it's impossible to look at, whereas this updating view gives you a real time progress for you.

## 2.1.3. Run the command again with the -ansi-log false option

If you want, you can run it again, and this time I'm going to use an additional Nextflow core argument with a single hyphen saying, "-ansi-log false". This uses the previous version of the Nextflow log output. And here you can see all the individual processes which have been launched.

It's up to you whether you do this or not. The output from Nextflow is exactly the same in both cases.

## 2.2. Ensure the output file names will be unique

Okay, let's have a look at the output files, then we'll go to results. But there's only a single output file. What's happened? We saw that the process had run lots of times. We can go into the work directory and see all of the different hashes, all the tasks were executed properly. But if you remember in our process here, we're saving everything to an output.txt file and then publishing that to this directory.

So the same file was created five times, and then it was overwritten five times. And we just have whichever task happen to execute last.

## 2.2.1. Construct a dynamic output file name

The way we fix this is by using a dynamic output file name. Here we already have a variable called greeting within the process, so we can use that in the output file name. I copy that and I do $greeting-output.txt.

I am going to surround this in quotes, just so that bash doesn't get confused by any spaces which might creep in here. And then I'm going to take the same file name and update the output here.

It's really important that the output matches this, because otherwise, this file won't be found and Nextflow will crash.

I'm going to make one more really important edit, which is I'm going to change these single quotes for double quotes. Note that the color of the code changed when I did that. This variable is only expanded if we use double quotes. If I use single quotes here, it's used as a literal value, and I'd get a single file called $greeting-output, which is not what I want.

## 2.2.2. Run the workflow

So let's put the double quotes back and give it a try.

I am just going to clear up my directory before I start, so it's easy to see the new files. I'm going to delete anything called .nextflow, work, and results.

And I'm going to run that Nextflow command again and let's see what files are created. So it runs the five processes there. If you were watching very closely, you might have seen that line update as it was running.

And now we can go into the results directory, and sure enough, we have five different outputs, and they're all prefixed with the different greeting.

If I open each of these, we'll see that they each contain the corresponding greeting. Fantastic. That's what we want.

## 3. Use an operator to transform the contents of a channel

Okay, so now we know what channels are and we know what channel factories are. What about operators? This is another term for part of the Nextflow language, which is a series of functions which allow us to operate on channels to do certain things to them. Nextflow, comes with a suite of operators, which allow us to manipulate channels in a variety of different ways.

## 3.1. Provide an array of values as input to the channel

Let's work through this with an example. Let's say that we want to take these input strings, but instead of just putting them directly into a channel factory, we want to define them as an array.

## 3.1.1. Set up the input variable

So I'm going to take these and do that as a new line above and say, greetings, array.

There we go. I'm going to take that array variable and put it into the channel.of, and hit save.

## 3.1.3. Run the workflow

Now, let's see what happens. Go back to my terminal. I'm just going to clear up all those temporary files again. And let's run the workflow.

Not good. Okay. It broke. That's okay. I expected it to break this time. Debugging what goes wrong when an Nextflow workflow fails is a key part of being an Nextflow developer. This will happen a lot and it's important to understand what the error message says and how to deal with it.

The Nextflow, error messages are actually quite structured. It tells us which process went wrong. It gives us an error message for a reason. It says what the command was that it tried to run within that particular task, what the exit status was, what the output was on where that task work directory was.

Note that I can option, click this in VS Code and it opens it in a sidebar so I can go straight there and view all of these hidden files, which we talked about in the previous chapter, including the .command.sh file. This you can see is the same as the commands which was executed here.

By looking at this file, we can get a feel for what might have gone wrong here instead of running a single task for each element in the array as it did last time, it just provided the entire array in one go as a string. So we need to unpack that array into individual values before we pass it into the channel. Let's go back and see if we can do that using an operator.

## 3.2. Use an operator to transform channel contents

In this case, we're not going to change the array before we pass it into the channel. We're going to adjust the channel so that it behaves in the way we expect. We're going to do that by using the flatten operator can do dot start typing and we can see that VS Code extension starts suggesting all the different operators we have available.

## 3.2.1. Add the flatten() operator

And I'm going to select flatten. Note that white space doesn't matter in this context for Nextflow. So you can put these operators on a new line if you want to. So I can drop this down here and indent it so it sits underneath ".of" and you'll see that people often chain lots of operators like this onto a channel and indent it in this way so that it is easier to read.

You can also see, like before I can hover over this and read what the flatten operator's doing, and also follow a link to the documentation if I want to.

So this operator is taking this channel, which has a single array within it, and separating out the array values.

## 3.2.2. Add view() to inspect channel contents

We can peek into the channels by using the special view operator, and I'm going to add a couple of them here. This is a bit like using print statements in other languages. So I'm going to do dot view and then I'm going to use these squiggly brackets.

This is called a closure. This basically gives additional code to the view operator, which it will execute on each item within the channel. In this case, I'm going to say greeting before flatten. Greeting.

I'm defining a variable here, which is just within the scope of this closure. So this variable is only used here and I could call it whatever I wanted. It doesn't really matter. I'm just using greeting to make it easy to read.

In some Nextflow pipelines, you might see people use a special implicit variable called "$it". Like that. This is a special variable within Nextflow code, which is a shorthand so that you don't have to do the little definition of a variable. However, over time we're thinking, this is not very clear to people who are new to Nextflow, and we discourage usage of "$it" now.

So I'm going to stick with the previous behavior of greeting and using it like this because that's more explicit and it's clearer on what's going on.

I'm going to then copy this line and do exactly the same thing again after the flatten arguments. The view operator's a bit special because it does something on the elements, but it also just continues passing them onto the next operator so we can chain it in the middle of a chain of operations like this, and it will print the status there and keep going. So hopefully this will show us what the channel looks like before and after the flatten operator.

## 3.2.3. Run the workflow

Let's try it out. Clean. Clean everything up in the workspace. Run the pipeline again.

Okay, so we can see it ran our five processes. Again, it didn't crash with an error, so that's definitely good. And now we have the before flatten and it sure enough we have our array and we have after flatten, printed five times once for each elements of the array. That's exactly what we were hoping for. So that's really good news. And that fits exactly with what we'd expect from the code.

We don't need these debug statements anymore, so I can either comment them out or delete them. I'm going to delete them just to keep my code nice and clean. Okay, great. This example is now working nicely and we can start to see how channels can do a bit more complicated logic.

## 4. Use an operator to parse input values from a CSV file

Now we're going to try and do this using a file with a series of inputs instead. This is a very common way to write Nextflow pipelines using a sample sheet or a CSV of metadata.

## 4.1. Modify the script to expect a CSV file as the source of greetings

If I go over to the sidebar, you can see greetings.csv in the example repository, and this is a very, very simple CSV file that just contains three lines with three different greetings. Let's see if we can use this CSV file within our workflow.

I'm now going to go back to using params like we did in chapter one, so that we can have a command line input.

I'm going to delete this greetings array.

## 4.1.1. Switch the input parameter to point to the CSV file

I'm going to set params greeting to the file name, which is greetings.csv, and I'm going to use this special variable to generate the channel. I'm going to put that in there, and the errors go away. Remember that this is setting this variable by default now. So if I run the pipeline without any arguments, it'll use greetings.csv, but I could do --greeting to overwrite this variable if I wanted to.

## 4.1.2. Switch to a channel factory designed to handle a file

Okay, we're passing a file now rather than a string or an array of strings, so we probably need a different channel factory.

We're going to get rid of "of" which we've been using so far, and instead use .fromPath. This does exactly what it sounds like. It creates a channel with paths instead of values, using a string file name or glob. I am also going to remove the flatten operator as we no longer need this, now that we're passing a file.

## 4.1.3. Run the workflow

I'm going to hit save, open up the terminal, run the workflow, and then see what happens.

Okay. It crashed again. Don't worry. I was expecting this one as well. Let's have a look at the error message and see if we can figure out what's going wrong. Here we can see the command executed, and a bit like before where we had the whole array printed. Now we have the file path being echoed into the command, rather than going through the contents of the file.

## 4.2. Use the splitCsv() operator to parse the file

So to use the contents of the file instead, we need another operator. The operator we're going to use for this one is called splitCsv. Makes sense, because it's a CSV file we're loading.

## 4.2.1. Apply splitCsv() to the channel

Ok, so splitCsv. Close bracket. We don't need any arguments here. And again, I'm going to use some view operators to give some insight as to what's going on here.

.view csv after splitCsv. Before split Cv.s

## 4.2.2. Run the workflow again

Okay, let's try running this and seeing what happens.

Okay, we've got a bit more output this time, but it still failed. We can look at the view statements, and here you can see before split CSV, and we have a file path as we saw in the previous error message. After split CSV, we now have three values corresponding to the three lines in the CSV file.

However, you can see that each of these values is surrounded by square brackets. So each one of those was an array in its own right, and that's given us the same area that we had before where it's trying to echo an array rather than just a single string.

If we think about a CSV file, this kind of makes sense. Typically, a CSV file will have rows and columns, so split CSV does two dimensional array. The first dimension of the array is each row, and then there's a second dimension, which is each column for each row.

So here we only have a single value on each line, so we have a single column, so we have a one element array for each line of the file.

That's fine. We just need another operator to collapse that array for each line of the parsed CSV file. Let's clean this up. Get rid of a terminal and see what we can do.

## 4.3. Use the map() operator to extract the greetings

Now we could use the flatten operator again, which we used before. We've seen how that can collapse an array into a series of values, which would work very well here. But I'm going to use the opportunity to demonstrate another operator, which is very common within workflows called the map operator.

## 4.3.1. Apply map() to the channel

I'm going to do dot map and I'm going to do item item[0].

If you write a lot of other code languages, you might be familiar with the map operator already. It takes an iterable, such as an array or a channel, and it does some operation on each value of that.

Here we're saying that we should define a variable called item within the scope of this closure, and then we want to return, just the first value in that array. So item index zero.

This is effectively flattening the array. You can see how we could extend this to be more complex, though: if our CSV file had six columns, but we're only interested in the fourth column, we could access a specific index here. Or do any other kind of operation on the value before we pass it onto downstream processing.

So the map operator is extremely flexible and very powerful for modifying channels in flight. Let's put in another view statement just so we can see what it's doing in our execution. Can adjudicate that line and move it down. And after map.

## 4.3.2. Run the workflow one more time

Let's bring up the terminal and try running the workflow.

Okay, no errors this time. That's a good sign. We can now go through all these different outputs from the view statements. Before split CSV, we had a single path. After split CSV, we had the single value arrays, and then after map, we have just the values without any array syntax. Let's go up to the results directory, and here are our output files behaving exactly as we wanted them to.

There's a little bonus here. You can actually see that the view operators are slightly mixed up in the order that they've done the output. This is because Nextflow is doing parallelization of these different tasks. So after it split the CSV, there are three elements in this channel, and it's handling the processing of those three elements in parallel automatically. That means that the order of the outputs is stochastic and can vary. In this case, it just happened that some of the view operators returned after the subsequent step had been completed, and so it came in this order.

If I run the same workflow again. Then sure enough, it's come in a different order and this time we've got the split CSVs and the maps in the order we'd expect.

So just bear in mind, you cannot rely on the order of outputs from a process task because Nextflow's handling this parallelization for you automatically. Nextflow does that for you with its data flow logic, and that's the real power of Nextflow.

Okay, this is probably one of the most important chapters of the whole training. Once you understand channels, channel factories and operators, you start to key into the strength of Nextflow and what makes it unique as a programming language. This functionality allows Nextflow to parallelize all of your workflows for you and generate extremely complex workflow logic with a very clean syntax and a push data flow model. It can be a bit of a strange concept at first, but once you get used to writing code like this, it will quickly feel natural and before you know it, you'll be writing fantastic workflows.

Have a break, cup of tea, walk around and let's move on to chapter three, where we start to extend these concepts into more complex workflows. See you in the next video.

[Next video transcript :octicons-arrow-right-24:](03_hello_workflow.md)
