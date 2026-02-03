# Part 4: Hello Modules - Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Xxp_menS0E8?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important notes"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../04_hello_modules.md).

    The section numbers shown in the transcript are provided for indicative purposes only and may not include all section numbers in the materials.

## Welcome

Hi, welcome to Part Four of the Hello Nextflow training course.

This chapter is called Hello Modules, and we'll be talking about how to modularize Nextflow code. What we're going to do is take our one workflow script and break it up into separate files.

This makes the code easier to navigate and maintain as your workflow gets bigger, and also makes it possible to share modules between pipelines so that if you have multiple pipelines using the same tool, you only need to write that process once.

A classic example of this is the nf-core modules repository, which has thousands of different tools in ready to go modules, which you can install and use in your workflow.

Nextflow can also work with sub workflows, which are like modules, but they have multiple processes. That's outside the scope of this training, but it works in basically the same way.

Okay. Let's take a look.

As usual, start by going to training.nextflow.io.

Go to "Hello Nextflow" in the sidebar, and we're doing part four: "Hello Modules".

I'm now going to jump into my GitHub Code Spaces environment and take a look at the file "hello- modules".

Just as before, we're starting off at the end point of the previous chapter, so this script should look familiar. We've got our three processes, say hello, convert to upper and collect greetings, and in a simple workflow, which runs these three commands and emits a message at the end. We have two parameters called greeting and the batch, which specifies the name, which is used for the collected output file at the end.

## 0. Warmup: Run hello-modules.nf

We can verify that this workflow still works as we expect by doing nextflow run hello, modules.

Great. It ran three tasks with each of this process, one collected task, and it told us that there are three greetings in this batch. If we go into results, we've got our different output files here, including the collected test, batch output.

## 1. Create a directory to store modules

Right. Let's do some modularization.

It is generally a good idea to put modules into a subfolder in your pipeline repository, just to keep things tidy. You can call this whatever you want, but by convention we usually call it modules.

So let's go ahead, go into a terminal and do make the modules. You can see it pop up in the sidebar and VS Code here.

## 2. Create a module for sayHello()

I'm then going to create a new file for my first module. You can do "touch" or "code" or you can do this in the side bar, it really doesn't matter. So I'm going to do code modules and I'm going to name it after the process. So sayHello.nf . NF is a traditional file extension for Nextflow files.

Going to hit save here and we can see our new module file turns up.

## 2.2. Move the sayHello process code to the module file

Right, next I'm going to take the module code from the workflow. I'm also going to take the hash bang here and copy that in first so that it's clearly a Nextflow file. And then I'm going to take this process and I'm going to cut. So I'm going to remove it from my main workflow script and I'm going to paste it into this new module.

That's all the content that this module file is going to contain. Just a single process, no workflow, no logic, just a process alone.

I can now close this file.

## 2.3. Add an import declaration before the workflow block

Now my workflow is missing that first process, so we need to bring it back by importing it. The syntax for this is very similar to other programming languages, so it may feel familiar. We do include squiggly brackets, the name of the process, say hello, and then from the file path modules, say hello, nf. Fantastic.

A couple of tricks here. The VS Code extension is clever about this. It recognizes this file path and you can hover over it and do follow link . Or I'm on Mac, I can do option click and it opens this file. So we can quickly jump to it.

This process name is now being used by the workflow down here, and we can do the same thing here. It shows us a little bit of information about that process, and again, I can hold option, click on it, and it will open it in the editor.

So it's a really fast way when you have lots of files for your different processes to quickly navigate around your code base in VS Code.

Okay. That's basically it for this chapter. Now we just do the same thing again for the other processes.

## 3. Modularize the convertToUpper() process

So let's create a new file here. Call it Convert to upper nf. Again, copy the hash bang. And then cut the process out.

Copy the process name there, include a new include statement with the new process name.

## 4. Modularize the collectGreetings() process

And then do the same for the third process. New file, connect. Greetings,

do the hash bang. Cut the process, paste the process, and do a new include statement.

Now you can see here I've got an error underline here saying invalid include source. And this is actually a genuine error that I made because I was moving a bit too quickly. If you look closely, you can see I missed the T and convert to upper

So VS Code very usefully has told me that I've made a mistake there. If I fix that file name, the error goes away. It's a good example of why the error checking within VS Code is so useful for writing Nextflow code. Otherwise I wouldn't have spotted that and I would've only found out much later when I tried to run the workflow.

Our main pipeline script is now looking much simpler. It doesn't have any processes in, we just have three include statements and our workflow. We haven't changed any of the logic of the workflow. We haven't changed any of the process code, so hopefully it should work in exactly the same way.

## 4.4. Run the workflow to verify that it does the same thing as before

Let's check. I'm going to open a terminal and I'm going to run exactly the same command as before.

Sure enough, it's run our processes, say hello, convert to upper collect greetings, and given us three greetings again.

So we've moved our code around, but we haven't changed anything about how the workflow executes and it's completely unchanged. The only difference is we now have cleaner code, easier to maintain, and easier to share with others.

And that's it. It was a short chapter. It's a simple concept, but it's very powerful and key to how we write more complex Nextflow workflows. So it's important that you understand it and get into the habit of using it.

In the next chapter, we're going to have a bit of a change of pace and stop thinking so much about the syntax of writing Nextflow code, and think a little bit about how we use software in the processes themselves. Join us in part five for Hello Containers.

[Next video transcript :octicons-arrow-right-24:](05_hello_containers.md)
