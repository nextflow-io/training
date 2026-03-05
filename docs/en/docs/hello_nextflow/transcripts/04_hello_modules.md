# Part 4: Hello Modules - Video Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/43Ot-f0iOME?si=0AWnXB7xqHAzJdJV&amp;list=PLPZ8WHdZGxmWKozQuzr27jyMGqp9kElVK&amp;cc_load_policy=1" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important notes"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../04_hello_modules.md).

    The section numbers shown in the transcript are provided for indicative purposes only and may not include all section numbers in the materials.

## Welcome

Hi, and welcome back to part four of Hello Nextflow. This section is all about modules, and it's quite a short section of the course. We're not actually really gonna write very much code, it's more about how we organize the code in our pipeline.

Up until now, we've been putting everything into a single file, which is fine, and that's actually how we used to build Nextflow pipelines back in the old days.

But as that pipeline grows, the script gets longer and longer and longer and harder and harder to navigate, to maintain, and also means that we can't really share any of the code.

Nextflow modules allow us to pull processes out of that main script and then import them. That means that the code is easier to navigate and it also means we can share that module code between different pipelines.

This little diagram on the, main page for the docs shows the concept nicely. Instead of one huge script, we are going to include these separate module files, from different module scripts, and it's all gonna get pulled into the workflow, but it's still gonna run in exactly the same way.

So let's jump into GitHub Codespaces and have a little look around. As before, I've cleaned up my workspace here a little bit. Removed the old Nextflow directories and the work directory and so on. But it doesn't matter if you still have those files around.

I'm gonna start working in the hello modules file, which is basically where we left it at the end of a previous chapter. We have our three processes in here. We have a couple of params, the workflow block, where we are running those three processes and stitching them together with channels. Then we publish the output channels and we have the output block saying how to publish those files.

## 1. Create a directory to store modules

Now, like I say, we're not really gonna write or edit very much code. We're just gonna move the code around that we have already. Nextflow module files typically have a single process in them, and by convention we normally keep them in a directory called modules. But you can call that whatever you want. But I'm gonna keep a modules directory in my repository here, and then I'm gonna create one file for each process. So I'm gonna say new file, sayHello.nf.

## 2. Create a module for sayHello()

Now I'm gonna take my process and I'm simply gonna select this code, cut it out of the main hello modules file and paste it in here.

Obviously that doesn't do anything by itself. Our main script still needs that process, so we need to pull it back in somehow. And we do that with the include statement.

So I type include and some squiggly brackets, and then I take the name of the process. And I say from, and then I give it a relative file path. So it says, starts with ./ because it's relative from where this script is saved. So it's modules sayHello.nf.

Notice that the VS code extension is quite helpful here. It tells us, if it can find this file and if it can find a process, which I'm naming. If I deliberately put in a typo here, it gives me an error straight away and it'll tell me that it can't find this process I'm trying to import. So just keep an eye on any errors that you find.

And that's really it. We have our process still here. There's no changes needed down here. The process has got the same name and it's executed in exactly the same way. It's just the actual code of the process is now in a separate file.

We can run the Nextflow workflow again, it is going to work in exactly the same way. And that's basically the rest of this chapter of the course is just moving these three processes out into their own files.

So let's do that now. I'm gonna quickly create a new module file for the second process: convertToUpper.nf. I'm gonna cut that code out, paste it in there. And then I'm gonna include that one. let's just, great.

And then I'm gonna create a new file for collectGreetings.nf. Cut that out.

Lots of cut, cutting and copying and pasting.

And now our main workflow script is suddenly looking much, much shorter, much more approachable and a lot easier to read.

And you can see how the project now starts to build up with our different files. We can dive into the detail in the places we want to. Navigate our way around to find specific steps in the pipeline much more easily, and get an overview of what the pipeline is doing quickly.

## Navigating modules with VS Code

Now, of course, the downside of doing this is that if you have a big pipeline, you'll have a lot of module files and they could be organized in multiple sub directories or all kinds of things. Now, again, one little tip here. The VS Code extension is pretty good at navigating your code base for you and also telling you about the code there.

You can see VS Code understands what this process is and gives me a little overview of it when I hover so I can see without having to go off and finding the source code, what the inputs and the outputs are, which is typically the most important thing when I'm using it in a workflow.

And also if I hold command, I'm on a Mac, and I click the process name, it opens the file directly straight away. Pulls it in. So I can just jump straight there without even thinking about what the actual file parts are. And that works anywhere, I can do that also, whether processes are being called. So that makes really fast.

## 4.4. Run the workflow

Okay, let's just check that the pipeline still does run as we expect it to. So bring up the terminal. Let's do "nextflow run hello modules", and see if it executes without any problems.

Hopefully the whole point of this is that the pipeline is basically unchanged, so you shouldn't, really see any changes from when we ran it before. The output here looks exactly the same, and you can see our results directory with all the same files, so that's great. No change is good.

## A note on nf-core/modules

Just before we wrap up, I want to quickly touch on the power of collaboration when it comes to modules. These files are sat in my repository, so it is not obvious straight away how we might collaborate on them. And there are many different ways that you can do this, but probably the biggest and best known example of this is nf-core.

If I go to the nf-core website, I go to resources, and modules. You can see that nf-core has a huge library of modules, nearly just under 1700 modules when I view this. And so I can type the name of any of my favorite tools, go and find if someone else has already written a module for it, and see this pre-written module process here with all the inputs, the outputs, the software containers, all this information, and you can see on the side here, how many different nf-core pipelines are all using this single shared process.

This is a bit of an extreme example, but you can see this is really reusing this code. And if I click through to the GitHub source for this, it's exactly the same as what we're doing. It's just a big process in a file.

Now on the nf-core side, we do some tricks to be able to share those files and bring them into different repositories. And if you want to know more about that, go and check out the course we have about using and building with nf-core specifically. But I wanted just to give you an idea of quite how powerful this concept of code reuse can be.

## Wrap up

Right, that's it for modules. I told you it was a short section of the course. Check out the quiz, make sure you understand it and make sure everything's still working properly. And I'll see you back in the next video, which is all about software containers. Thanks very much.

I.
