# Orientation - Video Transcript

<div class="video-wrapper">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/G3CV-FcV-rc?si=nyLvwhrSB2m1NPc5&amp;list=PLPZ8WHdZGxmXiHf8B26oB_fTfoKQdhlik" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen></iframe>
</div>

!!!note "Important note"

    This page shows the transcript only. For full step-by-step instructions, return to the [course material](../00_orientation.md).

## Welcome

Hi, welcome to Hello Nextflow. My name is Phil Ewels. I'm Product Manager for Open Source at Seqera, and I'm delighted to be here today to take you through this first Nextflow training course.

We're going to go through the basics of Nextflow, explaining how to write and run pipelines and configure them.

And you're going to build your own simple multi-step pipeline. We'll cover terminology like operators and channel factories, and by the end of the course, you'll be ready to get started building your own bioformatics pipelines.

If you have any questions, please reach out on community.seqera.io. We have a really active Nextflow community there's a section dedicated to training, so just let us know where you're stuck and someone will be able to help out.

Right. Let's get started.

## Training Website

All of the training material for the Nextflow courses lives at training.nextflow.io. You can go to it in your web browser. So boot that up now and we can have a look through.

I'll be running this with version 2.1.1. We push small updates and fixes here and there, so don't worry if it's a little bit off, but if a material's drifted too far, you can always use this version picker at the top to pick the exact version of the materials that I'm going to be talking through.

If you're more of a light mode person, you can change the theme for the website here.

See translations here, though at the time of recording, it's really only English, which covers this new material.

And also see all the source code for the training website and everything we'll be working with on GitHub.

The homepage here lists all of the different training material courses that we have. So I scroll down, we'll see Nextflow for newcomers with the Hello Nextflow course we'll be doing here. You can see all the other courses we also have, which work in a similar way.

## Environment Setup

I am actually going to start off by using this first one at the top, which is common for all of the training courses, and it's specifically about setting up our environment.

I click through, it takes me to this section, and we can see instructions for developing locally. If you want to use your own laptop with your own copy of VS Code and your own software installations, or what we expect most people to do, which is to use something called GitHub Codespaces.

Codespaces is a service provided by GitHub where they run a web server in the cloud, which you can connect to. That server has VS code installed, where you can run that in your web browser, or if you prefer, connect that to your local installation of VS code. All of the computation, all of the files, all of the editing happens remotely, which means that all the software you need comes pre-installed and is the same for everybody.

## Creating a GitHub Codespace

In order to create the codespace with everything we need, look for the buttons in the docs material, which say "Open in GitHub Codespaces". I'm going to click that now, open it in a new tab. And I'm presented it with this webpage. Now you can see this is pre-configured to set with an nextflow-io training.

I can just click create new codespace. But actually we recommend that we use a slightly bigger machine for the Nextflow training with four CPUs instead of two. You can change which version of a material it uses. So this is defaulting to 2.1.1 because that's the version of the docs I followed the link from. But I could also set it to a specific branch of the repository if I want to.

Now I'm going to click create codespace. And it's going to start setting up the environment for me.

## Codespace creation

Now, the first time you do this, it's going to take quite a long time, so now is a good time to go and have a cup of tea. Get yourself comfortable, chat to the person you're sat next to.

If you're interested, you can click building codespace down here to see the logs of the setup. And you can see here it's pulling a Docker image with everything I need and configuring of the environment.

Now, you only have to wait like this the first time you create a codespace. If you go to github.com/codespaces here, you'll see all the different Codespaces you have open. Here's the one I've just created. Next time you do this, you can go here and you can select the previous codespace and just jump straight back into it. And it's a much, much faster process to warm up that existing environment. That will also keep all of the changes that you've made to VS Code and to the files , so you won't lose your progress if you leave and come back.

You can click the three dots here to do other actions. For example, if you configured it with two CPUs and now you want four, you can change the machine type. Or if you wanna start from scratch and fresh, you can delete the codespace.

## Intro to VS Code

Okay, Codespaces is finished setting up my environment and now presented with VS Code in the web browser.

If you're used to VS code. This will feel very familiar if you've not used it before, it's pretty simple. There's a few different parts to the page you need to be aware of.

Over on the left here, we've got the sidebar. You can see the Explorer set up with all the different files in the GitHub repository from the training repo.

On these buttons down the left, can be different tooling. In the sidebar. I can search all the files in all the project. I can work with Git, can work with GitHub, all different things like that.

At the top here is the main menu. The file explorer is the one we'll have up most here, and you can right click any of these files and do the normal things you'd expect. You might need to click through some warnings like this where it like cut copy and you can download to your local machine as well.

When the codespace loads, it gives us a preview of the markdown file in this main area here. This is just the same as the one that renders on github.com. I can close that and if I double click that read Readme file, you'll see it opens it as code in the code editor and just as with any other file, we can edit this code directly.

Finally at the bottom here, we've got the terminal window. I was looking at the logs as it built, so that's what the current thing it's showing. I can also press this plus button to start a new terminal session. This isn't running on my machine. Remember, this is running in the cloud, and if I do tree three to depth of two, you'll see all the same files here, which were over on left.

## Showing just "hello-nextflow" files

This GitHub repository contains all the different training sets, not just the one we're doing. So if you like, you can focus on just the Hello Nextflow folder. One way to clean this up a little bit is to go to the menu file and then add folder to workspace.

We click that go to training. Hello nextflow, and click add. It will refresh your screen. And then over in the Explorer, we now have two different workspaces, the one we had before for training and one with just Hello Nextflow .

If you like, you can right click on training and click remove folder from workspace to get rid of it from the sidebar completely.

Now we've got just files for this particular training course in the side. I can hide that warning and now I can do the same thing in the terminal here and do CD for change directory. Hello, Nextflow. And again, we have the same files here, which are on the sidebar.

## Hello Nextflow: files

Looking at these files for the Hello Nextflow course.

We have a bunch of .nf files, which are for Nextflow, and there's one of these files for each of the chapters of the training course. We'll work through these files and modify them in the exercises.

We also have a nextflow.config file, which just has basic config settings for running Nextflow in this environment, which you don't really need to worry about at this point. A greetings.csv file, which we'll use for processing data, which will be introduced in the next part of this course, and a test-params.json file, which will be used in part six and you can ignore for now.

These Nextflow files are just the start of each exercise. If you want to see how they should look when they're finished, you can go into a solutions directory and there are the answers for each part of the training course, so you can see a working version of what you're aiming towards.

## Opening a terminal

If at any point you close the terminal and can't remember how to get back, don't worry about it. These buttons at the top, right open and close different panels in the workspace. So click this one for bottom panel and it will reappear. And just make sure you've got terminal selected here. You can also click this button here, the arrow on the right hand side of a terminal to make it full screen.

You'll see me doing that quite a lot because I have VS Code zoomed in so that you can read the text. Depending on your screen size, you may or may not need to do this. The same goes for minimizing the side panel.

Right. That's enough for environment. I think we're ready to get started. Join me back in the next video for chapter one.

[Next video transcript :octicons-arrow-right-24:](01_hello_world.md)
