# Booting a Jetstream Computer Instance for your use!

What we're going to do here is walk through starting up a running
computer (an "instance") on the Jetstream service. 
* Jetstream is run by NSF and provides elastic cloud computing services. 
* "Cloud" computing is a fancy word for being allowed to temporarily use someone else's computer somewhere else with full administrative privaleges to install whatever software we want. 
* Sometimes we need to use computing resources that are larger than our personal computers. Cloud computing lets us decide how much capacity we want to be using. 

If you would like to read more about cloud computing, see this [Carpentry Cloud Computing lesson](http://www.datacarpentry.org/cloud-genomics/01-why-cloud-computing/).

Below, we've provided screenshots of the whole process. You can click
on them to zoom in a bit.  The important areas to fill in are circled
in red or pointed out with arrows.

Some of the details may vary -- for example, if you have your own XSEDE
account, you may want to log in with that -- and the name of the operating
system or "Image" may also vary from "Ubuntu 18.04" or "DIBSI 2018" depending on the workshop.

-----

First, go to the Jetstream application at [https://use.jetstream-cloud.org/application](https://use.jetstream-cloud.org/application).

Now:

## Request to log in to the Jetstream Portal

Click the login link in the upper right.

[![login](images/login-1.thumb.png)](images/login-1.png)

## Use "XSEDE"

Choose "XSEDE" as your account provider (it should be the default) and click
on "Continue".
           
[![foo](images/login-2.thumb.png)](images/login-2.png)

## Fill in the username and password and click "Sign in"

Fill in the username and then the password (which we will tell you in class).

[![foo](images/login-3.thumb.png)](images/login-3.png)
           
## Select Projects and "Create New Project"

Now, this is something you only need to once if you have your own
account - but if you're using a shared account like tx160085, you will
need a way to keep your computers separate from everyone else's.

We'll do this with Projects, which give you a bit of a workspace in which
to keep things that belong to "you".

Click on "Projects" up along the top.

[![foo](images/login-5.thumb.png)](images/login-5.png)
           
## Name the project for yourself, click "create"

Enter your name into the Project Name, and something simple like "ANGUS"
into the description. Then click 'create'.

[![foo](images/login-6.thumb.png)](images/login-6.png)

## Boot an instance with a pre-built image 

Select [![foo](images/Images.thumb.png)](images/Images.png)

## Find the "DIBSI 2018 workshop image" image, click on it

Enter "DIBSI" into the search bar - make sure it's from
June 22nd, 2018 by Titus. This images is based on Ubuntu 18.04 devel and docker, with Rstudio and [bioconda](https://bioconda.github.io/) package manager added.

[![foo](images/Image_Search.thumb.png)](images/Image_Search.png)
           
Launch [![foo](images/DIBSI2018_launch.thumb.png)](images/DIBSI2018_launch.png)

## Name it something simple

Change the name after what we're doing - "Day1_workshop_tutorial", for example,
but it doesn't matter. Pull down the drop-down menu under 'Project' to select your name. Then make sure the appropriate Resources are selected. You probably won't have to change these. The 'Allocation Source' will already be selected. (This is our XSEDE allocation grant ID.) The 'm1.medium' instance size will already be chosen. This is the minimum instance size. A larger instance can be selected, depending on what we will be doing. The 'Provider' will be randomly chosen as either 'Jetstream - Indiana University' or 'Jetstream - TACC'.

[![foo](images/Launch_Instance.thumb.png)](images/Launch_Instance.png)

## Wait for it to become active

It will now be booting up! This will take 2-15 minutes, depending.
Just wait! Don't reload or anything.

[![foo](images/login-11.thumb.png)](images/login-11.png)
           
## Click on your new instance to get more information!

Now, you can either click "Open Web Shell", *or*, if you know how to use ssh,
you can ssh in as username (we will provide in the class) on the IP address of the machine - see
circled information below.  Note that you'll need to use the private key
file we sent around to everyone in last the pre-workshop e-mail if you decide to
use your system terminal. Here are the logging [instructions]((https://github.com/ngs-docs/angus/blob/update/schedule/jetstream/login.md)) using private-key

[![foo](images/login-12.thumb.png)](images/login-12.png)

## Miscellany

There's a possibility that you'll be confronted with this when you log in to jetstream:

[![foo](images/possible_instance_problem.thumb.png)](images/possible_instance_problem.png)

A refresh of the page should get you past it. Please try not to actually move any instances to
a new project; it's probably someone else's and it could confuse them :)

## Suspend your instance

You can save your workspace so you can return to your instance at a later time without losing any of your files or information stored in memory, similiar to putting your physical computer to sleep. At the Instance Details screen, select the "Suspend" button. 

[![foo](images/suspend-1.png)](images/suspend-1.png)

This will open up a dialogue window. Select the "Yes, suspend this instance" button.

[![foo](images/suspend-2.png)](images/suspend-2.png)

It may take Jetstream a few minutes to process, so wait until the progress bar says "Suspended."

### Resuming your instance

To *wake-up* your instance, select the "Resume" button.

[![foo](images/resume-1.png)](images/resume-1.png)

This will open up a dialogue window. Select the "Yes, resume this instance" button. 

[![foo](images/resume-2.png)](images/resume-2.png)

It may take Jetstream a few minutes to process, so wait until the progress bar says "Active." 

[![foo](images/resume-3.png)](images/resume-3.png)

## Shutting down your instance

You can shut down your workspace so you can return to your instance another day without losing any of your files, similiar to shutting down your physical computer. You will retain your files, but you will lose any information stored in memory, such as your history on the command line. At the Instance Details screen, select the "Stop" button. 

[![foo](images/stop-1.png)](images/stop-1.png)

This will open up a dialogue window. Select the "Yes, stop this instance" button.

[![foo](images/stop-2.png)](images/stop-2.png)

It may take Jetstream a few minutes to process, so wait until the progress bar says "Shutoff."

[![foo](images/stop-3.png)](images/stop-3.png)

[![foo](images/stop-4.png)](images/stop-4.png)

### Restarting your instance

To start your instance again, select the "Start" button.

[![foo](images/start-1.png)](images/start-1.png)

This will open up a dialogue window. Select the "Yes, start this instance" button. 

[![foo](images/start-2.png)](images/start-2.png)

It may take Jetstream a few minutes to process, so wait until the progress bar says "Active." 

[![foo](images/start-3.png)](images/start-3.png)

## Deleting your instance

To completely remove your instance, you can select the "delete" buttom from the instance details page. 

[![foo](images/delete-1.png)](images/delete-1.png)

This will open up a dialogue window. Select the "Yes, delete this instance" button.

[![foo](images/delete-2.png)](images/delete-2.png)

It may take Jetstream a few minutes to process your request. The instance should disappear from the project when it has been successfully deleted. 

[![foo](images/delete-3.png)](images/delete-3.png)

[![foo](images/delete-4.png)](images/delete-4.png)
