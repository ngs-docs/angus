# Jetstream: install bioconda on a blank Linux machine.

Learning objections:
* learn what bioconda is and how to install it
* understand basic `conda` commands
* learn how to list installed software packages 
* learn how to manage multiple installation environments

## Booting a "blank" machine

Follow [ANGUS instructions, with m1.medium](https://angus.readthedocs.io/en/2018/jetstream/boot.html), using "18.04 Ubuntu devel and docker" as the starting image you select – rather than "DIBSI 2018 workshop image".

Log in via the [Web shell or through ssh in your terminal](https://angus.readthedocs.io/en/2018/jetstream/boot.html#click-on-your-new-instance-to-get-more-information) if you are comfortable with that way now.

Note that neither RStudio nor conda are installed.

## What is bioconda?

See [the bioconda paper](https://www.biorxiv.org/content/early/2017/10/27/207092) and the [bioconda web site](http://bioconda.github.io).

Bioconda is a community-enabled repository of 3,000+ bioinformatics packages, installable via the `conda` package
manager.  It consists of a set of recipes, [like this one, for sourmash](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/sourmash/meta.yaml), that are maintained by the community.

It just works, and it's effin' magic!!

## What problems does conda (and therefore bioconda) solve?

Conda tracks installed packages and their versions.

Conda makes sure that different installed packages don't have
conflicting dependencies (we'll explain this below).

## Installing conda and enabling bioconda

Download and install conda:

```
curl -O -L https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Say "yes" to everything the installer asks, and accept default locations by pressing enter when it says "Miniconda3 will now be installed into this location". (If the prompt looks like this ">>>", then you are still within the installation process.)

When the installation is complete and the regular prompt returns, run the following command (or start a new terminal session) in order to activate the conda environment:

```
source ~/.bashrc
```

Next, enable various "channels" for software install, including bioconda:

```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
```

Try installing something:

```
conda install sourmash
```

and running it --
```
sourmash
```
will produce some output. (We'll tell you more about sourmash later.)

yay!

## Using conda

Conda is a "package manager" or software installer. See [the full list of commands](https://conda.io/docs/commands.html).

`conda install` to install a package.

`conda list` to list installed packages.

`conda search` to search packages. Note that you'll see one package for every version of the software and for every version of Python (e.g. `conda search sourmash`).

## Using bioconda

bioconda is a channel for conda, which just means that you
can "add" it to conda as a source of packages. That's what the `conda config` above does.

Note, Bioconda supports only 64-bit Linux and Mac OSX.

You can check out [the bioconda site](https://bioconda.github.io/).

### Finding bioconda packages

You can use `conda search`, or you can use google, or you can go visit [the list of recipes](https://bioconda.github.io/recipes.html#recipes).

### Freezing an environment

This will save the list of **conda-installed** software you have in a particular
environment to the file `packages.txt`:
```
conda list --export packages.txt
```
(it will not record the software versions for software not installed by conda.)

```
conda install --file=packages.txt
```
will install those packages in your local environment.

### Constructing and using multiple environments

A feature that we do not use much here, but that can be very
handy in some circumstances, is different environments.

"Environments" are multiple different collections of installed software. There are two reasons you might want to do this:

* first, you might want to try to exactly replicate a specific software install, so that you can replicate a paper or an old condition.
* second, you might be working with incompatible software, e.g. sometimes different software pipelines need different version of the same software. An example of this is older bioinformatics software that needs python2, while other software needs python3.

To create a new environment named `pony`, type:

```
conda create -n pony
```

Then to activate (switch to) that environment, type:

```
source activate pony
```

And now when you run `conda install`, it will install packages into this new environment, e.g.
```
conda install -y checkm-genome
```
(note here that checkm-genome *requires* python 2).

To list environments, type:
```
conda env list
```
and you will see that you have two environments, `base` and
`pony`, and `pony` has a `*` next to it because that's your
current environment.

And finally, to switch back to your base environment, do:

```
source activate base
```
and you'll be back in the original environment.

### Meditations on reproducibility and provenance

If you want to impress reviewers and also keep track of
what your software versions are, you can:

* manage all your software inside of conda
* use `conda list --export software.txt` to create a list of all your software and put it in your supplementary material.

This is also something that you can record for yourself, so that if you are trying to exactly reproduce 

### Using it on your own compute system (laptop or HPC)

conda works on Windows, Mac, and Linux.

bioconda works on Mac and Linux.

It does not require admin privileges to install, so you can
install it on your own local cluster quite easily.

## Bonus: installing RStudio Web on your running Linux box.

**Note:** this does require admin privileges, and you cannot run it on your local cluster. For your laptop, you can just install [the regular RStudio](https://www.rstudio.com).

----

Install necessary system software (gdebi and R):
```
sudo apt-get update && sudo apt-get -y install gdebi-core r-base
```

Now, download and install RStudio:
```
wget https://download2.rstudio.org/rstudio-server-1.1.453-amd64.deb
sudo gdebi -n rstudio-server-1.1.453-amd64.deb
```

At this point, RStudio will be running on port 8787, and you can [follow these instructions](https://angus.readthedocs.io/en/2018/visualizing-blast-scores-with-RStudio.html) to set your password and log into it.
