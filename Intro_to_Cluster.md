# Working on a computing cluster - DIBSI 2019
This tutorial is intended to be a general guideline towards researchers who are new to working on remote computers and fairly new to working on the command line.
Because cluster access is specific to the instituion granting you access, this tutorial has limited interactivity.

We will cover:
* Cluster structure overview
* Logging onto the cluster
* Working on the cluster
* Using job schedulers
* Finding resources on the cluster

 *This tutorial was written by Chissa Rivaldi and Shannon Joslin and was last updated July 2019*

---

## What is a "Cluster"
A cluster can be thought of as a group of computers which work together to allow you to log onto one computer (**head node**) and use resources to perform memory intensive functions from other connected computers.

<center><img src="https://i.imgur.com/2nl5zzP.png" width="100%"></a></center>
###### Image modified from (http://www.vrlab.umu.se/documentation/guides/beginner-guide)



## Logging on and using the cluster
1. **Open your terminal** (Windows - Bash, Git Bash, Cygwin, etc.)
2. Use the command `ssh` to **log on to a remote** computer.
You need to provide the command with an address to complete the connection:
    ```
    ssh <username>@<computer.address>.edu
    ```
The command you will use will look very similar to this. Sometimes there may be a few flags the sys admin will suggest you use when logging in––follow their instructions! The username and password will be set when you start an account. Sometimes your cluster admins will require a training session before your credentials will be granted access to the machines.

3. **Enter your password** and hit `enter`
This is the password you entered when you generated your rsa_key or the one you gave to the admins who set up your account. When you enter your password you'll probably not see your cursor move as you're typing - this is a little odd the first few times you do it, but rest assured your keys are being registered.

When you have successfully logged in, you'll probably see a lot of text that gives you some important information--READ THIS!--like when maintenence on the machines will take place, information about which machines are appropriate for what kind of work, and probably some information about space allocation. This will vary between institutions. It might look something like this:

<center><img src="https://i.imgur.com/JdDspJe.png" width="100%"></a></center>
 Also, notice when you log onto the cluster your prompt will change:

 <center><img src="https://i.imgur.com/igIqHiv.png" width="100%"></a></center>

When you log onto the cluster you automatically start interacting on the login/head node. **This is a computer that is shared by all the users** who have logged onto the cluster you are using. The login/head has a limited amount of resources (memory & RAM) dedicated to it. As such, make sure you only run commands that don't take up much memory like navigating around `cd`, looking at files with `less`, creating scripts or looking through directories `ls` on the head node.



4. **Start computing!**
Try some of the basic commands (such as `ls` or `pwd` or `cd`) as we've done previously in the command line!
Your cluster may run a different shell than `bash`, which we've been using throughout the workshop. To check which shell you are using type:
```
echo $0
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;If the output from your command is something other than bash (or /bin/bash), type `bash` to code in bash again!




5. **Warnings** for working on a remote computer

    - Do not use `sudo` on remote clusters- this is not your computer and the people who run this computer will send you an irritable email if you try to use sudo to force anything to install or delete, etc. Luckily, you almost never need to use sudo in a bioinformatics workflow.
    - Do not upload a lot of data to this computer unless you have a lot of space. Your space allocation will generally be fixed, and you will have issues of varying irritability if you run out. If your cluster has scratch space available to you, ask for an allocation of that to *usually temporarily* gain more space.
    - **Do not run memory-intensive programs on the login/head node**. Depending on your institution you may have _testing nodes_ and/or _interactive nodes_. You can make sure your code is working before you submit a job on these nodes. When you run big programs or use memory intensive commands (such as copying large files) on the **head node**, you consume A LOT of the resources and slow down the machines for everyone using them. So make sure you keep it simple on the head node!
    - Always back up your data and scripts! Occasionally clusters can have catastrophic failures where all or some data can be lost. Becuase of this chance, it is imperative to backup your data consistently!


6. Getting files on and off a remote space:
    - `scp` - Instructions -> https://kb.iu.edu/d/agye
    Note that some clusters require specific **transfer nodes** to transfer data to and from the cluster. To transfer from a specific port use:
    ```
    scp -p <port_number> <item(s)_to_copy> <destination>
    ```
    - `sftp` - Instructions -> https://hackmd.io/s/rybzCZasX#
        - bonus instructions for filezilla included in above link
    - `git` - clone/merge things into/from a remote repo. Not super efficient but works for most purposes.
    - **cyberduck** - User-friendly GUI. Slower than the other two, but useful when you're getting started: https://download.cnet.com/windows/cyberduck/3260-20_4-6244900-1.html

7. Modules
The admin who maintain the cluster will have some programs available for you to use. In order to use them you will need to load them into your cluster environment. The way these are organized may differ, but we can gain access to what is available using:
`module avail`
Usually this list is pretty long. We can use pipes with `grep` to speed up finding specific programs or modules we are interested in. For example:
`module avail | grep -i samtools`
This will only list modules that have samtools.
To **load** a module enter the following:
`module load samtools/1.9`
To **list** the modules you have loaded enter:
`module list`

| **Term** | **Definition** |
| -------- | --------|
| `module avail` | lists all available modules in `less    ` |
| `module list` | lists loaded modules |
| `module load <module>` | loads selected module into the shell environment |
| `module unload <module>` | unloads selected module |





# Job Scheduling
In order to carry out commands that are memory intensive we need to use auxillary computers that will not affect the login/head node. As a note, sometimes merely copying large files is memory intensive enough that we will need to use other computers! To request resources to run our scripts we use job schedulers. The purpose of a job scheduler is as follows:
* to allocate resources requested by the user to carry out scripts submitted to the scheduling system


There are a number of different flavors of job schedulers specific to the cluster you will use through your institution but they all have the following general structure:

![](https://i.imgur.com/9rSbIxR.png)



The commands you use to write and submit scripts on a cluster depending on which [scheduler](https://en.wikipedia.org/wiki/Job_scheduler) your institution uses to run the users' jobs. A couple of commons ones are UGE and SLURM, which we detail here.

We will go through two different kinds of job schedulers:
1. UGE (Univa grid engine)
2. slurm workload manager


## UGE  [(Grid Engine)](http://www.univa.com/products/)

This is one of two examples of job script schedulers.

### UGE setup for a job script
The first line directs the machine to run the rest of the lines using bash. This might differ depending on how the machine or scheduler has been set up.

```
#!/bin/bash
```

The next chunk of lines gives the machine specific information about how to run your job. Some of these parameters won't change from job to job (like your email), but they some do, like the resources to request if you're running a job that can be run in parallel. If you're not sure whether your job does this, you're most likely not ready to change this parameter without discussing it with someone more familiar with the cluster on which you're working.
The `-q` flag specifies where the scheduler should put your job. In my institution, we have a long queue for jobs that will take longer than 4 hours. There's also a short queue, and some more that are owned by research groups and are only available to users specified by that group.

```
#$ -M me@email.edu        # Email address for job notification
#$ -m abe                 # Send mail when job begins (b), ends (e) and aborts (a)
#$ -q long                # specify queue
#$ -N vsearch             # Specify job name
#$ -pe smp 1              # number of cores requested for parallel jobs processors
```

The modules you load will vary depending on which software you need to access to run your program. Their names will vary between institutions. If you're accessing something you've installed yourself on your space, you don't necessarily need to load any modules. If your job fails and the resulting error message indicates that a program is missing or otherwise unaccessible (and you usually have no issue using that particular program), the lack of a loaded module may be the reason behind the error.

```
module load xyz        # Required modules to load
```

### UGE job script continued
Now that your script has been setup and the machine has all the proper instructions, we can get started working with the data. The code you want to submit is generally no different than what you would put into a normal shell script. Here's a short example:

```
#cluster filtered reads at 97% - from previous output using the program Vsearch
vsearch --cluster_fast filtered.fasta --id 0.97 --centroids filtered.centroids
```

### UGE Whole file (for copy/paste convenience)
Example of what this all look like put together.

For UGE:
```
#!/bin/bash

#$ -M me@email.edu     # Email address for job notification
#$ -m abe                 # Send mail when job begins (b), ends (e) and aborts (a)
#$ -q long                # specify queue
#$ -N vsearch             # Specify job name
#$ -pe smp 1              # number of cores requested for parallel processors

module load xyz        # Required modules to load

#cluster filtered reads at 97% - from previous output using the program Vsearch
vsearch --cluster_fast filtered.fasta --id 0.97 --centroids filtered.centroids
```

### Interactive jobs
You can log directly into a node if you need to interact with your workflow. This command starts with `qrsh` and you can specify the amount of time and resources you need. You can run this command on the login node and your job will start immediately (if those resources are available).

```
#example for submission for an interactive job with one node and eight cores
loginnode$ qrsh -pe smp 8
```
### Useful Commands
`qsub` -- submits scripts to cluster
`qstat`-- checks the status of your job (qw - waiting in queue, r - running) use ```-u ,username>``` to only look at your jobs
`qdel jobID`-- abort job
Check the man pages (```man qsub```) to find out more about how to use these commands!


## SLURM [workload manager](https://slurm.schedmd.com/documentation.html)

slurm is an open source workload manager that is commonly used on compute cluster. It handles allocating the remote compute cluster's resources to users requesting resources

There are two main ways you can request resources using slurm:
#### 1.Run interactive sessions with `srun`

Interactive sessions allow you to work on computers that aren't the login/head node. Essentially you can do everything you've done at the command line interface on Jetstream on the compute cluster. This is really powerful for doing memory intensive commands that you may not need to keep track of. However, with this power comes a great danger as the commands you run will not be save in a script anywhere. So, if you wanted to go back and recreate an analysis, you won't know what you've run or with which versions of software.

To request and launch a basic interactive session that will last for two hours use the following:
```
srun --time=02:00:00 --pty /bin/bash
```

Pay close attention to the time you give to yourself using `srun`! Slurm will terminate the session immediately at the end of the allotted time. It, sadly, doesn't care if you are 99.99% of the way through your analysis :0]

Also, you can request more/different resources by using to following flags:
* `--mem=<number>Gb` request memory
* `-c <number>` request a certain number of CPUs


#### 2.Submit batch scripts with `sbatch`

Batch job scripts (also known as job scripts) are scripts that contain `#!/bin/bash` at the beginning of each script and are submitted to the slurm workload manager by using `sbatch`. When we submit a script to slurm it is considered a _job_ and gets a unique job ID assigned to it.

There are a few parameters which are required before the workload manager can & will accept the jobs we submit to run.

Required input:
* the **time** we request to run our job/submitted script. A job scheduler will not accept a job without a time parameter. We can request time with the following flag: `--time=01-02:03:04` Note: the format is dd-hh:mm:ss
* the **partition** we would like to use for our job––oftentimes this will also entail the _priority_ in which our job is submitted. We can request a partition by using the following flag: `-p <name_of_partition>`
* the **memory** required to run our job. We can request a specified amount of time with the following flag: `--mem=<number>Gb`

Optional input:
* we can have slurm **mail** us updates about our job, such as when it starts(`BEGIN`), ends(`END`), if it fails(`FAIL`) or all of the above (`ALL`). We can request slurm emails us with the following flag: `--mail-user=<your_email> --mail-type=ALL`
* we can also give jobs specific **names**. To name your job use: `-J <job_name>` Be careful, as there is a limit to the number of characters your job name can be.
* slurm automatically generates output scripts where all of the output from commands run from the script are printed to. These will take the form as `slurm12345.out` where 12345 is the unique identifying number slurm assigns to the file. We can change this to any output file name we want. To specify the name of your output file use `-o <file_name>.out`
* slurm can generate error files, where all of the errors from the script are printed to. We ask slurm to create err files and name them with `-e <file_name>.err`

If we were hard to ourselves we would write these out at the command line each time we submitted a job to slurm with `sbatch`. It would look something like this:
```
sbatch --time=01-02:03:04 -p <name_of_partition> --mem=<number>Gb --mail-user=<your_email> --mail-type=ALL -J <job_name> -o <file_name>.out -e <file_name>.err
```
We will ned to switch out the <text> with parameters specific to our preference, but hopefully you get the gist.

However, typing out the command above is quite lengthy and we want to be nice to our current and futureselves. So, we can put these run parameters into the beginning of our batch job scripts with the following **example `sbatch` script**:
```
#!/bin/bash
#
#SBATCH --mail-user=<your_email>    # Your email
#SBATCH --mail-type=ALL             # Options: ALL, NONE, BEGIN, END, FAIL, REQUEUE
#SBATCH -J fastQC                   # Job name (limit character usage)
#SBATCH -e fastQC.err               # batch script's standard error file
#SBATCH -o fastQC.out               # batch script's standard output file
#SBATCH --time=01-02:03:04          # time allocated to this job FORMAT: dd-hh:mm:ss
#SBATCH --mem=16Gb                  # memory allocated to each node
#SBATCH --partition <partition>     # partition to request resources from

set -e # exits job upon running into an error
set -x # print commands and their arguments as they are executed.


module load fastqc                  # Required modules to load (you 

fastqc <file>
```
Generally job scripts end in `.sh` and occasionally `.slurm`

### Other Useful Commands

**Check job status with `squeue`**
So we submit jobs to slurm but where do they go and what are they doing?? `squeue` allows you to check on the status of all job(s) submitted to slurm. You can check on the status of your job by entering `squeue -u <username>`
This will print out a list of all of the jobs you have currently in the queue or running with slurm.
You will see the

**Get info on resources availiable through slurm with `sinfo`**
How do we know what resources are available to us through slurm?? `sinfo` lets us look at the `nodes` and `partitions` that slurm has access to on the compute cluster. It will also tell us useful things like what the time limit for particular nodes are and what all the nodes are up to. Each cluster will have different outputs but some common ones are:
* `PARTITION` = the partition to which the nodes belong to
* `TIMELIMIT` = the timelimit for jobs submitted to a particular node
* `NODES` = the number of nodes on a particular partition at a particular state
* `STATE` = this lists the current state of the nodes
* `NODELIST` = this lists the node names that are currently at the same state on the same partition

**Cancel jobs with `scancel`**
We can cancel any and all jobs that we currently have submitted to slurm. 

To cancel a particular job, we will need the job number (remember you can find this by using `squeue -u <username>`). Then, we can cancel the job with 
```
scancel <job_number>
```

Additionally, we can cancel ALL of our jobs by typing:
```
scancel -u <username>
```



# Remote  Vocabulary & Thesaurus

| **Term** | **Definition** |
| ----------- | --------------|
|cluster| group of computers that is managed by institution |
| job | a program or series of programs found in a script submitted to the job scheduler |
| node | the functional unit of resources in a cluster|



| **Synonyms** | for all intensive purposes these are the same! |
| ---------- | ----------------- |
| login node | head node|
| cluster | remote cluster |
| batch job script | batch script | job script |
| job scheduler | workload manager |


# Other topics not covered here that you might come across

#### AFS -
This is a file system that some clusters use. Functionally, this will matter if you attempt to change permissions to a shared directory (e.g. your lab's working or storage space)
##### More details : https://web.stanford.edu/~consult/afsinfo/basic.shtml

#### Scratch space -
Extra storage space allocated by cluster admins. Should not be treated as permanant storage. Contact your admins to see what your options are here.

#### Software installation -
This varies quite a lot between and even within system. Some clusters allow you to install software on your own, some don't. Combined with the multiple ways you install different bits of software, this can get complicated. Best to ask the admins what their recommentation is, and follow from there. If you do have the ability to install software on your own space, you will (generally) have the easiest time using conda. Your cluster might already have this available for you to use, but if not it's a pretty quick install. If conda is not an option for you, the software installation documentation should be the first place you look (look for a file called INSTALL or README). When downloading source code and installing it, the command `./configure` allows you to set a prefix which is likely going to be your home directory. You command should look something like this:
```./configure prefix=/username/local/```
but will certainly change depening on your directory setup and preferences.
**Hard rule**: do not use `sudo` to force an install! Ever!
