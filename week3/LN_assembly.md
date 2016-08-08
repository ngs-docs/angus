# Assembly using velvet

Heavily based on [material](https://github.com/lexnederbragt/INF-BIOx121_fall2014_de_novo_assembly) developed for the *de novo* assembly part of the INF-BIOx121 course "High Throughput Sequencing technologies and bioinformatics analysis" at Univ. of Oslo Fall 2014. License: CC0.


### *De novo* assembly of Illumina reads using velvet

### Installing Velvet

Start an EC instance: m4.xlarge image with ubuntu 14.04. Update the base software and install some needed software:

```
sudo apt-get update && sudo apt-get install -y g++ gcc make git zlib1g-dev python
```

Download the source code and unpack it:

```
wget https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz
tar -xvzf velvet_1.2.10.tgz
```

We'll compile velvet allowing for kmers up to 127:

```
cd velvet_1.2.10
make 'MAXKMERLENGTH=127'
export PATH=${PWD}:$PATH
cd contrib/estimate-exp_cov/
export PATH=${PWD}:$PATH
```
Check that you actually can run velvet by typing

```
cd
velveth
```

You should see the velvet help text.

### Downloading reads

```
cd
mdir assembly_module
cd assembly_module
mkdir data
cd data
wget https://www.dropbox.com/s/kopguhd9z2ffbf6/MiSeq_Ecoli_MG1655_50x_R1.fastq
wget https://www.dropbox.com/s/i99h7dnaq61hrrc/MiSeq_Ecoli_MG1655_50x_R2.fastq
```


### Assembling short-reads with Velvet

We will use Velvet to assemble Illumina reads on their own. Velvet uses the *de Bruijn graph* approach. 

We will assemble *E. coli K12* strain MG1655 which was sequenced on an Illumina MiSeq. The instrument read 150 bases from each direction.

We wil first use paired end reads only: 

### Building the Velvet Index File

Velvet requires an index file to be built before the assembly takes place. We must choose a *k-* mer value for building the index. Longer *k-* mers result in a more stringent assembly, at the expense of coverage. There is no definitive value of *k* for any given project. However, there are several absolute rules:

* *k* must be less than the read length
* it should be an odd number. 

Firstly we are going to run Velvet in single-end mode, *ignoring the pairing information*. Later on we will incorporate this information.

First, we need to make sure we can use velvet:

Create the assembly folder in a suitable location:

```
cd 
mkdir velvet
cd velvet
```

#### A first assembly

Find a value of *k* (between 21 and 99) to start with, and record your choice in this google spreadsheet: <bit.ly/ngs2015Velvet>. Run `velveth` to build the hash index (see below).

Build the index as follows:

```
velveth ASM_NAME VALUE_OF_K \  
-short -separate -fastq \  
/home/ubuntu/data/MiSeq_Ecoli_MG1655_50x_R1.fastq \  
/home/ubuntu/data/MiSeq_Ecoli_MG1655_50x_R2.fastq  
```
**NOTES** 

* Change `ASM_NAME` to something else of your choosing
* Change `VALUE_OF_K` to the value you have picked
* The command is split over several lines by adding a space, and a `\` (backslash) to each line. This trick makes long commands more readable. If you want, you can write the whole command on one line instead.

After `velveth` is finished, look in the new folder that has the name you chose. You should see the following files:

```
Log
Roadmaps
Sequences
```


The '`Log`' file has a useful reminder of what commands you typed to get this assembly result, for reproducing results later on. '`Sequences`' contains the sequences we put in, and '`Roadmaps`' contains the index you just created.

Now we will run the assembly with default parameters:

```
velvetg ASM_NAME
```

Velvet will end with a text like this:

`Final graph has ... nodes and n50 of ..., max ..., total ..., using .../... reads`

The number of nodes represents the number of nodes in the graph, which (more or less) is the number of contigs. Velvet reports its N50 (as well as everything else) in 'kmer' space. The conversion to 'basespace' is as simple as adding k-1 to the reported length.

Look again at the folder `asm_name`, you should see the following extra files:

`contigs.fa`  
`Graph`  
`LastGraph`  
`PreGraph`  
`stats.txt`

The important files are:

`contigs.fa` - the assembly itself  
`Graph` - a textual representation of the contig graph  
`stats.txt` - a file containing statistics on each contig

**Questions**

* What k-mer did you use?
* What is the N50 of the assembly?
* What is the size of the largest contig?
* How many contigs are there in the `contigs.fa` file? Use `grep -c NODE contigs.fa`


Log your results in this google spreadsheet: `bit.ly/ngs2015Velvet`


**We will discuss the results together and determine *the optimal* k-mer for this dataset.**

*********************************************
**   Intermezzo: multiple choice question  **  
*********************************************

**Advanced tip:** You can also use VelvetOptimiser to automate this process of selecting appropriate *k*-mer values. VelvetOptimizer is included with the Velvet installation.

Now run `velveth` and `velvetg` for the kmer size determined by the whole class (k = 81). Use this kmer from now on!

#### Estimating and setting exp_cov

Much better assemblies are produced if Velvet understands the expected coverage for unique regions of your genome. This allows it to try and resolve repeats. The command `velvet-estimate-exp-cov.pl` is supplied with Velvet and will plot a histogram of k-mer frequency for each node in the graph, listing k-mer frequency, and the number of count of nodes with that frequency

`velvet-estimate-exp_cov.pl ASM_NAME/stats.txt`

The output shows:

* k-mer coverage
* count of contigs (nodes in the *de Bruijn* graph) with that coverage
* series of *'s making up a histogram 

The *peak value* in this histogram can be used as a guide to the best k-mer value for `exp_cov`.

**Question:**

* What do you think is the approximate expected k-mer coverage for your assembly?

Now run velvet again, supplying the value for `exp_cov` (k-mer coverage, *not* genome coverage) corresponding to your answer:

```
velvetg ASM_NAME -exp_cov PEAK_K_MER_COVERAGE
```
**Question:**

* What improvements do you see in the assembly by setting a value for `exp_cov`?

#### Setting *cov_cutoff*

You can also clean up the graph by removing low-frequency nodes from the *de Bruijn* graph using the `cov_cutoff` parameter. Low-frequency nodes can result from sequencing errors, or from parts of the genome with very little sequencing coverage. Removing them will often result in better assemblies, but setting the cut-off too high will also result in losing useful parts of the assembly. Using the histogram from previously, estimate a good value for `cov_cutoff`.

```
velvetg ASM_NAME -exp_cov YOUR_VALUE -cov_cutoff YOUR_VALUE  
```

Try some different values for `cov_cutoff`, keeping `exp_cov` the same and record your assembly results.

#### Asking velvet to determine the parameters

You can also ask Velvet to predict the values for you:

```
velvetg ASM_NAME -exp_cov auto -cov_cutoff auto
```

**Questions:**

* What values of *exp_cov* and *cov_cutoff* did Velvet choose?
* Check the output to the screen. Is this assembly better than your best one?


_____

Source of the reads: [http://www.illumina.com/science/data_library.ilmn](http://www.illumina.com/science/data_library.ilmn), random subsampling using seqtk [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)

____


For the rest of the module, see the [2014 course repository version](https://github.com/lexnederbragt/INF-BIOx121_fall2014_de_novo_assembly/blob/master/practicals/01_Assembly_using_velvet.md) (also check out [the other material](https://github.com/lexnederbragt/INF-BIOx121_fall2014_de_novo_assembly) developed for that course).
