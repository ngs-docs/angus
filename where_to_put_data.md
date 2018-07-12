Where to put your code and data for publication
=================

Let's say, hypothetically of course, you have a colleague who has a great script they made to generate figures and results for a paper:

[![](https://i.imgur.com/1EGQ5h3.png)](https://twitter.com/eliasoziolor/status/1014553520266022912)

See: Dr. Elias Oziolor's [markdown document](https://github.com/eoziolor/fgfh_post/blob/master/scripts/combined.Rmd) 

What would you tell your colleague?

In this lesson, you will be familiar with:
* available public options for placing raw data, data products and code for publication
* how the tools work
* how you can use/interact with the tools

We will explore these options and tools by taking the following steps during this lesson:
1. Download a few reads from SRA
2. Create a script that simply installs a conda environment with `fastqc` and `multiqc`, and then runs both
3. Upload SRA data and file products output to OSF
4. Upload script and jupyter notebook to GitHub
5. Link to Zenodo
6. Create a doi (but not actually create a doi unless you want to)
7. Connect GitHub to binder to have a running instance

# The modern paper

Blogs are cool, and can be a way to share results to a wide audience when results are obtained, before publication!

* [The top 10 reasons why blog posts are better than scientific papers, by Titus](http://ivory.idyll.org/blog/2017-top-ten-reasons-blog-posts.html)

Open data science and transparency is becoming a common practice. Frequently, I see colleagues sharing pre-prints and code and data before the peer-reviewed paper is released.

[![](https://i.imgur.com/TxfI0ni.png)](https://twitter.com/balkalarre/status/940107034619289606)

[![](https://i.imgur.com/NNuA4bJ.png)](https://twitter.com/micro_marian/status/940290769021227008)

This advances everyone's understanding of your cool science!

Here are some ways to share data:

# Raw Data

The [SRA](https://www.ncbi.nlm.nih.gov/sra) is the main repository for depositing raw sequencing data

**Question:** What's the best way to go through this during a lesson, discuss and show examples? Can't really go through submission process... 

* [Steps](https://www.ncbi.nlm.nih.gov/sra/docs/submitportal/):
    1. Create [Bioproject](https://submit.ncbi.nlm.nih.gov/subs/bioproject/): . Download the batch sample metadata template and add as much information as you can about the experimental units.
    2. Create [Biosamples](https://submit.ncbi.nlm.nih.gov/subs/biosample/) (eg., for popualtion genomics or transcriptomics): Download the sample metadata template and add info for each of your samples (one per biological tissue). There is a column in SRA metadata file where you identify the bioproject (add accession that you created in Step 1). 
    3. Create [SRA entry](https://submit.ncbi.nlm.nih.gov/subs/sra/). associate the SRA run entries with the experimental units that you specify in bioproject. That is, for some samples you may have more than one SRA entry, because you might have sequence reads spread across more than one lane of sequencing. After this is the point where you will upload your files once you have been manually granted access.

Additional references:
* [How to: Submit sequence data to NCBI](https://www.ncbi.nlm.nih.gov/guide/howto/submit-sequence-data/)
* [Sequence Read Archive Toolkit](https://www.ncbi.nlm.nih.gov/sites/books/NBK158899/)
* [Downloading data from the SRA](https://www.ncbi.nlm.nih.gov/books/NBK158899/)

# Getting raw data from the SRA

Install:
```
conda install -y sra-tools 
```

Let's also create a folder in our home directory, so that we keep things organized:

```
cd ~
mkdir openScienceTutorial
cd openScienceTutorial

```

Download [example set of reads](https://www.ncbi.nlm.nih.gov/sra/?term=SRR1300523):

and extract first 1000 paired reads (they call reads "spots"):

```
fastq-dump -X 1000 --split-files --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' SRR1300523
fastq-dump -X 1000 --split-files --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' SRR1300540
fastq-dump -X 1000 --split-files --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' SRR1300380
fastq-dump -X 1000 --split-files --defline-seq '@$ac.$si.$sg/$ri' --defline-qual '+' SRR1300313
```
(Because /1 and /2 keeping track of read pairs will not be included by default, see [issue](https://github.com/ncbi/sra-tools/issues/56) and [blog](https://standage.github.io/streaming-data-from-the-sra-with-fastq-dump.html))

Don't do this now! If you want the full set of reads (will take >5 min) 
```
fastq-dump SRR390728
```

Don't do this now, either! This is the same as doing it in two steps (will take > 5 min):

```
wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR390/SRR390728/SRR390728.sra
fastq-dump SRR390728.sra
```

Advanced challenge for later (Requires scripting, using bash, Python, or R): 
* Figure out how to get the`SraRunInfo.csv` from SRA for a large dataset, e.g. [719 Bioexperiments](https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=231566) in the Marine Microbial Eukaryotic Transcriptome Sequencing Project ([MMETSP](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA231566/)). 
* Write a script to loop through the SRR id for small subset, download first 5 samples in .csv and extract using commands above. 

# Data products

These are some options for sharing data products, such as transcriptomes, genomes and annotations.

The NCBI TSA (Transcriptome Shotgun Assembly Sequence Database) and Assembly are places to deposit data products, but you must be the owner of the original data. There is submission process. Read about the submission process here:

* [TSA](https://www.ncbi.nlm.nih.gov/genbank/tsa/)
* [Genome Assembly](https://www.ncbi.nlm.nih.gov/assembly)


## Public websites for sharing data products

There are several public websites that are available for sharing your data products. These are just several options. You might have know of other methods for sharing data. What are those?

**Exercise (think pair share)** 

Is data sharing a good idea? What are some of the venues where you prefer to share data with others?

These are some options we have tried, some features of each and some considerations we feel are important.

### Summary

Location | Funding | Interface, Download/Upload | Versioning | DOI | File size/number limits | Comments
--- | --- | --- | --- | --- | --- | --- |
[zenodo.org](https://zenodo.org/) | CERN (very sustainable) | website, manual down/upload | Yes | Yes | 50 GB/repository | Interface with GitHub
[figshare.com](https://figshare.com/) | Private | website, manual down/upload | Yes | Yes| Yes | Social media metrics
[osf.io](https://osf.io/) | Non-profit | website + commandline | No | Yes | 3 GB/file, unlimited files | 

### Figshare

* Features: https://figshare.com/features
* Requires [signing up for a Figshare account](https://figshare.com/). 
* Example Figshare dataset: highly used [dataset for Carpentry Ecology lessons](https://figshare.com/articles/Portal_Project_Teaching_Database/1314459) 
* Sabah Ul-Hassan's [Twitter decision-making flowchart](https://figshare.com/articles/Academic_Twitter_v3_0/6166361)
* Lisa Johnson's [MMETSP transcriptome re-assemblies](https://figshare.com/articles/Marine_Microbial_Eukaryotic_Transcriptome_Sequencing_Project_re-assemblies/3840153)
(although, I don't recommend maintainig a large multi-file dataset like this to Figshare)

### Zenodo
 
* Features: http://help.zenodo.org/
* Get Zenodo account: https://zenodo.org/ 
* a European general purpose open-access open-data repository optimized for sharing big data. 
* Integrated with GitHub to make code hosted in GitHub citable. 
* you can link to your GitHub or ORCiD account

Examples of repositories: 
* [Example 1](https://zenodo.org/record/1212585#.W0KDqxJKiL8)
* [Example 2](https://zenodo.org/record/257410#.W0KHcRJKiL8)

### OSF 

[Open science framework, operated by the non-profit, COS (Center for Open Science)](https://osf.io/4znzp/wiki/home/). Repositories for projects, include files and data products. Like GitHub, but only for file sharing. 5 GB/file size limit. Okay for raw fastq NGS data, although not if files are > 5GB.

Features:

* Great for classes, where you want to share many files with whole class
* Can incorporate download/upload in pipelines
* Can use this to upload data products before deleting Jetstream instances

Workflow for OSF client, something like this (have students make their own repository? Just listing an example project here, can make a new one for tutorial?):

1. Get [OSF account](https://osf.io/)
2. Make your own [repository](https://osf.io/b972e/) and set it to "public" (otherwise it won't allow download! Trust me. Took me a while to figure this out.)
3. Enable DOI to get citation, such as:

Johnson, L., & Psomopoulos, F. E. (2018, July 11). DIBSI2018. Retrieved from osf.io/gweqv

4. install [osf client](https://github.com/osfclient/osfclient), see [documentation](http://osfclient.readthedocs.io/en/stable/cli-usage.html)

```
pip install osfclient
```

5. Set username and password as variables:
```
export OSF_PASSWORD=
export OSF_USERNAME=
```

6. Download (clone) files from project, where address is https://osf.io/gweqv/ and Project = `gweqv`
```
osf -p gweqv clone
mv gweqv/osfstorage/scripts/ .
mv gweqv/osfstorage/Nematostella_annotation_files/ .
rm -rf gweqv
```
 
7. Switch to your own OSF project. (Substitute the `gweqv` below with your own project.

Upload one file:
```
osf -p c upload SRR1300540_2.fastq reads/SRR1300540_2.fastq
```

8. Upload many files (upload privaleges will only be enabled from those listed as "Contributors":
```
cd ~/openScienceTutorial
osf -p gweqv upload -r *.fastq reads/
```

Upload file products (remember to substitute your own project code in the command):
```
osf -p gweqv upload -r Nematostella_annotation_files/* annotation_files/

```


# Code - GitHub / Zenodo


Now that we have uploaded our input data (i.e. the two `SRA` files and the `dammit` output), we can add our scripts on a GitHub repository so that we can have the entire process available and linked.

Specifically, we will be uploading the script that does the QC for the `fastq` files as well as the `jupyter` notebook from dammit. Both files are also available here:

- [runQC.sh](https://osf.io/589xw/download)
- [jupyter notebook](https://osf.io/k2jhg/download)

```
cd ~/openScienceTutorial
wget -O runQC.sh https://osf.io/589xw/download
wget -O nema_annotation.ipynb https://osf.io/k2jhg/download
```

If all has been set up correctly, you should have the following structure in your `openScienceTutorial` folder:

```
dibtiger@js-170-21:~/openScienceTutorial$ ls -la

total 1056
drwxrwxr-x  2 dibtiger dibtiger   4096 Jul 11 13:12 .
drwx------ 27 dibtiger dibtiger   4096 Jul 11 13:08 ..
-rw-rw-r--  1 dibtiger dibtiger 128893 Jul 11 13:08 SRR1300313_1.fastq
-rw-rw-r--  1 dibtiger dibtiger 128893 Jul 11 13:08 SRR1300313_2.fastq
-rw-rw-r--  1 dibtiger dibtiger 128893 Jul 11 13:08 SRR1300380_1.fastq
-rw-rw-r--  1 dibtiger dibtiger 128893 Jul 11 13:08 SRR1300380_2.fastq
-rw-rw-r--  1 dibtiger dibtiger 128893 Jul 11 13:08 SRR1300523_1.fastq
-rw-rw-r--  1 dibtiger dibtiger 128893 Jul 11 13:08 SRR1300523_2.fastq
-rw-rw-r--  1 dibtiger dibtiger 128893 Jul 11 13:08 SRR1300540_1.fastq
-rw-rw-r--  1 dibtiger dibtiger 128893 Jul 11 13:08 SRR1300540_2.fastq
-rw-rw-r--  1 dibtiger dibtiger  11342 Jul 11 13:12 nema_annotation.ipynb
-rw-rw-r--  1 dibtiger dibtiger    163 Jul 11 13:12 runQC.sh
```


Run the `runQC.sh` script in the directory where you've downloaded the .fastq reads:

```
bash runQC.sh
```

**Independent Challenge!** 

Use git to push our code in the `sh` and `ipynb` files into a GitHub repository. You can use the same process as listed [here](http://angus.readthedocs.io/en/2018/rmarkdown_rnaseq.html#version-control-with-git-and-github).

1. Create GitHub repository in your account, e.g. http://www.github.com/<username>
2. Copy the URL with the green "Clone or download" button, e.g. https://github.com/<username>/<reponame>.git
3. Then copy the code `sh` and `ipynb` files into the directory:
```
git clone https://github.com/<username>/<reponame>.git
cd <reponame>
cp ~/openScienceTutorial/<filename1> .
cp ~/openScienceTutorial/<filename2> .
git add --all
git commit -m "initial commit"
git push origin master
```

3. Check https://github.com/<username>/<reponame> to see if your files appeared on GitHub

4. Make a change to your code file. Add a line or two. Repeat steps above to version control your code:

```
git add --all
git commit -m "changed commit"
git push origin master
```

4. Now, turn on feature which archives GitHub code in Zenodo, get a DOI and freezes code as a snapshot in time for publication.

Read about how to do this [here](https://guides.github.com/activities/citable-code/).
    
Read about why [GitHub + Zenodo makes your life great](http://ivory.idyll.org/blog/2016-using-zenodo-to-archive-github.html)!

# Binder

Now, we can link our GitHub repository to Binder. This will create an active version of the jupyter notebook in your GitHub repository:

[Binder](https://mybinder.org/)

**Challenge:**

Using the link above, how would you connect your GitHub repository to binder? 


# More information

* Digital Object Idntifiers (DOIs) and publication - [DOI](https://www.doi.org/faq.html), you can cite them.
* Licensing your code - There are different licenses. See Titus' blog: http://ivory.idyll.org/blog/2015-on-licensing-in-bioinformatics.html
