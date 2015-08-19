========================================
Genome comparison and phylogeny
========================================

This tutorial will introduce genome comparison techniques and some simple methods to compute phylogeny.
It will introduce the following pieces of software: 

`Mauve <http://gel.ahabs.wisc.edu/mauve>`__, for genome alignment
`PhyloSift <http://phylosift.wordpress.com>`__ & `FastTree <http://meta.microbesonline.org/fasttree/>`__, for phylogeny

We'll analyze some *E. coli* genome assemblies that were precomputed with techniques
previously introduced in the course.

Interactive visual genome comparison with Mauve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download a copy of the Mauve GUI installer for your platform:
`Mac OS X <http://gel.ahabs.wisc.edu/mauve/snapshots/2012/2012-02-03/MacOS/Mauve-snapshot_2012-02-03.dmg>`__
`Windows <http://gel.ahabs.wisc.edu/mauve/snapshots/2012/2012-03-03/windows/mauve_installer_20120303.exe>`__
`Linux <http://gel.ahabs.wisc.edu/mauve/snapshots/2012/2012-06-07/linux-x64/mauve_linux_snapshot_2012-06-07.tar.gz>`__

Download the following GenBank format genomes from NCBI:
`E. coli O157:H7 EDL933 <ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_O157_H7_EDL933_uid57831/NC_002655.gbk>`__
`E. coli CFT073 <ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_CFT073_uid57915/NC_004431.gbk>`__

An assembly of the Ribiero et al 2012 *E. coli* MiSeq data (accession SRR519926) made by the `A5-miseq pipeline <http://sourceforge.net/projects/ngopt/>`__ (`paper here <http://arxiv.org/abs/1401.5130>`__):

`E. coli A5-miseq assembly <http://edhar.genomecenter.ucdavis.edu/~koadman/ngs2014/ecoli_a5.final.scaffolds.fasta>`__

An assembly of the Chitsaz et al 2011 data made by the Velvet diginorm workflow:

`E. coli Velvet DN assembly <http://edhar.genomecenter.ucdavis.edu/~koadman/ngs2014/ecoli_dn_velvet25.fa>`__


Running a genome alignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Launch the Mauve software, either from the start menu on Windows, or by double-clicking the Mauve app on Mac.
On newer Mac OS it may be necessary to disable the nanny security features that prohibit opening downloaded software.
If you see error messages suggesting the disk image is damaged and can not be opened then you need to disable the security check.

Once Mauve is open, go to the File menu, select "Align with progressiveMauve..."
Drag & drop the NCBI genomes in first.
Then drag & drop in the draft *E. coli* assemblies

The assembly will take 20+ minutes to compute. Proceed to the next steps while this process runs.


Booting an Amazon AMI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Start up an Amazon computer (m1.large or m1.xlarge) using an ubuntu 64-bit linux instance 
(see :doc:`amazon/start-up-an-ec2-instance` and :doc:`amazon/starting-up-a-custom-ami`).

Log in `with Windows <amazon/log-in-with-ssh-win.html>`__ or
`from Mac OS X <amazon/log-in-with-ssh-mac.html>`__.

Ensure you have `dropbox installed <amazon/installing-dropbox.html>`__ on your virtual machine 
and mounted at /mnt/Dropbox


Logging in & updating the operating system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Copy and paste the following three commands
::

   sudo add-apt-repository ppa:webupd8team/java
   sudo apt-get update
   sudo apt-get install oracle-java6-installer
   sudo apt-get remove openjdk-7-jre openjdk-7-jre-headless

to update the computer with all the bundled software you'll need.
Mauve, when run in command-line mode, depends on the Oracle Java 6 runtime which the above commands install.
The above commands will also remove the openjdk if it's present on the system. If you would prefer to keep
this for some reason then look at the ``update-alternatives`` system.

Packages to install
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Install the latest `Mauve snapshot <http://gel.ahabs.wisc.edu/mauve/snapshots/>`__ to your home directory: ::

   cd
   curl -O http://gel.ahabs.wisc.edu/mauve/snapshots/2012/2012-06-07/linux-x64/mauve_linux_snapshot_2012-06-07.tar.gz
   tar xzf mauve_linux_snapshot_2012-06-07.tar.gz
   export PATH=$PATH:$HOME/mauve_snapshot_2012-06-07/linux-x64/

as well as the `PhyloSift <http://phylosift.wordpress.com>`__ software,
for marker-based analysis of genome & metagenome phylogeny: ::

   cd
   curl -O http://edhar.genomecenter.ucdavis.edu/~koadman/phylosift/releases/phylosift_v1.0.1.tar.bz2
   tar xjf phylosift_v1.0.1.tar.bz2


Getting the E. coli genome data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now, let's create a working directory::

   cd
   mkdir draft_genome
   cd draft_genome


Download some genome assemblies.  The first one is an assembly constructed with A5-miseq.
The second is an assembly constructed with Velvet. ::

   curl -O http://edhar.genomecenter.ucdavis.edu/~koadman/ngs2014/ecoli_a5.final.scaffolds.fasta
   curl -O http://edhar.genomecenter.ucdavis.edu/~koadman/ngs2014/ecoli_dn_velvet25.fa



What is the nearest reference genome?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this step, we will use the PhyloSift software to determine where our newly assembled genome
falls in the tree of life, and by extension what closely related genomes might be available for 
comparative genomics. In principle this would be a first step after sequencing a novel organism ::

   ~/phylosift_v1.0.1/bin/phylosift all --debug ecoli_dn_velvet25.fa
   
When that finishes (probably after a very long time), copy the results to dropbox for local viewing: ::

   cp -r PS_temp/ecoli_dn_velvet25.fa /mnt/Dropbox

once on your local computer, open the .html file in the ecoli_dn_velvet25.fa directory in your web browser.
It should be possible to launch the concat marker viewer in java from the link in the bottom left of the page.

Ordering the assembly contigs against a nearby reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the previous step we discovered that our genome was similar to E. coli.
We will now use an existing finished-quality E. coli genome as a reference for ordering
and orienting the contigs in the assembly. 
This is useful because the genome assembly process usually creates a large number of contigs or scaffolds rather
than complete reconstructions of the chromosome(s) and these sequences appear in an arbitrary order.
By ordering against a reference we can generate a candidate ordering which could be used for later manual closure efforts (e.g. via PCR) or other analyses. 
Let's use the *E. coli* K12 genome from NCBI as a reference: ::

   curl -O ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__DH10B_uid58979/NC_010473.fna

And now run the Mauve Contig Mover to order the contigs: ::

   java -Xmx500m -Djava.awt.headless=true -cp ~/mauve_linux_snapshot_2012-06-07/Mauve.jar  org.gel.mauve.contigs.ContigOrderer -output reorder_a5 -ref NC_010473.fna -draft ecoli_a5.final.scaffolds.fasta
   java -Xmx500m -Djava.awt.headless=true -cp ~/mauve_linux_snapshot_2012-06-07/Mauve.jar  org.gel.mauve.contigs.ContigOrderer -output reorder_a5 -ref NC_010473.fna -draft ecoli_dn_velvet25.fa

Let's copy the newly ordered genome to dropbox: ::

   cp reorder_a5/alignment3/ecoli_a5.final.scaffolds.fasta /mnt/Dropbox/

This could now be aligned with Mauve as demonstrated above to observe the improvement in contig order.

Making a phylogeny of many E. coli assemblies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For this component we will use a collection of related *E. coli* and *Shigella* genomes already on NCBI.
In practice, you might use your own collection of assemblies of these genomes.
Let's start out in a new directory and download these: ::

   mkdir ~/phylogeny ; cd ~/phylogeny
   
   # download a bunch of genomes from NCBI. Alternatively you can use the approach that 
   # Adina introduced in the previous lesson to programmatically download many files
   curl -O ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_O157_H7_EDL933_uid57831/NC_002655.fna
   curl -O ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_CFT073_uid57915/NC_004431.fna
   curl -O ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_K_12_substr__DH10B_uid58979/NC_010473.fna
   curl -O ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Escherichia_coli_O104_H4_2011C_3493_uid176127/NC_018658.fna
   curl -O ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Shigella_flexneri_2a_2457T_uid57991/NC_004741.fna
   curl -O ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Shigella_boydii_Sb227_uid58215/NC_007613.fna
   curl -O ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Shigella_dysenteriae_Sd197_uid58213/NC_007606.fna
   curl -O http://edhar.genomecenter.ucdavis.edu/~koadman/ngs2014/ecoli_a5.final.scaffolds.fasta
   curl -O http://edhar.genomecenter.ucdavis.edu/~koadman/ngs2014/ecoli_dn_velvet25.fa

   # find homologs of the elite marker genes
   find . -maxdepth 1 -name "*.fna" -exec ~/phylosift_v1.0.1/bin/phylosift search --isolate --besthit {} \;
   ~/phylosift_v1.0.1/bin/phylosift search --isolate --besthit ecoli_a5.final.scaffolds.fasta

   # align to the marker gene profile HMMs
   find . -maxdepth 1 -name "*.fna" -exec ~/phylosift_v1.0.1/bin/phylosift align --isolate --besthit {} \;
   ~/phylosift_v1.0.1/bin/phylosift align --isolate --besthit ecoli_a5.final.scaffolds.fasta

   # combine the aligned genes into a single file
   find . -type f -regex '.*alignDir/concat.codon.updated.1.fasta' -exec cat {} \; | sed -r 's/\.1\..*//'  > codon_alignment.fa

   # infer a phylogeny with FastTree
   ~/phylosift_v1.0.1/bin/FastTree -nt -gtr < codon_alignment.fa > codon_tree.tre

   # now copy the tree over to dropbox
   cp codon_tree.tre /mnt/Dropbox/

From tree file to figures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At last we have a phylogeny! The last steps are to view it, interpret it, and publish it.
There are many phylogeny viewer softwares, here we will use FigTree. You will need to 
`download and install <http://tree.bio.ed.ac.uk/software/figtree/>`__ FigTree to your computer.
Once installed, either launch by double-click (Mac) or via the start menu (Windows).
Now we can open the tree file codon_tree.tre from dropbox.

Once open, enable the node labels which show bootstrap confidence.
Optionally midpoint root the tree, adjust the line width, and export a PDF.


This document (c) copyleft 2014 Aaron Darling.
