# Running large and long command line jobs - using shmlast!

Our goal for this tutorial is for you become more familiar with running longer programs
on the command line. You'll be introduced to
[shmlast](http://joss.theoj.org/papers/3cde54de7dfbcada7c0fc04f569b36c7), which is implements
an algorithm for discovering potential orthologs between an RNA-seq assembly and a protein database.


## Installing shmlast

Install base packages:

```
sudo apt-get -y update && \
sudo apt-get install -y python3.5-dev python3.5-venv make \
    libc6-dev g++ zlib1g-dev last-align parallel
```

Then create a Python environment with `virtualenv`, which will isolate your python packages:
```
python3.5 -m venv ~/py3
. ~/py3/bin/activate
pip install -U pip
```

And now install shmlast 1.2:
```
pip install shmlast>=1.2
```

This downloads and installs the latest version of shmlast.

## Download some data

Next we need some data! Here we're going to grab one of the three
mouse RNA data sets,

```
curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.1.rna.fna.gz
```

and all 8 of the cow protein data sets.

```
for i in 1 2 3 4 5 6 7 8
do
   curl -O ftp://ftp.ncbi.nih.gov/refseq/B_taurus/mRNA_Prot/cow.$i.protein.faa.gz
done
```

shmlast wants one query database (here, we'll use mouse) and one
database to be searched (here, cow) - but first we have to combine
all of the databases into one:

```
gunzip -c cow.*.faa.gz > cow.faa
```

## Run shmlast!

Now run shmlast:
```
shmlast crbl -q mouse.1.rna.fna.gz -d cow.faa --n_threads=6
```
this will take 16 minutes (!!) and produce some large files.

----

Now, download this file & look at in a spreadsheet.  You should see:

```
E,EG2,E_scaled,ID,bitscore,q_aln_len,q_frame,q_len,q_name,q_start,q_strand,s_aln_len,s_len,s_name,s_start,s_strand,score
6.6e-24,9.8e-16,23.18045606445813,641897,109.65804469295703,89,1,390,"ref|NM_001013372.2| Mus musculus neural regeneration protein (Nrp), mRNA",64,+,89,389,ref|XP_005212262.1| PREDICTED: DNA oxidative demethylase ALKBH1 isoform X1 [Bos taurus],0,+,241.0
5.4e-194,4.4e-165,193.26760624017703,719314,605.7589445367834,313,0,331,"ref|NM_207235.1| Mus musculus olfactory receptor 358 (Olfr358), mRNA",0,+,313,313,ref|XP_607965.3| PREDICTED: olfactory receptor 1361 [Bos taurus],0,+,1365.0
2.8e-188,5e-160,187.5528419686578,423289,588.9868500580775,307,0,323,"ref|NM_146368.1| Mus musculus olfactory receptor 361 (Olfr361), mRNA",0,+,307,313,ref|XP_607965.3| PREDICTED: olfactory receptor 1361 [Bos taurus],0,+,1327.0
6.6e-183,5.6e-155,182.18045606445813,725159,572.2147555793716,307,0,318,"ref|NM_146622.1| Mus musculus olfactory receptor 360 (Olfr360), mRNA",0,+,307,313,ref|XP_607965.3| PREDICTED: olfactory receptor 1361 [Bos taurus],0,+,1289.0
5.4e-194,4.4e-165,193.26760624017703,719315,605.7589445367834,313,0,331,"ref|NM_207235.1| Mus musculus olfactory receptor 358 (Olfr358), mRNA",0,+,313,313,ref|XP_002691614.1| PREDICTED: olfactory receptor 1361 [Bos taurus],0,+,1365.0
2.8e-188,5e-160,187.5528419686578,423290,588.9868500580775,307,0,323,"ref|NM_146368.1| Mus musculus olfactory receptor 361 (Olfr361), mRNA",0,+,307,313,ref|XP_002691614.1| PREDICTED: olfactory receptor 1361 [Bos taurus],0,+,1327.0
6.6e-183,5.6e-155,182.18045606445813,725160,572.2147555793716,307,0,318,"ref|NM_146622.1| Mus musculus olfactory receptor 360 (Olfr360), mRNA",0,+,307,313,ref|XP_002691614.1| PREDICTED: olfactory receptor 1361 [Bos taurus],0,+,1289.0
4.8e-183,5.6e-155,182.3187587626244,373474,572.2147555793716,266,0,310,"ref|XR_001782298.1| PREDICTED: Mus musculus predicted gene 4786 (Gm4786), misc_RNA",29,+,266,266,ref|NP_001035610.1| 60S ribosomal protein L7a [Bos taurus],0,+,1289.0
3.2e-153,3.1e-138,152.4948500216801,643504,516.6020212552417,246,1,659,"ref|NR_003628.1| Mus musculus predicted gene 5766 (Gm5766), non-coding RNA",357,+,246,266,ref|NP_001035610.1| 60S ribosomal protein L7a [Bos taurus],0,+,1163.0
...
```

it is available at [_static/shmlast/mouse.1.rna.fna.gz.x.cow.faa.crbl.csv.gz](https://github.com/ngs-docs/angus/raw/17a0ba3b1d915de90a5b8bd1fbc1027eba47baf8/_static/shmlast/mouse.1.rna.fna.gz.x.cow.faa.crbl.csv.gz).

Also, plot distribution of e-value scores in `*.faa.crbl.csv`

[Viz blast scores with RStudio](visualizing-blast-scores-with-RStudio.html)
