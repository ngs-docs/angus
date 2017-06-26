# Running large and long command line jobs - using shmlast!

* should we introduce screen/tmux?

shmlast!

Install base packages

```
sudo apt-get -y update && \
sudo apt-get install -y python3.5-dev python3.5-venv make \
    libc6-dev g++ zlib1g-dev last-align parallel
```

Then create a Python virtualenv:
```
python3.5 -m venv ~/py3
. ~/py3/bin/activate
pip install -U pip
```

And now install shmlast 1.2:
```
pip install shmlast==1.2
```

More here.


```
curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.1.rna.fna.gz

for i in 1 2 3 4 5 6 7 8
do
   curl -O ftp://ftp.ncbi.nih.gov/refseq/B_taurus/mRNA_Prot/cow.$i.protein.faa.gz
done

gunzip -c cow.*.faa.gz > cow.faa

shmlast crbl -q mouse.1.rna.fna.gz -d cow.faa --n_threads=6
```

Now, download this file & look at in a spreadsheet.

Also, plot distribution of e-value scores in `*.faa.crbl.csv`

[Viz blast scores with RStudio](visualizing-blast-scores-with-RStudio.html)
