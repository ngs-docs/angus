# Annotating de novo transcriptomes with dammit

dammit!

[dammit](http://www.camillescott.org/dammit/index.html) is an annotation
pipeline written by [Camille
Scott](http://www.camillescott.org/). dammit runs a relatively standard annotation
protocol for transcriptomes: it begins by building gene models with [Transdecoder](http://transdecoder.github.io/),
and then
uses the following protein databases as evidence for annotation:
[Pfam-A](http://pfam.xfam.org/), [Rfam](http://rfam.xfam.org/),
[OrthoDB](http://www.orthodb.org/),
[uniref90](http://www.uniprot.org/help/uniref) (uniref is optional with
`--full`).

If a protein dataset is available, this can also be supplied to the
`dammit` pipeline with `--user-databases` as optional evidence for
annotation.

In addition, [BUSCO](http://busco.ezlab.org/) v3 is run, which will compare the gene content in your transcriptome
with a lineage-specific data set. The output is a proportion of your
transcriptome that matches with the data set, which can be used as an
estimate of the completeness of your transcriptome based on evolutionary
expectation ([Simho et al.
2015](http://bioinformatics.oxfordjournals.org/content/31/19/3210.full)).
There are several lineage-specific datasets available from the authors
of BUSCO. We will use the `metazoa` dataset for this transcriptome.

## Installation

Annotation necessarily requires a lot of software! dammit attempts to simplify this and
make it as reliable as possible, but we still have some dependencies..

    sudo apt-get -y install python3-dev hmmer unzip \
        infernal ncbi-blast+ liburi-escape-xs-perl emboss liburi-perl \
        build-essential libsm6 libxrender1 libfontconfig1 \
        parallel libx11-dev python-virtualenv

Create a python 3 environment for dammit:

    ipython3.5 -m venv ~/py3
    . ~/py3/bin/activate
    pip install -U pip

Install [shmlast](https://github.com/camillescott/shmlast) (we used this earlier this week!).

    pip install -y --file <(curl https://raw.githubusercontent.com/camillescott/shmlast/master/requirements.txt)
    pip install --upgrade pip
    pip install shmlast

Now install LAST:

    cd
    curl -LO http://last.cbrc.jp/last-658.zip
    unzip last-658.zip
    pushd last-658 && make && make install prefix=~ && popd

Install the proper version of GNU parallel:

    cd
    (wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - http://pi.dk/3) | bash
    sudo cp /home/ubuntu/bin/parallel /usr/bin/parallel

And then Transdecoder:

    cd
    curl -LO https://github.com/TransDecoder/TransDecoder/archive/2.0.1.tar.gz
    tar -xvzf 2.0.1.tar.gz
    cd TransDecoder-2.0.1; make

and BUSCO...

    cd
    git clone https://gitlab.com/ezlab/busco.git

Put everything in the path:

    echo export PATH=$HOME/last-658/src:$PATH >> /home/ubuntu/miniconda3/bin/activate
    echo export PATH=$HOME/last-658/scripts:$PATH >> /home/ubuntu/miniconda3/bin/activate
    echo export PATH=$HOME/busco:$PATH >> /home/ubuntu/miniconda3/bin/activate
    echo export PATH=$HOME/TransDecoder-2.0.1:$PATH >> /home/ubuntu/miniconda3/bin/activate

Install the proper version of matplotlib

    pip install https://pypi.python.org/packages/source/m/matplotlib/matplotlib-1.5.1.tar.gz

Finally, install dammit from the refactor/1.0 branch

    pip install https://github.com/camillescott/dammit/archive/refactor/1.0.zip

Install databases (this step alone takes \~15-20 min) \# Is there a
faster install? \# Don't need everything?

    dammit databases --install

By default, the metazoan busco group will be installed. For the
eukaryota database, use this:

    dammit databases --install --busco-group eukaryota

Make a directory for annotation and put files there

    cd ${PROJECT}
    mkdir -p annotation
    cd annotation
    ln -s ${PROJECT}/assembly/trinity_out_dir/Trinity.fasta .

Download a custom *Nematostella vectensis* protein database available
from JGI:

    curl -LO ftp://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000001593_45351.fasta.gz
    #curl -LO ftp://ftp.jgi-psf.org/pub/JGI_data/Nematostella_vectensis/v1.0/annotation/proteins.Nemve1FilteredModels1.fasta.gz
    gunzip proteins.Nemve1FilteredModels1.fasta.gz

Putnam NH, Srivastava M, Hellsten U, Dirks B, Chapman J, Salamov A,
Terry A, Shapiro H, Lindquist E, Kapitonov VV, Jurka J, Genikhovich G,
Grigoriev IV, Lucas SM, Steele RE, Finnerty JR, Technau U, Martindale
MQ, Rokhsar DS. (2007) Sea anemone genome reveals ancestral eumetazoan
gene repertoire and genomic organization. Science. 317, 86-94.

Run the `dammit` pipeline

    # after trinity
    deactivate

    source activate dammit

Run the command:

    dammit annotate Trinity.fasta --busco-group metazoa --user-databases proteins.Nemve1FilteredModels1.fasta --n_threads 2 | tee dammit_Trinity_outfile.log

If dammit runs successfully, there will be a directory
`Trinity.fasta.dammit` with \~dozen files inside, including
`Trinity.fasta.dammit.gff3`, `Trinity.fasta.dammit.fasta` and a data
frame matching new annotated contig id with the previous
assembler-generated contig id: `Trinity.fasta.dammit.namemap.csv`. If
the above `dammit` command is run again, there will be a message:
`**Pipeline is already completed!**`
