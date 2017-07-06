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
        parallel libx11-dev python3-virtualenv last-align

Create a python 3 environment for dammit:

    python3.5 -m venv ~/py3
    . ~/py3/bin/activate
    pip install -U pip

Install [shmlast](https://github.com/camillescott/shmlast) (we used this earlier this week!).

    pip install -r <(curl https://raw.githubusercontent.com/camillescott/shmlast/master/requirements.txt)
    pip install --upgrade pip
    pip install shmlast

Install the proper version of GNU parallel:

    cd
    (wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - http://pi.dk/3) | bash
    sudo cp $HOME/bin/parallel /usr/bin/parallel

and then BUSCO...

    cd
    git clone https://gitlab.com/ezlab/busco.git

    echo 'export PATH=$HOME/busco:$PATH' >> $HOME/.bashrc

Finally, install dammit from the refactor/1.0 branch:

    pip install https://github.com/camillescott/dammit/archive/refactor/1.0.zip

Now, we'll install the databases. In the interest of time, we're going to do a "quick"
run -- this will omit OrthoDB, uniref, Pfam, and Rafm (ie, all the homology searches other
than user-supplied databases). If you run without `--quick`, it will take a lot longer (about
a half hour), but you'll have access to the full annotation pipeline.

    dammit databases --install --busco-group eukaryota

We used the "eukaryota" BUSCO group. We can use any of the BUSCO databases, so long as we install
them with the `dammit databases` subcommand. You can see the whole list by running
`dammit databases -h`. You should try to match your species as closely as possible for the best
results.


Make a directory for annotation and put files there

    cd
    mkdir -p nema_annotation
    cd nema_annotation
    ln -s ../assembly/trinity_out_dir/Trinity.fasta .

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
