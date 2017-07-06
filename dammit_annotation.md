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
make it as reliable as possible, but we still have some dependencies.

    sudo apt-get -y install python3-dev hmmer unzip \
        infernal ncbi-blast+ liburi-escape-xs-perl emboss liburi-perl \
        build-essential libsm6 libxrender1 libfontconfig1 \
        parallel libx11-dev python3-venv last-align transdecoder

Create a python 3 environment for dammit:

    python3.5 -m venv ~/py3
    . ~/py3/bin/activate
    pip install -U pip

Install [shmlast](https://github.com/camillescott/shmlast) (we used this earlier this week!):

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
    pushd busco && python setup.py install && popd

    export PATH=$HOME/busco/scripts:$PATH
    echo 'export PATH=$HOME/busco/scripts:$PATH' >> $HOME/.bashrc

Finally, install dammit from the refactor/1.0 branch:

    pip install https://github.com/camillescott/dammit/archive/refactor/1.0.zip

## Database Preparation

dammit has two major subcommands: `dammit databases` and `dammit annotate`. `databases`
checks that the databases are installed and prepared, and if run with the `--install` flag,
will perform that installation and preparation. If you just run `dammit databases` on its
own, you should get a notification that some database tasks are not up-to-date -- we need
to install them!

In the interest of time, we're going to do a "quick"
run -- this will omit OrthoDB, uniref, Pfam, and Rfam (ie, all the homology searches other
than user-supplied databases). If you run without `--quick`, it will take a lot longer (about
a half hour), but you'll have access to the full annotation pipeline.

    dammit databases --install --quick --busco-group metazoa

We used the "metazoa" BUSCO group. We can use any of the BUSCO databases, so long as we install
them with the `dammit databases` subcommand. You can see the whole list by running
`dammit databases -h`. You should try to match your species as closely as possible for the best
results. If we want to install another, for example:

    dammit databases --install --quick --busco-group fungi

## Annotation

Keep things organized! Let's make a project directory:

    cd
    mkdir -p nema_annotation
    cd nema_annotation

You all ran Trinity earlier to generate an assembly, but just in case, we're going to download
a version of that assembly to annotate.

    curl -OL https://raw.githubusercontent.com/ngs-docs/angus/2017/_static/Trinity.fasta
    mv Trinity.fasta trinity.nema.fasta


Now we'll download a custom *Nematostella vectensis* protein database available
from JGI. Here, somebody has already created a proper database for us [1] (it has a reference proteome
available through uniprot). If your critter
is a non-model organism, you will
likely need to create your own with proteins from closely-related species. This will rely on your
knowledge of your system!

    curl -LO ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000001593_45351.fasta.gz
    gunzip -c UP000001593_45351.fasta.gz > nema.reference.prot.faa

Run the command:

    dammit annotate trinity.nema.fasta --quick --busco-group metazoa --user-databases nema.reference.prot.faa --n_threads 4

Note that `--quick` is used for both the `databases` subcommand and the `annotate` subcommand. If we
were to omit it at this point, it would warn us that the databases are not prepared and that we need
to install them -- because we skipped over the full database install earlier.

While dammit runs, it will print out which tasks its running to the terminal. dammit is
written with a library called [pydoit](www.pydoit.org), which is a python workflow library similar
to GNU Make. This not only helps organize the underlying workflow, but also means that if we
interrupt it, it will properly resume! 

After a successful run, you'll have a new directory called `trinity.nema.fasta.dammit`. If you
look inside, you'll see a lot of files:

    ls trinity.nema.fasta.dammit/
    annotate.doit.db                              trinity.nema.fasta.dammit.namemap.csv  trinity.nema.fasta.transdecoder.pep
    dammit.log                                    trinity.nema.fasta.dammit.stats.json   trinity.nema.fasta.x.nema.reference.prot.faa.crbl.csv
    run_trinity.nema.fasta.metazoa.busco.results  trinity.nema.fasta.transdecoder.bed    trinity.nema.fasta.x.nema.reference.prot.faa.crbl.gff3
    tmp                                           trinity.nema.fasta.transdecoder.cds    trinity.nema.fasta.x.nema.reference.prot.faa.crbl.model.csv
    trinity.nema.fasta                            trinity.nema.fasta.transdecoder_dir    trinity.nema.fasta.x.nema.reference.prot.faa.crbl.model.plot.pdf
    trinity.nema.fasta.dammit.fasta               trinity.nema.fasta.transdecoder.gff3
    trinity.nema.fasta.dammit.gff3                trinity.nema.fasta.transdecoder.mRNA

The most important files for you are `trinity.nema.fasta.dammit.fasta`,
`trinity.nema.fasta.dammit.gff3`, and `trinity.nema.fasta.dammit.stats.json`.

If the above `dammit` command is run again, there will be a message:
`**Pipeline is already completed!**`

## References

1. Putnam NH, Srivastava M, Hellsten U, Dirks B, Chapman J, Salamov A,
Terry A, Shapiro H, Lindquist E, Kapitonov VV, Jurka J, Genikhovich G,
Grigoriev IV, Lucas SM, Steele RE, Finnerty JR, Technau U, Martindale
MQ, Rokhsar DS. (2007) Sea anemone genome reveals ancestral eumetazoan
gene repertoire and genomic organization. Science. 317, 86-94.
