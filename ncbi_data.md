##SRA Toolkit

Log into Amazon with a small session and do [the normal linuxbrew install](http://angus.readthedocs.io/en/2016/linuxbrew_install.html).

To download SRA data on the CLI, we need the SRA Toolkit. Lets check installed it

````
brew install sratoolkit
```


There is [lots of documentation on how to use the toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/?view=toolkit_doc). Lets fetch some data

```
prefetch SRR1778760
```

This automatically downloads the data to `~/ncbi/public/sra/`. You can download multipe files

```
prefetch SRR1778760 SRR1778761
```

It knows what is already stored locally and doesn't get that one, just the new one. If you need the files to be stored somewhere other than `~/ncbi/public/sra/`, you can use the very weird interactive configuration menu:

```
vdb-config -i
```

We still have files in .sra format, which is an NCBI SRA thing (I.E. no other tools will deal with it). So we need to convert to a more standard format, usually fastq or fasta.

Convert to fastq in current directory 

```
fastq-dump SRR1778760
ls
```

Check out the file this yielded. This is actually still a bit odd - it arrived in one file. Often we prefer paired reads to be broken out into two files: a file with all the forward reads and a file with all the reverse reads.  The sra toolkit will let you do this:

```
fastq-dump --split-files SRR1778760
ls
```


##EUtils


The NCBI has a toolkit which they call Entrez Programming Utilities or **eutils** for short. You can read all about it in the [documentation](http://www.ncbi.nlm.nih.gov/books/NBK25501/). There are a lot of things you can do to interface with all of the different NCBI databases, including:
* ask for the total number of records in a database
* search with a text query
* search with a text query and respond with the number of matches
* provide a list of IDs and ask for information back in a certain format
* and more

Eutils will work to query these databases:
* BioProject
* BioSample
* Biosystems
* Books
* Conserved Domains
* dbGaP
* dbVar
* Epigenomics
* EST
* Gene
* Genome
* GEO Datasets
* GEO Profiles
* GSS
* HomoloGene
* MeSH
* NCBI C++ Toolkit
* NCBI Web Site
* NLM Catalog
* Nucleotide
* OMIA
* PopSet
* Probe
* Protein
* Protein Clusters
* PubChem BioAssay
* PubChem Compound
* PubChem Substance
* PubMed
* PubMed Central
* SNP
* SRA
* Structure
* Taxonomy
* UniGene
* UniSTS


> ## API = application programming interface
>
> In this case, we'll talk about an API in terms of the internet. A web API is a programmatic interface to a request-response message system. A user sends a request via the internet (via http), and the server responds with structured data. 
>
> Another definition: "“an interface through which you access someone else’s code or through which someone else’s code accesses yours – in effect the public methods and properties.”


Lets focus on the genome db for now. 

EUtil requests are basically about building a web URL. The base URL is always:

	http://eutils.ncbi.nlm.nih.gov/entrez/eutils/


Assuming you know the unique identifier of your genome of interest, you can start a download with this URL:

	http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000962&rettype=fasta&retmode=text

Let’s breakdown the command here:

Part | Explanation
-----|-----------
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi? | This is command telling your computer program (or your browser) to talk to the NCBI API tool efetch.
db=nuccore | This command tells the NCBI API that you’d like it to look in this particular database for some data. Other databases that the NCBI has available can be found here.
id=CP000962 | This command tells the NCBI API efetch the ID of the genome you want to find.
rettype=fasta&retmode=text | These two commands tells the NCBI how the data is returned. We asked for the FASTA sequence as a text file.

Lets try asking for the Genbank file instead:

	http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000962&rettype=gb&retmode=text

What changed?

[Here’s some elusive documentation on where to find these “return” objects](http://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly), ie which formats are available for each database (such as fasta format and genbank format for genome database entries).

How could we get this data into our newton account? With curl! (This is one of those instances where wget will work but it will save your data in a weird file name, so curl is better).

```
	curl "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CP000962&rettype=gb&retmode=text" > CP000962.gb
```
	
##Many nucleotide records
What if you want a whole bunch of records? Lets try a large set - go get a list of GI accessions from NCBI and scp them to your Amazon instance.

We know how to:
* Use the web portal and look up each FASTA
* Use the FTP site, find each genome, and download manually or with wget
* Build an eutils URL and use wget

These aren't great options if you have 300 records. However, this turns out to be very useful when you know how to write a bash script, or write some python code or R code or perl code or something that does loops. Here's is a bash script:

```
cat GIs.txt | while read line
do
  echo $line
  begin="http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="
  end="&rettype=fasta&retmode=text"
  url=$begin$line$end
  echo "${url}"
  curl $url >> GIs.fasta
done 
```

How would you alter this script to get genbank format instead of fasta format?

