SHELL := /bin/bash
SALMON := ~/SoftwareStaging/salmon/build/src/salmon
SEQ_DIR := /home/idworkin/DrosophilaSequences

dmel-all-transcript-r6.12.fasta.gz:
	wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/fasta/dmel-all-transcript-r6.12.fasta.gz

dmel-all-r6.12.gtf.gz:
	wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/gtf/dmel-all-r6.12.gtf.gz

# for transcript -> gene mappings
get_annotations: dmel-all-r6.12.gtf.gz
	gunzip -c dmel-all-r6.12.gtf.gz | grep "transcript_id" | cut -f 9 | sed "s/gene_id \"\([^\"]*\)\".*transcript_id \"\([^\"]*\)\".*/\2,\1/" | sort | uniq > txp_to_gene.csv

make_index: dmel-all-transcript-r6.12.fasta.gz
	$(SALMON) index -t <(gunzip -c dmel-all-transcript-r6.12.fasta.gz) -i salmon_index

quant_all: salmon_index
	ls $(SEQ_DIR)/*R1*.fastq.gz | sed 's/\(.*\)_R1.*/\1/' | parallel -j6 "$(SALMON) quant -i salmon_index -p 6 -l A -1 {}_R1_001.fastq.gz  -2 {}_R2_001.fastq.gz -o quants/{/}"

#grep "transcript_id" dmel-all-r6.12.gtf | cut -f 9 | sed "s/gene_id \"\([^\"]*\)\".*transcript_id \"\([^\"]*\)\".*/\2\t\1/" | sort | uniq > txp_to_gene.tsv
#"~/SoftwareStaging/salmon/build/src/salmo quant -i salmon_index -t 6 -l IU -1
