SHELL := /bin/bash
SALMON := ~/SoftwareStaging/salmon/build/src/salmon
SEQ_DIR := /home/idworkin/DrosophilaSequences
STAR := STAR-2.5.2a/bin/MacOSX_x86_64/STAR

get_ref: ref/dmel-all-chromosome-r6.12.fasta ref/dmel-all-r6.12.gtf
	echo "downloading and upacking reference sequences"

ref/dmel-all-chromosome-r6.12.fasta: 
	wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/fasta/dmel-all-chromosome-r6.12.fasta.gz
	mkdir -p ref
	mv dmel-all-chromosome-r6.12.fasta.gz ref
	cd ref && gunzip dmel-all-chromosome-r6.12.fasta.gz

dmel-all-transcript-r6.12.fasta.gz:
	wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/fasta/dmel-all-transcript-r6.12.fasta.gz

ref/dmel-all-r6.12.gtf:
	wget ftp://ftp.flybase.net/releases/FB2016_04/dmel_r6.12/gtf/dmel-all-r6.12.gtf.gz
	mkdir -p ref
	mv dmel-all-r6.12.gtf.gz ref
	cd ref && gunzip dmel-all-r6.12.gtf.gz

# for transcript -> gene mappings
get_annotations: dmel-all-r6.12.gtf.gz
	gunzip -c dmel-all-r6.12.gtf.gz | grep "transcript_id" | cut -f 9 | sed "s/gene_id \"\([^\"]*\)\".*transcript_id \"\([^\"]*\)\".*/\2,\1/" | sort | uniq > txp_to_gene.csv

make_star_index: star_index/SAindex

star_index/SAindex: ref/dmel-all-r6.12.gtf ref/dmel-all-chromosome-r6.12.fasta
	mkdir -p star_index
	$(STAR) --runThreadN 8 --runMode genomeGenerate\
		--genomeDir star_index \
		--genomeFastaFiles ref/dmel-all-chromosome-r6.12.fasta \
		--sjdbGTFfile ref/dmel-all-r6.12.gtf \
		--sjdbOverhang 99

quant_all: salmon_index
	ls $(SEQ_DIR)/*R1*.fastq.gz | sed 's/\(.*\)_R1.*/\1/' | parallel -j6 "$(SALMON) quant -i salmon_index -p 6 -l A -1 {}_R1_001.fastq.gz  -2 {}_R2_001.fastq.gz -o quants/{/}"

#grep "transcript_id" dmel-all-r6.12.gtf | cut -f 9 | sed "s/gene_id \"\([^\"]*\)\".*transcript_id \"\([^\"]*\)\".*/\2\t\1/" | sort | uniq > txp_to_gene.tsv
#"~/SoftwareStaging/salmon/build/src/salmo quant -i salmon_index -t 6 -l IU -1
