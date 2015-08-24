##De novo genome assembly with velvet

###Links to programs and data used during this tutorial

**Velvet installation file**  
https://www.ebi.ac.uk/~zerbino/velvet/velvet_1.2.10.tgz

**MiSeq paired end reads**  
https://www.dropbox.com/s/kopguhd9z2ffbf6/MiSeq_Ecoli_MG1655_50x_R1.fastq  
https://www.dropbox.com/s/i99h7dnaq61hrrc/MiSeq_Ecoli_MG1655_50x_R2.fastq  
Source: [http://www.illumina.com/science/data_library.ilmn](http://www.illumina.com/science/data_library.ilmn), random subsampling using seqtk [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)

**Nextera mate pair reads**  
https://www.dropbox.com/s/vd9toordtf181ki/Nextera_MP_R1_50x.fastq  
https://www.dropbox.com/s/9aziigagrwh0c6n/Nextera_MP_R2_50x.fastq  
Source: Illumina basespace [https://basespace.illumina.com/‎](https://basespace.illumina.com/‎), look for "Nextera Mate Pair (E. Coli)" [https://basespace.illumina.com/project/294296/Nextera-Mate-Pair-E-Coli](https://basespace.illumina.com/project/294296/Nextera-Mate-Pair-E-Coli)

**Assemblathon_stats script**  
https://raw.githubusercontent.com/lexnederbragt/sequencetools/master/assemblathon_stats.pl  
