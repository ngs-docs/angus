## Read group information

## a) Sample and library tags
*  SM = biological sample name. 
*  LB = name of DNA preparation library tube = {SM}.{library-specific identifier}  
(Important To identify PCR duplicates in MarkDuplicates step. Ignore in PCR free libraries)

Can be autmatically detected from current sample naming scheme:  
<Sample.ID>_<Index.Sequence>_<Lane.ID>_<R1or2>_<Set.number>.fastq

* SM = <Sample.ID>
* LB = <Sample.ID>_<Index.Sequence>   

## b) ID and PU (to enable merging replictes)
* ID = Read group identifier = {FLOWCELL_BARCODE}.{LANE}
* PU = Platform Unit = {FLOWCELL_BARCODE}.{LANE}.{library-specific identifier}. This is the most specific definition for a group of reads.

Also can be identified from the name of a sequence read in the Fastq file:  
@(instrument id):(run number):(flowcell ID):(lane):(tile):(x_pos):(y_pos)  (read):(is filtered):(control number):(index sequence)  
FLOWCELL_BARCODE =  @(instrument id):(run number):(flowcell ID)

## c) Others e.g. 
* PL = Platform (Valid values: ILLUMINA, SOLID, LS454, HELICOS and PACBIO)

> Notes:
> - One sample (SM) can have multiple libraries (e.g SE, PE50, and PE100) (LB), can run on multiple lanes and/or multiple flow cells (RGID), and can run on multiple platforms (PL).
> - One library can run on multiple lanes or multiple flow cells (PU).
> - If we have multiple samples for the same individual e.g. before and after treatment, each sample should have a different SM but unless you expect change of sequence, we can consider them multiple libraries of the same sample)
> - Multiple samples can share the same Read group ID (When manuals sya “must be unique”. They mean unique in a BAM file. So it is ok that multiple samples can share the same Read group ID)
> - if you have one library for each sample running on one lane of a sequencing machine then you can make SM=LB=RGID=PU

