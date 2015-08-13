for i in 1 2 3 4
do
curl -u single-cell:SanDiegoCA    http://bix.ucsd.edu/projects/singlecell/nbt_data/ecoli_mda_lane${i}.fastq.bz2 |bunzip2 - |head -400000 > ecoli_mda_lane${i}.fastq
bwa index sequence.fasta
bwa aln sequence.fasta ecoli_mda_lane${i}.fastq >ecoli_lane${i}.sai
bwa samse sequence.fasta ecoli_lane${i}.sai ecoli_mda_lane${i}.fastq > ecoli_lane${i}.sam
python ./sam-scan-errhist.py -o sam-scan-errhist_ecoli_lane${i}.out sequence.fasta ecoli_lane${i}.sam
done
