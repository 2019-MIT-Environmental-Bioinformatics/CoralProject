#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=Samtools
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbecker@whoi.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G
#SBATCH --time=02:00:00
#SBATCH --output=Samtools-bam-and-sort_%j.log
#export OMP_NUM_THREADS=1


for x in *bt2.sam
do
name=$(basename ${x} bt2.sam)
samtools view -bS $x > ${name}bt2.bam
samtools sort ${name}bt2.bam -o ${name}bt2.sorted.bam
done

# This script takes a bam file that was outputted from the bowtie2 alignment and view converts it to a BAM file.
# Then the sort function converts the BAM file to a sorted BAM file. 

# After this step --> 
# If I use express this is the usage: express [options]* <target_seqs.fasta> <aligned_reads.(sam/bam)>
# maybe use a wildcard, or give it a list
# I think the target sequences is our fasta file of our de novo transcriptome?

# Some options for generating counts:
# => eXpress
# => hw3 code
# => SAM_nameSorted_to_uniq_count_stats.pl perl script from Trinity?
