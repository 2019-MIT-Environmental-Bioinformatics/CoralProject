#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=trimming
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nhuntley@whoi.edu
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=5000
#SBATCH --time=10:00:00
#SBATCH --output=trimmed_%j.log
#export OMP_NUM_THREADS=1

for infile in *_R1_001.fastq.gz
do
   base=$(basename ${infile} _R1_001.fastq.gz)
   trimmomatic PE ${infile} ${base}_R2_001.fastq.gz \
                ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
                ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20
done

## no minlen specified 
## mp adaptor removal specified 

