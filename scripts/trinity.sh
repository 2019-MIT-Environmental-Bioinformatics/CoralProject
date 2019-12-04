#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=Trinity
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbecker@whoi.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --output=Trinitydenovo_%j.log
#export OMP_NUM_THREADS=1

singularity pull docker://trinityrnaseq/trinityrnaseq
# $1 = comma-separated list of R1 files
# $2 = comma-separated list of R2 files
# $3 = name of output directory Trinity will create to store results. This must include Trinity in the name, otherwise the job will terminate

Trinity --seqType fq --SS_lib_type RF --max_memory 124G --min_kmer_cov 1 --CPU 5 --left $1 --right $2 --output $3
