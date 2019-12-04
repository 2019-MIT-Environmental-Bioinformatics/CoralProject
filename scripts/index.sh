#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=IndexTranscriptome
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nhuntley@whoi.edu
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --output=Index_%j.log
#export OMP_NUM_THREADS=1

bowtie2-build Trinity.fasta denovoTrinityIndex
