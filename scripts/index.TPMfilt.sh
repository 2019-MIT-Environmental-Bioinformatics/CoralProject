#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=IndexLongestTPMfilt
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbecker@whoi.edu
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --output=Index_%j.log
#export OMP_NUM_THREADS=1

#You must have bowtie2 in your coral conda environment before running
#Description of script: bowtie2-build [denovoAssembly.fasta] [nameofIndexProducedFromScript]
#Adjust if name of assembly with only longest isoforms is there.

bowtie2-build Trinity.longest.TPMfilt.fasta denovoTrinity.longest.TPMfilt.Index
