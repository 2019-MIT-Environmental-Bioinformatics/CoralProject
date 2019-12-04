#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=index
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nhuntley@whoi.edu
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --mem=100G
#SBATCH --time=4:00:00
#SBATCH --output=salmonindex_%j.log
#export OMP_NUM_THREADS=2

## usage: sbatch salmon

salmon index -t Trinity.longest.fasta -i trinity-longest-index -k 25

##	-k : kmers 
##	-i : index 

