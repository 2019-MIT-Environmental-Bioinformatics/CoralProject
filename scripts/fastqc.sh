#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=fastqc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nhuntley@whoi.edu
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=5000
#SBATCH --time=2:00:00
#SBATCH --output=fastqc_%j.log
#export OMP_NUM_THREADS=1

## usage: sbatch fastqc.sh 

fastqc M*fastq.gz

## output looks like: M_faveolata__Mfav_DD_euk_1_2505__L1_TGACCA_L001_R1_001_fastqc.html    
## M_faveolata__Mfav_HH_euk_33_2503__L1_CGATGT_L001_R1_001_fastqc.html
## M_faveolata__Mfav_DD_euk_1_2505__L1_TGACCA_L001_R1_001_fastqc.zip     
## M_faveolata__Mfav_HH_euk_33_2503__L1_CGATGT_L001_R1_001_fastqc.zip
