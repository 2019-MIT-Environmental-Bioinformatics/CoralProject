#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=salmon
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nhuntley@whoi.edu
#SBATCH --ntasks=8
#SBATCH --mem=5000
#SBATCH --time=2:00:00
#SBATCH --output=salmon_%j.log
 

## Usage: sbatch salmon.sh 
	
salmon quant -i ./trinity-longest-index -l A -1 unfixrm_M_faveolata__Mfav_HH_euk_3_2502__L1_ATCACG_L001_1.trim.cor.fq unfixrm_M_faveolata__Mfav_HH_euk_33_2503__L1_CGATGT_L001_1.trim.cor.fq unfixrm_M_faveolata__Mfav_DD_euk_1_2505__L1_TGACCA_L001_1.trim.cor.fq unfixrm_M_faveolata__Mfav_DD_euk_31_2506__L1_ACAGTG_L001_1.trim.cor.fq -2 unfixrm_M_faveolata__Mfav_HH_euk_3_2502__L1_ATCACG_L001_2.trim.cor.fq unfixrm_M_faveolata__Mfav_HH_euk_33_2503__L1_CGATGT_L001_2.trim.cor.fq unfixrm_M_faveolata__Mfav_DD_euk_1_2505__L1_TGACCA_L001_2.trim.cor.fq unfixrm_M_faveolata__Mfav_DD_euk_31_2506__L1_ACAGTG_L001_2.trim.cor.fq -p 8 --validateMappings --output salmon/

##	-i: index directory
##	-l: libtype  type of sequencing library from which the reads come.  
## 	A: this tells the program to automatically decide the library type
## 	-1 , -2: -1 = forward reads, -2 = reverse reads
## 	-p:this is the number of threads 
## 	--validate mappings: salmon will use a more sensitive and accurate mapping algorithm. You can further add flags to specify control different parameters related to this. 
## 	--output: specify the name of the outfile, after this flag you specify the name. 
	
