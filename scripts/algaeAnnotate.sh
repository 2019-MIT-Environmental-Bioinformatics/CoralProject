#! /bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=Algaeannotate
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cbecker@whoi.edu
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=6
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --output=algaeannotate_%j.log
#export OMP_NUM_THREADS=4

#We are using BLASTN to map our de novo assembly against a symbiodinium database that is exactly the same as the paper one.

blastn -db algaeDB_all -query denovoAssembly_Trinity.fasta -out algaeAnnotate.tab -outfmt 6 -evalue 0.00001 -max_target_seqs 1 -num_threads 6

# - db : specifies our custom database 
# -query: this is the file being queried agains the custom database 
# -out: output file name 
# -outfmt: desired output format. In this case we chose 6 which produces a tab file 
# -evale: this is out max evalue 
# -max_target_seqs: this is the max number of sequences that will be kept. 
# -num_threads: number of threads this analysis should use. 
