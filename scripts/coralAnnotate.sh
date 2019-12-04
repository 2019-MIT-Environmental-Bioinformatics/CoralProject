#! /bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=Coralannotate
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nhuntley@whoi.edu
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=6
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --output=Coralannotate_%j.log
#export OMP_NUM_THREADS=4

#We are comparing our sequences to out custom coral database

blastn -db cnidarian_DB -query denovoAssembly_Trinity.fasta -out coralAnnotate.tab -outfmt 6 -evalue 0.00001 -max_target_seqs 1 -num_threads 6

#  -db : this is the custom cnidarian database that we created. Our de novo assembly will be compared to this 
# -query : this is out denovo assembly fasta file 
# -out : specifies the name of our output file 
# -outfmt : specifies the output format - we chose 6 to create a tab file 
# -evalue : specifies the max evalue 
# -max_target_seqs : this is the number of sequences that are reported  
# -num_threads: number of threads to use for this analysis 
