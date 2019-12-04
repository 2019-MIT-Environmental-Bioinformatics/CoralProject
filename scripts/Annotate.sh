#! /bin/bash
#SBATCH --partition=compute
#SBATCH --job-name=annotate
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nhuntley@whoi.edu
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=1
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=6
#SBATCH --mem=150G
#SBATCH --time=24:00:00
#SBATCH --output=annotate3_%j.log
#export OMP_NUM_THREADS=1

#We are annotating our de novo transcriptome against the SwissProt and TremBl protein databases.

blastx -query denovoAssembly_Trinity.fasta -db SwissProtDB -out denovo_annotateSwissProt.tab -outfmt 6 -evalue 0.00001 -max_target_seqs 1 -num_threads 6

# -query : this is our input denovo transcriptome. It must be in the same directory where this script is run.
# -db : specifies the database the program will use. This database must be in the same directory where the script is run, and do not add any file suffixes.
# -out : output file name 
# -outfmt : 6 specifies a tab file 
# -evalue : specify the max e-value 
# -max_target_seqs : this is the number of output annotations. We specified 1 based on talking with Carolyn, however normally the minimum that should be used is 5. 
# -num_threads : number of threads to use during analysis 
