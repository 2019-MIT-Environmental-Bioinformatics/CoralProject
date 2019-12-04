#! /bin/bash

#SBATCH --partition=compute
#SBATCH --job-name=BowtieLongest
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nhuntley@whoi.edu
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G
#SBATCH --time=12:00:00
#SBATCH --output=Bowtie2-LongestIsoform-align_%j.log
#export OMP_NUM_THREADS=1

#Run script with: sbatch bowtie.longest.sh.save NameOfDeNovoTranscriptomeIndex

# $1 is the name of the index file of the de novo transcriptome. Do not include extensions. Ex) denovoTrinityIndex or denovoTrinity.longest.Index (note the lack of extensions). 

for x in *1.trim.cor.fq
do
name=$(basename ${x} 1.trim.cor.fq)
echo ${name}
bowtie2 -p 16 --local -x $1 -q \
-1 ${name}1.trim.cor.fq \
-2 ${name}2.trim.cor.fq \
-S ${name}bt2.sam
done
