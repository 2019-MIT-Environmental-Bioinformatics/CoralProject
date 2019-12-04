setwd("~/Google Drive (cbecker@whoi.edu)/MIT_WHOI Coursework/Fall2019/Bioinformatics/FinalProject_Coral/")

#Bioinformatics Final Project

#Naomi Huntley and Cynthia Becker

#Goal of this script: Use the counts data from the alignment (Bowtie2) of sequences to the de novo assembled coral holobiont transcriptome (Trinity) that has been filtered to represent only the longest isoforms that are most highly expressed (TPM > 1) to generate differential expression stats in R.

#Attempting to reproduce the following statement from the paper, Daniels et al: 
# "imported into edgeR, and filtered according to gene loci >0-counts-per-million present in both samples. Subsequently, TMM library normalization and dispersion estimation were conducted.An exact test was performed on negative binomial fitted data to determine differentially expressed genes between HH and DD samples at an FDR cutoff of <0.1 (Benjamini and Hochberg, 1995). Expression data for eukaryotic genes are reported as log2 fold changes"

#Install the required packages for analysis

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("edgeR")
library(edgeR)

setwd("~/Google Drive (cbecker@whoi.edu)/MIT_WHOI Coursework/Fall2019/Bioinformatics/FinalProject_Coral/")

sessionInfo()
#other attached packages:
#[1] edgeR_3.28.0    limma_3.42.0    phyloseq_1.30.0

counts = as.matrix(read.csv("Mfav_counts_all_filttxm.tab", sep="\t", row.names="target_id"))
#check the dimensions of the counts matrix
dim(counts) #[1] 66044      5

counts = counts[,2:5]

treatment = matrix(c("DD", "HH",  "DD", "HH"), ncol = 1)
colnames(treatment) = "condition"
rownames(treatment) = c("DD_euk_1", "HH_euk_33", "DD_euk_31", "HH_euk_3")

#Rename the column names of the counts file so it is simpler to read

colnames(counts) <- rownames(treatment)


#1. Make the DGElist unit that edgeR works with. Combines counts and a group variable

DElist <- DGEList(counts = counts, group = treatment[,1], samples = treatment)


#2. Filtering "according to gene loci >0-counts-per-million present in both samples"

countsPerMillion <- cpm(DElist) #turn counts into counts per million
summary(countsPerMillion) #take a look at the summary statistics for counts per million for each sample

countCheck <- countsPerMillion > 0 #interested in only >0 cpm. This will make a matrix of TRUE and FALSES
head(countCheck)

keep <- which(rowSums(countCheck) >= 2) #group size is 2, so you want to make sure to keep contigs with >0 cpm in "both samples", which I interpret as the 2 samples in each group. Use the which, because it will select "which" data are TRUE.
DE.flt <- DElist[keep, , keep.lib.sizes=FALSE] #Also need to recalculate the library sizes
summary(cpm(DE.flt)) #check out resulting summary stats and compare them to the previous summary. It looks like the max cpm has decreased from 50s-70s down to 40s. There is now a median in some of the treatments. 

#Let's see what percentage of transcripts are left:
Filtstats <- data.frame("orig" = DElist$samples$lib.size, "postfilt" = DE.flt$samples$lib.size, "percentkept" = (DE.flt$samples$lib.size/DElist$samples$lib.size*100))
Filtstats$percentkept #Looks like only about 84-88% of reads were kept...[1] 88.79087 87.37074 84.68358 88.85941

#3. Normalization "TMM library normalization and dispersion estimation were conducted"

DE.flt.norm <- calcNormFactors(DE.flt, method = "TMM") #Normalize library within and between samples using TMM (trimmed mean of M-values) library normalization
Disp <- estimateCommonDisp(DE.flt.norm) #estimate common dispersion... 
Disp <- estimateTagwiseDisp(Disp) #and then tagwise dispersion, both using the qCML (quantile-adjusted conditional maximum likelihood) method.


#4. Differential expression. "An exact test was performed on negative binomial fitted data to determine differentially expressed genes between HH and DD samples at an FDR cutoff of <0.1"

DEgenes <- exactTest(Disp, pair = c("HH", "DD")) #Use the Disp DGElist object to then calculate the differentially expressed genes by using the exact test to test between negative-binomially distributed counts with dispersion estimates. I specified that DD is to be compared to the control, HH, which represents healthy coral.
topTags(DEgenes) #take a look at the top DE genes.

#Write DE genes table, but use the topTags because it can sort the table by p-value, or by logFC, and can also offer a false-discovery-rate p-value adjustment
dim(DEgenes) #Looks like there are 38,067 DE genes
DEtableP <- topTags(DEgenes, n = 56990, adjust.method = "BH", sort.by = "PValue", p.value = 1) 

write.table(DEtableP, file = "DEgenes_coral.tsv", sep = "\t")


#We have a problem....none of our DE genes are significant below an FDR cutoff of 0.1...WHYYYYYYYYYY. I shall go through the pipeline again and make sure I didn't mess anythng up.

#There are now significant PValues but the FDR cutoff p values is 1 for everything. I am not sure why. I may try and subset just the contigs with a p-value less than 0.1 as a substitute, since I am running into a problem here. 

DEtablePval.1 <- DEtableP$table[DEtableP$table[,"PValue"]<0.1,] #collect only the p values less than 0.1.
dim(DEtablePval.1) #[1] 1058    4

write.table(DEtablePval.1, file = "DEgenes_pval0.1_coral.tsv", sep = "\t")
