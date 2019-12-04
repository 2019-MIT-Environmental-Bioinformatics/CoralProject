# CoralProject
## Coral white plague disease holobiont metatranscriptome.
**Team:** Naomi Huntley and Cynthia Becker

**Proposed Paper:** Daniels CA, Baumgarten S, Yum LK, Michell CT, Bayer T, Arif C, Roder C, Weil E, Voolstra CR. 2015. Metatranscriptome analysis of the reef-building coral Orbicella faveolata indicates holobiont response to coral disease. Front Mar Sci. 2. doi:<10.3389/fmars.2015.00062>. [accessed 2019 Sep 15]. <http://journal.frontiersin.org/Article/10.3389/fmars.2015.00062/abstract>.    
        We propose reanalyzing the data from the paper titled, Metatranscriptome analysis of the reef-building coral Orbicella faveolata indicates holobiont response to coral disease, by Daniels et al. (2015). This paper compares differences in gene expression between two healthy and two white plague diseased coral hosts and their associated algae and bacterial symbionts. For our project we would focus on comparing four coral transcriptomes and four bacterial transcriptomes, each from two healthy and two white plague diseased Orbicella faveolata coral colonies. To cut down the analysis, we will omit the algal symbionts from our reanalysis. These data are available on NCBI (<https://www.ncbi.nlm.nih.gov/sra/SRX1096797[accn]>) and are approximately 30Gb total.     
        Our reanalysis will focus on recreating figure 1 and figure 3. Figure 1 is a comparison between differentially expressed genes associated with disease in the coral and bacterial “compartments” of the coral host. Figure 3 analyzes the bacteria that expressed the disease-associated genes. We found figure 3 challenging to interpret, and our reanalysis will also focus on finding a way to improve readability of the figure.
        After beginning the project, we ran into difficulty with MG-RAST, so we decided to move forward analyzing the algal and coral component of the holobiont, rather than the bacterial and coral components. As a result, we focused on recreating figure 1 and Data Sheet 4, which summarizes the differentially expressed genes that map to both the coral and algal compartments of the holobiont. The overarching goal was to reproduce an analysis that investigated how genes were differentially expressed based on the disease state of a coral host.    

## Folder Structure and outline of documents.
### /data/:Contains raw data. We did NOT push this to GitHub due to large file sizes
- /data/origfastq/ : original, untrimmed, unfiltered data files
- /data/TrimCorrectedFastq/ : contains all of the trimmed and corrected sequences. 
- /data/trimmed-corrected-filt/ : contains the trimmed, corrected, and filtered sequences
- /data/TrimmedFastq/ : contains the trimmed reads. 

### /envs/ : Contains a yaml file that can be used to generate the conda environments *coral* and *Trinity* used for data analysis. There are two separate environments because of package conflicts. 
- /envs/coralconda.yml
- /envs/trinityconda.yml     

### /jupyter-notebooks/ : Contains our jupyter-notebooks with pythons scripts.
- Combine_BLASTX-and-GO-Annotations.ipynb : python script to merbe the blastx and gene ontology annotations into one tab file          
- express_results_to_count_table_FilteredTXM.ipynb : python script to create a count table from the express results.  
- Final-Comparison.ipynb : Our major output comparing the results from our analysis to the paper results. 
- Split-algae-coral.ipynb : python script to assign contig to either coral or algae. 

### /scripts/ : Scripts used for analysis of project 			
- /scripts/algaeAnnotate.sh : compares sequences to custom algae database. 
    - Produces: algaeAnnotate.tab ( in the output folder) 
- /scripts/Annotate.sh : sbatch script for annotating de novo transcriptome against the SwissProt and TremBl protein databases. 
    - Produces:  denovo_annotateSwissProt.tab
- /scripts/bowtie.longest.TPMfilt.sh: Bowtie script used to map the longest isoforms that were filtered by TPM to our de novo transcriptome    
- /scripts/bowtie.sh : original bowtie script for unfiltered samples.                     
- /scripts/coralAnnotate.sh : compares sequences to our custom coral database. 
    - Produces: coralAnnotate.tab ( in the output folder) 
- /scripts/edgeR_coral_FinalProject.R : code to run edgeR to determine differentially expressed genes. 
- /scripts/eXpress.Longest.sh: This script takes each sorted BAM file from SAMtools and maps the reads back to the de novo transcriptome assembly from Trinity (example: Trinity.longest.TPMfilt.fasta)
- /scripts/eXpress.sh : This script takes each sorted BAM file from SAMtools and maps the reads back to the de novo transcriptome assembly from Trinity (Trinity.fasta). 
- /scripts/fastqc.sh : sbatch script to use fastqc to quality control reads. 
- /scripts/index.sh : Index for Bowtie 2 when using the original transcriptome
- /scripts/index.TPMfilt.sh : Index for Bowtie 2 when using the longest isoform only after TPM filter.     
- /scripts/Salmon-index.sh : Script to create the index needed for salmon.       
- /scripts/salmon.sh : Script to run salmon 
- /scripts/samtools.Longest.sh : This script takes a bam file that was outputted from the bowtie2 alignment and converts it to a BAM file. Then it sorts the BAM file.
- /scripts/Samtools.sh : This does the same as above, but was used on the unfiltered transcriptome. 
- /scripts/trim.sh : sbatch script to to trim the files with Trimmomatic. 
    - Produces : /output/untrimFastq/*un.trim.fastq.gz :files containing untrimmed and files containing trimmed reads in output/trimFastqc/*trim.fastq.gz and logs in /logs/trimmed_489211.log
-  /scripts/trinity.sh : code to generate a de novo transcriptome using trinity      

### /logs/ : Contains the output logs from data analysis runs 				
- /logs/Bowtie2-align_493063.log : this is the log file from the first time we used bowtie2 on the totally unfiltered sequences. 
- /logs/Bowtie2-LongestIsoform-align_495302.log : this is the log file from the second time we used bowtie2 on the longest isoform without TPM filtering. 
- /logs/Bowtie2-LongestIsoform-filt-align_495435.log : This is the log file from the second time we used bowtie2 to align the sequences that contained only the longest isoforms that were also TPM filtered. 
- /logs/eXpress_get_counts_LongestIso-495450.log :  this is the log from the express output.
- /logs/fastqc_489205.log: output log from Fast-QC analysis. 
- /logs/Index_denovoTrinity.longest.TPMfilt.log: Log from indexing the filtered TXM in preparation for Bowtie2. 
- /logs/Rcorrector*: This log document the use of Rcorrector to error-correct the trimmed fastq files. 
- /logs/salmonLogs/  : 
    - /salmon_495429.log : log from salmon quant step.
    - /salmonindex_495414.log : log from salmon index step. 
- /logs/transcriptlog/: Notes for commands we used and files we created as we went through the project. We treated it like a lab notebook for our computer pipeline.
- /logs/trimmed_489211.log: output from Trimmomatic analysis.
- /logs/Trinitydenovo_490484.log : Note this log file from the trinity de novo assembly is 70Mb. 
- /logs/unfixable/ : These logs are output log files generated from the `FilterUnfixableSeqs.txt` script and describe how many sequences were kept after filtering. 	   
	
### /output/ : contains the output of all analyses.
- /output/counts_longest_TPMfilt/: this is the count file for the longest isoforms that were TPM filtered. 
- /output/FinalDE_Table/
    - /algae-Copy1.tab : tab file containing contigs assigned to algae based on e-val and bitscore.           
    - /BLASTX_and_GO_merged.tab : This is a tab file containing the blastx and gene ontology annotations merged together into one file.  
    - /coral-Copy1.tab: tab file containing contigs assigned to coral based on e-val and bitscore.            
    - /DataSheet4_Danielsetal2015.txt : This datasheet is from the Daniels et al. paper and was our goal output for replication. 
    - /DataSheet4Replica-all-p-val.tab : This is our replica of datasheet 4 with no p-value cut off. 
    - /DataSheet4Replica.tab: This is our replica of datasheet 4 with a 0.1 p-value cut off. 
    - /DEgenes_coral_filttxm.tsv: This is our output file from edgeR listing all the DE genes recovered. This is from the filtered txm.
    - /DEgenes_pval0.1_coral_filttxm.tsv: This is the output file from edgeR listing all DE genes recovered with a p-value < 0.1. This is from the filtered txm.
- /output/Longest.TPMfilt.SAMBAMfiles/ : Bowtie2 and SAMtools output - SAM and BAM files of the TPM filtered sequences. 
- /output/QC : This folder contains all of the .html and .zip files for each of the samples following the quality control step using FastQC.
- /output/salmon/ : 
    - /quant.sf :  salmon pseduo-alignment output. 
    - /salmon_quant.isoform.counts.matrix :  salmon output converted to a matrix for TPM filtering. 
    - /salmon_quant.isoform.TPM.not_cross_norm : salmon output normalized for TPM filtering. 
- /output/swissprot_TrEMBL_DB/: contains the files to make the swissprot TrEMBLE database.   
- /output/transcriptome-database
    - /algaeAnnotate.tab : Trinity annotations against algae database
    - /algaefastas/ :folder that contains each of the algae genomes used to create a custom algae database.                 
    - /Algae.tab: tab file of trinity annotations that have been assigned to algae based on e-value and bitscore using Split-algae-coral.ipynb. 
    - /coralAnnotate.tab: Trinity annotations against algae database
    - /coralfastas/ :folder that contains each of the algae genomes used to create a custom coral database.        
    - /coral.tab : tab file of trinity annotations that have been assigned to coral based on e-value and bitscore using Split-algae-coral.ipynb. 
    - /Trinity.fasta : De novo transcriptome needed for splitting it into algae and coral parts
- /output/Transcriptomes/ : contains our trinity transcriptomes and the assembly statistics that correspond with each. 
- /output/trimFastqc : Folder that contains the trimmed fastq files for each sample.  
- /output/trinity.fasta-index: trinity output containing all isoforms.  
- /output/trinity-longest-index : trinity output against the longest isoforms only. 
- /output/TrinityIndex_TPMfilt_forBT2 : indexed transcriptome files for bowtie2
- /output/trinity_out_dir : trinity assembly output files
- /output/untrimFastq : orphaned reads that were not included in analysis 	

### Coral.gitignore
- This file contains paths to all files over 50Mb that are ignored, and not pushed to GitHub.

## Procedure and Instructions for Reproducing our Pipeline and Analysis

#### Raw sequence files were not uploaded properly to NCBI, so we had to reach out to the author to retrieve them. The files were downloaded to a local computer and copied onto the HPC using the `scp` command. All gzipped, raw files are located in the CoralProject folder (/env-bio/collaboration/CoralProject/) on the HPC at **/data/origfastq/**. In the following instructions, all files are in folders within the parent directory `/env-bio/collaboration/CoralProject/`.

### 1. Quality filtering and trimming

**FastQC v0.11.8.** Use the script, `/scripts/fastqc.sh` to perform standard quality control on the gzipped fastq files. To use the script, fastq files must be located within the same directory.  An example of our output fastqc files (.html and .zip) can be found in the `/output/QC/` folder. Feel free to look at the .html files to inform future trimming parameters.     
**Trimmomatic v0.39.** Use the script `/scripts/trim.sh` on the original .fastq.gz files to both merge and trim files. This script uses a package called Trimmomatic. Files without a mate will be placed into an *un.trim.fastq.gz file. Files for use downstream will have the suffix *trim.fastq.gz. The trimmomatic parameters used were general default parameters, which were outlined in the Daniels et al. paper at the beginning of the README.md. Trimmed fastq files are located at `/data/TrimmedFastq/`. 

### 2. Error correction

**Note:** The paper above, Daniels et al. (2015), uses the **ALLPATHS-LG** standalone error correction module, `ErrorCorrectReads.pl`, which is part of the package. Currently, this standalone module is no longer supported, so we used a different package for error correction, called **rCorrector v1.0.4**, which is available for a conda installation.     

**Error Correction.** Use the `/scripts/rCorrectScript.txt` to perform error correction. 
```
USAGE: sbatch rCorrectScript.txt <comma-sep list of Forward reads> <comma-sep list of Reverse reads>
```
The files began with a suffix: *trim.fastq.gz, and after running rCorrector with the above script, the output files will contain the suffix: *trim.cor.fq.gz. These files are located in: `/data/TrimCorrectedFastq/`. The “cor” indicates it is the corrected file.   
The output fastq files will have a different header that can mess up downstream assembly. You can adjust this by using a script developed by a Harvard Informatics group.     

Clone the repository of the Harvard Informatics group to get the python script: `git clone git@github.com:harvardinformatics/TranscriptomeAssemblyTools.git`. The script of interest is called: `FilterUncorrectabledPEfastq.py`. This is now in our `/scripts/` folder as well if you do not want to clone their repository.     
Execute this script using the batch submission script written: `/scripts/FilterUnfixableSeqs.txt`. The way to use this script follows. **Note** this script must be run once for each sample.
```
sbatch FilterUnfixableSeqs.txt <read1> <read2>
```
The output files will have the format: `unfixrm*[1 or 2].trim.cor.fq`. Ours are located in `/data/trimmed-corrected-filt/`.      
Conveniently, it will also output log files that provide summary statistics, and our logs are the following files located in `/logs/`:
rmunfixable_Mfav_DD_euk_1.log     
rmunfixable_Mfav_DD_euk_31.log     
rmunfixable_Mfav_HH_euk_33.log      
rmunfixable_Mfav_HH_euk_3.log     

### 3. Generate de novo Transcriptome with Trinity

Execute a de novo transcriptome assembly using **Trinity v2.8.5**. First, make sure Trinity is installed in the conda environment.     
Wrote a script to run trinity in `/scripts/trinity.sh`       
Execute with the following code. **Note:** Trinity can take comma-separated lists of the forward and reverse reads, which is what the long string of filenames are. This uses the .fastq files that have been trimmed and error-corrected. Our trimmed, error-corrected, and filtered fastq files used for the de novo Trinity assembly are in `/data/trimmed-corrected-filt/` folder. 
```
sbatch trinity.sh unfixrm_M_faveolata__Mfav_DD_euk_1_2505__L1_TGACCA_L001_1.trim.cor.fq,unfixrm_M_faveolata__Mfav_DD_euk_31_2506__L1_ACAGTG_L001_1.trim.cor.fq,unfixrm_M_faveolata__Mfav_HH_euk_3_2502__L1_ATCACG_L001_1.trim.cor.fq,unfixrm_M_faveolata__Mfav_HH_euk_33_2503__L1_CGATGT_L001_1.trim.cor.fq unfixrm_M_faveolata__Mfav_DD_euk_1_2505__L1_TGACCA_L001_2.trim.cor.fq,unfixrm_M_faveolata__Mfav_DD_euk_31_2506__L1_ACAGTG_L001_2.trim.cor.fq,unfixrm_M_faveolata__Mfav_HH_euk_3_2502__L1_ATCACG_L001_2.trim.cor.fq,unfixrm_M_faveolata__Mfav_HH_euk_33_2503__L1_CGATGT_L001_2.trim.cor.fq trinity_out_dir
```
The output directory is `/output/trinity_out_dir`. The transcriptome is `Trinity.fasta`, which has been saved to `/output/Transcriptomes/Trinity.fasta`.     
Generate some statistics on the transcriptome using a perl script that comes with Trinity. Execute the following command:
```
perl /vortexfs1/home/cbecker/.conda/envs/trinity/opt/TRINITY_HOME/util/TrinityStats.pl Trinity.fasta > Trinity.assemblystats.txt
```
Taking a look at the Trinity.assemblystats.txt will tell you about the N50, number of ‘genes’ and transcripts (often many transcripts or isoforms per ‘gene’). This document is now at `/output/Transcriptomes/Trinity.assemblystats.txt`. 

### 4. Filter the de novo Transcriptome

The Trinity.fasta transcriptome is very lengthy, with over 200,000 transcripts. By filtering this, we can get closer to the Daniels et al (2015) reported transcriptome length of ~67,000 genes. 

#### 4a. Select only the longest isoform

Use a perl script included with the Trinity package to eliminate all isoforms of a gene except the longest isoform. This will not use much compute power, so we ran this using `srun` rather than submitting a batch script. 
```
$ perl /vortexfs1/home/cbecker/.conda/envs/trinity/opt/TRINITY_HOME/util/misc/get_longest_isoform_seq_per_trinity_gene.pl Trinity.fasta > Trinity.longest.fasta
```
The output is a new transcriptome, which is saved at: `/output/Transcriptomes/Trinity.longest.fasta`    
We also generated summary statistics similarly to Step 3, which generated: `/output/Transcriptomes/Trinity.longest.assemblystats.txt`

#### 4b. Filter further based on expression values

Use **Salmon v1.0.0**, a pseudo-aligner to map reads back to the de novo Transcriptome (Trinity.longest.fasta) and keep genes in the transcriptome if they are mapped to at a TPM > 1. We used the pseudo-aligner, salmon, which is quick. First you must index the transcriptome for use with salmon, which is implemented in the `/scripts/salmon-index.sh` script. Next, generate a quantification of number of reads mapping to the transcriptome using the `/scripts/salmon.sh` script.    
Execute with:
```
sbatch salmon-index.sh     
```
The index command will output a folder, saved at: `/output/trinity-longest-index/` which contains indexing information that will be called by the salmon.sh script, so wait until the previous script is finished before running salmon:
```
sbatch salmon.sh
```
This will output a directory called `~/salmon/`. Within it is a `quant.sf` file that will be used for filtering the transcriptome.      

Use perl scripts included in the Trinity package to filter the transcriptome. Begin an interactive computing run with `srun`, then run the following. **Note** these scripts can be hard to find. Make sure the path works by running the scripts with no flags or data to verify the help screen pops up. Ex: `perl /vortexfs1/home/cbecker/.conda/envs/trinity/opt/TRINITY_HOME/util/abundance_estimates_to_matrix.pl` should open up the help file. If no directory is found, you may have to go back and find where the perl script is on your computer.    
```
perl /vortexfs1/home/cbecker/.conda/envs/trinity/opt/TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method salmon --cross_sample_norm TMM --out_prefix salmon_quant --gene_trans_map none ./salmon/quant.sf 
-reading file: ./salmon/quant.sf
```
The output files from this are saved at:    
1. /output/salmon/salmon_quant.isoform.counts.matrix     
2. /output/salmon/salmon_quant.isoform.TPM.not_cross_norm      

Next, actually filter the transcriptome at a threshold of 1 TPM. After this command runs, it will tell you how many transcripts were kept. Record this information.     
```
perl /vortexfs1/home/cbecker/.conda/envs/trinity/opt/TRINITY_HOME/util/filter_low_expr_transcripts.pl --matrix salmon_quant.isoform.TPM.not_cross_norm --transcripts Trinity.longest.fasta --min_expr_any 1 > Trinity.longest.TPMfilt.fasta
```
**Note** many sequences will be removed, and it will print information similar to:
```
-excluding TRINITY_DN57235_c1_g1_i1, max_expr: 0.246791 < 1
-excluding TRINITY_DN149853_c0_g1_i1, max_expr: 0.493582 < 1
-excluding TRINITY_DN121029_c0_g1_i1, max_expr: 0.370187 < 1
-excluding TRINITY_DN125156_c0_g1_i1, max_expr: 0.256289 < 1
```
When we ran this, it retained **66044** / 197483 = 33.44% of total transcripts.
This is close to the reported transcriptome length of Daniels et al (2015) = 67,593 genes.      

### 5. Map reads back to filtered transcriptome

We will use the package **Bowtie2 v2.3.5** to map reads from our filter and error-corrected fastq files for each sample back to the de novo transcriptome that has been filtered. This consists of a step where the de novo transcriptome is indexed (similar to what we did with salmon) and then reads are mapped. Indexing can be done by executing the script, `/scripts/index.TPMfilt.sh` :
```
sbatch index.TPMfilt.sh
```
The output of this command will create a list of files. In our case it generated 6 files:
1. denovoTrinity.longest.TPMfilt.Index.1.bt2     
2. denovoTrinity.longest.TPMfilt.Index.2.bt2     
3. denovoTrinity.longest.TPMfilt.Index.3.bt2      
4. denovoTrinity.longest.TPMfilt.Index.4.bt2     
5. denovoTrinity.longest.TPMfilt.Index.rev.1.bt2     
6. denovoTrinity.longest.TPMfilt.Index.rev.2.bt2       
These output files are in a folder called: `/output/TrinityIndex_TPMfilt_forBT2/`    

Then, execute the bowtie command that aligns reads back to the de novo transcriptome. Execute this with the following command:
```
sbatch bowtie.longest.TPMfilt.sh denovoTrinity.longest.TPMfilt.Index
```
Bowtie2 generates one `.sam` file for each sample with this script. These files are located in `/output/Longest.TPMfilt.SAMBAMfiles/`. Unfortunately, `.sam` files are not as useful for getting counts for differential expression analysis, so next, we use two tools to gather count files in preparation for differential expression analysis.          

### 6. SAM to BAM file conversion   

Bowtie2 generates SAM files, but we need BAM files for quantifying read counts, so we will use **SAMtools v1.9** to make the conversion and to sort the newly-created BAM files. Make sure SAMtools is added to the conda environment prior to running the scripts. The scripts needed for this is `/scripts/samtools.Longest.sh`, which can be executed with the command:
```
sbatch samtools.Longest.sh
```
The output files will have the suffixes `*.longfilt.bt2.bam` and `*longfilt.bt2.sorted.bam`. All output files are located in `/output/Longest.TPMfilt.SAMBAMfiles/`.

### 7. Generate count files and create a count matrix

We will use the package, **eXpress v1.5.1**, to generate count files for each of the sorted BAM files. Add eXpress to the conda environment if you haven’t, and run the following command to activate the eXpress script we wrote: `/scripts/eXpress.Longest.sh`
```
sbatch eXpress.Longest.sh
```
This will output 4 directories (one for each sample). For example, ours outputted the following directories:     

Within these directories is a `results.xprs` file. Copy over and rename the `results.xprs` file from each directory into one directory. For example, we made this directory: `/output/counts_longest_TPMfilt`. We populated this directory with one results file for each sample: 
Mfav_DD_1_counts.xprs      
Mfav_DD_31_counts.xprs      
Mfav_HH_33_counts.xprs      
Mfav_HH_3_counts.xprs     

### 8. Generate a count table with Python 3 in Jupyter notebooks
In the Jupyter notebook to work with the count files in `/output/counts_longest_TPMfilt/`: The notebook is located at `jupyter-notebooks/express_results_to_count_table_FilteredTXM.ipynb`

Run the code within that jupyter notebook from the /data/counts_longest_TPMfilt/ folder to generate the following count files in `/output/counts_longest_TPMfilt/`:      
Mfav_counts_all_filttxm.tab    
Mfav_counts_all_filttxm_decimals.tab     
The one made up of integers (not the decimals) tab file will be used in edgeR analyses.     

### 9. Calculate Differentially Expressed Genes using edgeR
We ran **edgeR v3.28.0** locally in RStudio running **R v.3.6.1**.     
Using `scp`, copy over the Mfav_counts_all_filttxm.tab file to a directory also containing the R script that can be run to generate edgeR statistics (/scripts/edgeR_filteredtxm_coral_FinalProject.R).  This R script should be able to be run on an HPC as well to generate the files.     
Running the R script with the input `Mfav_counts_all_filttxm.tab` file will generate the following tables, which are saved in the folder `/output/counts_longest_TPMfilt/` as:    
1. DEgenes_coral_filttxm.tsv     
2. DEgenes_pval0.1_coral_filttxm.tsv      
The `pval0.1` table includes genes that were differentially expressed at a p-value less than 0.1, and this is what will be needed for the final comparison. 

### 10. Annotate the de novo Transcriptome using SwissProt and TrEMBL

We are annotating the original de novo Transcriptome (Trinity.fasta) to identify putative functional genes. We will use **BLAST v2.9.0**, using the `blastx` command.        
The reviewed Uniprot and SwissProt databases were downloaded from the internet (https://www.uniprot.org/downloads) and unzipped. The unzipped file is located at `/output/swissprot_TrEMBL_DB/uniprot_sprot.fasta`. 

Use the following command from the BLAST program to generate a blast database:
```
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out SwissProtDB
```
The database this script created is located in the `/output/swissprot_TrEMBL_DB/` folder. 

Within the folder containing the SwissProtDB database files, run the script, `/scripts/Annotate.sh`, using the following code:
```
sbatch Annotate.sh
```
This code will generate a tab-delimited file containing all the annotations. Ours output file is located at: `/output/swissprot_TrEMBL_DB/denovo_annotateSwissProt.tab`    

### 11. Gene Ontology (GO) annotation
After retrieving annotations from SwissProt and TrEMBL, we retrieved the gene ontology annotations using an online interface: Uniprot Knowledgebase online, https://www.uniprot.org/uploadlists/

Copy the column of swissprot ID’s, or “entry names” from the Step 10 output file: `/output/swissprot_TrEMBL_DB/denovo_annotateSwissProt.tab` and placed it into a new file. For us, this was /output/swissprot_TrEMBL_DB/for-knowlKB5.txt     

Used `scp` command to copy `for-knowlKB5.txt` to local computer as a txt file. Uploaded file to Uniprot Knowledgebase online interface to get gene ontology annotations. 

Selected columns to include in the output file and under the gene ontology header, selected Gene ontology (biological process), Gene ontology (cellular component), Gene ontology (GO), Gene ontology (molecular function), and Gene ontology IDs. 

Selected the download option and chose to “download as a tab separated file”. This file was then uploaded to poseidon using `scp` command, and is now located at `/output/swissprot_TrEMBL_DB/uniprotKB.tab`     

The output was also downloaded as a .fasta file, and is located at: `/output/swissprot_TrEMBL_DB/uniprotKB.fasta`

### 12. Merge the SwissProt/Trembl annotations with GO annotations

The Jupyter Notebook, `/jupyter-notebooks/Combine_BLASTX-and-GO-Annotations.ipynb`, contains annotated code to merge the tab-delimited annotation outputs from Step 10 and 11. The code requires the following tables:    
1. `/output/swissprot_TrEMBL_DB/denovo_annotateSwissProt.tab`        
2. `/output/swissprot_TrEMBL_DB/uniprotKB.tab`         

It will generate or output the following table:     
`/output/swissprot_TrEMBL_DB/BLASTX_and_GO_merged.tab`

### 13. Generate Custom Databases for Coral and Symbiodiniaceae to sort the reference transcriptome

#### 13a. Generate a coral database
Verify Blast is installed into the active conda environment with `conda list`    
Retrieved the following fastq files manually from the internet, and use the `scp` command to copy them to a new folder. Then use the `gunzip *.gz` command to unzip the following files:     
The unzipped files are here:
1. *Orbicella faveolata*: /output/transcriptome-database/coralfastas/GCF_002042975.1_ofav_dov_v1_rna.fna
2. *Acropora hyacinthus*: /output/transcriptome-database/coralfastas/A_hyacinth_GDIF01.1.fsa_nt -out ahyadb.fna
3. *Hydra vulgaris*: /output/transcriptome-database/coralfastas/H.vulgaris_GCF_000004095.1_Hydra_RP_1.0_rna.fna
4. *Acropora digitifera*: /output/transcriptome-database/coralfastas/adi_transcriptome_assembly.v1.fa
5. *Exaiptasia pallida*: /output/transcriptome-database/coralfastas/E.pallida_GCF_001417965.1_Aiptasia_genome_1.1_rna.fna
6. *Stylophora pistillata*: /output/transcriptome-database/coralfastas/GCF_002571385.1_Stylophora_pistillata_v1_rna.fna
7. *Nematostella vectensis*: /output/transcriptome-database/coralfastas/N.vectensis_GCF_000209225.1_ASM20922v1_rna.fna

Used the following commands to create a database for each individual:
```
makeblastdb -dbtype nucl -in GCF_002042975.1_ofav_dov_v1_rna.fna -out ofavdb.fna -parse_seqids -title "ofavDB"
makeblastdb -dbtype nucl -in A_hyacinth_GDIF01.1.fsa_nt -out ahyadb.fna -parse_seqids -title "ahyaDB"
makeblastdb -dbtype nucl -in H.vulgaris_GCF_000004095.1_Hydra_RP_1.0_rna.fna -out HvulgDB.fna -parse_seqids -title "hvulDB"
makeblastdb -dbtype nucl -in adi_transcriptome_assembly.v1.fa -out adigiDB.fna -parse_seqids -title "adigiDB"
makeblastdb --dbtype nucl -in E.pallida_GCF_001417965.1_Aiptasia_genome_1.1_rna.fna -out ApallDB.fna -parse_seqids -title "ApallDB"
makeblastdb -dbtype nucl -in N.vectensis_GCF_000209225.1_ASM20922v1_rna.fna -out NvectDB.fna -parse_seqids -title "NvectDB"
makeblastdb -dbtype nucl -in GCF_002571385.1_Stylophora_pistillata_v1_rna.fna -out SpistDB.fna -parse_seqids -title "SpistDB"
```
After creating each database, use the following code to combine the above databases together. **Note**: ahyadb.fna did not work to join the other databases. 
```
blastdb_aliastool -dblist "ApallDB.fna HvulgDB.fna NvectDB.fna SpistDB.fna adigiDB.fna ofavdb.fna" -dbtype nucl -title "Cnidarian Database" -out cnidarian_DB
```

#### 13b. Generate a Symbiodinium database
Retrieved the files for 4 clades (now Genera) in the family Symbiodiniaceae from the following websites: http://medinalab.org/zoox/, https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=GAFO01, https://www.ncbi.nlm.nih.gov/Traces/wgs/?val=GAFP01     
After unzipping the files, they were placed in the following locations:  
1. *Clade A*: /output/transcriptome-database/algaefastas/kb8_assembly.fasta     
2. *Clade B*: /output/transcriptome-database/algaefastas/mf105_assembly.fasta     
3. *Clade C*: /output/transcriptome-database/algaefastas/GAFO01.1.fsa_nt      
4. *Clade D*: /output/transcriptome-database/algaefastas/GAFP01.1.fsa_nt 

The following code was used to make Blast databases from the fasta files. Then they were combined into one database.
```
makeblastdb -dbtype nucl -in kb8_assembly.fasta -out cladeA -parse_seqids -title "CladeA"
makeblastdb -dbtype nucl -in mf105_assembly.fasta -out cladeB -parse_seqids -title "CladeB"
makeblastdb -dbtype nucl -in GAFO01.1.fsa_nt -out cladeC -parse_seqids -title 
"CladeC"
makeblastdb -dbtype nucl -in GAFP01.1.fsa_nt -out cladeD -parse_seqids -title 
"CladeD"

blastdb_aliastool -dblist "cladeA cladeB cladeC cladeD" -dbtype nucl -title "All-Symbiodinium-database" -out algaeDB_all
```
The main output file for annotation is the /output/transcriptome-database/algaefastas/algaeDB_all.nal

#### 13c. Map reference transcriptome against custom databases

We queried our denovoAssembly_Trinity.fasta (copied to each `~/algaefastas/` and `~/coralfastas/` folder) file to the custom coral database, cnidarian_DB, using blastn (see coralAnnotate.sh) and to the custom algae database algae_DB, using blastn (see algaeAnnotate.sh).          
Queries were done using the BLASTN function written in `/scripts/coralAnnotate.sh` and `/scripts/algaeAnnotate.sh` scripts. Make sure to run these in the folder containing the de novo transcriptome and all database files. Execute with:
```
sbatch coralAnnotate.sh
sbatch algaeAnnotate.sh
```
This produced the files `/output/transcriptome-database/coralAnnotate.tab` and `/output/transcriptome-database/algaeAnnotate.tab`, respectively.

### 14. Assign contigs to either **Coral** or **Symbiodiniaceae**
The two files, `/output/transcriptome-database/coralAnnotate.tab` and `/output/transcriptome-database/algaeAnnotate.tab` now need to be sorted based on the e-value and bitscore. This was done using python in Jupyter Notebook, `/jupyter-notebooks/Split-algae-coral.ipynb`.     
E-values were compared between the coralAnnotate.tab and algaeAnnotate.tab files, and the file with the lower e-value was assigned to that organism -- ex contig A: coral e-value = 0.1e^-8  < contig A : algae e-value 0.2e^-4, then save contig to coral.tab file, save nothing to algae.tab file. If the e-values are equal, then the larger bitscore assigns the contig to that organism. If the bitscore and contig are equal, then they are assigned to neither (n = 15) organism and were not included in further data analysis.    
The output files were:    
1. `/output/transcriptome-database/algae.tab`      
2. `/output/transcriptome-database/coral.tab`

### 15. Combine major output files and compare to Data Sheet 4, from Daniels et al. (2015)

The jupyter notebook, /jupyter-notebooks/Final-Comparison.ipynb, contains code to combine the major files are needed for this analysis:         
1. /output/FinalDE_Table/BLASTX_and_GO_merged.tab 
2. /output/FinalDE_Table/algae-Copy1.tab
3. /output/FinalDE_Table/coral-Copy1.tab
4. /output/FinalDE_Table/DEgenes_pval0.1_coral_filttxm.tsv
5. /output/FinalDE_Table/DataSheet4_Danielsetal2015.txt

The main output file from the Jupyter Notebook is: `/output/FinalDE_Table/DataSheet4Replica.tab`
