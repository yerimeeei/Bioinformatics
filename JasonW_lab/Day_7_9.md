# Day 7 & 8 & 9 RNA sequencing of LINE 1 elements using GEO L1TD1 samples

### Duty for Day 6

- [x]  RPG Information session 2:30 - 4:15
- [x]  Try to run STAR on one FASTQ sample
- [x]  Try to run featureCounts of resulted BAM file
- [ ]  L1TD1 gene → for further studies, need to choose whether I want to do for other lines or only L1 → Study and **Read literature of L1 cancer basic** and make questions to ask.

## Links

- Batch Effect ([using combat in R](https://www.rdocumentation.org/packages/sva/versions/3.20.0/topics/ComBat))
- [GEO data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164541) that has normal, adenoma, and primary tumor
- Study [STAR](https://github.com/alexdobin/STAR) aligner (STAR → BAM → count ? (fc_output in R))
    - [STAR manual](https://gensoft.pasteur.fr/docs/STAR/2.7.3a/STARmanual.pdf)

# STAR alinger workflow

## [Command](https://m.blog.naver.com/naturelove87/221371566835) [Notes](https://github.com/stephenturner/oneliners/blob/master/README.md#awk--sed-for-bioinformatics)

## 2. Mapping reads to the genome

- FASTQ file

```bash
/storage2/temp/labs/jwlab/adenoma/fastq
```

- RUN STAR

```bash
/storage2/jwlab/yerim/STAR --runThreadN 5 
	--genomeDir /storage2/jwlab/yerim/genome_index/ 
	--readFilesIn /storage2/temp/labs/jwlab/adenoma/fastq/SRR3147729_1.fastq.gz 
		/storage2/temp/labs/jwlab/adenoma/fastq/SRR3147729_2.fastq.gz
	--readFilesCommand zcat
```

- resulted files
    - Aligned.out.sam
    - SJ.out.tab
    - Log.out
    - Log.final.out
    - Log.progress.out
- [Convert SAM to BAM](https://bio-info.tistory.com/26)

```bash
/storage2/jwlab/yerim/samtools view Aligned.out.sam > Aligned.out.bam
```

- [featureCounts](https://subread.sourceforge.net/SubreadUsersGuide.pdf)

```bash
/storage2/jwlab/yerim/featureCounts -p -O -T 5 
	-a /storage2/jwlab/sandy/data/gencode.v43.chr_patch_hapl_scaff.basic.annotation.gtf 
	-o featureCounts_output.txt /storage2/jwlab/yerim/rnaseq/results/Aligned.out.bam
```

<img width="564" alt="Screenshot_2023-06-13_at_10 21 29_PM" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/26d4ca82-eb99-48db-a5b5-cd974fccb32d">

### Lab Meeting

- ****[Pan-cancer analysis](https://www.nature.com/articles/s41588-023-01349-3) identifies tumor-specific antigens derived from transposable elements**

# Literature 3: LINE-1 in cancer basic

## **Literature Review**

### Structure and Retrotransposition Process of L1

- retrotransposition-competent human L1

<img width="685" alt="Screenshot_2023-06-13_at_11 24 59_AM" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/7ff5b917-afe3-4cf8-8e2a-b1360a0ab779">

- 5’ UTR, ORF1, ORF2, 3’ UTR (ending with poly (A) tail)
    - 5’ UTR: protect two internal promoters (sense, antisense)
        - sense promoter: binds RNA polymerase II and initiates L1 transcription from the 5’ end to 3’ end.
        - antisense promoter: give rise to chimeric RNAs transcribed partially from L1 5’ UTR and partially from neighboring sequences

# [GEO data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164541) for analysis of L1TD1 expression in colorectal cancer (CRC)

### To-do

- [x]  fastq file download
- [x]  BAM --> featureCounts
- [ ]  FeatureCounts values normalization
- [ ]  Normal, Adenoma, Primary Graph --> L1TD1 gene expression
- fast-dump

```bash
/storage2/jwlab/sandy/sratoolkit.3.0.1-ubuntu64/bin/fastq-dump
```

1. [Download](https://erilu.github.io/python-fastq-downloader/) sample data: download SRR_Acc_List and make download script

```bash
awk '{print"fastq-dump --split-3 --gzip "$0}' SRR_Acc_List.txt > srr_download.sh
```

- sir_download.sh

```bash
fastq-dump --split-3 --gzip SRR13401190
fastq-dump --split-3 --gzip SRR13401191
fastq-dump --split-3 --gzip SRR13401192
fastq-dump --split-3 --gzip SRR13401193
fastq-dump --split-3 --gzip SRR13401194
fastq-dump --split-3 --gzip SRR13401195
fastq-dump --split-3 --gzip SRR13401196
fastq-dump --split-3 --gzip SRR13401197
fastq-dump --split-3 --gzip SRR13401198
fastq-dump --split-3 --gzip SRR13401199
fastq-dump --split-3 --gzip SRR13401200
fastq-dump --split-3 --gzip SRR13401201
fastq-dump --split-3 --gzip SRR13401202
fastq-dump --split-3 --gzip SRR13401203
fastq-dump --split-3 --gzip SRR13401204
```

2. Start the job session to download fastq.gz file

```bash
parallel -j 5 < srr_download.sh
```

3. Run STAR

```bash
fq1=$1
fq2=$2
name=$3

mkdir /storage2/jwlab/yerim/rnaseq/BAM/$name
~/storage2/jwlab/yerim/STAR --runThreadN 5 
	--genomeDir /storage2/jwlab/yerim/genome_index/ 
	--readFilesIn $fq1 $fq2 --readFilesCommand zcat 
	--outFileNamePrefix /storage2/jwlab/yerim/rnaseq/BAM/$name/
	--outSAMtype BAM SortedByCoordinate --outSAMstrandField intronMotif --outSAMattributes All --outSAMattrIHstart 0 --outFilterMultimapNmax 100 --outFilterScoreMinOverLread 0.4 --outFilterMatchNminOverLread 0.4 --clip3pNbases 0 --winAnchorMultimapNmax 100 --alignEndsType EndToEnd --alignEndsProtrude 100 DiscordantPair --chimSegmentMin 250 --sjdbOverhang 249 --twopassMode Basic
```

- fastq_list.sh

```bash
bash run_star.sh /storage2/jwlab/yerim/rnaseq/fastq/SRR13401190_1.fastq.gz /storage2/jwlab/yerim/rnaseq/fastq/SRR13401190_2.fastq.gz SRR13401190
```

4. Start job session to run star behind the screen

```bash
parallel -j 5 < fastq_list.sh
```

5. Make fc.sh

```bash
name=$1
bam=$2

mkdir /storage2/jwlab/yerim/rnaseq/fc_output/$name
/storage2/jwlab/yerim/featureCounts -T 5 
	--verbose -O -M -a /storage2/jwlab/yerim/gencode.L1.gtf 
	-o /storage2/jwlab/yerim/rnaseq/fc_output/$name/fc_result.txt $bam

/storage2/jwlab/yerim/featureCounts -p -O -T 5
```

6. Make bam_list.sh and bam_fc_list.sh

```bash
awk '{print"/storage2/jwlab/yerim/rnaseq/BAM/"$0"/Aligned.sortedByCoord.out.bam"}' /storage2/jwlab/yerim/SRR_Acc_List.txt > bam_list.sh
```

```bash
awk '{print"bash fc.sh " $0 " /storage2/jwlab/yerim/rnaseq/BAM/" $0 "/Aligned.sortedByCoord.out.bam"}' /storage2/jwlab/yerim/SRR_Acc_List.txt > bam_fc_list.sh
```

- bam_fc_list.sh

```bash
bash fc.sh SRR13401190 /storage2/jwlab/yerim/rnaseq/BAM/SRR13401190/igv.bam
```

7. Run featureCounts

```bash
parallel -j 5 < bam_fc_list.sh
```

<img width="322" alt="Screenshot_2023-06-15_at_10 08 25_PM" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/4580a3e6-86d3-4d0c-a9b1-1ad0099aa2c8">

<img width="321" alt="Screenshot_2023-06-15_at_10 08 19_PM" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/e333189e-4b12-4b8b-b760-4a90d575477f">

## FeatureCounts value normalization

- [TPM, FPKM, or Normalized Counts](https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-021-02936-w)
- [RNA seq analysis](https://sbc.shef.ac.uk/workshops/2019-01-14-rna-seq-r/rna-seq-preprocessing.nb.html) in R
    - TPM R packages
    - [FPKM in R packages](https://cran.r-project.org/web/packages/countToFPKM/countToFPKM.pdf)

### TPM (Transcript per Million)

- Units: PM normalizes read counts by the total number of mapped reads in millions and the length of each transcript in kilobases, resulting in a **metric of transcripts per million.**
- Handling of transcription length: TPM normalization takes into account the length of each transcript, so that **longer transcripts are given less weight than shorter ones**.
- Handling of multi-mapping reads: TPM normalization assigns multi-mapping reads to **all transcripts** that they map to, **weighted by the length of each transcript.**
- Scaling: TPM normalization scales the expression values so that the sum of **TPM values across all genes is equal to a constant** (usually 1 million), allowing for **direct comparison of expression levels between samples**.

RPKM/FPKM methods account for, firstly, the **library size**, and secondly, the **gene lengths**. TPM also controls for both the library size and the gene lengths, however, with the TPM method, the read counts are first normalized by the gene length (per kilobase), and then gene-length normalized values are divided by the sum of the gene-length normalized values and multiplied by 10^6. Thus, the sum of normalized values for TPM will always be equal to 10^6 for each library, while the sum of RPKM/FPKM values do not sum to 10^6. Therefore, it is easier to interpret TPM values than RPKM/FPKM values.

### FPKM (Fragments Per kilobase of transcript per Million reads mapped)

- Units: normalizes read counts by the total number of mapped reads in millions and the length of each transcript in kilobases, resulting in a **metric of fragments per kilobase of transcript per million mapped reads.**
- Handling of transcription length: FPKM normalization also takes into account transcript length, but it **divides by the length of the transcript in kilobases to account for the fact that longer transcripts will have more fragments**, even if they are expressed at the same level as a shorter transcript.
- Handling of multi-mapping reads: FPKM normalization assigns multi-mapping reads **proportionally to the number of transcripts** they map to, **without weighting by transcript length**.
- Scaling: FPKM normalization **does not have a scaling factor**, so it is **not directly comparable between samples**.

E.g. fpkmheatmap

```r
# colorectal cancer
file_readcounts <- system.file("extdata", "RNA-seq.read.counts.csv",
                           package = "countToFPKM")
file_annotations <- system.file("extdata", "Biomart.annotations.hg38.txt",
                            package = "countToFPKM")
file_sample.metrics <- system.file("extdata", "RNA-seq.samples.metrics.txt",
                            package = "countToFPKM")

# Import the read count matrix data into R.
counts <- as.matrix(read.csv(file.readcounts))

# Import feature annotations.
# Assign feature length into a numeric vector.
gene_annotations <- read.table(file.annotations, sep="\t", header=TRUE)
featureLength <- gene.annotations$length
    
# Import sample metrics.
# Assign mean fragment length into a numeric vector.
samples_metrics <- read.table(file.sample.metrics, sep="\t", header=TRUE)
meanFragmentLength <- samples.metrics$meanFragmentLength
    
# Return FPKM into a numeric matrix.
fpkm_matrix <- fpkm (counts, featureLength, meanFragmentLength)

# Plot log10(FPKM+1) heatmap of top 30 highly variable features
    fpkmheatmap(fpkm_matrix, topvar=30, showfeaturenames=TRUE, return_log = TRUE)
```
