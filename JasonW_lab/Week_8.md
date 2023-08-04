# Week 8 L1TD1 GSEA, DNA Methylation, ATAC-seq

## DNA Methylation

- [GEO data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48684)
- [GEO dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534) sample
- ****[A cross-package Bioconductor workflow for analysing methylation array data](https://bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html)****
- DNA methylation is the most well-characterised epigenetic mark in humans. It is defined as the **addition of a methyl (CH3) group to DNA and in mammalian cells occurs primarily at the cytosine of cytosine-guanine dinucleotides (CpG)**.
- DNA methylation can modify the function of regulatory elements and gene expression and is therefore integral to normal human development and biological functioning. Perturbations to normal DNA methylation patterns can lead to dysregulation of cellular processes and are linked with disease. Widespread aberrations in DNA methylation are a well-established hallmark of many cancers [[1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5055731/#CR1)] and a growing body of literature shows a role for DNA methylation in the aetiology of other complex human diseases.

### CPG Island

**[Critical evaluation of the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5055731/)**

- **EPIC**
    - EPIC array is a significant improvement over the HM450 array, with increased genome coverage of regulatory regions and high reproducibility and reliability, providing a valuable tool for high-throughput human methylome analyses from diverse clinical samples.
- **HM450**

### Bisulfite sequencing

- **WGBS (whole-genome bisulphite sequencing)**
    - The current ‘gold standard’ technique for fine mapping of methylated cytosines is **whole-genome bisulphite sequencing (WGBS)** [[5](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5055731/#CR5)]. This is based on the **treatment of genomic DNA with sodium bisulphite, which converts unmethylated cytosines to uracils while leaving methylated cytosines unchanged, followed by whole-genome sequencing** [[6](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5055731/#CR6)]. WGBS has been successfully applied to a range of biological tissues and cell lines to provide a complete map of the ~28 million CpG sites in the human genome [[7](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5055731/#CR7)].
- **RRBS (reduced representation bisulphite sequencing)**
    
    ![methylation_plot](https://github.com/yerimeeei/Bioinformatics/assets/134043926/fd40bedd-eb7c-4bd6-b3ad-02f8c58166d2)

    ![methylation_plot_cpg_island](https://github.com/yerimeeei/Bioinformatics/assets/134043926/c3c4a335-a8cd-4f95-bf15-c0da874089c4)


**For all CpG probes boxplot**

```r
cpg_plot <- ggplot(beta_value, aes(x = sample, y = value, group = sample, fill = sample))
cpg_plot <- cpg_plot + geom_boxplot(width=0.05) + theme_classic() + scale_x_discrete(limits=c("normal colon", "adenoma", "CRC"))
cpg_plot <- cpg_plot + labs(title="DNA Methylation Normal colon vs Adenoma vs CRC", x="Sample Identities", y = "DNA methylation Beta Values") 
cpg_plot
```

![cpg_probes_boxplot](https://github.com/yerimeeei/Bioinformatics/assets/134043926/a880aceb-5543-4d55-bf1f-fb9ff1c13dfa)

- comment from Jason
    - I think to look at it properly and validate that there is no batch effect, we need to separate the CpG island probes and the intergenic CpGs (ideally in or near L1 elements). In cancer often CpG islands are hypermethylated and L1s are demethylated. If we see the opposite for the adenoma, then I will believe the data.

### **CpG**

- **CpG islands** are regions of DNA that are rich in CpG sites and are commonly found near the promoter regions of genes. **CpG island probes** are designed to detect DNA methylation at CpG islands.
- **Intergenic CpGs**, on the other hand, are CpG sites located in non-coding regions of the genome.

**STEP**

1. **island_cpg_probe_list (Linux)**

![Screenshot_2023-07-31_at_4 29 00_PM](https://github.com/yerimeeei/Bioinformatics/assets/134043926/09da8b72-242f-49f6-9363-a992dfc8b25a)

- [USCS Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1666944150_RS3FVaLyGZBjXljn1WNCNdFVt0h0&clade=mammal&org=Human&db=hg38&hgta_group=regulation&hgta_track=gtexGeneV8&hgta_table=0&hgta_regionType=genome&position=chr2%3A25%2C160%2C915-25%2C168%2C903&hgta_outputType=primaryTable&hgta_outFileName=cpg_island_hg38)

```bash
#head cpg_island_hg38
#tail -n+2 cpg_island_hg38 | head
#tail -n+2 cpg_island_hg38 | awk '{print$2"\t"$3"\t"$4}' > cpg_island_hg38.bed

wc -l HumanMethylation450_probe.txt
wc -l GPL13534_HumanMethylation450_15017482_v.1.1.csv
wc -l GSE48684_Matrix_signal_intensities_1.txt

tail -n+2 HumanMethylation450_probe.txt | head
tail -n+2 HumanMethylation450_probe.txt | cut -f 1,10 | head

awk -F'\t' '{if($26=="Island")print $0}' HumanMethylation450_probe.txt | head
awk -F'\t' '{if($26=="Island")print $0}' HumanMethylation450_probe.txt | cut -f 1 > island_cpg_probe_list

#rm cpg_island_hg38
#rm cpg_island_hg38.bed
```

2. **l1_probe_list (Linux)**

<img width="1454" alt="Screenshot_2023-07-31_at_4 25 37_PM" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/45adb495-cd4f-49b5-aa97-98be3ac1a1a2">

- [UCSC table browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1666944150_RS3FVaLyGZBjXljn1WNCNdFVt0h0&clade=mammal&org=Human&db=hg38&hgta_group=rep&hgta_track=cpgIslandExt&hgta_table=0&hgta_regionType=genome&position=chr2%3A25%2C160%2C915-25%2C168%2C903&hgta_outputType=primaryTable&hgta_outFileName=cpg_island_hg38)

```bash
cut -f 13 hg38_repeatmasker | sort | uniq
awk -F'\t' '{if($13=="L1")print $0}' hg38_repeatmasker | head
awk -F'\t' '{if($13=="L1")print $0}' hg38_repeatmasker > ../yerim/L1_list
awk '{print$6"\t"$7"\t"$8"\t"$11"\t"0"\t"$10}' L1_list | head
awk '{print$6"\t"$7"\t"$8"\t"$11"\t"0"\t"$10}' L1_list > L1_list.bed

tail -n+2 HumanMethylation450_probe.txt | awk -F'\t' '{print"chr"$12"\t"$34"\t"$35"\t"$1}' | head
tail -n+2 HumanMethylation450_probe.txt | awk -F'\t' '{if($26!="Island")print $0}' | awk -F'\t' '{print"chr"$12"\t"$34"\t"$35"\t"$1}' |  head
tail -n+2 HumanMethylation450_probe.txt | awk -F'\t' '{if($26!="Island")print $0}' | awk -F'\t' '{print"chr"$12"\t"$34"\t"$35"\t"$1}' > probe_coordinate.bed
#sort -k1,1 -k2,2n  probe_coordinate.bed | head

head probe_coordinate.bed
cut -f 1 probe_coordinate.bed | head
sort -k1,1 -k2,2n probe_coordinate.bed > probe_coordinate_sorted.bed
mv probe_coordinate_sorted.bed probe_coordinate.bed

conda install -c bioconda bedtools
bedtools intersect -wa -a probe_coordinate.bed -b ../L1_list.bed | head
#bedtools intersect -wa -a probe_coordinate.bed -wb -b ../L1_list.bed | head
bedtools intersect -wa -a probe_coordinate.bed -b ../L1_list.bed > l1_probe_list
```

- [bedtools explanation](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)
- [CpG island visualization](https://www.researchgate.net/figure/The-CpG-islands-and-CpG-island-shores_fig1_317256914)
3. **plot for each CpG island probes and the intergenic CpGs**

```r
#for CpG island probes
island_cpg_probe_list <- read.table("/storage2/jwlab/yerim/methylation/island_cpg_probe_list", quote="\"", comment.char="")
colnames(island_cpg_probe_list)[1] <- "ID_REF"
CpG_island <- merge(beta_value, island_cpg_probe_list, by = "ID_REF")
Cpg_island_plot <- ggplot(CpG_island, aes(x = sample, y = value, group = sample, fill = sample))
Cpg_island_plot <- Cpg_island_plot + geom_boxplot(width=0.05) + theme_classic() + scale_x_discrete(limits=c("normal colon", "adenoma", "CRC"))
Cpg_island_plot <- Cpg_island_plot + labs(title="DNA Methylation for CpG island probes Normal colon vs Adenoma vs CRC", x="Sample Identities", y = "DNA methylation Beta Values") 
Cpg_island_plot

#for intergenic CpGs
l1_probe_list <- read.delim("/storage2/jwlab/yerim/methylation/l1_probe_list", header=FALSE)
colnames(l1_probe_list)[4] <- "ID_REF"
CpG_intergenic <- merge(beta_value, l1_probe_list, by = "ID_REF")
library(ggbeeswarm)
Cpg_intergenic_plot <- ggplot(CpG_intergenic, aes(x = sample, y = value, group = sample, fill = sample))
Cpg_intergenic_plot <- Cpg_intergenic_plot + geom_boxplot(width=0.05) + theme_classic() + scale_x_discrete(limits=c("normal colon", "adenoma", "CRC"))
Cpg_intergenic_plot <- Cpg_intergenic_plot + labs(title="DNA Methylation for intergenic CpGs Normal colon vs Adenoma vs CRC", x="Sample Identities", y = "DNA methylation Beta Values") 
Cpg_intergenic_plot <- Cpg_intergenic_plot + geom_quasirandom(width=0.1) 
Cpg_intergenic_plot
```

![CpG_island_probes](https://github.com/yerimeeei/Bioinformatics/assets/134043926/33d6cc56-078c-449b-982c-147db7f9280e)

![intergenic_CpGs](https://github.com/yerimeeei/Bioinformatics/assets/134043926/17df4d3f-0987-42f1-8199-3d7ee06a6c92)

- ****[CpG island hypermethylation and tumor suppressor genes: a booming present, a brighter future](https://www.nature.com/articles/1205600)****
1. **cpg_island_probe_list (Linux)**

```bash
conda create -n bedtools
conda init bash
conda activate bedtools
conda install -c bioconda bedtools

tail -n+2 cpg_island_hg38 > cpg_island_list.bed
awk '{print substr($2,0,5)"\t"$3"\t"$4"\t"$5"\t"$6}' cpg_island_list.bed > island_list.bed
rm cpg_island_list.bed
sort -k1,1 -k2,2n island_list.bed > island_list_sorted.bed
mv island_list_sorted.bed island_list.bed

tail -n+2 HumanMethylation450_probe.txt | awk -F: 'match($1, "cg")' HumanMethylation450_probe.txt | awk -F'\t' '{print"chr"$12"\t"$34"\t"$35"\t"$1}' | head
tail -n+2 HumanMethylation450_probe.txt | awk -F: 'match($1, "cg")' HumanMethylation450_probe.txt | awk -F'\t' '{print"chr"$12"\t"$34"\t"$35"\t"$1}' > coordinate_probe.bed
sort -k1,1 -k2,2n coordinate_probe.bed > coordinate_probe_sorted.bed
mv coordinate_probe_sorted.bed coordinate_probe.bed

bedtools intersect -wa -a coordinate_probe.bed -b island_list.bed | head
bedtools intersect -wa -a coordinate_probe.bed -b island_list.bed > cpg_island_probe_list
```

2. **L1_probe_list (Linux)**

```bash
head hg38_repeatmasker
awk -F'\t' '{if($13=="L1")print $0}' hg38_repeatmasker | head
awk -F'\t' '{if($13=="L1")print $0}' hg38_repeatmasker > L1_list
awk '{print$6"\t"$7"\t"$8"\t"$11"\t"0"\t"$10}' L1_list | head
awk '{print$6"\t"$7"\t"$8"\t"$11"\t"0"\t"$10}' L1_list > L1_list.bed

bedtools intersect -wa -a coordinate_probe.bed -b L1_list.bed | head
bedtools intersect -wa -a coordinate_probe.bed -b L1_list.bed > L1_probe_list
```

3. **plot for each CpG island probes and the intergenic CpGs**

```bash
#for CpG island probes
island_cpg_probe_list <- read.table("/storage2/jwlab/yerim/methylation/cpg_island_probe_list", quote="\"", comment.char="")
colnames(island_cpg_probe_list)[4] <- "ID_REF"
CpG_island <- merge(beta_value, island_cpg_probe_list, by = "ID_REF")
Cpg_island_plot <- ggplot(CpG_island, aes(x = sample, y = value, group = sample, fill = sample))
Cpg_island_plot <- Cpg_island_plot + geom_boxplot(width=0.05) + theme_classic() + scale_x_discrete(limits=c("normal colon", "adenoma", "CRC"))
Cpg_island_plot <- Cpg_island_plot + labs(title="DNA Methylation for CpG island probes Normal colon vs Adenoma vs CRC", x="Sample Identities", y = "DNA methylation Beta Values") 
Cpg_island_plot

#for intergenic CpGs
l1_probe_list <- read.delim("/storage2/jwlab/yerim/methylation/L1_probe_list", header=FALSE)
colnames(l1_probe_list)[4] <- "ID_REF"
CpG_intergenic <- merge(beta_value, l1_probe_list, by = "ID_REF")
Cpg_intergenic_plot <- ggplot(CpG_intergenic, aes(x = sample, y = value, group = sample, fill = sample))
Cpg_intergenic_plot <- Cpg_intergenic_plot + geom_boxplot(width=0.05) + theme_classic() + scale_x_discrete(limits=c("normal colon", "adenoma", "CRC"))
Cpg_intergenic_plot <- Cpg_intergenic_plot + labs(title="DNA Methylation for intergenic CpGs Normal colon vs Adenoma vs CRC", x="Sample Identities", y = "DNA methylation Beta Values") 
Cpg_intergenic_plot
#library(ggbeeswarm)
#Cpg_intergenic_plot <- Cpg_intergenic_plot + geom_quasirandom(width=0.1) 
#Cpg_intergenic_plot
```

![CpG_island_probes_fixed](https://github.com/yerimeeei/Bioinformatics/assets/134043926/186a27d9-8dd0-4abc-af84-cb9c7d3f2588)

![intergenic_CpGs_fixed](https://github.com/yerimeeei/Bioinformatics/assets/134043926/0ae533a2-b439-4120-9cc7-b9e959a11231)

- We wanted to confirm if the L1TD1 site is demethylated. Because we checked colon cancer, the L1TD1 promoter site is methylated but demethylated in normal colon. To check if the resulting data of DNA methylation in the L1TD1 site is valid, we separated CpG island probes and the intergenic CpGs. Then the data above showed that the data is not biased because it isn’t showing that everything was demethylated (then the demethylation of L1TD1 might not be real). Since it isn't and most promoters are probably hypermethylated.
- **Thus, demethylation shown to be a possible reason why L1TD1 gene expression is high in adenoma samples.**

## GSEA

- [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4941619/)
- [GEO data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76987)

```r
mapIDs[mapIDs$mapIds.org.Hs.eg.db..keys...rownames.corr_test...keytype....ENSEMBL...=="NA",]
unique(mapIDs$mapIds.org.Hs.eg.db..keys...rownames.corr_test...keytype....ENSEMBL...)
mapIDs_df <- mapIDs %>% group_by(mapIds.org.Hs.eg.db..keys...rownames.corr_test...keytype....ENSEMBL...) %>% summarise(freq=n())
nrow(corrtest[corrtest$gene_symbol=="0",])
MANE_annotated_genes <- read.delim("/storage2/jwlab/yerim/MANE_annotated_genes", header=FALSE, comment.char="#")
head(substr(MANE_annotated_genes$V18, 1, 15))
MANE_annotated_genes$label <- substr(MANE_annotated_genes$V18, 1, 15)
exist <- rownames(corr_test)[corr_test$gene_symbol!="0"]
length(exist)
not_exist <- exist[exist %in% MANE_annotated_genes$label]
not_exist <- exist[!exist %in% MANE_annotated_genes$label]
```

![Screenshot_2023-07-30_at_1 32 58_PM](https://github.com/yerimeeei/Bioinformatics/assets/134043926/3db0c5b4-c1a3-4f0b-b0be-3fceefcf535f)

[Genome UCSC](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1665408986_grAteYqzn1HrbldjFQDHFLuCkDSL&clade=mammal&org=Human&db=hg38&hgta_group=genes&hgta_track=mane&hgta_table=0&hgta_regionType=range&position=chr2%3A25%2C160%2C915-25%2C168%2C903&hgta_outputType=primaryTable&hgta_outFileName=)

- [x]  find another method for mapping gene ID
    - [gconvert](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html)
    - select mapID
    - use map ID

<img width="1457" alt="Screenshot_2023-07-31_at_2 46 09_PM" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/42d712e9-a3ab-4481-b66f-6ffcafcd87fe">

## ATAC-seq

- single-cell Analysis [paper](https://www.nature.com/articles/s41588-022-01088-x)
- [GEO dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201349) ([ATAC-seq](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201336))
- **[ATAC-seq data analysis: from FASTQ to peaks](https://yiweiniu.github.io/blog/2019/03/ATAC-seq-data-analysis-from-FASTQ-to-peaks/)**
    
    ![Screenshot_2023-07-31_at_3 19 50_PM](https://github.com/yerimeeei/Bioinformatics/assets/134043926/f85645ae-a347-4c1e-a133-0416a0d6efa6)

- ****[Differential abundance testing with Milo](https://marionilab.github.io/miloR/articles/milo_demo.html)****
- find what is it about profileMatrix
- Milo Object
- peakAnnoEnrichment visualize
- hypergeometric enrichment test r
- [HuBMAP](https://portal.hubmapconsortium.org)

### ****Single-cell analyses define a continuum of cell state and composition changes in the malignant transformation of polyps to colorectal cancer****

**Samples**

- CRC
- Normal
- Polyp
- Unaffected

### HCT116

- [encode data](https://www.encodeproject.org/search/?type=Experiment&control_type!=*&status=released&perturbed=false)
- **a human colorectal carcinoma cell line initiated from an adult male**. The cells are adherent with an epithelial morphology. Following implantation into immunocompromised mice, the cells form primary tumors and distant metastases.

![Screen_Shot_2023-08-02_at_11 47 01_AM](https://github.com/yerimeeei/Bioinformatics/assets/134043926/81bbc1ef-9a52-4320-94ff-1c0b5cc5d1f4)

- ATAC-seq signal is lower for HCT116 in L1TD1
- From ATAC-seq visualization: normal에서 많이 open chromatin in L1TD1 → check correlation with potential TF (SMAD1-9)

### Correlation with potential TF (SMAD1-9) with L1TD1

- [GEO data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76987)

```r
#correlation with potential TF (SMAD1-9) in Normal, Adenomatous polyps, CRC
ColonCancerProcessed <- read_excel("/storage2/jwlab/yerim/colon/GSE76987_ColonCancerProcessed.xlsx")
RightColonProcessed <- read_excel("/storage2/jwlab/yerim/colon/GSE76987_RightColonProcessed.xlsx")
LeftColonProcessed <- read_excel("/storage2/jwlab/yerim/colon/GSE76987_LeftColonProcessed.xlsx")

ColonCancer_FPKM <- ColonCancerProcessed[-c(3:6)]
RightColon_FPKM <- RightColonProcessed[-c(3:59)]
LeftColon_FPKM <- LeftColonProcessed[-c(3:27)]

ColonCancer_FPKM <- melt(ColonCancer_FPKM)
RightColon_FPKM <- melt(RightColon_FPKM)
LeftColon_FPKM <- melt(LeftColon_FPKM)

merged_FPKM <- rbind(ColonCancer_FPKM, RightColon_FPKM)
colnames(LeftColon_FPKM)[1] <- "Ensembl_ID"
merged_FPKM <- rbind(merged_FPKM, LeftColon_FPKM)
merged_FPKM$Tissue <- sapply(strsplit(as.character(merged_FPKM$variable), "-"), "[", 1)
merged_FPKM$Tissue <- substr(merged_FPKM$Tissue, 6,10)
merged_FPKM$sample <- NA
merged_FPKM$sample[which(str_detect(merged_FPKM$Tissue, "^CA"))] <- "Colon adenocarcinoma"
merged_FPKM$sample[which(str_detect(merged_FPKM$Tissue, "^CL"))] <- "Control colon"
merged_FPKM$sample[which(str_detect(merged_FPKM$Tissue, "^HP"))] <- "Hyperplastic polyp"
merged_FPKM$sample[which(str_detect(merged_FPKM$Tissue, "^UL"))] <- "Uninvolved colon"
merged_FPKM$sample[which(str_detect(merged_FPKM$Tissue, "^AP"))] <- "Adenomatous polyp"
merged_FPKM$sample[which(str_detect(merged_FPKM$Tissue, "^CR"))] <- "Control colon"
merged_FPKM$sample[which(str_detect(merged_FPKM$Tissue, "^SSA/P"))] <- "Sessile serrated adenoma/polyp"
merged_FPKM$sample[which(str_detect(merged_FPKM$Tissue, "^UR"))] <- "Uninvolved colon"

#L1TD1
grep("ENSG00000240563", merged_FPKM$Ensembl_ID, value = TRUE)
#SMAD1
grep("ENSG00000170365", merged_FPKM$Ensembl_ID, value = TRUE)
#SMAD2
grep("ENSG00000175387", merged_FPKM$Ensembl_ID, value = TRUE)
#SMAD3
grep("ENSG00000166949", merged_FPKM$Ensembl_ID, value = TRUE)
#SMAD4
grep("ENSG00000141646", merged_FPKM$Ensembl_ID, value = TRUE)
#SMAD5
grep("ENSG00000113658", merged_FPKM$Ensembl_ID, value = TRUE)
#SMAD6
grep("ENSG00000137834", merged_FPKM$Ensembl_ID, value = TRUE)
#SMAD7
grep("ENSG00000101665", merged_FPKM$Ensembl_ID, value = TRUE)
#SMAD8 no SMAD8 ensembl ID
#SMAD9
grep("ENSG00000120693", merged_FPKM$Ensembl_ID, value = TRUE)

correlation_table <- merged_FPKM[merged_FPKM$Ensembl_ID =="ENSG00000240563" 
                                 | merged_FPKM$Ensembl_ID =="ENSG00000170365"
                                 | merged_FPKM$Ensembl_ID =="ENSG00000175387"
                                 | merged_FPKM$Ensembl_ID =="ENSG00000166949"
                                 | merged_FPKM$Ensembl_ID =="ENSG00000141646"
                                 | merged_FPKM$Ensembl_ID =="ENSG00000113658"
                                 | merged_FPKM$Ensembl_ID =="ENSG00000137834"
                                 | merged_FPKM$Ensembl_ID =="ENSG00000101665"
                                 | merged_FPKM$Ensembl_ID =="ENSG00000120693",]
L1TD1_SMAD1_9 <- dcast(correlation_table, variable + Tissue + sample ~ Ensembl_ID)
mapIDs <- data.frame(mapIds(org.Hs.eg.db, keys = colnames(L1TD1_SMAD1_9), keytype="ENSEMBL", column = "SYMBOL"))
colnames(L1TD1_SMAD1_9) <- mapIDs$mapIds.org.Hs.eg.db..keys...colnames.L1TD1_SMAD1_9...keytype....ENSEMBL...
colnames(L1TD1_SMAD1_9)[1] <- "variable"
colnames(L1TD1_SMAD1_9)[2] <- "Tissue"
colnames(L1TD1_SMAD1_9)[3] <- "sample"

#L1TD1 vs SMAD1
L1TD1_SMAD1_9_df <- dcast(correlation_table, Ensembl_ID ~variable, value.var = "value")
L1TD1_SMAD <- L1TD1_SMAD1_9_df[-1]
rownames(L1TD1_SMAD) <- L1TD1_SMAD1_9_df$Ensembl_ID
grep("ENSG00000240563", rownames(L1TD1_SMAD))
L1TD1_SMAD <- L1TD1_SMAD[c(9,1:8),]
mapIDs <- data.frame(mapIds(org.Hs.eg.db, keys = rownames(L1TD1_SMAD), keytype="ENSEMBL", column = "SYMBOL"))
rownames(L1TD1_SMAD) <- mapIDs$mapIds.org.Hs.eg.db..keys...rownames.L1TD1_SMAD...keytype....ENSEMBL...

##Control colon
normal_L1TD1_SMAD <- L1TD1_SMAD[c(15:24,62:71)]
normal_corr <- apply(normal_L1TD1_SMAD[2:nrow(normal_L1TD1_SMAD),], 1, function(i) cor(x = as.numeric(normal_L1TD1_SMAD[1,]), y = as.numeric(i), method = "spearman"))
normal_corr <- data.frame(normal_corr)

##Adenomatous Polyp
adenoma_L1TD1_SMAD <- L1TD1_SMAD[c(5:14)]
adenoma_corr <- apply(adenoma_L1TD1_SMAD[2:nrow(adenoma_L1TD1_SMAD),], 1, function(i) cor(x = as.numeric(adenoma_L1TD1_SMAD[1,]), y = as.numeric(i), method = "spearman"))
adenoma_corr <- data.frame(adenoma_corr)

##Colon adenocarcinoma
CRC_L1TD1_SMAD <- L1TD1_SMAD[c(1:4)]
CRC_corr <- apply(CRC_L1TD1_SMAD[2:nrow(CRC_L1TD1_SMAD),], 1, function(i) cor(x = as.numeric(CRC_L1TD1_SMAD[1,]), y = as.numeric(i), method = "spearman"))
CRC_corr <- data.frame(CRC_corr)

cor_L1TD1_SMAD <- cbind(normal_corr, adenoma_corr, CRC_corr)
write.csv(cor_L1TD1_SMAD, "/storage2/jwlab/yerim/cor_L1TD1_SMAD.csv", row.names=TRUE)
```

### Correlation of L1TD1 with All genes only in Adenomatous Polyps

```r
#correlation with all genes only in adenomatous polyps
adenomatous_polyps_genes <- genes[c(1:10)]
##log form correlation test
adenomatous_polyps_genes <- log2(adenomatous_polyps_genes)

adenomatous_polyps_corr <- apply(adenomatous_polyps_genes[2:nrow(adenomatous_polyps_genes),], 1, function(i) cor(x = as.numeric(adenomatous_polyps_genes[1,]), y = as.numeric(i), method = "spearman"))
adenomatous_polyps_corr <- data.frame(adenomatous_polyps_corr)

mapIDs <- data.frame(mapIds(org.Hs.eg.db, keys = rownames(adenomatous_polyps_corr), keytype="ENSEMBL", column = "SYMBOL"))
adenomatous_polyps_corr$gene_symbol <- mapIDs$mapIds.org.Hs.eg.db..keys...rownames.adenomatous_polyps_corr...
adenomatous_polyps_corr$rank <- rank(adenomatous_polyps_corr$adenomatous_polyps_corr)
adenomatous_polyps_corr <- adenomatous_polyps_corr[order(adenomatous_polyps_corr$rank, decreasing = TRUE),]
adenomatous_polyps_corr <- adenomatous_polyps_corr[-3]
adenomatous_polyps_corr <- adenomatous_polyps_corr %>% drop_na()
hist(adenomatous_polyps_corr$adenomatous_polyps_corr, breaks=100)

library(readxl)
MAGIC <- read_excel("/storage2/jwlab/yerim/MAGIC_l1td1_TF_enrichment_analysis.xlsx")
colnames(adenomatous_polyps_corr)[2] <- "TFs"
check_TF <- merge(MAGIC, adenomatous_polyps_corr, by = "TFs")
check_TF$rank <- rank(check_TF$`MaxValue'L1TD1'`)
check_TF <- check_TF[order(check_TF$rank, decreasing = TRUE),]
check_TF <- check_TF[-4]

check_TF <- check_TF[which(check_TF$adenomatous_polyps_corr > 0.5),]
sort_TSS_L1TD1 <- read.table("/storage2/jwlab/yerim/possum/sort_TSS_L1TD1.txt", quote="\"", comment.char="")
sort_TSS_L1TD1 <- sort_TSS_L1TD1[order(sort_TSS_L1TD1$V2, decreasing = TRUE),]
colnames(sort_TSS_L1TD1)[1] <- "TFs"
## possum에서 나온 결과랑 겹치는 친구들: EGR1, NFIA, SP1, ERF

library(tidyr)
colnames(adenomatous_polyps_corr)[2] <- "gene symbol"
write.csv(adenomatous_polyps_corr, "/storage2/jwlab/yerim/adenomatous_polyps_corr.csv", row.names=TRUE)
```

**MAGIC Tool**

- **[MAGIC TRANSCRIPTION FACTOR ENCODE ANALYSIS](https://www.mercuriolab.umassmed.edu/magic-transcription-factor-encode-analysis)**
- [paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007800)

### Permutation test

```r
##permutation test
temp <- adenomatous_polyps_genes[1,]
temp_col <- colnames(temp)
colnames(temp)[2:ncol(temp)] <- temp_col[1:(length(temp_col)-1)]
colnames(temp)[1] <- temp_col[length(temp_col)]
temp <- temp[,c(2:ncol(temp), 1)]
adenomatous_polyps_genes_permutation <- rbind(temp, adenomatous_polyps_genes[2:nrow(adenomatous_polyps_genes),])
adenomatous_polyps_genes_permutation <- log2(adenomatous_polyps_genes_permutation)

adenomatous_polyps_corr_permutation <- apply(adenomatous_polyps_genes_permutation[2:nrow(adenomatous_polyps_genes_permutation),], 1, function(i) cor(x = as.numeric(adenomatous_polyps_genes_permutation[1,]), y = as.numeric(i), method = "spearman"))
adenomatous_polyps_corr_permutation <- data.frame(adenomatous_polyps_corr_permutation)

mapIDs <- data.frame(mapIds(org.Hs.eg.db, keys = rownames(adenomatous_polyps_corr_permutation), keytype="ENSEMBL", column = "SYMBOL"))
adenomatous_polyps_corr_permutation$gene_symbol <- mapIDs$mapIds.org.Hs.eg.db..keys...rownames.adenomatous_polyps_corr_permutation...
adenomatous_polyps_corr_permutation$rank <- rank(adenomatous_polyps_corr_permutation$adenomatous_polyps_corr_permutation)
adenomatous_polyps_corr_permutation <- adenomatous_polyps_corr_permutation[order(adenomatous_polyps_corr_permutation$rank, decreasing = TRUE),]
adenomatous_polyps_corr_permutation <- adenomatous_polyps_corr_permutation[-3]
adenomatous_polyps_corr_permutation <- adenomatous_polyps_corr_permutation %>% drop_na()
hist(adenomatous_polyps_corr_permutation$adenomatous_polyps_corr_permutation, breaks=100)
```

- blue = correlation histogram plot (abs)
- lightpink = permutation histogram plot (abs)

![plot_zoom_png](https://github.com/yerimeeei/Bioinformatics/assets/134043926/da2bc373-073f-499e-9745-8504844515fc)

- **conclusion: cannot trust the correlation result!!**

## GSEA only for Adenomatous Polyps

<img width="1350" alt="GSEA_pathway_only_AP" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/3e0eab60-d1b3-411c-9e2a-6e5ae3fd90c1">

- **GSEA only for Adenomatous Polyps Interferon**

<img width="1455" alt="GSEA_pathway_only_AP_interferon" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/6d6201b0-9a42-48c7-9aeb-b9e6fe849bcf">
