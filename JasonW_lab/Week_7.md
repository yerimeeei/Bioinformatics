# Week 7 L1TD1 Gene Expression and Methylation

<img width="965" alt="genome_browser" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/9fb2837b-8b1a-4669-85de-01e51418ddc4">

**CpG Island**

- CpG islands (CGIs) are **regions of the genome that contain a large number of CpG dinucleotide repeats**. In mammalian genomes, CpG islands usually extend for 300–3000 base pairs. They are located within and close to sites of about 40% of mammalian gene promoters.
- CpG islands are DNA methylations regions in promoters known to **regulate gene expression through transcriptional silencing of the corresponding gene**. DNA methylation at CpG islands is crucial for gene expression and tissue-specific processes.
- [x]  L1TD1 vs gene --> spearman correlation
    - [x]  fgsea package gene enrichment profile
- [x]  [DNA methylation dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48684) --> beta value --> quasi_random plot among three groups (normal, adenoma, cancer)
- [ ]  ATAC seq 논문 읽어보고 ATAC-seq [Arch object](https://www.archrproject.com) [찾아보기](https://www.archrproject.com/bookdown/iterative-latent-semantic-indexing-lsi.html)

[fgsea package practice](https://github.com/ctlab/fgsea)

```r
library(data.table)
library(fgsea)
library(ggplot2)
data("examplePathways")
data("exampleRanks")
fgseaRes <- fgsea(pathways = examplePathways, stats = exampleRanks, eps = 0.0, minSize = 15, maxSize = 500)
plotEnrichment(examplePathways[["5990981_DNA_Replication"]], exampleRanks) + labs(title = "DNA Replication")
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, gseaParam = 0.5)
```

**Preparing for Correlation Test table**

```r
# L1TD1 vs gene 
ColonCancerProcessed <- read_excel("/storage2/jwlab/yerim/colon/GSE76987_ColonCancerProcessed.xlsx")
RightColonProcessed <- read_excel("/storage2/jwlab/yerim/colon/GSE76987_RightColonProcessed.xlsx")
LeftColonProcessed <- read_excel("/storage2/jwlab/yerim/colon/GSE76987_LeftColonProcessed.xlsx")

ColonCancer_FPKM <- ColonCancerProcessed[-c(3:6)]
RightColon_FPKM <- RightColonProcessed[-c(3:59)]
LeftColon_FPKM <- LeftColonProcessed[-c(3:27)]

ColonCancer_FPKM <- melt(ColonCancer_FPKM)
RightColon_FPKM <- melt(RightColon_FPKM)
LeftColon_FPKM <- melt(LeftColon_FPKM)

#length(RightColon_FPKM) <- length(ColonCancer_FPKM)
#length(LeftColon_FPKM) <- length(ColonCancer_FPKM)
#FPKM_merged <- cbind(ColonCancer_FPKM, RightColon_FPKM, LeftColon_FPKM)

#L1TD1 ENSG00000240563
#unique(ColonCancer_FPKM$Ensembl_ID)

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

merged_Polyps <- merged_FPKM[!grepl("Control colon|Uninvolved colon|Colon adenocarcinoma", merged_FPKM$sample),]

library(dplyr)
#merged_Polyps <- merged_Polyps[!grepl("ENSG00000240563", merged_Polyps$Ensembl_ID),]
#test <- merged_Polyps %>% filter(grepl(grep("41", table(merged_Polyps$Ensembl_ID)), merged_Polyps$Ensembl_ID))
#temp <- dcast(merged_Polyps, variable ~ Ensembl_ID, value.var = "value", fun.aggregate = list)
test <- dcast(merged_Polyps, Ensembl_ID ~variable, value.var = "value")
genes <- test[-1]
rownames(genes) <- test$Ensembl_ID
genes[is.na(genes)] <- 0
grep("ENSG00000240563", rownames(genes))
genes <- genes[c(19987,1:19986,19988:22960),]
corr_test <- apply(genes[2:nrow(genes),], 1, function(i) cor(x = as.numeric(genes[1,]), y = as.numeric(i), method = "spearman"))
corr_test <- data.frame(corr_test)
```

[An R package for Reactome Pathway Analysis](https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html)

****[DESeq results to pathways in 60 Seconds with the fgsea package](https://stephenturner.github.io/deseq-to-fgsea/)****

```r
library(fgsea)
library(org.Hs.eg.db)
mapIDs <- data.frame(mapIds(org.Hs.eg.db, keys = rownames(corr_test), keytype="ENSEMBL", column = "ENTREZID"))
corr_test$gene_symbol <- mapIDs$mapIds.org.Hs.eg.db..keys...rownames.corr_test...keytype....ENSEMBL...
corr_test <- corr_test[c(2,1)]
corr_test[is.na(corr_test)] <- 0
corr_test$rank <- rank(corr_test$corr_test)
corrtest <- corr_test[order(corr_test$rank),]
corrtest <- corrtest[,-3]

#density plot of correlation value
ggplot(corrtest, aes(x=corr_test)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.1,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")

#cutoff by density plot
corrtest_cutoff <- corrtest[which(abs(corrtest$corr_test) > 0.5),]
#rank <- deframe(corrtest)
rank_cutoff <- deframe(corrtest_cutoff)

library(tibble)
#ranks <- deframe(corr_test)
pathways <- gmtPathways("/storage2/jwlab/yerim/h.all.v2023.1.Hs.symbols.gmt")
curatedpathways <- gmtPathways("/storage2/jwlab/yerim/c2.all.v2023.1.Hs.symbols.gmt.txt")
entrez_pathways <- gmtPathways("/storage2/jwlab/yerim/c2.all.v2023.1.Hs.entrez.gmt.txt")

fgseaRes <- fgsea(pathways=entrez_pathways, stats=rank_cutoff, nperm=1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
# Pathways NES from GSEA
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA") + 
  theme_minimal()

# absolute
corr_test_abs <- corr_test
corr_test_abs$rank <- rank(abs(corr_test_abs$corr_test))
corrtest_abs <- corr_test_abs[order(corr_test_abs$rank),]
corrtest_abs <- corrtest_abs[,-3]
#cutoff by density plot
corrtest_cutoff_abs <- corrtest[which(abs(corrtest_abs$corr_test) > 0.5),]
rank_cutoff_abs <- deframe(corrtest_cutoff_abs)

fgseaResabs <- fgsea(pathways=entrez_pathways, stats=rank_cutoff_abs, nperm=1000)
fgseaResTidyabs <- fgseaResabs %>%
  as_tibble() %>%
  arrange(desc(NES))
fgseaResTidyabs %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
# Pathways NES from GSEA
ggplot(fgseaResTidyabs, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Pathways NES from GSEA") + 
  theme_minimal()
```

<img width="974" alt="Screenshot_2023-07-26_at_6 03 12_PM" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/7f9dc651-4af7-432d-b6f6-dc67a30dcc30">

not absolute

<img width="969" alt="Screenshot_2023-07-26_at_6 04 14_PM" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/36122051-7978-4ed2-ae7d-28edc1807081">

absolute

![GSEA](https://github.com/yerimeeei/Bioinformatics/assets/134043926/d153bf86-22d3-4ac2-a556-3c0e5096eff2)

**Significant pathway**

- Oxidative Phosphorylation
- MYC Targets V1
- Fatty Acid Metabolism
- Bile Acid Metabolism
- Adipogenesis

![GSEA_abs](https://github.com/yerimeeei/Bioinformatics/assets/134043926/84f9cea2-8394-49e6-ae57-db1e229672ea)

**Significant pathway**

- Oxidative Phosphorylation
- Fatty Acid Metabolism
- Adipogenesis
- Xenobiotic Metabolism
- MYC Targets V1

![GSEA_table_plot](https://github.com/yerimeeei/Bioinformatics/assets/134043926/48cef0ae-51f8-4e5d-89a3-98eb5b25428f)

![GSEA_abs_table_plot](https://github.com/yerimeeei/Bioinformatics/assets/134043926/1daaf600-d832-4dc7-bbe7-ce3ed195422f)

[GSEA analysis](https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html) RNA-seq analysis

[Pathway Dataset](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C2)

## DNA Methylation

### [**geo**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48684)

- normal individual 17 + normal colon sample from CRC patient 24 = 41
- adenoma patients 42
- CRC patients 64

**[A statistical model for the analysis of beta values in DNA methylation studies](https://www.notion.so/Week-7-L1TD1-Gene-Expression-and-Methylation-06c02dc525784d0fb84201afd6293541?pvs=21)**

extract signal intensities 2 (only adenoma samples)

```r
signal_intensities_2 <- GSE48684_Matrix_signal_intensities_2[-1]

library(dplyr)
signal_intensities_2 <- signal_intensities_2 %>%
  select(-ends_with("Pval"))
library(stringr)
signal_intensities_2 <- signal_intensities_2 %>%
  rename_with(~str_replace(.,"Signal_A","Unmethylated"))
signal_intensities_2 <- signal_intensities_2 %>%
  rename_with(~str_replace(.,"Signal_B","Methylated"))
signal_intensities_2$X16748.Beta <- signal_intensities_2$X16748.Methylated/(signal_intensities_2$X16748.Methylated+signal_intensities_2$X16748.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:2,74,3:73)]
signal_intensities_2$X16740.Beta <- signal_intensities_2$X16740.Methylated/(signal_intensities_2$X16740.Methylated+signal_intensities_2$X16740.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:6,75,7:74)]
signal_intensities_2$X16743.Beta <- signal_intensities_2$X16743.Methylated/(signal_intensities_2$X16743.Methylated+signal_intensities_2$X16743.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:10,76,11:75)]
#test <- signal_intensities_2 %>%
  #select(ends_with("Methylated"))/(select(ends_with("Methylated"))+select(ends_with("Unmethylated"))+100)
signal_intensities_2$X16725.Beta <- signal_intensities_2$X16725.Methylated/(signal_intensities_2$X16725.Methylated+signal_intensities_2$X16725.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:14,77,15:76)]
signal_intensities_2$X16730.Beta <- signal_intensities_2$X16730.Methylated/(signal_intensities_2$X16730.Methylated+signal_intensities_2$X16730.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:18,78,19:77)]
signal_intensities_2$X16749.Beta <- signal_intensities_2$X16749.Methylated/(signal_intensities_2$X16749.Methylated+signal_intensities_2$X16749.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:22,79,23:78)]
signal_intensities_2$X16733.Beta <- signal_intensities_2$X16733.Methylated/(signal_intensities_2$X16733.Methylated+signal_intensities_2$X16733.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:26,80,27:79)]
#test <- signal_intensities_2 %>%
#  mutate(across(matches("^Methylated$"), ~(.x / get(paste0("^Methylated$", substr(.col,2)))))+matches("^Unmethylated$")+100)
signal_intensities_2$X16736.Beta <- signal_intensities_2$X16736.Methylated/(signal_intensities_2$X16736.Methylated+signal_intensities_2$X16736.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:30,81,31:80)]
signal_intensities_2$X16747.Beta <- signal_intensities_2$X16747.Methylated/(signal_intensities_2$X16747.Methylated+signal_intensities_2$X16747.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:34,82,35:81)]
signal_intensities_2$X16734.Beta <- signal_intensities_2$X16734.Methylated/(signal_intensities_2$X16734.Methylated+signal_intensities_2$X16734.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:38,83,39:82)]
signal_intensities_2$X16732.Beta <- signal_intensities_2$X16732.Methylated/(signal_intensities_2$X16732.Methylated+signal_intensities_2$X16732.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:42,84,43:83)]
signal_intensities_2$X16723.Beta <- signal_intensities_2$X16723.Methylated/(signal_intensities_2$X16723.Methylated+signal_intensities_2$X16723.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:46,85,47:84)]
signal_intensities_2$X16739.Beta <- signal_intensities_2$X16739.Methylated/(signal_intensities_2$X16739.Methylated+signal_intensities_2$X16739.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:50,86,51:85)]
signal_intensities_2$X16742.Beta <- signal_intensities_2$X16742.Methylated/(signal_intensities_2$X16742.Methylated+signal_intensities_2$X16742.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:54,87,55:86)]
signal_intensities_2$X16745.Beta <- signal_intensities_2$X16745.Methylated/(signal_intensities_2$X16745.Methylated+signal_intensities_2$X16745.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:58,88,59:87)]
signal_intensities_2$X16750.Beta <- signal_intensities_2$X16750.Methylated/(signal_intensities_2$X16750.Methylated+signal_intensities_2$X16750.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:62,89,63:88)]
signal_intensities_2$TY1226B.Beta <- signal_intensities_2$TY1226B.Methylated/(signal_intensities_2$TY1226B.Methylated+signal_intensities_2$TY1226B.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:66,90,67:89)]
signal_intensities_2$TY0652B.Beta <- signal_intensities_2$TY0652B.Methylated/(signal_intensities_2$TY0652B.Methylated+signal_intensities_2$TY0652B.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:70,91,71:90)]
signal_intensities_2$TY1374A.Beta <- signal_intensities_2$TY1374A.Methylated/(signal_intensities_2$TY1374A.Methylated+signal_intensities_2$TY1374A.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:74,92,75:91)]
signal_intensities_2$TY0905A.Beta <- signal_intensities_2$TY0905A.Methylated/(signal_intensities_2$TY0905A.Methylated+signal_intensities_2$TY0905A.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:78,93,79:92)]
signal_intensities_2$X1902.Beta <- signal_intensities_2$X1902.Methylated/(signal_intensities_2$X1902.Methylated+signal_intensities_2$X1902.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:82,94,83:93)]
signal_intensities_2$X1186A.Beta <- signal_intensities_2$X1186A.Methylated/(signal_intensities_2$X1186A.Methylated+signal_intensities_2$X1186A.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:86,95,87:94)]
signal_intensities_2$X7497.Beta <- signal_intensities_2$X7497.Methylated/(signal_intensities_2$X7497.Methylated+signal_intensities_2$X7497.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:90,96,91:95)]
signal_intensities_2$X2026.Beta <- signal_intensities_2$X2026.Methylated/(signal_intensities_2$X2026.Methylated+signal_intensities_2$X2026.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:94,97,95:96)]
```

[signal intensities 1 (normal, adenoma, CRC)](https://www.notion.so/signal-intensities-1-normal-adenoma-CRC-642af0b1a0f0475ba99882be971e9e0a?pvs=21)

- cannot figure out which sample is normal, adenoma, or colorectal cancer

```r
# codes I used but not useful anymore
#library(readxl)
#description <- read_excel("/storage/jwlab/yerim/methylation/GPL13534_450K_Manifest_header_Descriptions.xlsx")
GSE48684_Matrix_signal_intensities_1 <- read.delim("/storage2/jwlab/yerim/GSE48684_Matrix_signal_intensities_1.txt")
GSE48684_Matrix_signal_intensities_2 <- read.delim("/storage2/jwlab/yerim/GSE48684_Matrix_signal_intensities_2.txt", comment.char="#")
signal_intensities_1 <- GSE48684_Matrix_signal_intensities_1[-1]
signal_intensities_2 <- GSE48684_Matrix_signal_intensities_2[-1]

#signal_intensities_2 data sorting
library(dplyr)
signal_intensities_2 <- signal_intensities_2 %>%
  select(-ends_with("Pval"))
library(stringr)
signal_intensities_2 <- signal_intensities_2 %>%
  rename_with(~str_replace(.,"Signal_A","Unmethylated"))
signal_intensities_2 <- signal_intensities_2 %>%
  rename_with(~str_replace(.,"Signal_B","Methylated"))
signal_intensities_2$X16748.Beta <- signal_intensities_2$X16748.Methylated/(signal_intensities_2$X16748.Methylated+signal_intensities_2$X16748.Unmethylated+100)
signal_intensities_2 <- signal_intensities_2[c(1:2,74,3:73)]
#test <- signal_intensities_2 %>%
  #select(ends_with("Methylated"))/(select(ends_with("Methylated"))+select(ends_with("Unmethylated"))+100)
#test <- signal_intensities_2 %>%
#  mutate(across(matches("^Methylated$"), ~(.x / get(paste0("^Methylated$", substr(.col,2)))))+matches("^Unmethylated$")+100)

#signal_intensities_1 data sorting
signal_intensities_1 <- GSE48684_Matrix_signal_intensities_1[-1]
signal_intensities_1 <- signal_intensities_1[-2]
signal_intensities_1 <- signal_intensities_1[-2]
library(dplyr)
signal_intensities_1 <- signal_intensities_1 %>%
  select(-ends_with("Pval"))
library(stringr)
signal_intensities_1 <- signal_intensities_1 %>%
  rename_with(~str_replace(.,"Signal_A","Unmethylated"))
signal_intensities_1 <- signal_intensities_1 %>%
  rename_with(~str_replace(.,"Signal_B","Methylated"))

signal_intensities_1$X3325.normal.Beta <- signal_intensities_1$X3325.normal.Methylated/(signal_intensities_1$X3325.normal.Methylated+signal_intensities_1$X3325.normal.Unmethylated+100)
signal_intensities_1 <- signal_intensities_1[c(1,248,2:247)]

#extract only Beta values
signal_intensities1 <- signal_intensities_1 %>%
  select(-ends_with("Unmethylated"))
signal_intensities1 <- signal_intensities1 %>%
  select(-ends_with("Methylated"))

signal_intensities2 <- signal_intensities_2 %>%
  select(-ends_with("Unmethylated"))
signal_intensities2 <- signal_intensities2 %>%
  select(-ends_with("Methylated"))
signal_intensities2 <- signal_intensities2 %>%
  select(-ends_with("AVG_Beta"))
signal_intensities2 <- signal_intensities2 %>%
  rename_with(~str_replace(.,"Beta","adenoma.Beta"))

#merge signal intensities table
signal_intensities <- merge(signal_intensities1, signal_intensities2, by = "TargetID")
signal_intensities <- melt(signal_intensities)
colnames(signal_intensities)[3] <- "Beta_value"
#signal_intensities$sample <- NA
#signal_intensities$sample[which(str_detect(signal_intensities$variable, "^$normal$"))] <- "Normal"
#signal_intensities$sample[which(str_detect(signal_intensities$variable, "^$adenoma$"))] <- "Adenoma"
#signal_intensities$sample[which(str_detect(signal_intensities$variable, "^$cancer$"))] <- "CRC"

signal_intensities <- signal_intensities %>%
  mutate(samples = str_sub(variable, start=1, end = -6))
```

- download matrix series file instead!!!

```r
GSE48684_series_matrix_beta <- read.delim("/storage2/jwlab/yerim/methylation/GSE48684_series_matrix_beta.txt")
View(GSE48684_series_matrix_beta)
beta_matrix <- GSE48684_series_matrix_beta[-485578,]
beta_value <- melt(beta_matrix)

beta_value$sample <- NA
beta_value$sample[which(str_detect(beta_value$variable, "^GSM1183439"))] <- "normal colon"
beta_value$sample[which(str_detect(beta_value$variable, "^GSM1183440"))] <- "normal colon"
beta_value$sample[which(str_detect(beta_value$variable, "^GSM1183444"))] <- "normal colon"
beta_value$sample[which(str_detect(beta_value$variable, "^GSM1183445"))] <- "adenoma"
beta_value$sample[which(str_detect(beta_value$variable, "^GSM1183446"))] <- "adenoma"
beta_value$sample[which(str_detect(beta_value$variable, "^GSM1183447"))] <- "adenoma"
beta_value$sample[which(str_detect(beta_value$variable, "^GSM1183453"))] <- "CRC"
beta_value$sample[which(str_detect(beta_value$variable, "^GSM1183454"))] <- "CRC"
beta_value$sample[which(str_detect(beta_value$variable, "^GSM1183455"))] <- "CRC"

#Cut cpg
HumanMethylation450_probe <- read.delim("/storage2/jwlab/yerim/methylation/HumanMethylation450_probe.txt")
cpgprobe <- HumanMethylation450_probe[c(1,15,16,22,25,30,34,35)]
cpgprobe <- cpgprobe[cpgprobe$Chromosome_36 == "1",]
cpgprobe$Coordinate_36 <- as.numeric(cpgprobe$Coordinate_36)

cut <- cpgprobe[str_detect(cpgprobe$UCSC_RefGene_Name, "^L1TD1"),]
beta_cut1 <- beta_value[beta_value$ID_REF == "cg03731268",]
#beta_cut2 <- beta_value[beta_value$ID_REF == "cg03987985",]
beta_cut3 <- beta_value[beta_value$ID_REF == "cg06190807",]
#beta_cut4 <- beta_value[beta_value$ID_REF == "cg08008065",]
#beta_cut5 <- beta_value[beta_value$ID_REF == "cg11146821",]
beta_cut6 <- beta_value[beta_value$ID_REF == "cg12640000",]
beta_cut7 <- beta_value[beta_value$ID_REF == "cg14254748",]
beta_cut8 <- beta_value[beta_value$ID_REF == "cg14850026",]
beta_cut9 <- beta_value[beta_value$ID_REF == "cg19594218",]
beta_cut10 <- beta_value[beta_value$ID_REF == "cg23049458",]
beta_cut11 <- beta_value[beta_value$ID_REF == "cg27300647",]
#beta_cut12 <- beta_value[beta_value$ID_REF == "cg27424906",]
beta_cut <- rbind(beta_cut1,beta_cut3,beta_cut6,
                  beta_cut7,beta_cut8,beta_cut9,beta_cut10,beta_cut11)

#cpgprobe_cut <- cpgprobe[cpgprobe$RANGE_START > 62194134,]
#cpgprobe_cut <- cpgprobe_cut[cpgprobe_cut$RANGE_START < 62212537,]

#cpgprobe_cut <- cpgprobe[cpgprobe$Coordinate_36 > 62194134,]
#cpgprobe_cut <- cpgprobe_cut[cpgprobe_cut$Coordinate_36 < 62212537,]

#beta_value_cut <- beta_value[beta_value$ID_REF == cpgprobe_cut$ID,]
#beta_value_cut1 <- beta_value[beta_value$ID_REF == "cg00819857",]
#beta_value_cut2 <- beta_value[beta_value$ID_REF == "cg03536320",]
#beta_value_cut3 <- beta_value[beta_value$ID_REF == "cg03605365",]
#beta_value_cut4 <- beta_value[beta_value$ID_REF == "cg04914562",]
#beta_value_cut5 <- beta_value[beta_value$ID_REF == "cg08894594",]
#beta_value_cut6 <- beta_value[beta_value$ID_REF == "cg10139846",]
#beta_value_cut7 <- beta_value[beta_value$ID_REF == "cg10316635",]
#beta_value_cut8 <- beta_value[beta_value$ID_REF == "cg10704177",]
#beta_value_cut9 <- beta_value[beta_value$ID_REF == "cg12128740",]
#beta_value_cut10 <- beta_value[beta_value$ID_REF == "cg13919148",]
#beta_value_cut <- rbind(beta_value_cut1,beta_value_cut2,beta_value_cut3,beta_value_cut4,beta_value_cut5,
#                        beta_value_cut6,beta_value_cut7,beta_value_cut8,beta_value_cut9,beta_value_cut10)
```

### Quasirandom Plot

- [ggbeeswarm package examples](https://cran.r-project.org/web/packages/ggbeeswarm/vignettes/usageExamples.pdf)

```r
library(ggplot2)
library(ggbeeswarm)
methylation_plot <- ggplot(beta_value, aes(x=sample, y=value, group=sample, fill=sample)) 
methylation_plot <- methylation_plot + geom_quasirandom(width=0.1) 
methylation_plot <- methylation_plot + labs(title="DNA Methylation Normal colon vs Adenoma vs CRC", x="Sample Identities", y = "Average Beta Value") 
methylation_plot <- methylation_plot + geom_boxplot(width=0.1)+ theme_classic() + scale_x_discrete(limits=c("normal colon", "adenoma", "CRC"))
methylation_plot

#for each cpg island
methylation_plot <- methylation_plot + facet_wrap(~ID_REF)
methylation_plot
```

![methylation_plot](https://github.com/yerimeeei/Bioinformatics/assets/134043926/78fe85bc-d5ca-4a70-8928-d7b272a9dc92)

![methylation_plot_cpg_island](https://github.com/yerimeeei/Bioinformatics/assets/134043926/7ad761f1-12bb-4a17-8810-1dee4ca16c97)

Website that can visualize genome information

- [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway?hgsid=1665308484_lsRL1MuCON4lIjVUp0ClMiOSOqsa)
- ENCODE (can search [by experiment](https://www.encodeproject.org))

### CPG Island

**[Critical evaluation of the Illumina MethylationEPIC BeadChip microarray for whole-genome DNA methylation profiling](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5055731/)**

- EPIC
- HM450

### Bisulfite sequencing

- WGBS
- RRBS
