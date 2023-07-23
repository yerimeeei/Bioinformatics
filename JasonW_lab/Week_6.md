# Week 6 L1TD1 Gene Expression Single Cell Analysis

- [paper](https://www.nature.com/articles/s41588-022-01088-x)
- [github code](https://github.com/winstonbecker/scCRC_continuum/blob/main/README.md) form paper
- Seurat [scRNA-seq integration](https://satijalab.org/seurat/articles/integration_introduction.html) (not detailed)
- Seurat [quality control](https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html)
- Seurat deatailed [introduction to scRNA-seq](https://hbctraining.github.io/scRNA-seq/lessons/09_merged_SC_marker_identification.html)
- Seurat [clustering tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#run-non-linear-dimensional-reduction-umaptsne)

## Motif Analysis

```r
"L1TD1"    "FGFRL1"    "ART3"    "ENC1"      "MGC32805"  "SLC12A2"   "MSX2"     
"MET"      "NEBL"      "APIP"    "TCN1"      "SLC22A11"  "SPTBN2"    "SMAD9"    
"NETO2"    "RNF43.1"   "AXIN2"   "LINC00668" "SCML2P1"   "ZNRF3"
```

### Conclusion

**Differential expressed gene exclusively in Adenoma**

- SMAD9 belongs to SMAD family of TFs.
- Binding of other SMAD factors in pluripotent stem cell: SMAD2/3 might bind to L1TD1 promoter (in the motif list, SMAD2 exist).

**L1TD1 motif analysis**

- TF with highest binding affinity: ETV5::FOXI1 (family: ets related::FOX)

**Relationship between SMAD and FOX**

- SMAD2/3 in FOXH1 (TGF-b signaling) [paper](https://pubmed.ncbi.nlm.nih.gov/31582430/)

## Single Cell Analysis

- hypothesis

**Find differentially expressed markers**

```r
# find markers between specific clusters
findMarkers()
# find markers in a cluster compared to all other clusters
# useful when there is only 1 condition
findAllMarkers()
# find markers conserved across conditions
# useful when compare multiple conditions
findConservedMarkers()

condition <- DimPlot(seurat_phase, reduction = "umap", group.by = "sample")
clusters <- DimPlot(seurat_phase, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)
condition|clusters
```

![plot_zoom_png](https://github.com/yerimeeei/Bioinformatics/assets/134043926/d1ca0042-cc80-450b-8998-9800128b5c6d)

**To match with celltypes table (epithelial, immune, stromal)**

```r
seurat_phase$samplename <- substr(seurat_phase$seq_folder, 12,21)
metadata_phase <- seurat_phase@meta.data
seurat_phase$cells <- sub("^[^_]*_","", rownames(metadata_phase))
seurat_phase$cells <- substr(seurat_phase$cells, 1,16)
seurat_phase$name <- paste(seurat_phase$samplename, seurat_phase$cells, sep = "_")
View(seurat_phase@meta.data)
```

- integrate data
- assign number

```r
unique(merged_celltypes$CellType)
 [1] "Stem"                           "Immature Goblet"                "CyclingTA"                     
 [4] "TA2"                            "TA1"                            "Goblet"                        
 [7] "Tuft"                           "Best4+ Enterocytes"             "Enterocytes"                   
[10] "Enterocyte Progenitors"         "Immature Enterocytes"           "Enteroendocrine"               
[13] "Memory B"                       "CD4+"                           "GC"                            
[16] "Naive B"                        "CD8+"                           "Naive T"                       
[19] "Macrophages"                    "Mast"                           "Plasma"                        
[22] "Tregs"                          "NK"                             "DC"                            
[25] "ILCs"                           "Villus Fibroblasts WNT5B+"      "Pericytes"                     
[28] "Endothelial"                    "Crypt Fibroblasts 1"            "Crypt Fibroblasts 3"           
[31] "Myofibroblasts/Smooth Muscle 3" "Crypt Fibroblasts 4"            "Lymphatic endothelial cells"   
[34] "Adipocytes"                     "Glia"                           "Myofibroblasts/Smooth Muscle 2"
[37] "Myofibroblasts/Smooth Muscle 1" "Crypt Fibroblasts 2"            "Neurons"                       
[40] "Unknown"                        "Cancer Associated Fibroblasts"

grep("0", Idents(seurat_phase), value = TRUE)
0 == CyclingTA or TA2 or Stem
1 == TA2 or Stem
2 == Immature Goblet
3 == CyclingTA
4 == Enterocytes (or Immature Enterocytes or Enterocyte Progenitors)
5 == TA1 or TA2
6 == CyclingTA or Stem or TA2
7 == CD8+ or CD4+ NK or Naive T or Tregs
8 == Crypt Fibroblasts
9 == Memory B or Naive B
10 == TA2 or CyclingTA or Immature Goblet
11 == Myofibroblasts/Smooth Muscle 1 or 2
12 == TA2 or CyclingTA or Stem
13 == Enterocytes
14 == CyclingTA
15 == Macrophages
16 == TA2 or CyclingTA or Stem
17 == Stem or TA2
18 == Best4+ Enterocytes
19 == Unknown
20 == Endothelial
21 == Tuft
22 == Plasma
23 == GC
24 == GC (or Goblet)
25 == Crypt Fibroblasts 2
26 == Enteroendocrine
27 == Glia
28 == CRC2 CyclingTA
29 == CD8+
cellname <- RenameIdents(object = seurat_phase, "0" = "CyclingTA or TA2 or Stem",
                         "1" = "TA2 or Stem",
                         "2" = "Immature Goblet",
                         "3" = "CyclingTA",
                         "4" = "Enterocytes or Immature Enterocytes",
                         "5" = "TA1 or TA2",
                         "6" = "CyclingTA or Stem or TA2",
                         "7" = "CD8+ or CD4+ NK or Naive T or Tregs",
                         "8" = "Crypt Fibroblasts",
                         "9" = "Memory B or Naive B",
                         "10" = "TA2 or CyclingTA or Immature Goblet",
                         "11" = "Myofibroblasts/Smooth Muscle",
                         "12" = "TA2 or CyclingTA or Stem",
                         "13" = "Enterocytes",
                         "14" = "CyclingTA",
                         "15" = "Macrophages",
                         "16" = "TA2 or CyclingTA or Stem",
                         "17" = "Stem or TA2",
                         "18" = "Best4+ Enterocytes",
                         "19" = "Unknown",
                         "20" = "Endothelial",
                         "21" = "Tuft",
                         "22" = "Plasma",
                         "23" = "GC",
                         "24" = "GC or Goblet",
                         "25" = "Crypt Fibroblasts 2",
                         "26" = "Enteroendocrine",
                         "27" = "Glia",
                         "28" = "CRC2 CyclingTA",
                         "29" = "CD8+")
DimPlot(object = cellname, reduction = "umap", 
        label = TRUE,label.size = 3,repel = TRUE)
```

![celltype_plot](https://github.com/yerimeeei/Bioinformatics/assets/134043926/f35213e1-2bea-4fa2-b98b-f45631b1445f)

![L1TD1_in_each_sample](https://github.com/yerimeeei/Bioinformatics/assets/134043926/c8ca6589-c1a2-481c-9bc0-aac098781823)

- Polyp: A polyp is a growth inside of your body. **Most aren't cancerous (benign), but a polyp contains abnormal cells or cells that may become abnormal (malignant)**. A polyp is usually a flat bump or shaped like a mushroom. Cancerous polyps can develop in many places in your body, such as your colon or uterus.
- L1TD1, SMAD2/3/9 expression and by disease state

![L1TD1_SMAD239](https://github.com/yerimeeei/Bioinformatics/assets/134043926/538a9302-21b4-471a-84d4-d7e09ed392bb)

![L1TD1_SMAD239 1](https://github.com/yerimeeei/Bioinformatics/assets/134043926/a25c5ae2-fb97-4f9e-a12c-96f9421d4e9d)

### TODO

1. check [SMAD1-9](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10038740/) in bulk RNAseq (see expression if pattern similar with L1TD1/highly expressed in adenoma tissue)
2. SMAD1-9 sample type별로 expression pattern 확인 → L1TD1과 비교하기

![0e657d39-ae9d-4fe7-9091-1833f909e28b](https://github.com/yerimeeei/Bioinformatics/assets/134043926/1406bc3a-1618-4f52-958d-010abd7c9e3b)

![aee32074-0d94-4e69-863f-84dbec1ed0ce](https://github.com/yerimeeei/Bioinformatics/assets/134043926/c8928a13-2915-4784-8fb5-f221902b40bf)

![3c2f36ef-004b-40da-99b0-d2514268a70c](https://github.com/yerimeeei/Bioinformatics/assets/134043926/6bba7b8e-1f39-474e-b64c-031845dc5a75)

TODO

1. SMAD1,4,7,8 (5,6 expression [known](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10038740/) to be low) [SMAD](https://www.notion.so/SMAD-0578f1f8ce814f72b12f823e42d8c54b?pvs=21) 

In bulk RNA-seq data, SMAD5,6,7,8,9 are highly expressed in adenoma but not in colorectal cancer. For the UMAP of single cell analysis, the expression pattern of SMAD9 was most similar to L1TD1, but not others. I also looked up if SMAD family are differentially expressed genes, only SMAD9 was exclusively highly expressed in adenoma but not others.

1. find adenoma랑 normal있는 cohort (bulk RNA-seq)

[**GEO data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207398) of Analysis of the effects of carcinogen treatment on gene expression in colon organoids**

### [GEO data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76987) of RNA sequencing analysis of gene expression in serrated colon polyps, uninvolved colon and control colon

- [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4941619/)

```bash
vi getGEOSuppFiles.R
	#if (!require("BiocManager", quietly = TRUE))
	#    install.packages("BiocManager")

	#BiocManager::install("GEOquery")
	library("GEOquery")
	getGEOSuppFiles("GSE76987", makeDirectory=FALSE, 
									baseDir = "/storage2/jwlab/yerim/colon", fetch_files=TRUE)
Rscript getGEOSuppFiles.R
```

L1TD1 (Ensembl ID: ENSG00000240563)

```r
install.packages("readxl")
library("readxl")
# xlsx files
ColonCancerProcessed <- read_excel("/storage2/jwlab/yerim/colon/GSE76987_ColonCancerProcessed.xlsx")
RightColonProcessed <- read_excel("/storage2/jwlab/yerim/colon/GSE76987_RightColonProcessed.xlsx")
LeftColonProcessed <- read_excel("/storage2/jwlab/yerim/colon/GSE76987_LeftColonProcessed.xlsx")

grep("ENSG00000240563", ColonCancerProcessed$Ensembl_ID, value = TRUE)
grep("ENSG00000240563", LeftColonProcessed$Ensemble_ID, value = TRUE)
grep("ENSG00000240563", RightColonProcessed$Ensembl_ID, value = TRUE)

Colon_L1TD1 <- ColonCancerProcessed[ColonCancerProcessed$Ensembl_ID == "ENSG00000240563",]
Left_L1TD1 <- LeftColonProcessed[LeftColonProcessed$Ensemble_ID == "ENSG00000240563",]
Right_L1TD1 <- RightColonProcessed[RightColonProcessed$Ensembl_ID == "ENSG00000240563",]

Colon_L1TD1 <- melt(Colon_L1TD1)
Colon_L1TD1 <- Colon_L1TD1[-c(1,2,3,4),]
Colon_L1TD1$Tissue <- substr(Colon_L1TD1$variable, 6,7)

Left_L1TD1 <- melt(Left_L1TD1)
Left_L1TD1 <- Left_L1TD1[-c(1:25),]
Left_L1TD1$Tissue <- substr(Left_L1TD1$variable, 6,7)

Right_L1TD1 <- melt(Right_L1TD1)
Right_L1TD1 <- Right_L1TD1[-c(1:57),]
Right_L1TD1$Tissue <- sapply(strsplit(as.character(Right_L1TD1$variable), "-"), "[", 1)
Right_L1TD1$Tissue <- substr(Right_L1TD1$Tissue, 6,10)

colnames(Left_L1TD1)[1] <- "Ensembl_ID"
merged_L1TD1 <- rbind(Colon_L1TD1, Left_L1TD1, Right_L1TD1)

library(ggplot2)
vplot_L1TD1 <- ggplot(merged_L1TD1, aes(x=Tissue, y=value, group=Tissue, fill=Tissue)) 
vplot_L1TD1 <- vplot_L1TD1 + geom_violin(trim=FALSE, fill="gray") 
vplot_L1TD1 <- vplot_L1TD1 + labs(title="L1TD1 expression in serrated colon polyps, uninvolved colon and control colon", x="Sample Identities", y = "log2(FPKM+1)") 
vplot_L1TD1 <- vplot_L1TD1 + geom_boxplot(width=0.1)+ theme_classic() + scale_x_discrete(limits=c("CA","CL","HP","UL","AP","CR","SSA/P","UR"))
vplot_L1TD1
```

![vplot](https://github.com/yerimeeei/Bioinformatics/assets/134043926/47a782cd-0db8-40bb-9ab4-8f2396b1894a)

- 86 colon RNA sequencing datasets
    - 21 sessile serrated adenomas/polyps
    - 10 hyperplastic polyps
    - 10 adenomatous polyps
    - 21 uninvolved colon
    - 20 control colon
    - 4 colon cancer

```r
merged_L1TD1$sample <- NA
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^CA"))] <- "Colon adenocarcinoma"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^CL"))] <- "Control left colon"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^HP"))] <- "Hyperplastic polyp"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^UL"))] <- "Uninvolved left colon"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^AP"))] <- "Adenomatous polyp"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^CR"))] <- "Control right colon"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^SSA/P"))] <- "Sessile serrated adenoma/polyp"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^UR"))] <- "Uninvolved right colon"

merged_L1TD1$sample <- NA
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^CA"))] <- "Colon adenocarcinoma"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^CL"))] <- "Control colon"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^HP"))] <- "Hyperplastic polyp"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^UL"))] <- "Uninvolved colon"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^AP"))] <- "Adenomatous polyp"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^CR"))] <- "Control colon"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^SSA/P"))] <- "Sessile serrated adenoma/polyp"
merged_L1TD1$sample[which(str_detect(merged_L1TD1$Tissue, "^UR"))] <- "Uninvolved colon"

library(ggplot2)
vplot_L1TD1 <- ggplot(merged_L1TD1, aes(x=sample, y=value, group=sample, fill=sample)) 
vplot_L1TD1 <- vplot_L1TD1 + geom_violin(trim=FALSE, fill="gray") 
vplot_L1TD1 <- vplot_L1TD1 + labs(title="L1TD1 expression in serrated colon polyps, uninvolved colon and control colon", x="Sample Identities", y = "log2(FPKM+1)") 
vplot_L1TD1 <- vplot_L1TD1 + geom_boxplot(width=0.1)+ theme_classic() + scale_x_discrete(limits=c("Control colon","Uninvolved colon","Adenomatous polyp","Hyperplastic polyp","Sessile serrated adenoma/polyp","Colon adenocarcinoma"))
vplot_L1TD1
```

![plot_zoom_png 1](https://github.com/yerimeeei/Bioinformatics/assets/134043926/bbf45fdb-e572-4de6-a338-2c5c748c41eb)

**Serrated colon polyps**

- classified into three groups: Hyperplastic polyps (HPs), sessile serrated adenoma/polyps (SSA/Ps), and traditional serrated adenomas (TSAs)
- SSA/Ps often have basilar crypt dilation, which may present as an L-shaped or inverted T-shaped morphology. HPs lack these specific features
- However, differentiating SSA/Ps from HPs by colonoscopy or histopathology remains difficult due to overlapping morphological and pathological features
- HPs
    - Some types of polyps (called adenomas) have the potential to become cancerous, while **others (hyperplastic or inflammatory polyps) have virtually no chance of becoming cancerous**.

```r
cor.test(merged_L1TD1$value, merged_SMAD9$value, method = "spearm", alternative = "g")
HP_L1TD1 <- merged_L1TD1[15:24,]
AP_L1TD1 <- merged_L1TD1[30:39,]
SSA_P_L1TD1 <- merged_L1TD1[50:70,]
Polyp_L1TD1 <- rbind(HP_L1TD1,AP_L1TD1,SSA_P_L1TD1)

HP_SMAD9 <- merged_SMAD9[15:24,]
AP_SMAD9 <- merged_SMAD9[30:39,]
SSA_P_SMAD9 <- merged_SMAD9[50:70,]
Polyp_SMAD9 <- rbind(HP_SMAD9,AP_SMAD9,SSA_P_SMAD9)

# combined polp samples correlation between L1TD1 and SMAD9
cor.test(Polyp_L1TD1$value, Polyp_SMAD9$value, method = "spearm", alternative = "g")

	Spearman's rank correlation rho

data:  Polyp_L1TD1$value and Polyp_SMAD9$value
S = 7282, p-value = 0.009628
alternative hypothesis: true rho is greater than 0
sample estimates:
      rho 
0.3656794

# Hyperplastic polyp
cor.test(HP_L1TD1$value, HP_SMAD9$value, method = "spearm", alternative = "g")

	Spearman's rank correlation rho

data:  HP_L1TD1$value and HP_SMAD9$value
S = 106, p-value = 0.1564
alternative hypothesis: true rho is greater than 0
sample estimates:
      rho 
0.3575758

# Adenomatous polyp
cor.test(AP_L1TD1$value, AP_SMAD9$value, method = "spearm", alternative = "g")

	Spearman's rank correlation rho

data:  AP_L1TD1$value and AP_SMAD9$value
S = 140, p-value = 0.3409
alternative hypothesis: true rho is greater than 0
sample estimates:
      rho 
0.1515152

# Sessile serrated adenoma/polyp
cor.test(SSA_P_L1TD1$value, SSA_P_SMAD9$value, method = "spearm", alternative = "g")

	Spearman's rank correlation rho

data:  SSA_P_L1TD1$value and SSA_P_SMAD9$value
S = 828, p-value = 0.01813
alternative hypothesis: true rho is greater than 0
sample estimates:
      rho 
0.4623377
```

| Entire polyp (n=41) | Hyperplastic polyp (n=10) | Adenomatous polyp (n=10) | Sessile serrated adenoma/polyp (n=21) |
| --- | --- | --- | --- |
| sample estimates:
      rho 
0.3656794 | sample estimates:
      rho 
0.3575758 | sample estimates:
      rho 
0.1515152 | sample estimates:
      rho 
0.4623377 |

```r
AP <- merge(x = AP_L1TD1, y = AP_SMAD9, by = "variable")
library(ggpubr)
ggscatter(AP, x = "value.x", y = "value.y", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "L1TD1", ylab = "SMAD9")
```

![scatter](https://github.com/yerimeeei/Bioinformatics/assets/134043926/fef8ebb6-b032-48e7-b032-c35f8d40cfea)

when n = 40 (Polyp)

```r
cor.test(Polyp$value.x, Polyp$value.y, method = "spearm", alternative = "g")

	Spearman's rank correlation rho

data:  Polyp$value.x and Polyp$value.y
S = 7282, p-value = 0.02341
alternative hypothesis: true rho is greater than 0
sample estimates:
      rho 
0.3168856
```

![scatter_polyp](https://github.com/yerimeeei/Bioinformatics/assets/134043926/adf351b7-2a86-41bd-81f2-aac0148394d8)

- [ ]  co expression in Immature Goblet in scRNA-seq (+ can also look up for ATAC-seq)
- [ ]  immunogenicity (CD4 & CD8 T cell)
- [ ]  Normal (open, demethylated, but not expressed), Adenoma (open, expressed), CRC (closed, methylated)
- [ ]  functional pathway (if happening in specific cell type like Immature Goblet?)
- [ ]  CPG Island?
