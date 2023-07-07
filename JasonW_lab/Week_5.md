# Week 5 L1TD1 Gene Expression Single Cell Analysis

![TF_binding_affinity](https://github.com/yerimeeei/Bioinformatics/assets/134043926/b4812cf6-9601-4162-84ba-99aeea52d367)

## L1TD1 Transcription Factor Motif

| name | class | family |  |
| --- | --- | --- | --- |
| ETV5::FOXI1 | Tryptophan cluster factors::Fork head/winged helix factors | Ets-related::FOX | 10.8 |
| HOXC11 | Homeo domain factors | HOX | 10.1 |
| HMBOX1 | Homeo domain factors | POU domain factors | 9.42 |
| IKZF1 | C2H2 zinc finger factors | Factors with multiple dispersed zinc fingers | 9.32 |
| HOXC9 | Homeo domain factors | HOX | 9.27 |
| POU2F1 | Homeo domain factors | POU domain factors | 8.89 |
| ZNF331 | C2H2 zinc finger factors | More than 3 adjacent zinc fingers | 8.80 |
| ERF::FOXO1 | Tryptophan cluster factors::Fork head/winged helix factors | Ets-related::FOX | 8.65 |
| GCM2 | GCM domain factors  | GCM | 8.48 |
| NR2F1 | Nuclear receptors with C4 zinc fingers  | RXR-related receptors (NR2)  | 8.44 |
| FOXD2 | Fork head/winged helix factors | FOX | 8.32 |
| ZNF214 | C2H2 zinc finger factors | More than 3 adjacent zinc fingers | 8.29 |
| HOXD9 | Homeo domain factors | HOX | 8.26 |

### Which on Differentially Expressed Gene is a Transcript Factor

from the list:

"L1TD1"    "FGFRL1"    "ART3"    "ENC1"      "MGC32805"  "SLC12A2"   "MSX2"     
"MET"      "NEBL"      "APIP"    "TCN1"      "SLC22A11"  "SPTBN2"    "SMAD9"    
"NETO2"    "RNF43.1"   "AXIN2"   "LINC00668" "SCML2P1"   "ZNRF3"

| HGNC Symbol | DBD Family | Is TF? | TF Assessment | Binding Activity | Motif Status | Notes | Comments |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MSX2 | Homeodomain | Yes | Known motif | 1 Monomer or homomultimer | High-throughput in vitro |  | Class: Homeo domain factors Family: NK |
| MET | Unknown | No | Unlikely to be sequence specific TF | 4 Not a DNA binding protein | No motif |  | Its a cell surface receptor and kinase - included only because Vasquerizas 2009 included it with an x. |
| SMAD9 | SMAD | Yes | Known motif | 2 Obligate heteromer | In vivo/Misc source | Only known motifs are from Transfac or HocoMoco - origin is uncertain | Class: SMAD/NF-1 DNA-binding domain factors Family: SMAD factors |

Transcription factor binding motifs (TFBMs) are **genomic sequences that specifically bind to transcription factors**. The consensus sequence of a TFBM is variable, and there are a number of possible bases at certain positions in the motif, whereas other positions have a fixed base.

### Match motif (gene list only highly expressed in adenoma and L1TD1 promoter/gene)

# Single Cell Analysis

- [singel cell analysis CRC](https://www.nature.com/articles/s41588-022-01088-x)
1. Download GEO file
- [european nucleotide archive](https://www.ebi.ac.uk/ena/browser/view/PRJNA830912) (ENA)
- [download GEOSuppFiles](https://rdrr.io/bioc/GEOquery/man/getGEOSuppFiles.html)

```bash
vi yerim.R
	#if (!require("BiocManager", quietly = TRUE))
	#    install.packages("BiocManager")

	#BiocManager::install("GEOquery")
	library("GEOquery")
	getGEOSuppFiles("GSE201348", makeDirectory=FALSE, baseDir = "/storage2/jwlab/yerim/singlecell", fetch_files=TRUE)
Rscript yerim.R
tar -xvf GSE201348_RAW.tar
# gunzip GSM6*.txt.gz

# move each barcode, gene, and matrix file to according subdirectory
mkdir GSE201348_RAW
find /storage2/jwlab/yerim/singlecell/*_barcodes.tsv.gz | head
		#/storage2/jwlab/yerim/singlecell/GSM6061645_A001-C-007_barcodes.tsv.gz
		#/storage2/jwlab/yerim/singlecell/GSM6061646_A001-C-014_barcodes.tsv.gz
		#/storage2/jwlab/yerim/singlecell/GSM6061647_A001-C-023_barcodes.tsv.gz
find /storage2/jwlab/yerim/singlecell/*_barcodes.tsv.gz | awk -F'_barcodes' '{print$1}' | head
		#/storage2/jwlab/yerim/singlecell/GSM6061645_A001-C-007
		#/storage2/jwlab/yerim/singlecell/GSM6061646_A001-C-014
		#/storage2/jwlab/yerim/singlecell/GSM6061647_A001-C-023
find /storage2/jwlab/yerim/singlecell/*_barcodes.tsv.gz | awk -F'_barcodes' '{print$1}' | awk -F'/' '{print$6}' | head
		#GSM6061645_A001-C-007
		#GSM6061646_A001-C-014
		#GSM6061647_A001-C-023
find /storage2/jwlab/yerim/singlecell/*_barcodes.tsv.gz | awk -F'_barcodes' '{print$1}' | awk -F'/' '{print$6}' | awk '{print"mkdir /storage2/jwlab/yerim/singlecell/GSE201348_raw/"$0}' > temp.sh
		#mkdir /storage2/jwlab/yerim/singlecell/GSE201348_RAW/GSM6061645_A001-C-007
		#mkdir /storage2/jwlab/yerim/singlecell/GSE201348_RAW/GSM6061646_A001-C-014
bash temp.sh
# move files into each directory
find /storage2/jwlab/yerim/singlecell/*_barcodes.tsv.gz | awk -F'_barcodes' '{print$1}' | awk -F'/' '{print$6}' | awk '{print"mv /storage2/jwlab/yerim/singlecell/"$0"* /storage2/jwlab/yerim/singlecell/GSE201348_raw/"$0}' > temp.sh
bash temp.sh
```

- Rename files to barcodes.tsv, features.tsv, and matrix.mtx to do Read10X in R

```bash
find /storage2/jwlab/yerim/singlecell/GSE201348_RAW/GSM6061*/*_barcodes.tsv | awk '{print"mv "$0 "/storage2/jwlab/yerim/singlecell/GSE201348_RAW/"}'| awk -F'/' '{print$7}' > list.sh

# list of names of barcode files
find /storage2/jwlab/yerim/singlecell/GSE201348_RAW/GSM6061*/*_barcodes.tsv

# rename to barcodes.tsv, features.tsv, and matrix.mtx
awk '{print"mv /storage2/jwlab/yerim/singlecell/GSE201348_RAW/"$0"/"$0"_barcodes.tsv /storage2/jwlab/yerim/singlecell/GSE201348_RAW/"$0"/barcodes.tsv"}' /storage2/jwlab/yerim/singlecell/list.sh > rename.sh
parallel < rename.sh
awk '{print"mv /storage2/jwlab/yerim/singlecell/GSE201348_RAW/"$0"/"$0"_features.tsv /storage2/jwlab/yerim/singlecell/GSE201348_RAW/"$0"/features.tsv"}' /storage2/jwlab/yerim/singlecell/list.sh > rename.sh
bash rename.sh

# for unzipped files
awk '{print"mv /storage2/jwlab/yerim/singlecell/GSE201348_raw/"$0"/"$0"_barcodes.tsv.gz /storage2/jwlab/yerim/singlecell/GSE201348_raw/"$0"/barcodes.tsv.gz"}' /storage2/jwlab/yerim/singlecell/list.sh > rename.sh
parallel < rename.sh
awk '{print"mv /storage2/jwlab/yerim/singlecell/GSE201348_raw/"$0"/"$0"_features.tsv.gz /storage2/jwlab/yerim/singlecell/GSE201348_raw/"$0"/features.tsv.gz"}' /storage2/jwlab/yerim/singlecell/list.sh > rename.sh
bash rename.sh
awk '{print"mv /storage2/jwlab/yerim/singlecell/GSE201348_raw/"$0"/"$0"_matrix.mtx.gz /storage2/jwlab/yerim/singlecell/GSE201348_raw/"$0"/matrix.mtx.gz"}' /storage2/jwlab/yerim/singlecell/list.sh > rename.sh
bash rename.sh
```

- In R

```r
install.packages('Seurat')
library(Seurat)

data_dir <- '/storage2/jwlab/yerim/singlecell/GSE201348_raw/GSM6061645_A001-C-007'
list.files(data_dir)
raw_data <- Read10X(
  data.dir = data_dir,
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE)
```

## Seurat (R package) Tutorial

- [load in data](https://satijalab.org/seurat/reference/read10x) from 10X
- [guided clustering tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
- [introduction to scRNAseq](https://hbctraining.github.io/scRNA-seq/lessons/03_SC_quality_control-setup.html)

```r
# Load dataset and Initialize the Seurat object with the raw (non-normalized data)
for (file in c("GSM6061645_A001-C-007","GSM6061646_A001-C-014","GSM6061647_A001-C-023",
               "GSM6061648_A001-C-104","GSM6061649_A001-C-108","GSM6061650_A001-C-119",
               "GSM6061651_A001-C-124","GSM6061652_A001-C-202","GSM6061653_A001-C-203",
               "GSM6061654_A001-C-207","GSM6061655_A001-C-223","GSM6061656_A002-C-010", 
               "GSM6061657_A002-C-010-R0","GSM6061658_A002-C-016","GSM6061659_A002-C-021",
               "GSM6061660_A002-C-024","GSM6061661_A002-C-106","GSM6061662_A002-C-114",
               "GSM6061663_A002-C-116","GSM6061664_A002-C-121","GSM6061665_A002-C-121-R0",
               "GSM6061666_A002-C-201","GSM6061667_A002-C-203","GSM6061668_A002-C-204",
               "GSM6061669_A002-C-205","GSM6061670_A002-C-212","GSM6061671_A008-E-008",
               "GSM6061672_A008-E-015","GSM6061673_A010-E-018","GSM6061674_A010-E-023",
               "GSM6061675_A014-C-001","GSM6061676_A014-C-008","GSM6061677_A014-C-040",
               "GSM6061678_A014-C-043","GSM6061679_A014-C-052","GSM6061680_A014-C-054",
               "GSM6061681_A014-C-101","GSM6061682_A014-C-108","GSM6061683_A014-C-111",
               "GSM6061684_A014-C-114","GSM6061685_A014-C-201","GSM6061686_A015-C-001", 
               "GSM6061687_A015-C-002","GSM6061688_A015-C-005","GSM6061689_A015-C-006",
               "GSM6061690_A015-C-008","GSM6061691_A015-C-010","GSM6061692_A015-C-104",
               "GSM6061693_A015-C-106","GSM6061694_A015-C-109","GSM6061695_A015-C-202",
               "GSM6061696_A015-C-203","GSM6061697_A015-C-204","GSM6061698_A015-C-208",
               "GSM6061699_A018-E-013","GSM6061700_A018-E-020","GSM6061701_A022-E-022",
               "GSM6061702_CRC1_8810","GSM6061703_CRC2_15564","GSM6061704_CRC3_11773",
               "GSM6061705_F007","GSM6061706_F034","GSM6061707_F072B",  
               "GSM6061708_F091","GSM6061709_B001-A-301","GSM6061710_B001-A-401",
               "GSM6061711_B001-A-406","GSM6061712_B001-A-501","GSM6061713_B004-A-004",
               "GSM6061714_B004-A-008","GSM6061715_B004-A-104","GSM6061716_B004-A-204")){
  seurat_data <- Read10X(data.dir = paste0("/storage2/jwlab/yerim/singlecell/GSE201348_raw/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

tail(unique(colnames(`GSM6061645_A001-C-007`)))
# [1] "TTTCATGCACTGGAAG-1" "TTTCCTCCAAATCGTC-1" "TTTGACTGTTCGAACT-1" "TTTGATCGTTTCCCAC-1"
# [5] "TTTGTTGCAAGTCCAT-1" "TTTGTTGGTCCATAGT-1"
head(unique(row.names(`GSM6061645_A001-C-007`)))
# [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3"  "AL627309.2" 
grep("L1TD1", row.names(`GSM6061645_A001-C-007`), value = TRUE)
# [1] "L1TD1"

# subset 'L1TD1'
L1TD1_645 <- subset(`GSM6061645_A001-C-007`, features = c("L1TD1"))
View(L1TD1_645@meta.data)
L1TD1_645 = subset(x = L1TD1_645, subset = nCount_RNA > 0)
```

## ****Single-cell analyses define a continuum of cell state and composition changes in the malignant transformation of polyps to colorectal cancer****

The study aimed to understand the **changes in cell composition and state during the transformation of healthy colon to precancerous adenomas to colorectal cancer (CRC)**. To achieve this, single-cell chromatin accessibility profiles and single-cell transcriptomes were generated from 1,000 to 10,000 cells per sample for 48 polyps, 27 normal tissues, and 6 CRCs collected from patients with or without germline APC mutations.

The study found that a large fraction of polyp and CRC cells exhibit a stem-like phenotype, and a continuum of epigenetic and transcriptional changesoccur in these stem-like cells as they progress from homeostasis to CRC. 

- Advanced polyps contain increasing numbers of stem-like cells, regulatory T cells, and a subtype of pre-cancer-associated fibroblasts.
- In the cancerous state, the study observed T cell exhaustion, RUNX1-regulated cancer-associated fibroblasts, and increasing accessibility associated with HNF4A motifs in epithelia.
- DNA methylation changes in sporadic CRC were strongly anti-correlated with accessibility changes along this continuum, further identifying regulatory markers for molecular staging of polyps.

### Introduction

- need to identify genes and pathways that drive the formation of invasive cancers, with previous studies focusing on bulk profiling of advanced stage tumors and ignoring premalignant lesions.
- The article proposes that **colorectal cancer (CRC) is an ideal** system to study because it **follows a stereotyped progression from normal to atypical to carcinoma, including the formation of precancerous polyps that can subsequently give rise to CRCs.** An estimated 80-90% of colorectal tumors are initiated by loss of APC, resulting in β-catenin stabilization and increased WNT signaling leading to intestinal hyperplasia. Subsequent mutations in other cancer driver genes result in the transformation to carcinoma. Patients with familial adenomatous polyposis (FAP), who have germline mutations in APC, are a suitable population in which to study the natural progression of polyposis, as they typically develop hundreds of polyps by early adulthood, providing numerous polyps of varied molecular ages and stages of progression, all arising in the same germline background.
    - FAP (Familial adenomatous polyposis) have germline mutations in APC (Adenomatous polyposis coli): **a gene that suppresses tumor growth**. If the APC gene is defective, it makes the gene unstable and more susceptible to additional changes that may lead to colon and rectal cancers
- The article describes a study that aimed to **identify the regulatory and transcriptomic changes that occur on the phenotypic continuum from healthy colon to invasive carcinoma.** The study used single-nuclei transcriptomes and epigenomes of healthy colon, polyps, and CRCs, with many polyps obtained from patients with FAP who underwent surgical colectomies.
- The analysis of these datasets **identified immune, stromal, and epithelial cell types, and large shifts in fibroblast subpopulations that occur along the transition from normal colon to CRC.** The study also found a subpopulation of exhausted T cells present only in CRC tissue and a much larger fraction of cells exhibiting a stem-like state within polyps and CRCs.
- The study identified an epigenetic and transcriptional continuum from normal colon to CRC characterized by sequential opening and closing of chromatin and upregulation and downregulation of genes associated with the cancer state. The study identified regulatory elementsand transcription factors associated with different stages of transformation from normal colon to carcinoma, including early increases in accessibility of regions containing TCF and LEF motifs and loss of accessibility in regions containing KLF motifs. In the final stage of this pathway, malignant transformation, the study observed increased accessibility in regions containing HNF4A motifs.
- The study also showed that accessibility changes in polyps are strongly anti-correlated with DNA methylation changes in sporadic CRC and identified a subset of these regions that change their accessibility state early in the malignant continuum, suggesting **potential strategies for detection of premalignant polyps**.

## scRNA-seq

- [github introduction](https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html)

```r
#code not working
#list2 <- lapply (X = list2, FUN = function(x){
#  x <- ScaleData(x,verbose = FALSE)
#  x <- RunPCA(x,verbose = FALSE)
#  if (cluster){
#    x <- FindNeighbors(x, dims = 1:20)
#    x <- FindClusters(x, resolution = 0.5)
#  }
#  x <- RunUMAP(x, dims = 1:20)
#  return(x)
#})
```

Metadata of merged Seurat includes

- `orig.ident`: this often contains the sample identity if known, but will default to `project` as we had assigned it
- `nCount_RNA`: number of UMIs per cell
- `nFeature_RNA`: number of genes detected per cell

We need to calculate some additional metrics for plotting:

- **number of genes detected per UMI:** this metric with give us an idea of the complexity of our dataset (more genes detected per UMI, more complex our data)
- **mitochondrial ratio:** this metric will give us a percentage of cell reads originating from the mitochondrial genes

```r
merged_seurat <- merge(`GSM6061645_A001-C-007`, c(`GSM6061646_A001-C-014`,`GSM6061647_A001-C-023`,
                                                  `GSM6061648_A001-C-104`,`GSM6061649_A001-C-108`,`GSM6061650_A001-C-119`,
                                                  `GSM6061651_A001-C-124`,`GSM6061652_A001-C-202`,`GSM6061653_A001-C-203`,
                                                  `GSM6061654_A001-C-207`,`GSM6061655_A001-C-223`,`GSM6061656_A002-C-010`, 
                                                  `GSM6061657_A002-C-010-R0`,`GSM6061658_A002-C-016`,`GSM6061659_A002-C-021`,
                                                  `GSM6061660_A002-C-024`,`GSM6061661_A002-C-106`,`GSM6061662_A002-C-114`,
                                                  `GSM6061663_A002-C-116`,`GSM6061664_A002-C-121`,`GSM6061665_A002-C-121-R0`,
                                                  `GSM6061666_A002-C-201`,`GSM6061667_A002-C-203`,`GSM6061668_A002-C-204`,
                                                  `GSM6061669_A002-C-205`,`GSM6061670_A002-C-212`,`GSM6061671_A008-E-008`,
                                                  `GSM6061672_A008-E-015`,`GSM6061673_A010-E-018`,`GSM6061674_A010-E-023`,
                                                  `GSM6061675_A014-C-001`,`GSM6061676_A014-C-008`,`GSM6061677_A014-C-040`,
                                                  `GSM6061678_A014-C-043`,`GSM6061679_A014-C-052`,`GSM6061680_A014-C-054`,
                                                  `GSM6061681_A014-C-101`,`GSM6061682_A014-C-108`,`GSM6061683_A014-C-111`,
                                                  `GSM6061684_A014-C-114`,`GSM6061685_A014-C-201`,`GSM6061686_A015-C-001`, 
                                                  `GSM6061687_A015-C-002`,`GSM6061688_A015-C-005`,`GSM6061689_A015-C-006`,
                                                  `GSM6061690_A015-C-008`,`GSM6061691_A015-C-010`,`GSM6061692_A015-C-104`,
                                                  `GSM6061693_A015-C-106`,`GSM6061694_A015-C-109`,`GSM6061695_A015-C-202`,
                                                  `GSM6061696_A015-C-203`,`GSM6061697_A015-C-204`,`GSM6061698_A015-C-208`,
                                                  `GSM6061699_A018-E-013`,`GSM6061700_A018-E-020`,`GSM6061701_A022-E-022`,
                                                  GSM6061702_CRC1_8810,GSM6061703_CRC2_15564,GSM6061704_CRC3_11773,GSM6061705_F007,
                                                  GSM6061706_F034,GSM6061707_F072B,GSM6061708_F091,`GSM6061709_B001-A-301`,
                                                  `GSM6061710_B001-A-401`,`GSM6061711_B001-A-406`,`GSM6061712_B001-A-501`,
                                                  `GSM6061713_B004-A-004`,`GSM6061714_B004-A-008`,`GSM6061715_B004-A-104`,
                                                  `GSM6061716_B004-A-204`), add.cell.ids = c("CRC", "Polyp","Unaffected", "Polyp",
                                                                                             "Polyp", "Polyp","Unaffected","Polyp","Polyp",
                                                                                             "Polyp","Unaffected","Polyp","Polyp","Polyp","Polyp",
                                                                                             "Unaffected","Polyp","Polyp","Polyp","Unaffected","Unaffected",
                                                                                             "Polyp","Polyp","Polyp","Polyp","Unaffected","Polyp","Polyp","Polyp",
                                                                                             "Polyp","Polyp","Polyp","Polyp","Polyp","Unaffected","Unaffected",
                                                                                             "Polyp","Polyp","Unaffected","Unaffected","Unaffected","CRC",
                                                                                             "Polyp","Polyp","Polyp","Unaffected","Unaffected","Polyp","Polyp",
                                                                                             "Unaffected","Polyp","Polyp","Polyp","Unaffeccted","Polyp","Polyp","Polyp",
                                                                                             "CRC","CRC","CRC","Polyp","Polyp","Polyp","Polyp","Normal","Normal",
                                                                                             "Normal","Normal","Normal","Normal","Normal","Normal"))
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^CRC_"))] <- "CRC"
metadata$sample[which(str_detect(metadata$cells, "^Polyp_"))] <- "Polyp"
metadata$sample[which(str_detect(metadata$cells, "^Normal_"))] <- "Normal"
metadata$sample[which(str_detect(metadata$cells, "^Unaffected_"))] <- "Unaffected"
#add metadata back to Seurat object
merged_seurat@meta.data <- metadata
```

![1aa7d706-10ab-4604-af86-7c673adafd97](https://github.com/yerimeeei/Bioinformatics/assets/134043926/c9e208ee-b8c8-4232-bd95-bd06204f2783)

```r
# Visualize the number of cell counts per sample
metadata %>% 
  	ggplot(aes(x=sample, fill=sample)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")
```

![ef98ef80-b5af-48a5-8022-35f350b75b75](https://github.com/yerimeeei/Bioinformatics/assets/134043926/cdbcbfaf-01c3-473d-8960-19eced5d1ba2)

```r
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~sample)
```

Filtering

```r
#Cell-level filtering
# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                          subset= (nUMI >= 500) & 
                            (nGene >= 250) & 
                            (log10GenesPerUMI > 0.85)&
                            (mitoRatio < 0.20))

#Gene-level filtering
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
metadata_clean <- filtered_seurat@meta.data
```

Cell cycle scoring

```r
seurat_phase <- NormalizeData(filtered_seurat)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = cc.genes$g2m.genes, 
                                 s.features = cc.genes$s.genes)
View(seurat_phase@meta.data)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)
PC_ 1 
Positive:  CALD1, CARMN, SYNPO2, MYH11, NEGR1, TNC, MSRB3, PDZRN4, PCDH7, COL4A2 
	   ACTG2, ADAMTS9-AS2, FBXL7, MEIS2, CACNA2D1, DMD, RBPMS, TNS1, NEXN, CCBE1 
	   ARHGAP6, MEIS1, COL6A1, LMOD1, FERMT2, FILIP1, COL4A1, ARHGAP10, FHL1, CACNA1C 
Negative:  XACT, LINC00511, CFTR, MECOM, SATB2, MYH14, PPARG, SELENBP1, ADAMTSL1, B4GALNT3 
	   RBFOX1, PLS1, SHROOM3, SEMA5A, TMPRSS2, BCAS1, HHLA2, AC073050.1, RHPN2, LMO7 
	   BTNL8, SLC4A4, AL589669.1, LINC01876, OLFM4, AC019330.1, TFCP2L1, PRR5L, MUC12, NR5A2 
PC_ 2 
Positive:  TCF7L2, SHROOM3, PDE3A, PCSK5, PRKG1, HHLA2, MYH14, MICAL2, PPARG, FLNB 
	   PLS1, TMPRSS2, DST, NEAT1, LINC00511, SCNN1B, PDE9A, B4GALNT3, NEDD4L, SATB2 
	   MS4A12, SLC4A4, BTNL8, GCNT3, CEACAM1, SVIL, GDA, MAST2, C1orf21, MYLK 
Negative:  ARHGAP15, PTPRC, CHST11, RIPOR2, DOCK2, ANKRD44, IKZF1, RHOH, RCSD1, PRKCB 
	   BACH2, BLK, BANK1, FLI1, IKZF3, CD53, PPP1R16B, TMEM131L, CCND3, FCRL1 
	   DOCK8, AFF3, APBB1IP, MARCH1, PAX5, EBF1, KCNQ5, INPP5D, MS4A1, CARD11 
PC_ 3 
Positive:  MS4A12, PCSK5, FKBP5, PAG1, HHLA2, PDE3A, SCNN1B, GCNT3, PLAC8, TRPM6 
	   CLCA4, TMIGD1, DHRS9, SLC17A4, GDA, MXI1, CLIC5, PTPRR, CEACAM1, EPB41L3 
	   PRSS12, UGP2, CPNE8, SGK1, PDE9A, RNF152, HPGD, CEACAM7, HDAC9, FOXO1 
Negative:  DACH1, NKD1, CEMIP, OLFM4, ZNRF3, GRIN2B, LGR5, PGGHG, SNTB1, SAMD5 
	   MT-CO2, GRM8, LINC01811, APCDD1, MYRIP, NRXN3, MT-CO3, CENPF, CADPS, PRUNE2 
	   HSP90AA1, HSPD1, MKI67, APOLD1, AL033504.1, XACT, AC009975.1, AC099792.1, CLCA1, LINC01594 
PC_ 4 
Positive:  RBFOX3, PDZRN4, CCBE1, GPM6A, MYH11, AC131025.2, ACTG2, CARMN, GJC1, EPHA7 
	   JPH2, MYOCD, LDB3, TNC, CPXM2, LMOD1, LINC01798, LINC00578, EPHA6, DES 
	   ATRNL1, AL031663.3, CNN1, LMO3, PRKCB, SYNPO2, SYNM, REEP1, RBM20, PCA3 
Negative:  GPC6, LAMA2, FBN1, TSHZ2, PREX2, BICC1, PRR16, LDB2, ZNF521, SVEP1 
	   NR2F1-AS1, GAS7, CFH, VCAN, ADGRL3, PAPPA, MIR99AHG, LIFR, PIEZO2, SGIP1 
	   DCN, ABCA8, EPHA3, NRP1, PDE7B, KCNT2, MSC-AS1, TFPI, UNC5C, LHFPL6 
PC_ 5 
Positive:  MGAM2, AOAH, PTPRR, AC011287.1, MT-CO2, MT-CO3, GDA, SLCO2A1, TRIM31, SLC26A3 
	   CEACAM1, EMP1, DST, CLCA4, SLC30A10, SLC15A1, ABCB1, PLAUR, MT-ND3, FYB1 
	   SLC9A3, LINC01811, ABCG2, MTRNR2L12, NKD1, HKDC1, TMIGD1, LEF1, MT-ND4, THEMIS 
Negative:  EFNA5, RBFOX1, BRINP3, CNTN4, SATB2, DPP10, AC019330.1, EDIL3, FOXP2, IQCM 
	   MECOM, ADAMTSL1, LINC01876, PDE3A, CACNB2, TOX, PPARGC1A, EYA2, LINC01687, HTR4 
	   HMGCS2, NEDD4L, PRKG1, SLC4A4, PLXDC2, PTGDR, CFTR, AC113383.1, NR5A2, AL139383.1
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
PCAPlot(seurat_phase,
        +         split.by = "sample")
seurat_phase <- FindNeighbors(seurat_phase, dims = 1:30)
seurat_phase <- FindClusters(seurat_phase, resolution = 0.5)
seurat_phase <- RunUMAP(seurat_phase, dims = 1:40,reduction = "pca")
DimPlot(seurat_phase)
```

![0fca79c6-e94e-4fb5-8109-d2d1e051f8e0](https://github.com/yerimeeei/Bioinformatics/assets/134043926/5dedcf43-5956-4201-b49b-e8d3950f8024)

![22865a1f-83c6-4476-89cb-f2fc02f79ab2](https://github.com/yerimeeei/Bioinformatics/assets/134043926/c6a220e9-2b05-46a5-8c31-4b9f5990b51a)

![050ebd02-3dbe-416a-8c30-6acc262f1bda](https://github.com/yerimeeei/Bioinformatics/assets/134043926/d73f7fb4-44dc-4869-9f5b-6af52d6e0dfc)

To check

```r
"L1TD1"     "INSL5"     "CCN1"      "GSTM5"     "DPT"    "PRELP"     "GREM2"    
"MEIS1"     "ECRG4"     "CYP27A1"   "TSC22D3"   "PDZD4"  "ZNRF3"     "FBLN1"
"DES"       "SCG2"      "TDGF1"     "CCDC80"    "BCHE"   "FGFRL1"    "ART3"     
"BMP3"      "SFRP2"     "HAND2"     "FERMT1"    "CCN5"   "HSPB6"     "ADAM33"
"C7"        "ENC1"      "MGC32805"  "SLC12A2"   "SLIT3"  "MSX2"      "FOXQ1"   
"LY6G6F"    "VIP"       "GTF2IRD1"  "TM6SF2"    "CD22"   "CNN1"      "NOTUM" 
"ACHE"      "MET"       "AKR1B1"    "FABP4"     "PCAT2"  "IQANK1"    "CCL21"
"OGN"       "HMCN2"     "ITIH5"     "SCML2P1"   "PALM"   "LINC00668" "FADS6"
"NEBL"      "APIP"      "TCN1"      "SLC22A11"  "SPTBN2" "MMP7"      "CRYAB"
"APOA4"     "GRIN2B"    "PRPH"      "AXIN2"     "NGFR"   "RNF43.1"   "PER1"
"KRT80"     "PTPRR"     "SMAD9"     "SLC15A1"   "F10"    "NDN"       "CEMIP"
"NETO2"     "NKD1"      "CDH3"      "SLC9A3R1"
```
