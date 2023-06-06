# Day 2 L1TD1 gene in TCGA samples

# To Start

- Read 4 literatures
- L1TD1 summary of TCGA samples comparing between normal vs tumor (11A normal patient, 01A tumor patient, 06A metastasis)
- Normal (TCGA) vs adenoma (GEO~~) vs tumor (TCGA)

# To Start using raw data

1. Download GDC Colon cancer data from [https://xenabrowser.net/](https://xenabrowser.net/)
2. Get the [HTSeq - FPKM-UQ](https://xenabrowser.net/datapages/?dataset=TCGA-COAD.htseq_fpkm-uq.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) and Phenotype data.
3. Extract expression for L1TD1 (ENSG00000240563.1)
4. Divide samples into normal and tumour. Then for tumour divide into pathological Stage
5. Make boxplot or violin plot comparing expression of L1TD1 for normal, and tumour stages

# Day 2 06-06-2023

### Literature 2: L1TD1 - a prognostic marker for colon cancer

### **Literature Review**

### Methods

1. **Microarray data sets**
2. **Gene expression analysis**
3. **Gene list descriptions**
    - Interaction partners
        - 311 interaction partners of L1TD1 determined
            - 306 identified **by performing a mass spectrometry analysis** on co-immunoprecipitated proteins with 2 different anti-L1TD1 antibodies.
            - 5 proteins (NANOG, OCT4 (POU5F1), SOX2, DNMT3B, and TRIM28), which were challenging to detect using mass spectrometry, interactions were shown using **immunoprecipitation** and **Western Blotting**.
        - 285/311 corresponded to genes that had probes associated to them in the microarray platforms used in this study.
    - Top 20 Interaction partners and Top 20 co-expressed genes with L1TD1 in colon cancer determined on the basis of their co-expression with L1TD1 in the seminoma and stem cell data sets.
    - **Question**
        - **What are the implications that comparing L1TD1 with Interaction partners and co-expressed genes?**
        - **Co-expression of L1TD1 with its currently known interaction partners was investigated to understand distinctive role of L1TD1 in different cancers.**
4. **Survival analysis of microarray data**
    - Disease-free survival (DFS) analysis: **Kaplan-Meier method as implented in R package “survival”**
    - Survival curves plotting: **R package “survminer”**
    - Log-rank test used to compare survival rates between L1TD1 high and L1TD1 low groups.
5. **Association between L1TD1 expression and clinicopathological variables**
    - clinicopathological variables: cancer stage, prior therapy received by patients, tumor location, chromosomal instability, CpG island methylation status, DNA mismatch repair proficiency, mutation status of BRAF, KRAS, and TP53.
        - **Question: what are mutation status of BRAF, KRAS, and TP53?**
            - **BRAF** is a gene that encodes a protein BRAF kinase, which is involved in the regulation of cell growth and division. When the BRAF gene is mutated, it can lead to the activation of this protein and contribute to the **development of cancer**.
            - **KRAS** is a gene that encodes a protein KRAS, which is involved in the regulation of cell growth and division. When the KRAS gene is **mutated**, it can lead to the activation of this protein and contribute to the **progression of cancer**.
            - **TP53** is a **tumor suppressor gene** that codes for a protein p53, which plays a key role in regulating cell division and preventing the formation of cancerous cells. When TP53 is mutated, it can lead to the loss of p53 function, allowing cells to form tumors.
    - Wilcoxon rank sum test for variables with only 2 categories, Kruskal-Wallis test for variables with more than 2 categories, and Pearson correlation for association with age.
6. **Analysis of TCGA Colon adenocarcinoma RNA-seq dat set**
    - RNA-seq data from TCGA-COAD qcruied from gdc → FPMK-UQ normalized RNA-seq counts used to validate correlation analyses → fitted 2 Gaussian distributions and evaluated 2 different thresholds → 2 thresholds used to perform survival analysis

### Results

1. **High expression of L1TD1 associates with longer DFS**
    - Across **colon cancer** microarray data sets: 26.7% of **patients to have high L1TD1** expression, which was **lower than observed in seminoma and stem cell.**
    - Kaplan-Meier analysis: colon cancer samples with **high L1TD1 expression had longer DFS** compared to no/low L1TD1 expression.
    - **L1TD1 expression higher from early cancer stages** compared to later stages (but differences between later stages not statistically significant).
    - L1TD1 expression high for samples with mutated KRAS, WT TP53, and negative chromosomal instability marker.
    - **Significant associations between L1TD1 expression and tumor location** (tumor differentiation status), while no statistically significant associations between L1TD1 expression and age, sex, prior therapy, BRAF mutation status, CpG island methylation status, or DNA mismatch repair proficiency observed.
    
    **L1TD1 expression higher from early cancer stages and for samples with mutated KRAS, WT TP53, and negative chromosomal instability marker. L1TD1 expression has significant associations with tumor location other than clinicopathological variables.**
    
    - **Question**
        - **what are samples with negative chromosomal instability marker?**
            - **Chromosomal instability (CIN)** is a phenomenon in which cells accumulate abnormal numbers of chromosomes or chromosomal rearrangements, which can contribute to the **development of cancer**. Therefore, a **negative CIN marker** would indicate that there is **no evidence of chromosomal instability** in the sample being tested.
        - **Mutated KRAS** can lead to the activation of this protein and contribute to the progression of cancer while **WT TP53** is a tumor suppressor gene that codes for a protein p53, which plays a key role in regulating cell division, thus preventing the formation of tumor. **Why L1TD1 expression high for both mutated KRAS and WT TP53? Why L1TD1 expression high for samples with muated KRAS, which can lead to tumor cell?**
2. **Interactome of L1TD1 is not co-expressed in colon cancer**
    - potential role of previously identified interaction partners of L1TD1, Spearman rank correlation matrices calculated.
    - Result: High positive correlatin observed among L1TD1 and its top 20 interaction partners in seminoma and stem cell data sets was **absent in al 3 colon cancer data sets**. However, interaction partners did not consistenly improve the predictive prognostic power obtained with L1TD1 alone.
3. **Genes co-expressed with L1TD1 in colon cancer**
    - **None of top 20 co-expressed genes outperformed L1TD1 as independent prognostic marker** for colon cancer in all 3 data sets.
    - **5 genes had statistically significant impact on surviva**l in at least 2/3 colon cancer data sets.
        - SPINK4, RETNLB, ASRGL1, CLCA1, FCGBP

**None of the top 20 co-expressed genes were found to be a better independent prognostic marker than L1TD1. 5 genes had impact on survival for some colon data sets.**

4. **Validation in TCGA Colon adenocarcinoma RNA-seq data set**
    - Analysis of TCGA-COAD for validate findings
    - Correlation analysis: indicated lack of correlation between L1TD1 and top 20 interacting partners, confirmed significant correlation between **L1TD1 and genes that were co-expressed with L1TD1 in the colon cancer microarray data sets.**

### Discussion

- the study found evidence of L1TD1 being postiive prognostic marker for colon cancer.
    - positive prognostic marker: characteristic that is associated with a better prognosis or outcome in patients with a particular disease.
- High L1TD1 expression
    - Healthy tissues: embryonic stem cells, brain, and colon.
    - Besides: seminoma, embryonic carcinomas, medulloblastoma, and **colon adenocarcinoma**.
- Hypothesis: high expression of L1TD1 in colon cancer might be associated with prognosis.
- Result from study suggested in colon cancer, high expression of L1TD1 is linked to better prognosis.
    - expression of L1TD1 was associated with samples of low clinical cancer stage → perheps be a reason for its prognostic significance.
- Co-expression of L1TD1 with its currently known interaction partners was investigated to understand distinctive role of L1TD1 in different cancers.
    - outcome that L1TD1 was **not co-expressed with its interaction partners** in **colon cancer** points to potential participation of L1TD1’s interaction partners in the **contrasting prognostic outcome**.
        - recent study shown **association of high L1TD1 expression with poor clinical outcome** and significant **co-expression between L1TD1 and its interaction partner OCT4**.
    - Finding suggest: co-expression of L1TD1 with its interaction partners might be **required for clearly showing an aggressive and detrimental phenotype**.
        - aggressive phenotype refers to a set of observable traits in an organism that are associated with aggression.
        - detrimental phenotype refers to a set of characteristics that are harmful to the organism's survival and reproductive success

## Using Xenabrowser

- [heat map of TCGA colon cancer samples of L1TD1 gene expression by tumor stages](https://xenabrowser.net/heatmap/)
    
    <img width="474" alt="heat_map_of_TCGA_colon_cancer_samples_of_L1TD1_gene_expression_by_tumor_stages" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/46d4a809-10c8-491a-b2aa-f2b9483b4aae">

- Kaplan Meier plot

<img width="850" alt="Kaplan_Meier_gene_expression_RNAseq_-_HTSeq_-_FPKM-UQ" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/a2e7e302-68cc-45d2-9076-6c654c1790e8">

- Downloaded the GDC Colon cancer data to get HTSeq - FPKM-UQ and Phenotype data.

### **Expression DIY using GEPIA in different tumor stages of L1TD1 gene expression in [COAD and READ](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations)**

<img width="466" alt="L1TD1_stageplot" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/58816510-085b-420d-9c54-621216de07a2">

## Using R

- [x]  Extract expression for L1TD1 (ENSG00000240563.1)
- [x]  Divide samples into normal and tumour. Then for tumour divide into pathological Stage
- [ ]  Make boxplot or violin plot comparing expression of L1TD1 for normal, and tumour stages
1. **Extract expression for L1TD1**

```r
library(stringr)
TCGA.COAD.GDC_phenotype <- read.delim("~/Desktop/intern2023/TCGA-COAD.GDC_phenotype.tsv")
TCGA.COAD.htseq_fpkm.uq <- read.delim("~/Desktop/intern2023/TCGA-COAD.htseq_fpkm-uq.tsv")

## extract expression for L1TD1 (ENSG00000240563.1)
L1TD1_exp <- TCGA.COAD.htseq_fpkm.uq[TCGA.COAD.htseq_fpkm.uq$Ensembl_ID ==
																		 "ENSG00000240563.1",]

## change - to . of submitter_id.samples name in phenotype file to match
TCGA.COAD.GDC_phenotype$submitter_id.samples <- gsub("-",".",
										TCGA.COAD.GDC_phenotype$submitter_id.samples, fixed=TRUE)
```

2. **Extract expression for “A”**

```r
## Only extracting "A" 
# colnames(TCGA.COAD.htseq_fpkm.uq)[str_detect(colnames(TCGA.COAD.htseq_fpkm.uq), "A")]
View(TCGA.COAD.GDC_phenotype[,c(1,116)])
L1TD1_get <- TCGA.COAD.GDC_phenotype[,c(1,116)]
L1TD1_get$extract <- substr(L1TD1_get$submitter_id.samples, 14,16)
L1TD1_get <- L1TD1_get[str_detect(L1TD1_get$extract, "A"),]
```

3. **Divide samples into normal (11A) and tumor.**

```r
## Extract L1TD1 Solid Tissue Normal patient
L1TD1_nor <- L1TD1_get[str_detect(L1TD1_get$extract, "11A"),]
View (L1TD1_nor)
```

4. **Divide tumor samples into pathological stages: primary tumor (01A), metastasis (06A) and recurrent tumor (02A).**

```r
## L1TD1 Primary Tumor patient
L1TD1_tum <- L1TD1_get[str_detect(L1TD1_get$extract, "01A"),]
View(L1TD1_tum)

## L1TD1 Metastatic patient
L1TD1_met <- L1TD1_get[str_detect(L1TD1_get$extract, "06A"),]
View(L1TD1_met)

## To view unique values in 'extract' column
unique(L1TD1_get$extract)

## L1TD1 Recurrent Tumor patient
L1TD1_rec <- L1TD1_get[str_detect(L1TD1_get$extract, "02A"),]
View(L1TD1_rec)
```

- **Question**
    - **Can I make boxplot or violin plot from multiple different tables in R?**
    - **How to match phenotype data with HTSeq - FPKM-UQ?**
5. **Make boxplot or violin plot comparing expression of L1TD1 for normal, and tumour stages.**

### Duty for Tmr

- Read literature to understand and make questions to ask.
- Try using R to match phenotype data with HTSeq and make boxplot or violin plot.
    - study R commands used for the code.
- Ask question to Sandy & Dr. JW
    - regarding L1TD1 gene expression
    - regarding using the R (trial & progress).
