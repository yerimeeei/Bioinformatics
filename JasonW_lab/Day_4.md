# Day 4 L1TD1 gene in TCGA colorectal samples

# To Start

- Read 4 literatures
- Normal (TCGA) vs adenoma (GEO~~) vs tumor (TCGA)
    - Using raw data from GDC Colon cancer from [https://xenabrowser.net/](https://xenabrowser.net/).
    - Get the [HTSeq - FPKM-UQ](https://xenabrowser.net/datapages/?dataset=TCGA-COAD.htseq_fpkm-uq.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) and Phenotype data.
        - **HTSeq-FPKM-UQ**: method for quantifying gene expression levels from RNA sequencing (RNA-seq) data. RNA-Seq-based expression normalization method across different samples.
        - Resulting HTSeq-FPKM-UQ values provide a **measure of gene expression that is normalized for differences in library size, gene length, and the influence of highly expressed genes.**
    - Extract expression for L1TD1 (ENSG00000240563.1)
    - Divide samples into normal and tumour. Then for tumour divide into pathological Stage
    - Make boxplot or violin plot comparing expression of L1TD1 for normal, and tumour stages

### Today

- [ ]  Read literature of L1 cancer basic in the morning and bring questions to ask.
- [x]  Meeting with Sandy & Dr. JW 10am.
- [x]  Using R to plot t-test graph.
- [x]  Draw violin plot for different pathological stages of tumor.
- [ ]  Find a paper or dataset that has all normal, tumor, and adenoma data samples to analyze (normalization)

# L1TD1 gene expression TCGA data Analysis

## Analysis Using R ([COAD](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations) & READ)

### T-Test Chart between normal and tumor samples (COAD & READ)

<img width="1100" alt="Violin_plot_of_L1TD1_expression_in_both_COAD_and_READ_(Normal_vs_Tumor)_with_T_test" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/a5cb170e-3e74-4ec0-b6a2-73f5c8a6d0b4">

<img width="517" alt="Violin_plot_comparing_expression_of_L1TD1_in_READ_with_T_test" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/de8a1ea8-9a6e-4532-b45f-8c4db2023fa9">

## Analysis by pathological stages of cancer Using R ([COAD](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations) & READ)

- Cancer Staging: way of describing where the cancer is located, if or where it has spread, and whether it is affecting other parts of the body.
- **TNM staging system**
    - **Tumor (T):** “T” plus a letter or number (0 to 4) is used to describe the size and location of the tumor. Tumor size is measured in centimeters (cm).
        
        
        | TX | The primary tumor cannot be evaluated. |
        | --- | --- |
        | T0 | There is no evidence of a tumor. |
        | T1 | The tumor is 2 centimeters (cm) or smaller and limited to the thyroid. |
        | T1a | The tumor is 1 cm or smaller. |
        | T1b | The tumor is larger than 1 cm but less than 2 cm. |
        | T2 | The tumor is larger than 2 cm but smaller than 4 cm and is limited to the thyroid. |
        | T3 | The tumor is larger than 4 cm, but the tumor does not extend beyond the thyroid gland. |
        | T4 | The tumor is any size and has extended beyond the thyroid. |
        | T4a | The tumor has spread beyond the thyroid to nearby soft tissues, the larynx, trachea, esophagus, or recurrent laryngeal nerve. |
        | T4b | The tumor has spread beyond the regions in T4a (above). |
        | Tis | cancer cells are only growing in the layer of cells where they started, without growing into deeper layers |
    - **Node (N):** “N” in the TNM staging system stands for lymph nodes. Careful evaluation of lymph nodes is an important part of staging thyroid cancer. There are many regional lymph nodes located in the head and neck area. Lymph nodes in other parts of the body are called distant lymph nodes.
        
        
        | NX | The regional lymph nodes cannot be evaluated. |
        | --- | --- |
        | N0 | There is no evidence of cancer in the regional lymph nodes. |
        | N1 | Cancer has spread to the lymph nodes. |
        | N1a | Cancer has spread to the lymph nodes around the thyroid (called the central compartment; the pretracheal, paratracheal, and prelaryngeal lymph nodes). |
        | N1b | Cancer has spread beyond the central compartment, including unilateral cervical (lymph nodes on 1 side of the neck), bilateral cervical (lymph nodes on both sides of the neck), contralateral cervical (the opposite side of the tumor), or mediastinal (the chest) lymph nodes. |
    - **Metastasis (M):** “M” in the TNM system describes whether cancer has spread to other parts of the body, called metastasis.
        
        
        | MX | Distant metastasis cannot be evaluated. |
        | --- | --- |
        | M0 | Cancer has not spread to other parts of the body. |
        | M1 | Cancer has spread to other parts of the body. |
    
    <img width="513" alt="Violin_plot_comparing_expression_of_L1TD1_in_COAD_with_T_test" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/c787c7dd-636e-4f38-935e-4a861a5b716d">

    <img width="1438" alt="Violin_plot_comparing_expression_of_L1TD1_in_READ_pathological_T_stages" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/aa1b9217-32b5-42b1-9101-49d3156d5613">


**Summary**

1. Violin plot of L1TD1 expression in both COAD and READ (Normal vs Tumor)

```r
library(gridExtra)
mer_L1TD1_vplot <- grid.arrange(L1TD1_vplot_C, L1TD1_vplot_R, ncol = 2)
```

<img width="1133" alt="Violin_plot_of_L1TD1_expression_in_both_COAD_and_READ_(Normal_vs_Tumor)" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/427c5f94-a4dd-4fd0-8835-55e826f8ec6a">

1. Violin plot of L1TD1 expression in both COAD and READ (Normal vs Tumor) with T-test

```r
library(gridExtra)
mer_L1TD1_vplot_T <- grid.arrange(L1TD1_vplot_C_T, L1TD1_vplot_R_T, ncol = 2)
```

<img width="1100" alt="Violin_plot_of_L1TD1_expression_in_both_COAD_and_READ_(Normal_vs_Tumor)_with_T_test" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/7765d56c-efbc-43e1-9f64-4bb01762e47d">

1. Violin plot of L1TD1 expression in both COAD and READ (pathological T stages)

```r
library(gridExtra)
mer_L1TD1_vplot_patT <- grid.arrange(L1TD1_vplot_patT_C, L1TD1_vplot_patT_R, 
												ncol = 2)
```

<img width="1436" alt="Violin_plot_of_L1TD1_expression_in_both_COAD_and_READ_(pathological_T_stages)" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/e0e19746-4146-42b2-9af4-6eaaaf75e472">

- Sample number

| Type | COAD | READ |
| --- | --- | --- |
| “” | 2 | 5 |
| T1 | 11 | 11 |
| T2 | 92 | 33 |
| T3 | 379 | 129 |
| T4 | 38 | 5 |
| T4a | 20 | 8 |
| T4b | 12 | 1 |
| Tis | 1 |  |

### Summary (to Dr.JW and Sandy)

1. COAD, READ (01A, 11A, 02A, 06A) violin plots

```r
library(gridExtra)
mer_L1TD1_vplot <- grid.arrange(L1TD1_vplot_C, L1TD1_vplot_R, ncol = 2)
```

<img width="1133" alt="Violin_plot_of_L1TD1_expression_in_both_COAD_and_READ_(Normal_vs_Tumor)" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/2efa95ba-91f3-421d-8bb3-8c82c9a818e5">

1. COAD + READ combined (01A, 11A, 02A, 06A) violin plot
    
    <img width="580" alt="COAD__READ_combined_(01A_11A_02A_06A)_violin_plot" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/3b0e4cff-20bd-4c9a-b6ac-c5c511b1de92">

2. COAD, READ pathological stages violin plots + bar charts showing proportion

```r
mer_L1TD1_patT <- grid.arrange(L1TD1_vplot_patT_C, L1TD1_vplot_patT_R, 
									L1TD1_patT_C_bar, L1TD1_patT_R_bar, nrow = 2, ncol = 2)
```

<img width="1243" alt="COAD_READ_pathological_stages_violin_plots__bar_charts_showing_proportion" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/253347ea-35bb-4a46-bacf-7227a75cc3e6">

1. COAD+READ combined pathological stages violin plot + bar chart showing proportion

```r
mer_CRC_L1TD1_patT <- grid.arrange(CRC_L1TD1_NT_vplot, CRC_L1TD1_patT_bar, 
											ncol = 2)
```

<img width="1198" alt="COADREAD_combined_pathological_stages_violin_plot__bar_prop" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/e63c0ad8-91e1-490a-ad3f-e87873b02f26">

1. COAD+READ combined pathological stages (combined T4, T4a, T4b) violin plot

```r
CRC_L1TD1_patT_vplot_T4
```

<img width="930" alt="COADREAD_combined_pathological_stages_(combined_T4)_violin_plot" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/9b41b31d-afd8-42ab-a086-5e6fa9c1adcf">

# TODO

- [ ]  Read literature of L1 cancer basic in the morning and bring questions to ask.
1. Normal, Adenoma, Primary Tumor, Meta CRC RNA-seq (transcriptome) 에 해당되는 GEO/study 찾아보기
2. STAR aligner go through ([https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html))
3. FASTQ, BAM이 무엇인지 공부해보기
4. Try to run STAR, featureCounts on one FASTQ sample

```r
/storage2/temp/labs/jwlab/adenoma/fastq
```

- Batch Effect ([using combat in R](https://www.rdocumentation.org/packages/sva/versions/3.20.0/topics/ComBat))
- [GEO data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164541) that has normal, adenoma, and primary tumor
- Installing [SRA](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) Toolkit
- Study [STAR](https://github.com/alexdobin/STAR) aligner (STAR → BAM → count ? (fc_output in R))
    - [STAR manual](https://gensoft.pasteur.fr/docs/STAR/2.7.3a/STARmanual.pdf)

# Literature 3: LINE-1 in cancer basic

## **Literature Review**

### Background

- Interspersed repetitive sequences originated from retrotransposons.
    - **Restrotransposons** that have remained unequivocally active: L1, Alu, and SVA.
        - Retrotransposons are a type of transposable element, which is a segment of DNA that can move from one location to another within a genome. Retrotransposons are unique in that they use a "copy and paste" mechanism to move.
    - L1: autonomous (capable of self-propagating through RNA intermediates).
    - Alu, SVA: nonautonomous → rely on L1 mobilization.
- **L1**
    - majority **lost retrotransposition competency** due to 5’ truncations, inverted rearrangements, or point mutations (occur during reverse transcription or subsequent chromosomal replication of the inserted element).
        - hot L1s: comprising majority of above activity (~5-10% of active elements).
        - active L1s utilize “copy-and-paste” mechanism to insert themselves throughout the genome (potential disruptive effects on neighboring genes or regulatory sequences).
            
            → active L1s keep **reshaping human genome** and become a **source of endogenous mutagenesis** (cause individual genome variation) and can **participate in the pathogenesis** (occurrence) **of many genetic disease** (include cancer).
            
- Cancer: disease resulting from accumulated genetic mutations (L1 can be one contributor).
- In this paper:
    - discuss putative multilayered functions of L1s in cancer
    - potential of L1s for clinical implications

### Duty for Tmr

- Read literature of L1 cancer basic and make questions to ask.
- Using R to plot N and M stages of combined COAD + READ.
- find some datas from GEO or study
- Study STAR aligner, FASTQ, BAM and run STAR, featureCounts on one FASTQ sample
