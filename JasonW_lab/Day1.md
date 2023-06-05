# Day 1 L1TD1 gene in TCGA samples

# To Start

1. Read 4 literatures
2. Using DrBioRight look at L1TD1 summary of TCGA samples in general
3. Using DrBioRight look at L1TD1 summary of TCGA samples comparing between normal vs tumor (11A normal patient, 01A tumor patient)
4. Normal (TCGA) vs adenoma (GEO~~) vs tumor (TCGA)

# Day 1 05-06-2023

### Literature 1: **DrBioRight**

- has immense potential to increase efficiency and reproducibility of omics data analysis.
- Two subsystems: 1) user-friendly web interface. 2) backend compute server.
1. Mocualized framework
    - has flexible modualized framework, based on which a new computational analysis can be added with steps: 1) adding necessary modules. 2) training the modules using natural human languages.
2. NGS analysis
    - supports bioinformatics analysis from raw next-generation sequencing reads.
    - analysis include: quality control, read mapping, gene expression quantification, differential expression analysis, and gene set enrichment analysis.
3. Reproducibility analysis 
    - enables users to conveniently check the reproducibility of publilshed results.
- **5 Key features of next-generation data analytics**
    1. *Natural Language Understanding (NLU)*
    2. *Artificial Intelligence (AI)*
    - data-driven predictive models: identify appropriate datasets and algorithms, select informative visualization.
    1. Transparency
    - reproducibility: next-generation analytics should be able to generate detailed analysis, ensuring they are transparent and that obtained results are reproducible.
    - important to provide functionalities that allow users to check the reproducibility of omics results from studies.
    1. Mobile and Social Media Friendliness
    - Mobile: no restriction on place and time
    - Social: online chat interface
    1. Crowdsourcing
    - build open-development user center
    - build data-sharing system

### Using DrBioRight

- Summary about **tcga** **pancan** **L1TD1** **mrnaseq**
    - **Question arise**
        1. **How to interpret the plot that DrBioRight provided?**
        2. **How to know the code DrBioRight used to generate scatter box plot? (When I want to make my own code to analyze)**

<img width="1441" alt="Screenshot_2023-06-05_at_3 02 46_PM" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/91603101-4167-4aab-a48a-413b040a38d2">

![tcga_pancan_L1TD1_mrnaseq](https://github.com/yerimeeei/Bioinformatics/assets/134043926/a57fdd9b-7ec8-4f5c-9a01-def76c920ab6)

### Literature 2: L1TD1 - a prognostic marker for colon cancer

### **Background**

- **L1TD1 (LINE-1 type transposase domain containing 1)**
    - encodes for a protein that contains a domain similar to the transposase domain of LINE-1 retrotransposons, which are mobile genetic elements that can move around in the genome.
    - potential tumor suppressor functions that may have implications for cancer diagnosis and treatment.
    - A **transposase domain** is a protein domain that is responsible for the **movement of transposable elements** (TEs) within a genome.
        - **TEs** are **DNA sequences that can move** from one location to another within a genome, and they are found in all organisms, ranging from bacteria to humans.
- **Colorectal Cancer (colon cancer)**
    - type of cancer that forms in the colon or rectum.
    - colon: large intestine, which is a long, muscular tube that absorbs water and nutrients from the food we eat.
    - rectum: lower part of the colon, which connects to the anus and helps control bowel movements.
    - Figure: Arise of tumor-initiating cells from aberrant colon crypt and subsequent transition of early Figure 1. Arise of tumor-initiating cells from aberrant colon crypt and subsequent transition of early
    adenoma to metastatic cancer.
    
    ![colorectal_adenomas](https://github.com/yerimeeei/Bioinformatics/assets/134043926/8d72ca30-a790-42bd-989b-6f2c7f6965f3)

- finding suggest: **L1TD1 may have tumor suppressor functions in colon cancer** and could be a potential target for future therapies.
    - Research has shown that L1TD1 expression is downregulated in colon cancer tissues compared to normal colon tissues. Loss of L1TD1 expression has been associated with more advanced stages of colon cancer and poorer patient outcomes.
    - One study found that L1TD1 can **inhibit colon cancer cell proliferation** and **induce apoptosis** (programmed cell death) in vitro and in vivo.
    - Another study suggested that L1TD1 can **suppress the growth and metastasis of colon cancer cells** by **inhibiting the activity of a protein** called **AKT**, which is involved in cell survival and growth.

### Using DrBioRight

- [survival analysis of **L1TD1** **gene expression** in **tcga** **colon** patients](https://drbioright.org/reports/647d934304d0b297134fcf45/report.html)

<img width="1445" alt="survival_analysis_of_L1TD1_gene_expression_in_tcga_colon_patients" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/711e6612-b0b5-4bfb-bdce-21c9dcad42a3">

- **Pipes in R** (code given by DrBioRight has %>% pipe)
- tried code of ****Summary of survival type and patient number****

```r
#params$report_data$survtype <- as.character(params$report_data$survtype)
#params$result[, "Survival type"] <- as.character(params$result[, "Survival type"])
> out_data <- params$report_data %>%
#mutate(Disease = toupper(disease)) %>%
+ mutate(Disease = as.character(toupper(disease))) %>%
#mutate(Survtype = survtype) %>%
+ mutate(Survtype = as.character(survtype)) %>%
+ group_by(Disease, Survtype) %>%
+     summarise(`Total patients` = n(),
								`Missing survival time` = sum(is.na(time), na.rm = T), 
								`Missing survival status` = sum(is.na(status), na.rm = T),
								`Survival time equals zero`= sum(time==0, na.rm = T),
								`Patients with valid survival data` = sum(!is.na(time)
								 & !is.na(status) & time>0, na.rm = T)) %>%
								left_join(params$result %>%
+ left_join(params$result %>%
+     select(`Survival type`, Disease, `Num of samples`) %>%
+     mutate(`Survival type` = as.character(`Survival type`), 
		Disease = as.character(Disease)), 
		by=c("Survtype"="Survival type", "Disease"="Disease")) %>%
+ rename(`Patients with valid survival and valid biomarker data` = 
					`Num of samples`) %>%
+ mutate(`Patients with valid survival and valid biomarker data` = 
					as.numeric(as.character(`Patients with valid survival and 
					valid biomarker data`)))

# my_function <- function(data, column) {
#   print(column)
#   for (each_col in column) {
#     print(each_col)
#     quo_column <- enquo(each_col)
#     new_col1 <- paste0("Patients with ", quo_column, " data")[2]
#     new_col2 <- paste0("Patients with valid ", quo_column, " data")[2]
#     print(new_col1)
#     print(data %>%
#       mutate(Disease = toupper(disease)) %>%
#       mutate(Survtype = survtype) %>%
#       group_by(Disease, Survtype) %>%
#       summarise(!!new_col1 := sum(!is.na((!!quo_column))))
#     )
#   }
# }
# 
# my_function(report_data, c(feature_name))

+ datatable(out_data)
```

`Error: unexpected symbol in:"    mutate(`Patients with valid survival 
	and valid biomarker data` = as.numeric(as.character(`Patients with 
	valid survival and valid biomarker data`)))datatable"`

- **Question**
    - **I want to perform the L1TD1 gene expression profile of the TCGA colon tumor patients and compare it with normal patients using DrBioRight. Do you know how?**
- **Lackage of DrBioRight**
    - **only get it to summarise across tumour types and differential analysis, but it refuses to give me a boxplot comparing the expression of L1TD1 by tumour stage**
    - To resolve: Use GEPIA2 to start analysis with.

### Using GEPIA2

- [Differential Expression Analysis of COAD](http://gepia2.cancer-pku.cn/#degenes)

![chr_dist_COAD_vDbaW](https://github.com/yerimeeei/Bioinformatics/assets/134043926/0d1a9fc8-aa14-4f51-94cb-16fb0fbe5a50)

- **[Expression DIY of L1TD1 gene in COAD](http://gepia2.cancer-pku.cn/#analysis) matching with TCGA normal data**

![L1TD1_profile_XkINC](https://github.com/yerimeeei/Bioinformatics/assets/134043926/52b57ed3-4dc8-45b0-8278-7f592f70b0cf)

- **[Expression DIY of L1TD1 gene in COAD](http://gepia2.cancer-pku.cn/#analysis) matching with TCGA normal data and GTEx data**

![L1TD1_profile_oSvkq](https://github.com/yerimeeei/Bioinformatics/assets/134043926/e27d8b3b-69cc-4cd6-a9fd-a0890a328f47)

- **Question: what is difference dataset in TCGA normal data and GTEx data?**

## **Literature Review**

### **Abstract**

- **Background**
    - prognostic markers (in this case L1TD1) to a particular cancer type (colorectal cancer) can assist in the evaluation of survival probability of patients → help clinicians to access the available treatment modalities.
- **Methods**
    - Data Analysis from 3 independent colon cancer microarray gene expression data sets (N=1052).
    - Survival Analysis: performed for 3 data sets, arranged by the expression level of L1TD1.
    - Correlation Analysis: performed to investigate the role of the interactome of L1TD1 in colon cancer patients.
- **Results**
    - L1TD1: novel positive prognostic marker for colon cancer; Increased expression of L1TD1 associated with longer disease-free survival in all 3 datasets.
    - Contrast to previous study on medulloblstoma (cancerous brain tumor), which L1TD1 linked with poor prognosis.
        - medulloblastoma L1TD1 was co-expressed with its interaction partners.
    - Analysis in colon cancer in this study revealed lack of co-expression of L1TD1 with its interaction partners.
- **Conclusions**
    - Increased expression of L1TD1 as a prognostic marker predicting longer disease-free survival in colon cancer patients.

### Background

- **L1TD1**: RNA-binding protein required for self-renewal of undifferentiated embryonic stem cells.
    - L1TD1 protein: shown to form core interaction network with the canonical pluripotency factors OCT4, NANOG, LIN28, and SOX2 in human embryonic stem cells (hESCs) → **L1TD1 depletion** resulted in downregulation of the pulripotenxy markers *OCT4, NANOG, and LIN28* in hESCs.
    - L1TD1 previously shown to be essential for self-renewal of embryonal carcinoma cells and to support growth of seminoma cells.
- Study on **L1TD1 immunoexpression** in colon adenocarcinoma tissue to assess the prognostic significance of L1TD1 in colon cancer.

### Methods

1. **Microarray data sets**
    - Raw microarray data sets from Gene Expression Omnibus (GEO)
        - Total of 1052 samples, but 124 samples excluded.
        - Addition of 2 seminoma and 1 stem cell gene expression microarray data sets were analyzed.
        - Stem Cell Data Sets composed of samples from 10 hESCs, 49 induced iPSCs, 5 cancer cell lines, and 6 non-cancerous somatic cell lines.
2. Gene expression analysis
    - Normalize CEL file containing the probe intensity measurements of the Affymetrix probes.
        - using Universal exPression Code (UPC) normalization method from the Bioconductor package “SCAN.UPC”
        - using Robust Multiarray Average (RMA) normalization method from the Bioconductor package “affy”
    - UPC normalization method
        - provides score, which represents the probability of a particular gene being expressed in a particular sample, was used to categorize the samples in all data sets based on their L1TD1 expression status.
        - UPC threshold of 0.6 determined
    - RMA
        - provides normalized $log_2$ intensity values.
        - RMA normalized gene expression values were used to calculate pairwise correlations between genes.
        - to correct for multiple testing
            - false discovery rate (FDR) was controlled using the Benjamini-Hochberg procedure.
    - primary probe: probe “219955_at” was chosen for quntification of L1TD1
        - as it was present in both of the Affymetrix platforms in this study.
3. Gene list descriptions
