# Day 3 L1TD1 gene in TCGA colorectal samples

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

### Duty for today

- [x]  Finish reading L1TD1 literature and make questions.
- [ ]  Read line-1 in cancer basic paper to understand and make questions.
- [x]  Try using R to match phenotype data with HTSeq and make boxplot or violin plot.
    - study R commands used for the code.
- [x]  Ask question to Sandy & Dr. JW
    - regarding L1TD1 gene expression
    - regarding using the R (trial & progress).

- Read the basic paper on LINE-1 in cancer to better understand the role of LINE-1 in cancer progression.
- Ask questions about L1TD1 gene expression and using R for data analysis.
- **Regarding the co-expression of L1TD1 with its interaction partners: would it also be important to investigate their expression patterns in different types of cancer?**

# Literature 2: L1TD1 - a prognostic marker for colon cancer

## **Literature Review**

### Discussion

- Co-expression of L1TD1 with its currently known interaction partners was investigated to understand distinctive role of L1TD1 in different cancers.
    - outcome that L1TD1 was **not co-expressed with its interaction partners** in **colon cancer** points to potential participation of L1TD1’s interaction partners in the **contrasting prognostic outcome**.
    - Question
        - The study found that high expression of L1TD1 gene linked to higher survival for colon cancer patients. Is the interaction partners previously known to poor prognosis in other type of cancers?
            - In **liver cancer,** L1TD1 has been found to interact with the tumor suppressor protein p53, and its overexpression has been associated with poor prognosis and **reduced overall survival**.
            - In other types of cancer, such as **lung cancer and ovarian cancer**, L1TD1 expression has been correlated with **better prognosis**, indicating that the role of L1TD1 in cancer progression may be context-dependent and may vary depending on the specific cancer type.
        - **Why L1TD1 in cancer progression is dependent and express different progression in different specific types of cancer?**
        
        ***** (thus co-expression of L1TD1 with its interaction partners might be investigated?)**
        
    - Finding suggest: co-expression of L1TD1 with its interaction partners might be **required for clearly showing an aggressive and detrimental phenotype**.
- **Co-expressed Genes in Colon Cancer**

| Gene | Function |
| --- | --- |
| CLCA1 | Tumor suppressor protein that regulates differentiation and proliferation. |
|  | Low expression Associated with tumorigenesis, metastasis, and chromosomal instability, as well as poor prognosis. |
| KLF4 | Target of the tumor suppressor gene APC. Overexpression reduces cell migration and invasion. |
| GMDS | Shown to have exon deletions linked to progression of COAD and READ. |
| MUC2 | High expression associated with longer DFS. |
| PCCA | Frameshift mutation result in premature protein synthesis termination of it (also in gastric cancer). |
| SERPINA | controversial result |
- Several of co-expressed genes link to other cancers

| Gene | Function | Cancer |
| --- | --- | --- |
| FCGBP | Downregulation is associated with decreased overall survival. | Gallbladder adenocarcinoma |
|  | Progression | Prostate cancer (e.g. https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations) |
| ST6GALNAC1 | Upregulation associated with good prognosis. | breast cancer |
|  | siRNA-mediated silencing shown to lead to reduced growth, migration, and invasion. | Gastric cancer cells in vitro |
| KIAA1324, LINC00261, ITLN1 | Tumor suppressors with decreased expression associated with poor prognosis. | Gastric cancer |
| ASRGL1 | Low expression as a marker for poor prognosis. | endometrial carcinoma |
| SLC27A2 | Reduced levels associated with poor survival. | lung cancer |
| SLIT and SLITRK6 | Bladder tumor antigen that are under investigation as a target for antibody-drug conjugate therapy. |  |
| HEPACAM2 | paralog of HEPACAM, Act as a tumor suppressor by promoting differentiation. (not well-studied yet) |  |

### Conclusion

- This study: L1TD1 as a marker for good prognosis in colon cancer.
- Further study: need for validation of L1TD1 as a potential prognostic marker in larger cohorts of colon cancer.
- Also emphasize potential merits of investigating co-expressed genes to markers of interest.

# L1TD1 gene expression Raw data Analysis

## Analysis Using R ([COAD](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations))

- [x]  Extract expression for L1TD1 (ENSG00000240563.1)
- [x]  Divide samples into normal and tumour. Then for tumour divide into pathological Stage
- [x]  Make boxplot or violin plot comparing expression of L1TD1 for normal, and tumour stages
- [x]  T-test between Normal and Primary Tumor samples
- [ ]  Graph t-test results
1. **Extract expression for L1TD1**
2. **Extract expression for “A”**
3. **Divide samples into normal (11A) and tumor.**
4. **Divide tumor samples into pathological stages: primary tumor (01A), metastasis (06A) and recurrent tumor (02A).**
- **Question**
    - **Can I make boxplot or violin plot from multiple different tables in R?**
        - One way to do this is by combining the data from the different tables into a single data frame or matrix. → then use the `ggplot2` package to create boxplots or violin plots using the combined data.
    - **How to match phenotype data with HTSeq - FPKM-UQ?**
        - I tried to merge two tables, phenotype table that are divided into tumor stagesand values from HTSeq - FPKM-UQ table, but it didn’t work. (hmm..)

→ Try another method.

- divide by tumor stages from HTSeq-FPKM-UQ table. → not worked.

**Solution**

- using melt in reshape library to sort data table

```r
## using melt in reshape library to sort data table
library(reshape)
L1TD1_exp_df <- melt(L1TD1_exp)
```

<img width="875" alt="day3_r_progress_1" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/fe268d38-4d42-4bf5-ba5c-a45377970385">

<img width="875" alt="day3_r_progress_2" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/fea2be35-cd6f-412e-b486-e8401e212ba6">

- merge table by "submitter_id.samples", prior to this step, change column names in L1TD1_exp_df

```r
colnames(L1TD1_exp_df)[2] <- "submitter_id.samples"
merged_L1TD1 <- merge(L1TD1_exp_df, L1TD1_get, by = "submitter_id.samples")
```

5. **Make boxplot or violin plot comparing expression of L1TD1 for normal, and tumour stages.**

```r
library(ggplot2)

# Basic violin plot
L1TD1_vplot <- ggplot(merged_L1TD1, aes(x=extract, y=value)) 
							+ geom_violin()

# Set trim argument to FALSE to don't trim the tail
L1TD1_vplot <- ggplot(merged_L1TD1, aes(x=extract, y=value)) 
							+ geom_violin(trim=FALSE)
# violin plot with mean points
L1TD1_vplot + stat_summary(fun.y = mean, geom = "point", shape=23, size=2)
# violin plot with median points
L1TD1_vplot + stat_summary(fun.y=median, geom="point", size=2, color="red")
# violin plot with median and quartile
L1TD1_vplot + geom_boxplot(width=0.1)
# violin plot with mean and standard deviation
L1TD1_vplot + stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.2 )
L1TD1_vplot + stat_summary(fun.data=mean_sdl, geom="pointrange", color="red")
# Change violin plot line colors by groups
ggplot(merged_L1TD1, aes(x=extract, y=value, color=extract)) 
+ geom_violin(trim=FALSE)

# FINAL violin plot
L1TD1_vplot <- ggplot(merged_L1TD1, aes(x=extract, y=value, group=extract, 
						fill=extract)) + geom_violin(trim=FALSE, fill="gray") 
						+ labs(title="L1TD1 expression in COAD normal vs tumor",
						x="Sample Types", y = "HTseq Value") 
						+ geom_boxplot(width=0.1)+ theme_classic()
```

- [custom ggplot2 colors](http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually)

### Violin plot with boxplot comparing expression of L1TD1 in COAD for normal, and tumour stages

<img width="560" alt="Violin_plot_with_boxplot_comparing_expression_of_L1TD1_in_COAD_for_normal_and_tumour_stages" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/c4064cb0-64e6-433d-bd44-82e66a9a474e">

### T-test between normal and tumor samples using R (COAD)

```r
## T-test between normal and tumor samples
t.test(merged_L1TD1$value[merged_L1TD1$extract == "01A"], 
			 merged_L1TD1$value[merged_L1TD1$extract == "11A"])

			Welch Two Sample t-test

data:  merged_L1TD1$value[merged_L1TD1$extract == "01A"] and 
			 merged_L1TD1$value[merged_L1TD1$extract == "11A"]
t = -4.1385, df = 72.514, p-value = 9.311e-05
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval: 
-2.0390113 -0.7133698
sample estimates:
mean of x mean of y  
14.60799  15.98418
```

- Interpret result of t-test
    - **Welch Two Sample t-test**: This line indicates the type of t-test that was performed.
    - **data**: This line indicates the names of the two samples that were compared in the t-test.
    - **t**: This line shows the value of the t-statistic for the test. The larger the absolute value of t, the more evidence there is for a difference between the two sample means.
    - **df**: This line shows the degrees of freedom for the t-test. This value is used to calculate the p-value.
    - **p-value**: This line shows the probability of observing a difference as large as the one observed in the sample, assuming that there is no difference between the population means. If the p-value is less than the significance level(usually 0.05), then we reject the null hypothesis and conclude that there is evidence for a difference between the population means.
    - **alternative hypothesis**: This line shows the alternative hypothesis for the t-test. In this case, the alternative hypothesis is that the true difference in means is not equal to zero.
    - **95 percent confidence interval**: This line shows the range of values within which the true difference in population means is likely to fall with 95% confidence. In this case, the true difference in means is likely to be between -2.0390113 and -0.7133698.
    - **sample estimates**: This line shows the mean of each sample. In this case, the mean of sample x is 14.60799, and the mean of sample y is 15.98418.
- **Interpret t-test data result of COAD**
    - t = -4.1385 (large absolute value of t, more evidence there is difference from normal (11A) and tumor (01A) sample means).
    - p-value = 9.311e-05 (less than 0.05 → reject null and conclude evidence for difference between normal and tumor sample means)
    - 95% confidence interval: true difference in means is likely to be between -2.0390113 and -0.7133698.
    - sample means: tumor(01A) mean: 14.60799, normal(11A) mean: 15.98418.

## Analysis Using R ([READ](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations))

```r
## import datasets of phenotype and HTseq_FPKM.UQ
## extract L1TD1 gene expression
L1TD1_exp_R <- TCGA.READ.htseq_fpkm.uq[TCGA.READ.htseq_fpkm.uq$Ensembl_ID 
								== "ENSG00000240563.1",]

## reshape HTseq tablelibrary(reshape)
L1TD1_exp_R_df <- melt(L1TD1_exp_R)

## change - to . of ID samples in phenotype table
TCGA.READ.GDC_phenotype$submitter_id.samples <- gsub("-",".",
							TCGA.READ.GDC_phenotype$submitter_id.samples, fixed=TRUE)

## Only extract sample types and extract tumor stages of only A samples
L1TD1_get_R <- TCGA.READ.GDC_phenotype[,c(1,114)]
L1TD1_get_R$extract <- substr(L1TD1_get_R$submitter_id.samples, 14,16)
L1TD1_get_R <- L1TD1_get_R[str_detect(L1TD1_get_R$extract, "A"),]

## change column names of HTseq table and merge with phenotype table
colnames(L1TD1_exp_R_df)[2] <- "submitter_id.samples"
merged_L1TD1_R <- merge(L1TD1_exp_R_df, L1TD1_get_R, by = "submitter_id.samples")

## violin plot with boxplot
library(ggplot2)
L1TD1_vplot_R <- ggplot(merged_L1TD1_R, aes(x=extract, y=value)) + geom_violin()
## Final violin plot
L1TD1_vplot_R <- ggplot(merged_L1TD1_R, aes(x=extract, y=value, group=extract, fill=extract)) 
					+ geom_violin(trim=FALSE, fill="gray") + labs(title=
					"L1TD1 expression in READ normal vs tumor",
					x="Sample Types", y = "HTseq Value") + geom_boxplot
					(width=0.1)+ theme_classic()
```

### [Violin plot](http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization) with boxplot comparing expression of L1TD1 in READ for normal, and tumour stages

<img width="560" alt="Violin_plot_with_boxplot_comparing_expression_of_L1TD1_in_READ_for_normal_and_tumour_stages" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/d6c3e3a2-7407-4a44-a496-3bb299ccc3f1">

- Question
    - I have made a violin plot of the expression of L1TD1 for the normal and primary tumors for both COAD and READ samples. But the there are only 1 or none of data for the metastatic sample and recurrent Tumor sample. So.. I am curious if I also have to look up for adenoma data!

### T-test between normal and tumor samples using R (READ)

```r
## T-test between normal and tumor samples
t.test(merged_L1TD1_R$value[merged_L1TD1_R$extract == "01A"], 
			 merged_L1TD1_R$value[merged_L1TD1_R$extract == "11A"])

			Welch Two Sample t-test
data:  merged_L1TD1_R$value[merged_L1TD1_R$extract == "01A"] and 
			 merged_L1TD1_R$value[merged_L1TD1_R$extract == "11A"]
t = -0.64425, df = 13.043, p-value = 0.5306
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval: 
-1.864651  1.007773
sample estimates:
mean of x mean of y  
13.76239  14.19083
```

- **Interpret [t-test](http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/) data result of READ**
    - t = -0.64425 (small absolute value of t, less evidence there is difference from normal (11A) and tumor (01A) sample means).
    - p-value = 0.5306 (higher than 0.05 → accept null and conclude no evidence for difference between normal and tumor sample means)
    - 95% confidence interval: true difference in means is likely to be between -1.864651 and 1.007773
    - sample means: tumor(01A) mean: 13.76239, normal(11A) mean: 14.19083.
- **Question: WHY statistically not proven that there is difference between normal and tumor samples?**

# Literature 3: LINE-1 in cancer basic

## **Literature Review**

### Abstract

- L1 retrotransposons: jumping genes that utilize a “copy-and-paste” mechanism to spread themselves throughout the genome.
    - by RNA intermediates, which process termed retrotransposition.
- L1s active in germ line and during embryogenesis. But yet **epigenetically suppressed in somatic cells.**
- In cancer cells: L1 abnormally activated → have role in genome instability (one of hallmarks of cancer).
    - methylation states and retrotransposition activities are associated with
    - fluctuate during cancer initiation and progression
    
    → representing promising diagnostic biomarkers and therapeutic targets.
    
- During Tumorigenesis (initial formation of tumor), L1s exert both retrotransposition-dependent and retrotransposition-independent functions.

| Function | Retroposition-dependent | Retroposition-independent |
| --- | --- | --- |
| Result | alterations in target gene expression or chromosomal rearrangement, or drive Alu and SVA (possible events in tumorigenesis) | potentially exert 1) epigenetic regulation by generating endo-siRNAs, forming chimeric L1 transcripts. |
|  |  | 2) changing the expression of adjacent genes by providing novel splicing sites or alternative promoters. |
- L1 encoded proteins, ORF1p and ORF2p, may have pro-oncogenic potential.
    - e.g. by activatign oncogenic transcriptional factors or sequestering onsocuppressors.
- In this study:
    - introduce components and mechanisms of L1 retrotransposition, discuss landscape, possible functions, and regulation of L1 activity in cancer.
    - seek their potential as diagnostic biomarkers and therapeutic targets.

### Duty for Tmr

- Read literature of L1 cancer basic in the morning and bring questions to ask.
- Meeting with Sandy & Dr. JW 10am.
- Using R to plot t-test graph.
