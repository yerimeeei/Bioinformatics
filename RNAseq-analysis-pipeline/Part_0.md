# Part 0: TCGA,cBIoportal,GEPIA2 application

## The Cancer Genome Atlas Program (TCGA)

[https://portal.gdc.cancer.gov/](https://portal.gdc.cancer.gov/)

> The Cancer Genome Atlas (TCGA), a landmark cancer genomics program, molecularly characterized over 20,000 primary cancer and matched normal samples spanning 33 cancer types. This joint effort between the National Cancer Institute and the National Human Genome Research Institute began in 2006, bringing together researchers from diverse disciplines and multiple institutions.Over the next dozen years, TCGA generated over 2.5 petabytes of genomic, epigenomic, transcriptomic, and proteomic data. The data, which has already lead to improvements in our ability to diagnose, treat, and prevent cancer, will remain publicly available for anyone in the research community to use.
> 

![Untitled](https://github.com/yerimeeei/Bioinformatics/assets/134043926/73cf6a82-2337-4bb4-93d7-f89bbd3446fe)

![Untitled 1](https://github.com/yerimeeei/Bioinformatics/assets/134043926/4f86d522-664f-4a3b-bf04-2501ffad1f82)

---

Q1: How many cases are involved in TCGA database? How many projects?

**86,962 cases, 78 projects**

Q2: How many primary site involved in TCGA database, why the cases size are different amoung various primary sites?

**68 primary sites involved. The availability and quality of the samples can vary between different cancer types. Also, the amount of funding and resources available for each primary site may have varied, which could influence the number of cases that were included in the study.**

Q3: what are the difference between projects,cases and files?  (click "Repository" to see more)

- **TCGA projects refer to specific cancer types that have been studied.**
- **TCGA cases refer to individual patients within those cancer types.**
- **TCGA files refer to specific data files associated with those patients.**

---

**Task1: Find projects about breast cancer RNA-seq cohort.**

"project" —— select "breast" in *Primary Site* and "RNA-seq" in *Experimental Strategy*

In the end there are 6 projects (1700+ cases). ******************************Now 9 projects (2023-06).******************************

![Untitled 2](https://github.com/yerimeeei/Bioinformatics/assets/134043926/47ac5a37-5f97-4266-8053-ca44e818aa2b)

**Task2: Find/Download a WGS mutation information (vcf files) from brain cancer case(age from 20 to 40).**

"Repository" —— "Cases" —— select "brain" in *Primary Site* —— set age in *Age at Diagnosis*

then "Files" —— select "VCF" in *Data Format*，click one samples into detail pages and click "Download".

![Untitled 3](https://github.com/yerimeeei/Bioinformatics/assets/134043926/5dea971c-1d7d-4f9f-9452-bf9c76be5f23)

![Untitled 4](https://github.com/yerimeeei/Bioinformatics/assets/134043926/73194c68-f93c-4a8d-bad0-9f5aae4d3643)

TCGA as a cancer samples data base, you can easily get a opened and small file with "Download" button.

For large file and batch download, you can "Add to Cart" first, then download Manifest file and use [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) in server.

![Untitled 5](https://github.com/yerimeeei/Bioinformatics/assets/134043926/867a7213-ccf1-4ec3-89b8-66159975ff2d)

Since most RNA-seq data in TCGA are controlled, if you want to do basic the gene analysis among TCGA samples, online tools GEPIA2 can help you!

## Gene Expression Profiling Interactive Analysis(GEPIA2)

[http://gepia2.cancer-pku.cn/#index](http://gepia2.cancer-pku.cn/#index)

![Untitled 6](https://github.com/yerimeeei/Bioinformatics/assets/134043926/dbe2624e-3ada-4532-8133-7336893299e0)

**Task1: Find differential expressed genes in ACC.**

"Expression Analysis" —— "Differential Genes" —— select Dataset(cancer name) —— set Log2FC and qvalue cutoff (use default here) —— "List", then you will get the DE genes in ACC.

![Untitled 7](https://github.com/yerimeeei/Bioinformatics/assets/134043926/f394f431-e299-474e-94dd-49fcf998fe82)

**Task2: Check gene IGF2 expression in tumor samples and normal samples in ACC and BLCA.**

"Expression Analysis" —— "Expression DIY" —— "Box Plot" —— add Dataset(cancer name)

![Untitled 8](https://github.com/yerimeeei/Bioinformatics/assets/134043926/3af08627-b463-4ab8-875e-24891b74ddc9)

![Untitled 9](https://github.com/yerimeeei/Bioinformatics/assets/134043926/2f72573a-15a4-4e9e-8541-b0519be7ee29)

---

**Challenge1: Try to do a survival analysis, see if the IGF2 gene expression shows difference survival percentage in ACC.**

<img width="504" alt="Untitled 10" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/861b303c-3796-4faf-84b2-a00529c79942">

**Challenge2：GEPIA2 is a helpful tool in TCGA RNA-seq basic analysis. Could you please come up with some questions related gene expression in Cancer samples, and try to answer them with GEPIA2?**

- Q1: What is the correlation between the expression levels of two genes, GRB7 and EGFR, in a specific cancer type, BLCA Tumor?

<img width="505" alt="Untitled 11" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/b0229e7b-1e1b-4996-8f49-37e54c7609b7">

---

Above analysis mostly based on gene expression. For mutations in genes, tool cbioportal is helpful.

## cbioportal for cancer genome

[https://www.cbioportal.org/](https://www.cbioportal.org/)

![Untitled 12](https://github.com/yerimeeei/Bioinformatics/assets/134043926/38a5512b-8a48-4737-8dbf-e686e92014d7)

**Task1：check TP53 gene mutations happened in all cancer.**

["Quick search" —— select "TP53" in —— "Mutations", you can search keyword to selected cancer type.](https://www.cbioportal.org/results/mutations?case_set_id=all&gene_list=TP53&cancer_study_list=5c8a7d55e4b046111fee2296)

![Untitled 13](https://github.com/yerimeeei/Bioinformatics/assets/134043926/67e2fab9-be63-484d-b47c-47f07a0e76e8)

**Task2：What is the most common hotspot point mutation in KRAS in lung cancer?.**

"Query" —— select "Lung" —— select study *Lung Adenocarcinoma(TCGA,Nature 2014)* —— "Query By Gene"

[unselect "Putative CNA from GISTIC" —— type KRAS in *Enter Genes* —— "Submit Query"  —— "Mutations"  , G12C](https://www.cbioportal.org/results/mutations?cancer_study_list=luad_tcga_pub&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations%2Cstructural_variants&case_set_id=luad_tcga_pub_sequenced&gene_list=KRAS&geneset_list=%20&tab_index=tab_visualize&Action=Submit)

![Untitled 14](https://github.com/yerimeeei/Bioinformatics/assets/134043926/7ea97f3f-8cdf-40f0-9a4c-f07b30d26edf)

![Untitled 15](https://github.com/yerimeeei/Bioinformatics/assets/134043926/822f7bb1-605e-477e-8794-9a9d5a5113c8)

![Untitled 16](https://github.com/yerimeeei/Bioinformatics/assets/134043926/f4e036ef-0aab-4f81-8e94-f000c15dbf29)

---

Q1: What is the meaning of "Missense"?

**A missense mutation is a DNA change that results in different amino acids being encoded at a particular position in the resulting protein.**

Q2: How many samples are supporting this hostpot point mutation? how to get detail information of those samples?

******************************************************************************************************************************************There are 76 mutations. Click the sample ID to get datil information.******************************************************************************************************************************************

---

**Task3：use cbioportal to look at difference in survival for lung cancer patients with TP53 and/or KRAS mutations .**

"Query" —— select "Lung" —— select study *Lung Adenocarcinoma(TCGA,Nature 2014)* —— "Query By Gene"

[unselect "Putative CNA from GISTIC" —— type KRAS and TP53 in *Enter Genes* —— "Submit Query"  —— "Survival"  , Not statistical significant in survival analysis](https://www.cbioportal.org/results/comparison?cancer_study_list=luad_tcga_pub&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations%2Cstructural_variants&case_set_id=luad_tcga_pub_sequenced&gene_list=KRAS%2520TP53&geneset_list=%20&tab_index=tab_visualize&Action=Submit&comparison_subtab=survival)

![Untitled 17](https://github.com/yerimeeei/Bioinformatics/assets/134043926/fc9959f5-f2f9-41f0-b073-413dc27b9ec6)

**How about divided into 3 group "TP53 only", "TP53 + KRAS", and "KRAS only".**

"overlap" —— select "TP53(75)" and "KRAS(107)" group, we can see there are some overlap samples(25)

group samples：click "82" —— "Creat Group From Seleted Diagram Areas" —— "Submit" (name the other 2 groups in the same way)

[select 3 new groups only and click "Survival"](https://www.cbioportal.org/results/comparison?cancer_study_list=luad_tcga_pub&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations%2Cstructural_variants&case_set_id=luad_tcga_pub_sequenced&gene_list=KRAS%2520TP53&geneset_list=%20&tab_index=tab_visualize&Action=Submit&comparison_subtab=survival&comparison_selectedGroups=%5B%22TP53%20only%22%2C%22TP53%20%2B%20KRAS%22%2C%22KRAS%20only%22%5D&comparison_createdGroupsSessionId=647b65f19f013a0b001fa61a)

![Untitled 18](https://github.com/yerimeeei/Bioinformatics/assets/134043926/33d9e14c-f736-4134-93be-59e53a5c9f49)

![Untitled 19](https://github.com/yerimeeei/Bioinformatics/assets/134043926/d8e31c03-d793-444a-8a55-ed585ce2f440)

---

**[Challenge1: What is the most significantly higher expressed gene in ERCC2 mutant versus WT bladder cancer?](https://www.cbioportal.org/results/comparison?cancer_study_list=blca_tcga_pan_can_atlas_2018&Z_SCORE_THRESHOLD=2.0&RPPA_SCORE_THRESHOLD=2.0&profileFilter=mutations%2Cstructural_variants&case_set_id=blca_tcga_pan_can_atlas_2018_sequenced&gene_list=ERCC2&geneset_list=%20&tab_index=tab_visualize&Action=Submit&comparison_selectedGroups=%5B%22WT%22%2C%22ERCC2%22%2C%22Unaltered%20group%22%5D&comparison_createdGroupsSessionId=647b6ce3c7548976dbf00d7e&comparison_subtab=mrna)**

SGO1

---
