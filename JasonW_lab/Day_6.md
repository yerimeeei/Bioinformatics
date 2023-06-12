# Day 6 L1TD1 gene expression for colorectal cancer

### Duty for Today

- [ ]  Read literature of L1 cancer basic and make questions to ask.
- [ ]  L1TD1 gene → for further studies, need to choose whether I want to do for other lines or only L1 → Study and Read about LINE
- [x]  Study STAR, FASTQ, BAM, featureCounts
- [x]  log in to SERVER

## Links

- Batch Effect ([using combat in R](https://www.rdocumentation.org/packages/sva/versions/3.20.0/topics/ComBat))
- [GEO data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164541) that has normal, adenoma, and primary tumor
- Study [STAR](https://github.com/alexdobin/STAR) aligner (STAR → BAM → count ? (fc_output in R))
    - [STAR manual](https://gensoft.pasteur.fr/docs/STAR/2.7.3a/STARmanual.pdf)
- [Youtube](https://www.youtube.com/watch?v=WLvST3OvAJ0) Video

### Questions

- server에서 STAR을 어떻게 run하는지
- how to do **[COMPILING](https://github.com/alexdobin/STAR) FROM SOURCE of STAR? ([Error](https://cosmy.tistory.com/entry/Mac-OS에서-R-패키지-설치시-Unsupported-option-fopenmp-문제-해결법))**
- GEO data를 서버로 어떻게 다운로드 받는지

# STAR alinger

Spliced Transcripts Alignment to a Reference (STAR) is **a fast RNA-seq read mapper, with support for splice-junction and fusion read detection**. STAR aligns reads by finding the Maximal Mappable Prefix (MMP) hits between reads (or read pairs) and the genome, using a Suffix Array index.

### Mapping reads to the genome

In this step user supplies the genome files generated in the 1st step, as well as the RNA-seq reads (sequences) in the form of **FASTA or FASTQ files**. **STAR maps the reads to the genome, and writes several output files, such as alignments (SAM/BAM)**, mapping summary statistics, splice junctions, unmapped reads, signal (wiggle) tracks etc. Output files are described in Section 4. Output files. Mapping is controlled by a variety of input parameters (options) that are described in brief in Section 3. Running mapping jobs, and in more detail in Section 14. Description of all options.

**STAR command line has the following format:**

```bash
 STAR --option1-name option1-value(s)--option2-name option2-value(s) ...
```

If an option can accept multiple values, they are separated by spaces, and in a few cases - by commas.

**FeatureCounts**

- featureCounts is a tool to quantify RNA-seq and gDNA-seq data as counts. It is also suitable for single-cell RNA-seq (scRNA-seq) data. It supports multi-threading. The featureCounts is part of the Subread package

## **Introduction to [RNA-Seq using high-performance computing](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html)**

**Read Alinger**

![Untitled](https://github.com/yerimeeei/Bioinformatics/assets/134043926/ca16714b-e219-4db8-910a-5131c33936ef)

### STAR Alinger

To determine where on the human genome our reads originated from, we will align our reads to the reference genome using [STAR](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) (Spliced Transcripts Alignment to a Reference).

**STRATEGY: two-step process 1) Seed searching 2) Clustering, stitching, and scoring**

1. **Seed searching**
- STAR will search for the longest sequence that exactly matches one or more locations on the reference genome. These longest matching sequences are called the **Maximal Mappable Prefixes (MMPs):**

![Untitled 1](https://github.com/yerimeeei/Bioinformatics/assets/134043926/0110af45-78da-40c5-92b6-2a48bf91a8bd)

- The **different parts of the read that are mapped separately are called ‘seeds’.** So the first MMP that is mapped to the genome is called *seed1*. STAR will then search again for only the unmapped portion of the read to find the next longest sequence that exactly matches the reference genome, or the next MMP, which will be *seed2*.

![Untitled 2](https://github.com/yerimeeei/Bioinformatics/assets/134043926/8e3d9562-1187-4cf5-aec8-d2926c1a2f1c)

- STAR uses an uncompressed suffix array (SA) to efficiently search for the MMPs, this allows for quick searching against even the largest reference genomes.
- **If STAR does not find an exact matching sequence** for each part of the read due to mismatches or indels, the previous MMPs will be extended.

![Untitled 3](https://github.com/yerimeeei/Bioinformatics/assets/134043926/e974eeae-7da9-4a76-9609-be2faf1d4009)

- **If extension does not give a good alignment**, then the poor quality or adapter sequence (or other contaminating sequence) will be soft clipped.

![Untitled 4](https://github.com/yerimeeei/Bioinformatics/assets/134043926/5c1b08aa-dcf6-4944-95ba-9df18d9180f0)

2. **Clustering, stitching, and scoring**

The separate seeds are stitched together to create a complete read by first clustering the seeds together based on proximity to a set of ‘anchor’ seeds, or seeds that are not multi-mapping.

Then the seeds are stitched together based on the best alignment for the read (scoring based on mismatches, indels, gaps, etc.).

![Untitled 5](https://github.com/yerimeeei/Bioinformatics/assets/134043926/fbe05eab-1dc6-42a6-8abe-69983e91147f)

**RUNNING STAR**

# STAR Workflow

## 1. Generating Genome Indexes

In this step user supplied the reference genome sequences (FASTA files) and annota- tions (GTF file), from which STAR generate genome indexes that are utilized in the 2nd (mapping) step. The genome indexes are saved to disk and need only be generated once for each genome/annotation combination. 

- path to genome FASTA file (Genome Directory)

```bash
/storage/jwlab/gencode.v43.chr_patch_hapl_scaff.basic.annotation.gtfsandy/kallisto_analysis/kallisto_index/GRCh38.primary_assembly.genome.fa
```

- path to annotations.gtf file

```bash
/storage2/jwlab/sandy/data/gencode.v43.chr_patch_hapl_scaff.basic.annotation.gtf
```

- **To use STAR command**

```bash
/storage2/jwlab/yerim/STAR
```

- Genome Index

```bash
/storage2/jwlab/yerim/STAR --runThreadN 5 --runMode genomeGenerate 
	--genomeDir genome_index 
	--genomeFastaFiles /storage/jwlab/sandy/kallisto_analysis/kallisto_index/GRCh38.primary_assembly.genome.fa 
	--sjdbGTFfile /storage2/jwlab/sandy/data/gencode.v43.chr_patch_hap1_scaff.basic.annotation.gtf 
	--sjdbOverhang 99
```

- screen to keep running job

To use the `screen` command, follow these steps:

1. Open a terminal window.
2. Type `screen` and press Enter to start a new screen session.
3. Run the command or program you want to keep running in the background.
4. **To detach from the screen session and return to your normal terminal, press the key combination `Ctrl-a` followed by `d`.**
5. You can now close your terminal window or log out of your SSH session without stopping the program running in the screen session.
6. **When you want to resume the screen session, open a new terminal window and type `screen -r` and press Enter. This will reattach you to your previous screen session and you can continue working on the program you left running.**

Note: If you have multiple screen sessions running, you can use `screen -ls` to list them and `screen -r <session-id>` to reattach to a specific session.

## 2. Mapping reads to the genome
