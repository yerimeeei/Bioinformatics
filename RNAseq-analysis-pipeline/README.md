# RNAseq-analysis-pipeline

# bioinformatic-basic-practices

## Cancer Genomics related

You can search and get to know some background knowledge before starting this practices,

1. Try to understand and explain the following concepts:
- Central Dogma
- Exon: sections of DNA that code for proteins
- Intron: noncoding sections of an RNA transcript.
- Gene
- Gene expression (it is a relative value or absolute value?)
    - Relative value
        - determines fold changes in expression between two samples.
        - measures the relative change in mRNA expression lebels.
    - Absolute value
        - determines expression levels in absolute numbers of copies.
        - relates the PCR signal to input copy number using a calibration curve.
- SNP (Single Nucleotide Polymorphisms)
- germline mutation & somatic mutation
    - germline mutation
        - changes to DNA that you inherit from the egg and sperm cells during conception.
        - occur in gametes and can be passed onto offspring (every cell in the entire organism will be affected).
    - somatic mutation
        - changes to DNA that happen after conception to cells other than the egg and sperm.
        - occur in a single body cell and cannot be inherited (only tissues derived from mutated cell are affected).

2. Explain the principle of high throughput sequencing（like illumina) in your own words.
    - High throughput sequencing (Next-generation sequencing)
        - allows for parallel sequencing of fragments of DNA.
        - Steps: Illumina sequencing process begins with the preparation of a DNA library, which involves fragmenting the DNA into small pieces, adding adapters to the ends of the fragments, and amplifying them by PCR. The resulting library contains millions of DNA fragments, each with a unique adapter sequence. Library is then loaded onto an Illumina sequencer, which uses a process called bridge amplification to create clusters of DNA fragments on a flow cell. Each cluster contains multiple copies of the same DNA fragment, and each cluster is spatially separated from the other clusters on the flow cell. During sequencing, fluorescently labeled nucleotides are added to the flow cell, and the nucleotides are incorporated into the growing DNA strands by DNA polymerase. As each nucleotide is added, a camera takes a picture of the flow cell to capture the fluorescent signal, and the identity of the nucleotide is determined based on the color of the signal.
        - e.g. RNA sequencing
            - provides relatively higher coverage and higher throughput and also generates additional data that can help discover novel transcipts, understand allele specific information, and identify alternatively spliced genes.
            - Steps: isolation of RNA of interest → conversion of RNA to cDNA (ensure molecule’s stability, easy handling) → adapters attached to DNA fragments to enable sequencing → sequencing data aligned and used to generate corresponding RNA sequence map.
            
3. What is the format for fastq file and gtf file?
    - FASTQ file
        - text-based file format used for storing both a biological sequence (usually nucleotide sequence) and its corresponding quality scores. It is a commonly used format for storing the output of high-throughput sequencing technologies such as Illumina, Ion Torrent, and PacBio. The format consists of four lines per sequence, with each line starting with a specific character:
            1. Line 1 starts with the "@" character and contains a sequence identifier and an optional description.
            2. Line 2 contains the actual nucleotide sequence.
            3. Line 3 starts with the "+" character and may contain an optional description (which is often identical to the description in line 1).
            4. Line 4 contains the quality scores for the nucleotide sequence in line 2.
        - example of a FASTQ entry:
        
        ```bash
        @SEQ_ID
        GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
        +
        !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
        ```
        
        - Quality Score
            - measure of the confidence that a particular nucleotide was sequenced correctly. It represents the quality of the base call for each corresponding nucleotide in line 2. The quality score is typically represented as a single ASCII character, with each character corresponding to a Phred score, which is a logarithmic measure of the probability that a nucleotide is incorrect. The Phred score is calculated using the formula:
            
            $$
            Q = -10 * log(P)
            $$
            
            - Q: Phred score
            - P: probability that the nucleotide is incorrect
            - [ASCII character Illumina](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm)
    - GTF file
        - tab-delimited file format that is used to describe gene and transcript models in genome annotation. It contains information about gene features such as exons, introns, and UTRs, as well as information about the gene's location on the genome, its name, and other attributes. The GTF format consists of one line per feature, with the following columns:
            1. Chromosome name
            2. Source of the feature (e.g., "ensembl", "refseq")
            3. Feature type (e.g., "gene", "exon", "CDS")
            4. Start position of the feature on the chromosome
            5. End position of the feature on the chromosome
            6. Score (usually set to ".")
            7. Strand (+ or -)
            8. Frame (0, 1, or 2 for CDS features, or "." for non-CDS features)
            9. Attributes (a semicolon-separated list of key-value pairs that describe additional information about the feature, such as gene ID, transcript ID, and gene name)
        - example of GTF file

```bash
chr1	ENSEMBL	exon	11869	12227	.	+	.	
gene_id "ENSG00000223972"; transcript_id "ENST00000456328"; 
exon_number "1"; gene_name "DDX11L1"; 
gene_biotype "transcribed_unprocessed_pseudogene"; 
transcript_name "DDX11L1-202"; exon_id "ENSE00002234944";
```

4. Introduce the TCGA database in your own words.
    - TCGA: The Cancer Genome Atlas
        - publicly available database that contains a vast amount of molecular and clinical data on various types of cancer

## Do analysis with database and online tools

Part 0: TCGA,cBIoportal,GEPIA2 application

## Coding practices

Part I: LINUX pactices & simple project.

Part II: RNA-seq analysis pipeline

Part III: R + Differential expression analysis basic plots
