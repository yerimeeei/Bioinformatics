# Day 5 L1TD1 gene expression for colorectal cancer

# TODO Today

- [x]  Using R to plot N and M stages of combined COAD + READ.
- [x]  4966 briefing session 11am-12pm
- [ ]  Read literature of L1 cancer basic and make questions to ask.
- [ ]  find some datasets from GEO or study
- [ ]  Study STAR aligner, FASTQ, BAM and run STAR, featureCounts on one FASTQ sample

## Links

- Batch Effect ([using combat in R](https://www.rdocumentation.org/packages/sva/versions/3.20.0/topics/ComBat))
- [GEO data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164541) that has normal, adenoma, and primary tumor
- Installing [SRA](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) Toolkit
- Study [STAR](https://github.com/alexdobin/STAR) aligner (STAR → BAM → count ? (fc_output in R))
    - [STAR manual](https://gensoft.pasteur.fr/docs/STAR/2.7.3a/STARmanual.pdf)
- **FASTA, FASTQ, BAM, SAM [**file types](https://hhj6212.github.io/biology/tech/2020/08/26/Bioinformatics-fileformats.html)
- [Youtube](https://www.youtube.com/watch?v=WLvST3OvAJ0) Video

# Summary of N and M staging using R

1. **COAD + READ pathological N staging**

<img width="918" alt="COAD__READ_pathological_N_staging" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/77e034a7-9c80-4096-bd22-dbbbdd25cc90">

- COAD cohort

<img width="844" alt="COAD_pathological_N_staging" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/7017e8ce-33a0-45d3-8ff1-8cd7f1d13e06">

- value distribution for each COAD and READ

<img width="885" alt="dis_N1_COAD" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/8271cc76-3ad3-4854-b1f3-252892e3885c">

<img width="885" alt="dis_N1_READ" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/863bd913-0208-4f16-93c6-0cc522c38f2b">

2. **COAD + READ pathological M staging**

<img width="900" alt="COAD__READ_pathological_M_staging" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/994f0b45-17f2-4105-9c00-67a927db2af5">

- COAD cohort

<img width="967" alt="COAD_pathological_M_staging" src="https://github.com/yerimeeei/Bioinformatics/assets/134043926/f2aa7cec-861d-47e4-99d7-22e37c43e5f6">

# Study STAR aligner, FASTQ, BAM and run STAR, featureCounts on one FASTQ sample

- **FASTA, FASTQ, BAM, SAM [**file types](https://hhj6212.github.io/biology/tech/2020/08/26/Bioinformatics-fileformats.html)

| FASTA | 특정 분자의 서열을 나타내는 데에 사용됩니다.주로 Genome 의 각 chromosome 마다의 서열을 저장하는 데 쓰이죠.혹은 gene/transcript/protein 각각의 염기서열 및 아미노산 서열을 저장할 때 사용되는 형식 |
| --- | | --- |
| FASTQ | NGS (Next generation sequencing data) 의 결과를 저장.NGS 실험을 진행하면, 그 결과로 cDNA library 서열을 읽어서 데이터로 얻을 수 있습니다. 즉, 각 cDNA library 의 염기서열을 알 수 있는거죠. 이 서열 하나를 ‘read’ 라고 하는데요. fastq 파일은 여러 read 의 정보를 한 파일에 저장. 보통 이 파일을 reference genome 에 align 한 뒤에 활용. |
|  | 1. sequence ID: @ 로 시작하며, 해당 서열의 이름을 나타냅니다. 2. sequence: 실제로 읽은 염기서열 정보입니다. 3. description: ‘+’ 글자로 시작하는데, + 하나만 있기도 하고 sequence ID를 넣거나 설명 포함 4. quality: 각 염기서열이 얼마나 정확히 읽혔는지를 나타냄. Phred quality score 라는 표현법을 사용. |
| BAM | binary alignment map 이라는 형식. 이 파일은 위에서 설명한 fastq 파일을 reference genome에 align 했을 때 만들어지는 파일. 즉, 각 cDNA library 조각이 reference genome 의 어느부분에서 나왔구나~ 하는 정보를 담았다는 거죠. fastq 에서는 각 read 의 염기서열과 그 품질을 알 수 있다면, BAM 파일은 염기서열과 reference 에서의 위치정보 를 알 수 있습니다. |
| SAM | BAM 파일을 ‘볼 수 있게’ 해놓은 파일이 있습니다. 바로 SAM 파일입니다.BAM 파일과 SAM 파일은 동일한 정보를 가지고 있고, 서로 변환이 가능.  |
|  | SAM 파일은 header 부분과 alignment 부분으로 이루어져 있습니다.- header: 파일에 대한 설명을 주는 부분입니다. @ 로 시작하는 라인들입니다.- alignment: 각 read에 대한 alignment 정보를 제공하는 부분입니다. 필수적인 11개의 컬럼으로 이루어져 있고, 추가로 몇 개의 컬럼이 더 있을 수도 있음. |
| SAM alignment column | 1. QNAME: read 이름 2. FLAG: 2진수로 된 read alignment 에 대한 설명 3. RNAME: refernce sequence 의 이름 4. POS: reference sequence 에서 align 된 위치 5. MAPQ: mapping quality. 즉 얼마나 정확히 align 되었는지. 6. CIGAR string: alignment 정보를 표현한 문자열. Match, Gap 등의 설명을 각 염기마다 표현. 7. RNEXT: 다음 read 의 reference sequence 이름. 주로 paired end read 에 대한 분석을 위해 사용. 8. PNEXT: 다음 read 의 align 된 위치. 주로 paired end read 에 대한 분석을 위해 사용. 9. TLEN: Template length. paired-end read 둘의 left-end 부터 right-end 까지의 길이. 10. SEQ: segment sequence. 염기 서열을 나타냄. 11. QUAL: Phread quality score |
- [Youtube](https://www.youtube.com/watch?v=WLvST3OvAJ0) Video
- Installing [SRA](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) Toolkit
- Study [STAR](https://github.com/alexdobin/STAR) aligner (STAR → BAM → count ? (fc_output in R))
    - [STAR manual](https://gensoft.pasteur.fr/docs/STAR/2.7.3a/STARmanual.pdf)

## After Meeting 3:30pm

- [x]  send COAD violin plot except blank data to Dr.JW → indicate there’s no problem
- [x]  Make T, N, M violin plots for CRC combining a,b staging
    - T, N, M 각각 T1(T1a, T1b) 다 합치기 (CRC (COAD+READ) 합쳐서 그래프 완성)
- [ ]  L1TD1 gene → for further studies, need to choose whether I want to do for other lines or only L1 → Study and Read about LINE
- [ ]  Install STAR, study FASTQ, BAM, featureCounts in Terminal

# Violin plots of L1TD1 expression in CRC

1. T Staging

![T_Staging](https://github.com/yerimeeei/Bioinformatics/assets/134043926/f927d973-0f8e-4212-98cb-19f3e0f14205)

2. N Staging

![N_Staging](https://github.com/yerimeeei/Bioinformatics/assets/134043926/a138f388-cb6d-4415-a7c2-845859699f86)

3. M Staging

![M_Staging](https://github.com/yerimeeei/Bioinformatics/assets/134043926/ba512d37-13aa-4f06-b7dd-5c4a0d8ffbc5)

### Duty for Next week

- Read literature of L1 cancer basic and make questions to ask.
- find some datas from GEO or study
- L1TD1 gene → for further studies, need to choose whether I want to do for other lines or only L1 → Study and Read about LINE
- Install STAR, study FASTQ, BAM, featureCounts in Terminal
