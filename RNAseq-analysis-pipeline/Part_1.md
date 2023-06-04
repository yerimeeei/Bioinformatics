# Part I: LINUX pactices & simple project.

## Basic linux practices

### Introduction

**UNIX** is an operating system which was first developed in the 1960s, and has been under constant development ever since.

**Linux** is an entire family of open-source Unix operating systems, that are based on the Linux Kernel. This includes all of the most popular Linux based systems like Ubuntu, Fedora, Mint, Debian, and others. More accurately, they’re called distributions or distros.

The **shell** acts as an interface between the user and the kernel. Most Linux distributions use a graphic user interface (GUI) as their shell, mainly to provide ease of use for their users，like virtual machine with linux system.  But here we highly recommend using a command-line interface (CLI) because it’s more powerful and effective. Tasks that require a multi-step process through GUI can be done in a matter of seconds by typing commands into the CLI.  The CLI can be easily get access to by Terminal in Mac and Putty/Xshell (install required) in Windows.

### Practices tutorial

- [UNIX Tutorial for Beginners](http://www.ee.surrey.ac.uk/Teaching/Unix/)
- UNIX Introduction
    - UNIX operating system: the kernel, the shell, and the programs.
        - The kernel: allocates time and memory to programs & handles the filestore and communications in reponse to system calls.
        - The shell: interface between the user and the kernel. command line interpreter (CLI) that interprets the commands the user types in and arranges for them to be carried out.
- Command
    
    
    | Command | Meaning |
    | --- | --- |
    | ls | list files and directories |
    | ls -a | list all files and directories |
    | mkdir | make a directory |
    | cd directory | change to named directory |
    | cd | change to home-directory |
    | cd ~ | change to home-directory |
    | cd .. | change to parent directory |
    | pwd | display the path of the current directory |
    | cp file1 file2 | copy file1 and call it file2 |
    | mv file1 file2 | move or rename file1 to file2 |
    | rm file | remove a file |
    | rmdir directory | remove a directory |
    | cat file | display a file |
    | less file | display a file a page at a time |
    | head file | display the first few lines of a file |
    | tail file | display the last few lines of a file |
    | grep ‘keyword’ file | search a file for keywords |
    | wc file | count number of lines/words/characters in file |
    | command > file | redirect standard output to a file |
    | command >> file | append standard output to a file |
    | command < file | redirect standard input from a file |
    | command1 | command 2 | pipe the output of command1 to the input of command2 |
    | cat file1 file2 > file0 | concatenate file1 and file2 to file0 |
    | sort | sort data |
    | who | list users currently logged in |
    | * | match any number of characters |
    | ? | match one character |
    | man command | read the online manual page for a command |
    | whatis command | brief description of a command |
    | apropos keyword | match commands with keyword in their man pages |
    | ls -lag | list access rights for all files |
    | chmod [options] file | change access rights for named file |
    | command & | run command in background |
    | ^C | kill the job running in the foreground |
    | ^Z | suspend the job running in the foreground |
    | bg | background the suspended job |
    | jobs | list current jobs |
    | fg %1 | foreground job number 1 |
    | kill %1 | kill job number 1 |
    | ps | list current processes |
    | kill 26152 | kill process number 26152 |

---

*Q1: Define relative path and absolute path*

- **Relative path: location of a file relative to the current (working) directory.**
- **Absolute path: location from the root directory.**

*Q2：What ERRORs did you meet and how with deal with those errors?*

******No ERRORs yet******

*Q3: if you are using remote server, how to transfer files between server and your local PC? (hint: scp for Mac, xftp for Windows)*

- **Using a SCP command on Mac**
    1. Open a terminal window on your local Mac computer.
    2. Type the following command: `scp [options] [source] [destination]`
    3. Replace `[options]` with any additional options you want to use, such as -`r` for recursive copying of directories.
    4. Replace `[source]` with the path to the file or directory you want to copy from the remote server.
    5. Replace `[destination]` with the path to the directory where you want to copy the files on your local computer.
    6. Press Enter to run the command.
    
    For example, if you want to copy a file named `example.txt` from the remote server at IP address 192.168.0.1 to your local Downloads folder, you would type the following command:
    
    `scp username@192.168.0.1:/path/to/example.txt ~/Downloads`
    

---

### Linux_practice

(provide hints and answer from the example file only)

1. Download a gtf file for plactise; or you can use any other gtf file you have.

```bash
wget <http://hgdownload.soe.ucsc.edu/goldenPath/felCat9/bigZips/genes/felCat9.refGene.gtf.gz>
gunzip -d felCat9.refGene.gtf.gz
```

2. How many rows in that gtf file?

```bash
wc -l felCat9.refGene.gtf
8487
```

3. How to extract the first 10 rows of the file? (try "head")

```bash
head felCat9.refGene.gtf
```

4. Save all "CDS" transcript information to a new file called "CDS.txt". (try "grep" and ">" )

```bash
grep -i CDS felCat9.refGene.gtf > CDS.txt
```

5. print out the first column (ChromosomeID) of file CDS.txt. (try "cut")

```bash
cut -f1 CDS.txt
```

6. search keyword "gene" in gtf file , why the result from "grep -c" and "grep -o ** | wc -l" are different?

```bash
grep -c "gene" felCat9.refGene.gtf
8487

grep -o "gene" felCat9.refGene.gtf | wc -l
16974
```

- **grep cases**
    - **-i: to ignore upper/lower case distinctions**
    - **-v: display those lines that do NOT match**
    - **-n: precede each matching line with the line number**
    - **-c: print only the total count of matched lines**
    - **-o: output only the matched parts of the pattern or expression, rather than the entire line that continas the match**
7. How many kind of transcript(column 3) listed in the gtf file? and what's the number size of each kind of transcript respectively? (combine "cut" "sort" and "uniq" together)

```bash
cut -f3 felCat9.refGene.gtf | sort | uniq | wc -l
7
cut -f3 felCat9.refGene.gtf | sort | uniq -c
274 3UTR
 302 5UTR
3269 CDS
3347 exon
 433 start_codon
 425 stop_codon
 437 transcript
```

8. use vim editor to add a header at the first line  of gtf file. ("vim" and "I" to get access to the editing mode."Esc" to the control mode,":wq" to save and quit)

```bash
vim felCat9.refGene.gtf
```

- i: enter insert mode
- Esc: exit insert mode
- :wq: save the changes and exit Vim
9. Jump to line 1000 of gtf file in vim editor

```bash
vim felCat9.refGene.gtf
:1000
```

10. replace all "refGene" with "CatRef". (try "vim" or "sed")

```bash
vim felCat9.refGene.gtf
:%s/refGene/CatRef/g
:wq
```

## ***Project: Analysis on a patient cohort.

### understand the data (biomed1 server ~/data/fastq/Genie_9.0/)

---

What information are involved in each column/row in every file?

How to define "patient" and "sample"?

How to merge sample and patient information?

How to search/locate a target SNP?

What is the difference between germline mutation and somatic mutation?

---

**Documents**

*introduction:data_guide.pdf;*

*Sample information: data_clinical_patient.txt; data_clinical_sample.txt;*    [clinical format](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#clinical-data)

*Mutation information: data_mutations_extended.txt;*  [MAF format](https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#mutation-data)

### Task 1: complete the following clinical form

summary the Primary & Metastasis tumor samples

| Cancer Type | total | female | male |
| --- | --- | --- | --- |
| LUAD | 12533 | 7545 | 4978 |
| IDC | 8026 | 7958 | 63 |
| COAD | 6268 | 2973 | 3289 |
| PRAD | 3781 |  |  |
| PAAD | 3573 |  |  |

---

points for discussion：

1. How to get all answers(from multi-cancer types) in one command line.
2. why the total number of female and male not equal to the total size? how to deal with those extra samples?

**challenge 1: Does LUAD patient has higher probability to get metastasis than other cancer type?**

### Task 1 Answer

```bash
##calculate the number of each cancer typ,remember to keep tumor samples only.  (SAMPLE_TYPE)
awk '{if ($5=="Metastasis"||$5=="Primary") print $0}' data_clinical_sample.txt  >tumor_sample_tmp.txt
cut -f 4 tumor_sample_tmp.txt|sort|uniq -c|sort -rnk1|head

##merge sample information with patient information
awk 'NR==FNR{a[$1]=$2}NR>FNR{print $1,$2,$3,$4,a[$1]}'  data_clinical_patient.txt tumor_sample_tmp.txt > Gender_tmp.txt
grep "LUAD" Gender_tmp.txt |cut -d " "  -f 5 |sort |uniq -c

##count the gender in one command line, CancernameList_tmp.txt file to list all CancerType that is gonna be analysis.
#cat CancernameList_tmp.txt
#LUAD
#IDC
#COAD
#PRAD
#PAAD
#head -1 data_clinical_patient.txt
#Patient Identifier	Sex	Primary Race	Ethnicity Category	Center
#head -1 Primary_tumor_tmp.txt
#PATIENT_ID	SAMPLE_ID	AGE_AT_SEQ_REPORT	ONCOTREE_CODE	SAMPLE_TYPE	SEQ_ASSAY_ID	CANCER_TYPE	CANCER_TYPE_DETAILESAMPLE_TYPE_DETAILED
for i in `cat CancernameList_tmp.txt`;do echo "$i"; awk 'NR==FNR{a[$1]=$2}NR>FNR{print $1,$2,$3,$4,a[$1]}'  data_clinical_patient.txt tumor_sample_tmp.txt |grep "$i\\s" |cut -d " " -f 5|sort|uniq -c;done
```

---

### Task 2: SNP analysis on paired samples

ref: [https://ascopubs-org.eproxy.lib.hku.hk/doi/full/10.1200/PO.19.00394](https://ascopubs-org.eproxy.lib.hku.hk/doi/full/10.1200/PO.19.00394)

- What information can be used to locate one SNP?
- Find paired sample (Primary and Metastasis sample from the same patient)
- Do the metastasis sample always have more mutation than paired primary sample?

---

**challenge 2: Can you find some mutated gene signatures in metastasis cancer ?**

---

### Task 2 Answer

```bash
##keep "SNP" data and needed information only
grep "SNP" data_mutations_extended.txt |cut -f 5,6,12,13,14,17 > SNP_tmp.txt
wc -l SNP_tmp.txt
(727146)

#merge sample information with mutation information.
awk -v OFS="(control+V+I)" 'NR==FNR{a[$2]=$1"(control+V+I)"$5}NR>FNR{print $0,a[$6]}'  data_clinical_sample.txt SNP_tmp.txt > SNP_Sampletype_tmp.txt

#extract mutation in Primary Samples and in Metastasis Samples
grep "Primary" SNP_Sampletype_tmp.txt >SNP_Primary_tmp.txt
grep "Metastasis" SNP_Sampletype_tmp.txt >SNP_Metastasis_tmp.txt

#find paired-samples
cut -f 7 SNP_Primary_tmp.txt |sort |uniq >Primary_patient_tmp.txt
cut -f 7 SNP_Metastasis_tmp.txt |sort |uniq >Metastasis_patient_tmp.txt
cat Primary_patient_tmp.txt Metastasis_patient_tmp.txt |sort|uniq -d >paired_patient_tmp.txt
for i in `cat paired_patient_tmp.txt`;do echo "$i";grep "$i" SNP_Sampletype_tmp.txt|cut -f 8|sort|uniq -c ;done|head
#GENIE-COLU-00191
#     18 Metastasis
#      2 Primary
#GENIE-COLU-00294
#      3 Metastasis
#      3 Primary
#...
```
