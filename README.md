# <p align="center">🤩 HackBio2021-Team Epigenomics1 🤩
 ![Visitor](https://visitor-badge.laobi.icu/badge?page_id=DSJawanth.HACKBIO_2021_TEAM_EPIGENOMICS1)
  ![HackBio](https://media-exp1.licdn.com/dms/image/C561BAQHKcVQGbcedOA/company-background_10000/0/1598491473588?e=2159024400&v=beta&t=rxECjvQ_YSc28Dn0n9YOtDoFFmvXjatRiqc__C2mpU0) <br>
  Image source: HackBio (LinkedIn) 
  
  <p align="center"> <img src="https://i.redd.it/sc0r7yw63b231.gif"> 
  
 This is ***Team EPIGENOMICS1 *** for Hackbio 2021 virtual internship.
 
- 🌱 About: Hackbio Internship is a 5-weeks virtual research internship that is practice oriented and focused on equipping scientists globally with advanced bioinformatics and      computational biology skills. 
- :desktop_computer: Hackbio official website: https://thehackbio.com/
- 📫 Contact: contact@hackbio.com
- :man_technologist: We are a diverse team of awesome 17 members 👩‍💻.
<!---
![Top Langs](https://github-readme-stats.vercel.app/api/top-langs/?username=DSJawanth&layout=compact)
--->

<br>  

##  :scroll: HEADING OVER TO THE TASK :scroll:
## <p align="center"> STAGE-2-PROBE

👉 SEARCH FOR DESIRED TUTORIALS OF BIOSTACK (EPIGENOMICS) ANALYSIS.<br>
👉 REPRODUCE THE ANALYSIS.<br>
👉 CREATE A COMPREHENSIVE MARKDOWN.<br>
👉 ADD STEP BY STEP OF ANALYSIS TO GITHUB LINK.<br>

#Dont be so crazy to look up our workflow 
BEFORE LETS KNOW WHAT IS MY TEAM BIOSTACK IS ABOUT

##  <p align="center"> **ATAC-Seq Analysis**
 
 ### **INTRODUCTION**

ATAC-Seq, short for Assay for Transposase-Accessible Chromatin using Sequencing, is a technique for assessing chromatin accessibility and how it affects gene expression. The chromatin consists of several nucleosome units. A nucleosome is a complex of DNA and histone proteins that keeps DNA compact through packaging. For transcription to occur, DNA must be loosened from the nucleosomes to make it accessible to transcription factors. This accessibility of DNA mediated by chromatin is a key component of epigenetics.

The ATAC-Seq analysis described below is based on a [tutorial](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html#trimming-reads) on the Galaxy Project and uses purified CD4+ T cells (GM12878).

#### Aim:
 - To compare the predicted open chromatin regions to the known binding sites of CTCF, a DNA-binding protein implicated in 3D structure: CTCF. CTCF is known to bind to thousands of sites in the genome and thus it can be used as a positive control for assessing if the ATAC-Seq experiment is good quality. Good ATAC-Seq data would have accessible regions both within and outside of TSS, for example, at some CTCF binding sites. For that reason, we will download

#### Datasets
 - Data is gotten from the study of Buenrostro et al. 2013. The data from the original dataset is downsized to 200,000 randomly selected reads and about 200,000 reads pairs that will map to chromosome 22. Binding sites of CTCF identified by ChIP in the same cell line from ENCODE (ENCSR000AKB, dataset ENCFF933NTR are also used.

 **Downloading Dataset: **

 ***Download files from the [here](https://zenodo.org/record/3862793#.YRze2XUvNH4) directly for galaxy  or use the wget command along with link in the linux***
 
 1. [https://zenodo.org/record/3862793/files/ENCFF933NTR.bed.gz] <br>
 2. [https://zenodo.org/record/3862793/files/SRR891268\_chr22\_enriched\_R1.fastq.gz]<br>
 3. [https://zenodo.org/record/3862793/files/SRR891268\_chr22\_enriched\_R2.fastq.gz] <br>
 
 ## **WORKFLOW/METHODOLOGY**
 
 
 ###GRAPHICAL WORKFLOW DESIGN
 
 <p align="center"> <img src="https://hackbiointernship2021.slack.com/files/U029LFE0FSA/F02BS25EW2H/workflow.png">  #### result image

 
 
 
### STEP 1 :- PREPROCESSING 
 
 <details><summary><h3>A)Data Upload for galaxy</h3></summary><br>
  
***Create a new history***
  
1. Add or import DATASETS (2 fastq.gz files and 1 bed.gz file) for analysis (via link or from a data library)
2. Add tags to R1 and R2 files. To add tags:
    - Click on the dataset
    - Click on the tag icon
    - Add a tag starting with **#: FOR  #SRR891268\_R1** FILE to the R1 file and  FOR **#SRR891268\_R2** FILE to the R2 file.
    - Check that tag appears below the dataset name

***Check datatype of files and edit if necessary***

1. Click on the pencil icon for the dataset to edit its attributes
2. In the central panel, click on the Datatypes tab on the top
3. Select correct datatype (fastqsanger.gz for the FASTQ files, and encodepeak for the bed.gz file)
4. Click the Change datatype button
 </details>

 #### B Obtain Annotation for hg38 genes
 <details>
 <summary>FOR GALAXY IMPLEMEMTATION</summary>
 
1. Select the **USCS Main table browser tool** with the following parameters

- &quot;clade&quot;: Mammal
- &quot;genome&quot;: Human
- &quot;assembly&quot;: Dec. 2013 (GRCh38/hg38)
- &quot;group&quot;: Genes and Gene Prediction
- &quot;track&quot;: All GENCODE V37
- &quot;table&quot;: Basic
- &quot;region&quot;: position chr22
- &quot;output format&quot;: all fields from selected table
- &quot;Send output to&quot;: Galaxy

2. Click get output  #### result image


3. Click Send Query to Galaxy

4. Select the **Cut columns from a table tool** with the following parameters

1. &quot;Cut columns&quot;: c3,c5,c6,c13,c12,c4
2. &quot;Delimited by&quot;: Tab
3. param-file &quot;From&quot;: UCSC Main on Human: wgEncodeGencodeBasicV37 (chr22:1-50,818,468)

_Rename the dataset as chr22 genes_

1. Click on the pencil icon for the dataset to edit its attributes
2. In the central panel, change the name field
3. Click the Save button

_Change the datatype to a BED format_

1. Click on the pencil icon for the dataset to edit its attributes
2. In the central panel, click on the Datatypes tab on the top
3. Select  **bed**
4. Click the Change datatype button

Click on the eye icon to check changes effected. There should be matching column names in each column of the dataset.
  </details>
 
  <details><summary>FOR LINUX IMPLEMEMTATION</summary>
   Go to [http://genome.ucsc.edu/cgi-bin/hgTables](http://genome.ucsc.edu/cgi-bin/hgTables) and set the parameters as-

  - &quot;_ **clade&quot;** _ **:**  **Mammal**
  - &quot;_ **genome&quot;** _ **:**  **Human**
  - &quot;_ **assembly&quot;** _ **:** **Dec. 2013 (GRCh38/hg38)**
  - &quot;_ **group&quot;** _ **:**  **Genes and Gene Prediction**
  - &quot;_ **track&quot;** _ **:**  **All GENCODE V37**
  - &quot;_ **table&quot;** _ **:**  **Basic**
  - &quot;_ **region&quot;** _ **:**  **position**** chr22**
  - &quot;_ **output format&quot;** _ **:**  **all fields from selected table**
  - &quot;_ **output filename:&quot; chr22** _
  - &quot;_ **file type returned:&quot; gzipped compressed** _

And then select **Get output**

Thus, chr22.gz file will be downloaded.

- **Converting chr22 file into a bed file: 
   
1. Unzip the downloaded chr22.gz using 
   ````gunzip chr22.gz ````
   command
2. ````awk -F &quot;\t&quot; &#39;OFS=&quot;\t&quot; {print $3, $5, $6, $13, $12, $4 \&gt; (&quot;chr22.bed&quot;)}&#39; chr22 ````
   (to get only expected columns into a newly created chr22.bed file)
   
3. Output should be as follows-

![](RackMultipart20210820-4-1srrkc3_html_56d6bf868d184b2e.png) 
   #### result image  </details>

  
<details><summary>DATA PREPROCESSING FOR LINUX</summary>
  ***You can unzip the sequence files with gunzip***

````$ gunzip SRR891268_chr22_enriched_R1.fastq.gz ````<bR> 
````$ gunzip SRR891268_chr22_enriched_R2.fastq.gz ````
 </details>
 
#### C)QUALITY CONTROL
 
 <details>
<summary>Galaxy Implementation</summary>
<br>
Select the **FastQC tool** with the following parameters
1. &quot;Short read data from your current history&quot;: Choose here either only the SRR891268\_R1 file with param-file or use param-files; use Multiple datasets to choose both SRR891268\_R1 and SRR891268\_R2.
2. Inspect the web page output of FastQC tool for the SRR891268\_R1 sample. Check what adapters are found at the end of the reads.
</details>  
 
 <details >
<summary>Linux Implementation</summary>
<br>

- Download the FastQC module
Note: FASTQC requires java and javac installed for implementation and you need to run the fastqc file from the folder (using the relative/absolute links to the sequence reads)<bR> 
```$ sudo apt install default-jre```<br>
```$ sudo apt install default-jdk```<bR> 
Make the “fastqc” an executable file<bR>
```python $ chmod 755 fastqc```<bR> 
- Run the fastqc on all sequenced reads from its folder<bR> 
```python  $ fastqc SRR891268_chr22_enriched_R1.fastq```<bR> 
```SRR891268_chr22_enriched_R2.fastq ```<bR> 
The report for each file is generated as an html file and a zip file containing more files that can be customised for reports. Look into the html files.
</details>

#### image 
 
#### TRIMMING READ
<details>
<summary>Galaxy Implementation</summary>
<br>

 - Select the **Cutadapt tool** with the following parameters

1. &quot;Single-end or Paired-end reads?&quot;: Paired-end

2. param-file &quot;FASTQ/A file #1&quot;: select SRR891268\_R1
3. param-file &quot;FASTQ/A file #2&quot;: select SRR891268\_R2
4. In &quot;Read 1 Options&quot;:

In &quot;3&#39; (End) Adapters&quot;:

param-repeat &quot;Insert 3&#39; (End) Adapters&quot;

&quot;Source&quot;: Enter custom sequence

&quot;Enter custom 3&#39; adapter name (Optional if Multiple output is &#39;No&#39;)&quot;: Nextera R1

&quot;Enter custom 3&#39; adapter sequence&quot;: CTGTCTCTTATACACATCTCCGAGCCCACGAGAC

1. In &quot;Read 2 Options&quot;:

In &quot;3&#39; (End) Adapters&quot;:

param-repeat &quot;Insert 3&#39; (End) Adapters&quot;

&quot;Source&quot;: Enter custom sequence

&quot;Enter custom 3&#39; adapter name (Optional)&quot;: Nextera R2

&quot;Enter custom 3&#39; adapter sequence&quot;: CTGTCTCTTATACACATCTGACGCTGCCGACGA

1. In &quot;Filter Options&quot;:

&quot;Minimum length&quot;: 20

1. In &quot;Read Modification Options&quot;:

&quot;Quality cutoff&quot;: 20

1. In &quot;Output Options&quot;:

&quot;Report&quot;: Yes

1. Click on the galaxy-eye (eye) icon of the report and read the first lines.
  - Check Adapter Removal

Select the **Fast QC tool** with the following parameters

1. &quot;Short read data from your current history&quot;: select the output of Cutadapt param files; use; Multiple datasets to choose both Read 1 Output and Read 2 Output.
2. Click on the galaxy-eye (eye) icon of the report and read the first lines.
  </details>
 
<details >
<summary>Linux Implementation</summary>
<br> 
 
- ##### Adapter Trimming 

The fastqc report indicates the presence of an overrepresented sequence and fastqc identifies it as &quot;Nextera Transposase Sequence &#39;&#39;. This sequence is similar to but longer than the one given in the tutorial.

```SRR891268\_chr22\_enriched\_R1 = CTGTCTCTTATACACATCTCCGAGCCCACGAGACTAAGGCGAATCTCGTA (fastqc)``` <br>

```SRR891268\_chr22\_enriched\_R1 = CTGTCTCTTATACACATCTCCGAGCCCACGAGAC (Galaxy tutorial)```<br>

```SRR891268\_chr22\_enriched\_R2 = CTGTCTCTTATACACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGT (fastqc)```<br>


```SRR891268\_chr22\_enriched\_R2 = CTGTCTCTTATACACATCTGACGCTGCCGACGA (Galaxy tutorial)```<br>


 - ##### Adapter Trimming with Cutadapt 

Install cutadapt running-

```$ sudo apt install cutadapt```

For paired end trimming-
 
 ```$ cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA --minimum-length 20 -q 20 -o trimmed\_1.fastq -p trimmed\_2.fastq SRR891268\_chr22\_enriched\_R1.fastq SRR891268\_chr22\_enriched\_R2.fastq```
</details>

### STEP2 :- MAPPING 
  
<details>
<summary>Galaxy Implementation</summary>
<br>
 
 ### **Mapping reads to reference genome**

 - Select the **Bowtie2**   **tool** with the following parameters:

1. &quot;Is this single or paired library&quot;: Paired-end
2. param-file &quot;FASTQ/A file #1&quot;: select the output of Cutadapt tool &quot;Read 1 Output&quot;
3. param-file &quot;FASTQ/A file #2&quot;: select the output of Cutadapt tool &quot;Read 2 Output&quot;
4. &quot;Do you want to set paired-end options?&quot;: Yes

   - &quot;Set the maximum fragment length for valid paired-end alignments&quot;: 1000
&quot;Allow mate dovetailing&quot;: Yes

1. &quot;Will you select a reference genome from your history or use a built-in index?&quot;: Use a built-in genome index
2. &quot;Select reference genome&quot;: Human (Homo sapiens): hg38 Canonical
3. &quot;Set read groups information?&quot;: Do not set
4. &quot;Select analysis mode&quot;: 1: Default setting only

   - &quot;Do you want to use presets?&quot;: Very sensitive end-to-end (--very-sensitive)

1. &quot;Do you want to tweak SAM/BAM Options?&quot;: No
2. &quot;Save the bowtie2 mapping statistics to the history&quot;: Yes
3.Click on the galaxy-eye (eye) icon of the mapping stats.
  </details>

<details >
<summary>Linux Implementation</summary>
<br>
 
- Mapping and Alignment 

  Pulling the sequence for chromosome 22 for indexing and mapping
```$ wget --timestamping &#39;ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz&#39; -O chr22.fa.gz```

- For mapping to chr22-

1. install bowtie2
2. Create index for Chromosome 22:``` bowtie2-build chr22.fa.gz indexed\_chr22```
3. Start mapping for the parameters specified by Galaxy: ```bowtie2 --very-sensitive --maxins 1000 --dovetail -x indexed\_chr22 -1 trimmed\_1.fastq -2 trimmed\_2.fastq -S Aligned\_output.sam```

Output should be as follows-
 
 
 </details>
  
### STEP 3 :- Filtering Mapped Reads 
  
#### A) Filter Uninformative Reads
  
  <details>
<summary>Galaxy Implementation</summary>
<br>
   
 ****_Filtering of uninformative mapped reads_****

Select the **Filter BAM datasets on a variety of attributes**   **tool** with the following parameters

param-file &quot;BAM dataset(s) to filter&quot;: Select the output of Bowtie2 tool &quot;alignments&quot;

In &quot;Condition&quot;:

1. param-repeat &quot;Insert Condition&quot;
2. In &quot;Filter&quot;:
3. param-repeat &quot;Insert Filter&quot;

&quot;Select BAM property to filter on&quot;: mapQuality

&quot;Filter on read mapping quality (phred scale)&quot;: \&gt;=30

1. param-repeat &quot;Insert Filter&quot;

&quot;Select BAM property to filter on&quot;: isProperPair

&quot;Select properly paired reads&quot;: Yes

1. param-repeat &quot;Insert Filter&quot;

&quot;Select BAM property to filter on&quot;: reference

&quot;Filter on the reference name for the read&quot;: !chrM

&quot;Would you like to set rules?&quot;: No

Click on the input and the output BAM files of the filtering step. Check the size of the files.

   
 </details >
  
<details >
<summary>Linux Implementation</summary>
<br>
 
****_Filtering of uninformative mapped reads_****

1. Install samtools
2. ```samtools view -q 30 -f 0x2 -b -h Aligned\_output.sam \&gt; Filtered\_output.bam```

This will filter out uninformative reads (Mapping quality \&gt;= 30 &amp; Properly Paired)

 </details> 
  
#### B) Filter duplicate reads
  
  <details>
<summary>Galaxy Implementation</summary>
<br>
   
***_Remove duplicates_***

Select the **MarkDuplicates**   **tool** with the following parameters

1. param-file &quot;Select SAM/BAM dataset or dataset collection&quot;: Select the output of Filter tool &quot;BAM&quot;
2. &quot;If true do not write duplicates to the output file instead of writing them with appropriate flags set&quot;: Yes

Click on the eye icon of the MarkDuplicate metrics.
  
  </details>
  
 <details >
<summary>Linux Implementation</summary>
<br>
  
 **_Mark Duplicate Reads_**

- Download picard.jar in your working folder from [here](https://github.com/broadinstitute/picard/releases/download/2.26.0/picard.jar)
- From that directory, run ````java -jar picard.jar -h```` to check whether it works (you can skip this step)
- For sorting the output file from last step use-```` samtools sort -T temp -O bam -o filtered\_output\_sorted.bam Filtered\_output.bam````
- Finally, run ````java -jar picard.jar MarkDuplicates I=filtered\_output\_sorted.bam O=marked\_dup.bam M=marked\_dup.metrics.txt```` for marking duplicates
- If you can have a look into the metrics in the metrics.txt file
  
    </details>
  
  
  #### C) Check Insert Sizes 
  
  <details>
<summary>Galaxy Implementation</summary>
<br>
   
****_Plot the distribution of fragment sizes_****

Select **Paired-end histogram tool** with the following parameters

1. param-file &quot;BAM file&quot;: Select the output of MarkDuplicates tool &quot;BAM output&quot;
2. &quot;Lower bp limit (optional)&quot;: 0
3. &quot;Upper bp limit (optional)&quot;: 1000

Click on the galaxy-eye (eye) icon of the lower one of the 2 outputs (the png file).

</details>
 
  <details >
<summary>Linux Implementation</summary>
<br>

***_Check Insert Sizes_***
   
Check Insert Size tells us the size of the DNA fragment the read pairs came from. For this step we have to make a plot of the frequencies of the reads in the bam file to observe the peaks around where there are likely Tn5 transposase activities into nucleosome-free regions.

````$ sudo apt install r-base````<br>

````$ java -jar picard.jar CollectInsertSizeMetrics I=marked\_dup.bam O=chart.txt H=insertSizePlot.pdf M=0.5````

![](RackMultipart20210820-4-1srrkc3_html_c751a794d88d647c.png)

Two peaks can be observed around the 200bp and 400bp from the plot

  </details>

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
