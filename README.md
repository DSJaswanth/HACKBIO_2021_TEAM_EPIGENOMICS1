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
 
 <p align="center"> <img src="https://hackbiointernship2021.slack.com/files/U029LFE0FSA/F02BS25EW2H/workflow.png">
 
 
 
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

 ### B Obtain Annotation for hg38 genes
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

2. Click get output

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

- **Converting chr22 file into a bed file: (added by @Nishat)**

1. Unzip the downloaded chr22.gz using 
   ````gunzip chr22.gz ````
   command
2. ````awk -F &quot;\t&quot; &#39;OFS=&quot;\t&quot; {print $3, $5, $6, $13, $12, $4 \&gt; (&quot;chr22.bed&quot;)}&#39; chr22 ````
   (to get only expected columns into a newly created chr22.bed file)
   
3. Output should be as follows-

![](RackMultipart20210820-4-1srrkc3_html_56d6bf868d184b2e.png)


  

 
 
 
 
 
 
 
 
 
 
 



         
 


