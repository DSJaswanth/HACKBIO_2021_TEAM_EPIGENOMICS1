**Quality Control**

<details >
<summary>Galaxy Implementation</summary>
<br>
-Select the **FastQC tool** with the following parameters
1. &quot;Short read data from your current history&quot;: Choose here either only the SRR891268\_R1 file with param-file or use param-files; use Multiple datasets to choose both SRR891268\_R1 and SRR891268\_R2.
2. Inspect the web page output of FastQC tool for the SRR891268\_R1 sample. Check what adapters are found at the end of the reads.
</details>  
<details >
<summary>Linux Implementation</summary>
<br>
  
-Download the FastQC module
Note: FASTQC requires java and javac installed for implementation and you need to run the fastqc file from the folder (using the relative/absolute links to the sequence reads)
$ sudo apt install default-jre
$ sudo apt install default-jdk
-Make the “fastqc” an executable file
```python
$ chmod 755 fastqc
```
-Run the fastqc on all sequenced reads from its folder
```python  
$ fastqc SRR891268_chr22_enriched_R1.fastq SRR891268_chr22_enriched_R2.fastq  
``` 
The report for each file is generated as an html file and a zip file containing more files that can be customised for reports. Look into the html files
</details>

<details >
<summary>Fig of Results</summary>
<br>  
<img src=“https://user-images.githubusercontent.com/81503326/130264630-4905fa8f-3f34-4dfc-be63-dd63eac3e49c.PNG”>

</details>   

**Trimming Reads**
  
<details >
<summary>Galaxy Implementation</summary>
<br>
  
-Select the **Cutadapt tool** with the following parameters

1. &quot;Single-end or Paired-end reads?&quot;: Paired-end

1. param-file &quot;FASTQ/A file #1&quot;: select SRR891268\_R1
2. param-file &quot;FASTQ/A file #2&quot;: select SRR891268\_R2
3. In &quot;Read 1 Options&quot;:

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

-*Check Adapter Removal*

Select the **Fast QC tool** with the following parameters

1. &quot;Short read data from your current history&quot;: select the output of Cutadapt param files; use; Multiple datasets to choose both Read 1 Output and Read 2 Output.
2. Click on the galaxy-eye (eye) icon of the report and read the first lines.
</details>  
<details >
<summary>Linux Implementation</summary>
<br>
-The fastqc report  indicates the presence of an overrepresented sequence and fastqc identifies it as “Nextera Transposase Sequence ''. This sequence is similar to but longer than the one given in the tutorial.

SRR891268_chr22_enriched_R1 = CTGTCTCTTATACACATCTCCGAGCCCACGAGACTAAGGCGAATCTCGTA (fastqc)
SRR891268_chr22_enriched_R1 = CTGTCTCTTATACACATCTCCGAGCCCACGAGAC (Galaxy tutorial)
SRR891268_chr22_enriched_R2 = CTGTCTCTTATACACATCTGACGCTGCCGACGAGTGTAGATCTCGGTGGT (fastqc)
SRR891268_chr22_enriched_R2 = CTGTCTCTTATACACATCTGACGCTGCCGACGA (Galaxy tutorial)
  
**Install cutadapt running**
 ```python 
-$ sudo apt install cutadapt
```
**For paired end trimming**
 ```python 
-$ cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA --minimum-length 20 -q 20 -o trimmed_1.fastq -p trimmed_2.fastq SRR891268_chr22_enriched_R1.fastq SRR891268_chr22_enriched_R2.fastq
 ```
</details> 
<details >
<summary>Fig of Results</summary>
<br>
</details> 
  
**Mapping**
  
<details >
<summary>Galaxy Implementation</summary>
<br>
-### **Mapping reads to reference genome**

Select the **Bowtie2**   **tool** with the following parameters:

1. &quot;Is this single or paired library&quot;: Paired-end

1. param-file &quot;FASTQ/A file #1&quot;: select the output of Cutadapt tool &quot;Read 1 Output&quot;
2. param-file &quot;FASTQ/A file #2&quot;: select the output of Cutadapt tool &quot;Read 2 Output&quot;
3. &quot;Do you want to set paired-end options?&quot;: Yes

&quot;Set the maximum fragment length for valid paired-end alignments&quot;: 1000

&quot;Allow mate dovetailing&quot;: Yes

1. &quot;Will you select a reference genome from your history or use a built-in index?&quot;: Use a built-in genome index

1. &quot;Select reference genome&quot;: Human (Homo sapiens): hg38 Canonical
2. &quot;Set read groups information?&quot;: Do not set
3. &quot;Select analysis mode&quot;: 1: Default setting only

&quot;Do you want to use presets?&quot;: Very sensitive end-to-end (--very-sensitive)

1. &quot;Do you want to tweak SAM/BAM Options?&quot;: No
2. &quot;Save the bowtie2 mapping statistics to the history&quot;: Yes

1. Click on the galaxy-eye (eye) icon of the mapping stats.
</details>  
<details >
<summary>Linux Implementation</summary>
<br>
  
**Pulling the sequence for chromosome 22 for indexing and mapping**
```python   
-$ wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz' -O chr22.fa.gz
   For mapping to chr22-
```
**install bowtie2**
  
-Create index for Chromosome 22:  
```python  bowtie2-build chr22.fa.gz indexed_chr22 ```
-Start mapping for the parameters specified by Galaxy: 
```pythonbowtie2 --very-sensitive --maxins 1000 --dovetail -x indexed_chr22 -1 trimmed_1.fastq -2 trimmed_2.fastq -S Aligned_output.sam ```
</details>
  
<details >
<summary>Fig of Results</summary>
<br>
</details>    

  

  
  

























