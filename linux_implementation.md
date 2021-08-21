# ATAC-seq Data Analysis-Linux Implementation

The genome of eukaryotes such has human is packed and organized in place with the assembling of smaller structures called nucleosome. This structure is a complex of 8 histone protein with dna tightly wound around it, ~147bp. The nature of packing (open/loose or tight) of these nucleosomes complex plays a major role in Transcription and it determines the accessibility of transcription factors and protein(RNA polymerase) to bind to the transcription site. Several factors are responsible for the accessibility of the DNA including histone modifications, position of nucleosomes etc., these regulates the activation and deactivation of genes.

Assay for Transposase-Accessible Chromatin using sequencing (ATAC-Seq) is a method to investigate the accessibility of chromatin and thus a method to determine regulatory mechanisms of gene expression. It can help in identifying promoter regions and potential enhancers and silencers sites. ATAC-Seq has become popular for identifying accessible regions of the genome as it&#39;s easier, faster and requires less cells than alternative techniques, such as FAIRE-Seq and DNase-Seq.

To find accessible (open) chromatin regions with ATAC-Seq, the genome is treated with a hyperactive derivative of the **Tn5 transposase**. A transposase can bind to a transposable element (jumping genes) within the genome. During ATAC-Seq, the modified Tn5 inserts DNA sequences corresponding to truncated Nextera adapters into open regions of the genome and concurrently, the DNA is sheared by the transposase activity. The read library is then prepared for sequencing, including PCR amplification with full Nextera adapters and purification steps. Paired-end reads are recommended for ATAC-Seq.

This tutorial uses data from a study by Buenrostro et al. 2013 ([https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html#Buenrostro2013](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html#Buenrostro2013)). The data is from a human cell line of purified CD4+ T cells, called GM12878. The original dataset had 2 x 200 million reads and would be too big to process in a training session, it was down sampled the original dataset to 200,000 randomly selected reads. We also added about 200,000 reads pairs that will map to chromosome 22 to have a good profile on this chromosome, similar to what you might get with a typical ATAC-Seq sample (2 x 20 million reads in original FASTQ). In addition, we want to compare the predicted open chromatin regions to the known binding sites of CTCF, a DNA-binding protein implicated in 3D structure: CTCF. CTCF is known to bind to thousands of sites in the genome and thus it can be used as a positive control for assessing if the ATAC-Seq experiment is good quality.

Good ATAC-Seq data would have accessible regions both within and outside of Transcription Start Site (TSS), for example, at some CTCF binding sites. For that reason, we will download binding sites of CTCF identified by ChIP in the same cell line from ENCODE (ENCSR000AKB, dataset ENCFF933NTR).

## Tutorial Sections

### Preprocessing

1. Get Data
2. Quality Control
3. Trimming Reads

### Mapping

1. Mapping Reads to Reference Genome

### Filtering Mapped Reads

1. Filter Uninformative Reads
2. Filter Duplicate Reads
3. Check Insert Sizes

### Peak calling

1. Call Peaks

### Visualisation of Coverage

1. Prepare the Datasets
2. Create heatmap of coverage at TSS with deepTools
3. Visualise Regions with pyGenomeTracks

## Conclusion

## Tools Used

- FastQC
- Cutadapt
- Bowtie2
- Samtools
- Picard
- Bedtools
- MACS2
- Deeptool
- pyGenome Tracks

# Lets get started

## Getting Tutorial dataset

The data that will be used from this tutorial will be imported from Zenodo. Get files using wget command

$ wget [https://zenodo.org/record/3862793/files/ENCFF933NTR.bed.gz](https://zenodo.org/record/3862793/files/ENCFF933NTR.bed.gz)

$ wget [https://zenodo.org/record/3862793/files/SRR891268\_chr22\_enriched\_R1.fastq.gz](https://zenodo.org/record/3862793/files/SRR891268_chr22_enriched_R1.fastq.gz)

$ wget [https://zenodo.org/record/3862793/files/SRR891268\_chr22\_enriched\_R2.fastq.gz](https://zenodo.org/record/3862793/files/SRR891268_chr22_enriched_R2.fastq.gz)

## Obtain Annotation for hg38 genes

Go to [http://genome.ucsc.edu/cgi-bin/hgTables](http://genome.ucsc.edu/cgi-bin/hgTables) and set the parameters as follows-

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

And then select **Get output.** This downloads a zipped file, chr22.gz

# Cutting columns from table and converting to .bed file

$ gunzip chr22.gz (this uncompressed the zipped file)

To filter out columns 3, 5, 6, 13, 12 and 4 and save in a new file

$ awk -F &quot;\t&quot; &#39;BEGIN{OFS = &quot;\t&quot;}{print $3,$5,$6,$13,$12,$4}&#39; chr22 \&gt; chr22.bed

Output should look like this-

![colunms of 22.bed](C:\Users\KEHINDE\Desktop\task3\awk(1).png "chr22.bed")



## Quality Control

The first step is to check the quality of the reads and the presence of the Nextera adapters. When we perform ATAC-Seq, we can get DNA fragments of about 40 bp if two adjacent Tn5 transposases cut the DNA Adey et al. 2010. This can be smaller than the sequencing length so we expect to have Nextera adapters at the end of those reads.

To carry out quality check on read files, first unzip files as follows:

$ gunzip SRR891268\_chr22\_enriched\_R1.fastq.gz

$ gunzip SRR891268\_chr22\_enriched\_R2.fastq.gz

FastQC tool will be used to perform quality checks on read files. Download the FastQC module.

_Note: FASTQC requires java and javac installed for implementation and you need to run the fasta file from the folder (using the relative/absolute links to the sequence reads)_

$ sudo apt install default-jre

$ sudo apt install default-jdk

$ chmod 755 fastqc (this makes the &quot;fastqc&quot; an executable file)

To perform quality check on both read files (forward and reverse)

$ fastqc SRR891268\_chr22\_enriched\_R1.fastq SRR891268\_chr22\_enriched\_R2.fastq

This generates a HTML report for each file which looks like this, indicating which parameters are good (great tick), bad (red star) or warning (orange).

![QC_report](C:\Users\KEHINDE\Desktop\task3\QC(2).png "QC report")

## Trimming

The fastqc report indicates the presence of an overrepresented sequence and fastqc identifies it as &quot;Nextera Transposase Sequence &#39;&#39; and also adapter content sequence. This sequence is a 3`-adapter sequence which has to be trimmed out.

Read 1 3&#39; adapter sequence&quot;: CTGTCTCTTATACACATCTCCGAGCCCACGAGAC

Read 2 3&#39; adapter sequence&quot;: CTGTCTCTTATACACATCTGACGCTGCCGACGA

Trimming (paired end) is done with Cutadapt tool:

$ cutadapt -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -q 20 --minimum-length 20 -o SRR891268\_chr22\_enriched\_R1\_out.fastq -p SRR891268\_chr22\_enriched\_R2\_out.fastq SRR891268\_chr22\_enriched\_R1.fastq

SRR891268\_chr22\_enriched\_R2.fastq

This command cut out these sequences at the 3` from the forward and reverse reads respectively, filters out reads with quality score less than 20 and length less than 20 long

![filtering](C:\Users\KEHINDE\Desktop\task3\filter(3).png "filtered reads")


**Screenshot after trimming and Rerunning fastqc on the output file**

![QC report](C:\Users\KEHINDE\Desktop\task3\QC4.png "QC report after trimming")


# Mapping

## Mapping Reads to Reference Genome

Now we will map the trimmed reads to the human reference genome using **Bowtie2**. We will extend the maximum fragment length (distance between read pairs) from 500 to 1000 because we know some valid read pairs are from this fragment length, then use --very-sensitive parameter to have more chance to get the best match even if it takes a bit longer to run. We will run the end-to-end mode because we trimmed the adapters so we expect the whole read to map, no clipping of ends is needed. Reference genome chosen is Human chr22 .

Firstly, we have to download our reference genome, chr22 which would be used for mapping.

$ wget --timestamping &#39;ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz&#39; -O chr22.fa.gz

Then we create an index file for our reference genome

$ bowtie2-build chr22.fa.gz indexed\_chr22

Mapping our reads (forward and reverse) to reference genome

$ bowtie2 --very-sensitive --maxins 1000 --dovetail -x indexed\_chr22 -1 SRR891268\_chr22\_enriched\_R1\_out.fastq -2 SRR891268\_chr22\_enriched\_R2\_out

.fastq -S SRR891268\_chr22\_enriched\_out.sam

Output should be as follows-

![alignment results](C:\Users\KEHINDE\Desktop\task3\align5.png "reads alignment")


# Filtering Mapped Reads

## Filter Uninformative Reads

We apply some filters to the reads after the mapping. Here remove reads with low mapping quality and reads that are not properly paired **.**

Firstly, we convert our sam file to bam

$ samtools view -bSo SRR891268\_chr22\_enriched\_out.bam SRR891268\_chr22\_enriched\_out.sam

Then we filter

$ samtools view -q 30 -f 0x2 -b -h SRR891268\_chr22\_enriched\_out.bam \&gt; SRR891268\_chr22\_enriched\_out.filt.bam

This will filter out uninformative reads (Mapping quality \&gt;= 30 &amp; Properly Paired)

## Filter Duplicate Reads

Because of the PCR amplification, there might be read duplicates (different reads mapping to exactly the same genomic region) from overamplification of some regions. As the Tn5 insertion is random within an accessible region, we do not expect to see fragments with the same coordinates. We consider such fragments to be PCR duplicates. We will remove them with Picard MarkDuplicates Module.

Download picard.jar in your working folder from [here](https://github.com/broadinstitute/picard/releases/download/2.26.0/picard.jar). From that directory,

$ run java -jar picard.jar -h

This confirms the picard tool is fuctioning. Next we&#39;ll sort the filtered bam file;

$ samtools sort -T temp -O bam -o SRR891268\_chr22\_enriched\_out.sorted.bam SRR891268\_chr22\_enriched\_out.filt.bam

Then, mark out duplicates reads:

$ java -jar picard.jar MarkDuplicates I=SRR891268\_chr22\_enriched\_out.sorted.bam O=SRR891268\_chr22\_enriched\_out.markdup.bam M=SRR891268\_chr22\_enriched\_out.metrics

This is what the output metrics file looks like.

![Mark Duplicates](C:\Users\KEHINDE\Desktop\task3\markdup6.png "Mark Duplicates")


Next, index markdup file

$ samtools index SRR891268\_chr22\_enriched\_out.markdup.bam

## Check Insert Sizes

Check Insert Size tells us the size of the DNA fragment the read pairs they came from. For this step we have to make a plot of the frequencies of the reads in the bam file to observe the peaks around where there are likely Tn5 transposase activities into nucleosome-free regions.

Installing an R environment;

$ sudo apt install r-base

$ java -jar picard.jar CollectInsertSizeMetrics I=SRR891268\_chr22\_enriched\_out.markdup.bam O=chart.txt H=insertSizePlot.pdf M=0.5

Output is this:

Two peaks can be observed around the 200bp and 400bp from the plot

![reads insert size](C:\Users\KEHINDE\Desktop\task3\insertsize7.png "insert size"))


# Peak calling

## Call Peaks

First, convert the markdup bam file to bed file:

$ bedtools bamtobed -i SRR891268\_chr22\_enriched\_out.markdup.bam \&gt; SRR891268\_chr22\_enriched\_out.markdup.bed

Using macs2 tool, perform peak calling:

$ macs2 callpeak -t SRR891268\_chr22\_enriched\_out.markdup.bed -n macs\_output -g 50818468 --nomodel --shift -100 --extsize 200 --keep-dup all --call-summits –bdg

# Visualisation of Coverage

## Extract CTCF peaks on chr22 in intergenic regions

filter dataset with condition c1==chr22, this command replaces chr22 with c1 in the file

$ sed &#39;s/chr22/c1/&#39; ENCFF933NTR.bed \&gt; ENCFF933NTR\_filt.bed (replaces chr22 with c1)

Extract colunms with c1 into a new file

$ grep c1 ENCFF933NTR\_filt.bed \&gt; ENCFF933NTR\_chr22.bed

Change c1 back to chr22

$ sed &#39;s/c1/chr22/&#39; ENCFF933NTR\_chr22.bed \&gt; ENCFF933NTR\_CHR22genes.bed

Finding overlapping intervals

$ bedtools intersect -v -a ENCFF933NTR\_CHR22genes.bed -b chr22\_genes.bed \&gt; intergenic\_CTCF\_peaks\_chr22

## Convert bedgraph from MACS2 to bigwig

Execute the following commands to convert the output bedGraph file from macs2 to bigwig (refer to this link for better understand the commands [https://www.biostars.org/p/176875/](https://www.biostars.org/p/176875/) )

$ awk &#39;NR!=1&#39; macs\_output\_treat\_pileup.bdg \&gt; macs.deheader.bedGraph

$ sort -k1,1 -k2,2n macs.deheader.bedGraph \&gt; macs.sorted.bedGraph

$ touch chrom22.sizes

$ nano hg19.chrom.sizes → write only one line (tab delimited) in this file chr22 51304566

$ awk &#39;{print $1,$2,$3,$4}&#39; macs.sorted.bedGraph \&gt; macs.sorted.4.bedGraph

$ bedGraphToBigWig macs.sorted.4.bedGraph hg19.chrom.sizes macs.bw

## Create heatmap of coverage at TSS with deepTools

To check the coverage on specific regions, you can compute a heatmap. We will use the deepTools plotHeatmap

## Using computeMatrix generate the matrix;

• Remove the first header line from chr22.bed file

• Then run

$ computeMatrix reference-point --referencePoint TSS -R chr22.bed -S macs.bw --missingDataAsZero -o output\_from\_computeMatrix.gz

We will now generate a heatmap. Each line will be a transcript. The coverage is summarized with a color code from red (no coverage) to blue (maximum coverage). All TSS will be aligned in the middle of the figure and only the 2 kb around the TSS will be displayed. Another plot, on top of the heatmap, will show the mean signal at the TSS.

plotHeatmap will generate the plot using the output of computeMatrix

$ plotHeatmap -m output\_from\_computeMatrix.gz -out plotHeatMap.png

![Heatmap of coverage at TSS](C:\Users\KEHINDE\Desktop\task3\heat8.png "Heatmap at TSS")


Repeat the previous two steps for plotting **CTCF peaks of chr22 in intergenic regions** with slight moderation:

$ computeMatrix reference-point --referencePoint center -R intergenic\_ctcf\_peaks\_chr22 -S macs.bw --missingDataAsZero -o peak\_output\_from\_computeMatrix.gz

$ plotHeatmap -m peak\_output\_from\_computeMatrix.gz -out intragenic\_plotHeatMap.png

![Heatmap for CTCF peaks](C:\Users\KEHINDE\Desktop\task3\heat9.png "Heatmap")


## Visualise Regions with pyGenomeTracks

In order to visualize a specific region (e.g. the gene RAC2), we can either use a genome browser like IGV or UCSC browser, or use pyGenomeTracks to make publishable figures.

Set up the config.ini file with the following parameters-

_**[test bedgraph]**_

_ **file = macs.bw** _

_ **color = blue** _

_ **height = 5** _

_**title = Coverage from MACS2 (extended +/-100bp)**_

_ **min\_value = 0** _

_**[spacer]**_

_ **height = 0.5** _

_**[narrow]**_

_ **file = intergenic\_ctcf\_peaks\_chr22.encodepeak** _

_ **line\_width = 2** _

_**title = Peaks from MACS2 (extended +/-100bp)**_

_ **type = box** _

_ **color = red** _

_ **show\_labels = false** _

_ **file\_type = narrow\_peak** _

_**[spacer]**_

_ **height = 0.5** _

_**[genes 0]**_

_ **file = chr22.bed** _

_ **height = 7** _

_ **title = genes** _

_ **height = 5** _

_ **color = #ffbbff** _

_**[spacer]**_

_ **height = 0.5** _

_**[narrow 1]**_

_ **file = ENCFF933NTR\_sorted.bed** _

_ **color = #A020F0** _

_ **line\_width = 2** _

_ **title = CTCF peaks** _

_ **type = box** _

_ **show\_labels = false** _

_**[x-axis]**_

Then, sort ENCFF933NTR.bed file-

$ sort -k 1,1 -k2,2n ENCFF933NTR.bed \&gt; ENCFF933NTR\_sorted.bed

Next, install pyGenomeTracks

$ conda -install -c bioconda pyGenomeTracks

To visualize regions

$ pyGenomeTracks --tracks config.ini --region chr22:37,193,000-37,252,000 -o Genome\_track\_plot.png

![tracks visualiazation](C:\Users\KEHINDE\Desktop\task3\10).png "Tracks")
