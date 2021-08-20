<details>
  <summary><h1 style="display:inline-block">ATAC SEQUCENCE BASED ON GALAXY </h1></summary>


# <p align="center"> TEAM EPIGENOMICS1
# <p align="center">**HACKBIO INTERNSHIP**

##  <p align="center"> **ATAC-Seq Analysis**

### **INTRODUCTION**

ATAC-Seq, short for Assay for Transposase-Accessible Chromatin using Sequencing, is a technique for assessing chromatin accessibility and how it affects gene expression. The chromatin consists of several nucleosome units. A nucleosome is a complex of DNA and histone proteins that keeps DNA compact through packaging. For transcription to occur, DNA must be loosened from the nucleosomes to make it accessible to transcription factors. This accessibility of DNA mediated by chromatin is a key component of epigenetics.

The ATAC-Seq analysis described below is based on a [tutorial](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html#trimming-reads) on the Galaxy Project and uses purified CD4+ T cells (GM12878).

Aim; To compare the predicted open chromatin regions to the known binding sites of CTCF, a DNA-binding protein implicated in 3D structure: CTCF. CTCF is known to bind to thousands of sites in the genome and thus it can be used as a positive control for assessing if the ATAC-Seq experiment is good quality. Good ATAC-Seq data would have accessible regions both within and outside of TSS, for example, at some CTCF binding sites. For that reason, we will download

Datasets; Data is gotten from the study of Buenrostro et al. 2013. The data from the original dataset is downsized to 200,000 randomly selected reads and about 200,000 reads pairs that will map to chromosome 22. Binding sites of CTCF identified by ChIP in the same cell line from ENCODE (ENCSR000AKB, dataset ENCFF933NTR are also used.

**WORKFLOW**

**Preprocessing**

**Data Upload**

_Create a new history_

1. Add or import files (2 fastq.gz files and 1 bed.gz file) for analysis (via link or from a data library)
2. Add tags to R1 and R2 files. To add tags:

1. Click on the dataset
2. Click on the tag icon
3. Add a tag starting with **#: #SRR891268\_R1** to the R1 file and **#SRR891268\_R2** to the R2 file.
4. Check that tag appears below the dataset name

_Check datatype of files and edit if necessary_.

1. Click on the pencil icon for the dataset to edit its attributes
2. In the central panel, click on the Datatypes tab on the top
3. Select correct datatype (fastqsanger.gz for the FASTQ files, and encodepeak for the bed.gz file)
4. Click the Change datatype button

**Obtain Annotation for hg38 genes**

1. Select the **USCS Main table browser tool** with the following parameters

1. &quot;clade&quot;: Mammal
2. &quot;genome&quot;: Human
3. &quot;assembly&quot;: Dec. 2013 (GRCh38/hg38)
4. &quot;group&quot;: Genes and Gene Prediction
5. &quot;track&quot;: All GENCODE V37
6. &quot;table&quot;: Basic
7. &quot;region&quot;: position chr22
8. &quot;output format&quot;: all fields from selected table
9. &quot;Send output to&quot;: Galaxy

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

**Quality Control**

Select the **FastQC tool** with the following parameters

1. &quot;Short read data from your current history&quot;: Choose here either only the SRR891268\_R1 file with param-file or use param-files; use Multiple datasets to choose both SRR891268\_R1 and SRR891268\_R2.
2. Inspect the web page output of FastQC tool for the SRR891268\_R1 sample. Check what adapters are found at the end of the reads.

**Trimming Reads**

Select the **Cutadapt tool** with the following parameters

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

**Check Adapter Removal**

Select the **Fast QC tool** with the following parameters

1. &quot;Short read data from your current history&quot;: select the output of Cutadapt param files; use; Multiple datasets to choose both Read 1 Output and Read 2 Output.
2. Click on the galaxy-eye (eye) icon of the report and read the first lines.

**Mapping**

### **Mapping reads to reference genome**

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

**Filtering Mapped Reads**

_Filtering of uninformative mapped reads_

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

**Filter duplicate reads**

_Remove duplicates_

Select the **MarkDuplicates**   **tool** with the following parameters

1. param-file &quot;Select SAM/BAM dataset or dataset collection&quot;: Select the output of Filter tool &quot;BAM&quot;
2. &quot;If true do not write duplicates to the output file instead of writing them with appropriate flags set&quot;: Yes

Click on the eye icon of the MarkDuplicate metrics.

**Check Insert Sizes**

_Plot the distribution of fragment sizes_

Select **Paired-end histogram tool** with the following parameters

1. param-file &quot;BAM file&quot;: Select the output of MarkDuplicates tool &quot;BAM output&quot;
2. &quot;Lower bp limit (optional)&quot;: 0
3. &quot;Upper bp limit (optional)&quot;: 1000

Click on the galaxy-eye (eye) icon of the lower one of the 2 outputs (the png file).

**Peak calling**

**MACS2** is used to find regions corresponding to potential open chromatin regions and to identify regions where reads have piled up (peaks) greater than the background read coverage.

1. **Convert BAM to BED**

Convert BAM file (output of MarkDuplicates) into BED format by **bedtools BAM to BED converter**.

1. **MACS2 callpeak**

**MACS2 callpeak** with the following parameters:

- Are you pooling Treatment Files?: No
  - Select the output of  **bedtools BAM to BED**  converter tool
- Do you have a Control File?: No
- Format of Input Files: Single-end BED
- Effective genome size: _H. sapiens_ (2.7e9)
- Build Model: Do not build the shifting model (--nomodel)
  - Set extension size: 200
  - Set shift size: -100. It needs to be - half the extension size to be centered on the 5&#39;.
- In Additional Outputs:
  - Check Peaks as tabular file (compatible with MultiQC)
  - Check Peak summits
  - Check Scores in bedGraph files
- In  Advanced Options:
  - Composite broad regions: No broad regions
    - Use a more sophisticated signal processing approach to find subpeak summits in each enriched peak region: Yes
  - How many duplicate tags at the exact same location are allowed?: all

**Visualisation of Coverage**

1. **Prepare the Datasets**

1. **Extract CTCF peaks on chr22 in intergenic regions**

 In order to get the list of intergenic CTCF peaks of chr22, select the peaks on chr22 and then exclude the one which overlap with genes.

- **Filter data on any column using simple expressions**   with the following parameters:
  - Filter : Select the first dataset: ENCFF933NTR.bed.gz
  - With following condition: c1==&#39;chr22&#39;
- **bedtools Intersect intervals find overlapping intervals in various ways** with the following parameters:
  - File A to intersect with B: Select the output of  **Filter**  data on any column using simple expressions tool
  - Combined or separate output files: One output file per &#39;input B&#39; file
    - File B to intersect with A: Select the dataset chr22 genes
  - What should be written to the output file?: Write the original entry in A for each overlap (-wa)
  - Required overlap: Default: 1bp
  - Report only those alignments that  **do not**  overlap with file(s) B: Yes
- Rename the datasets intergenic CTCF peaks chr22.

1. **Convert bedgraph from MACS2 to bigwig**

- **Wig/BedGraph-to-bigWig** with the following parameters:
  - Convert: Select the output of  **MACS2**  tool (Bedgraph Treatment).
  - Converter settings to use: Default
- Rename the datasets MACS2 bigwig.

1. **Create heatmap of coverage at TSS with deepTools**

1. **Generate computeMatrix**

- **computeMatrix**   with the following parameters:
  - In  Select regions:
    - Insert Select regions
      - Regions to plot: Select the dataset chr22 genes
  - Sample order matters: No
    - Score file: Select the output of  **Wig/BedGraph-to-bigWig**  tool that should be named MACS2 bigwig.
  - computeMatrix has two main output options: reference-point
  - The reference point for the plotting: beginning of region (e.g. TSS)
  - Show advanced output settings: no
  - Show advanced options: yes
    - Convert missing values to 0?: Yes

1. **Plot with plotHeatmap**

- **plotHeatmap**   with the following parameters:
  - Matrix file from the computeMatrix tool: Select the output of  **computeMatrix**  tool.
  - Show advanced output settings: no
  - Show advanced options: no

![](RackMultipart20210820-4-1t5hz5a_html_a652aa7889402a30.png)

**Figure: plotHeatmap output**

Repeat the procedure for CTCF peaks of chr22 in intergenic regions:

**Generate the matrix**

- **computeMatrix** with the following parameters:
  - In Select regions:
    - Insert Select regions
      - Regions to plot: Select the dataset intergenic CTCF peaks chr22
  - Sample order matters: No
    - Score file: Select the output of  **Wig/BedGraph-to-bigWig**  tool that should be named MACS2 bigwig.
  - Would you like custom sample labels?: No, use sample names in the history
  - computeMatrix has two main output options: reference-point
    - The reference point for the plotting: center of region
  - Show advanced output settings: no
  - Show advanced options: yes
    - Convert missing values to 0?: Yes

- **plotHeatmap**   with the following parameters:
  - Matrix file from the computeMatrix tool: Select the output of  **computeMatrix**  tool.
  - Show advanced output settings: no
  - Show advanced options: yes
    - In Colormap to use for each sample:
      - Insert Colormap to use for each sample
        1. Color map to use for the heatmap: your choice
    - The x-axis label: distance from peak center (bp)
    - The y-axis label for the top panel: CTCF peaks
    - Reference point label: peak center
    - Labels for the regions plotted in the heatmap: CTCF\_peaks
    - Did you compute the matrix with more than one groups of regions?: Yes, I used multiple groups of regions

![](RackMultipart20210820-4-1t5hz5a_html_4f5763259562217e.png)

**Figure:  plotHeatmap output on CTCF**

1. **Visualise Regions with pyGenomeTracks**

- **pyGenomeTracks**  Tool with the following parameters:
  - Region of the genome to limit the operation: chr22:37,193,000-37,252,000
  - In Include tracks in your plot:
    - Insert Include tracks in your plot
      - Choose style of the track: Bigwig track
        1. Plot title: Coverage from MACS2 (extended +/-100bp)
        2. Track file(s) bigwig format: Select the output of Wig/BedGraph-to-bigWig tool called MACS2 bigwig.
        3. Color of track: Select the color of your choice
        4. Minimum value: 0
        5. height: 5
        6. Show visualization of data range: Yes
    - Insert Include tracks in your plot
      - Choose style of the track: NarrowPeak track
        1. Plot title: Peaks from MACS2 (extended +/-100bp)
        2. Track file(s) encodepeak or bed format: Select the output of MACS2 tool (narrow Peaks).
        3. Color of track: Select the color of your choice
        4. display to use: box: Draw a box
        5. Plot labels (name, p-val, q-val): No
    - Insert Include tracks in your plot
      - Choose style of the track: Gene track / Bed track
        1. Plot title: Genes
        2. Track file(s) bed or gtf format: chr22 genes
        3. Color of track: Select the color of your choice
        4. height: 5
        5. Plot labels: yes
          1. Put all labels inside the plotted region: Yes
          2. Allow to put labels in the right margin: Yes
    - Insert Include tracks in your plot
      - Choose style of the track: NarrowPeak track
        1. Plot title: CTCF peaks
        2. Track file(s) encodepeak or bed format: Select the first dataset: ENCFF933NTR.bed.gz
        3. Color of track: Select the color of your choice
        4. display to use: box: Draw a box
        5. Plot labels (name, p-val, q-val): No
    - param-repeat Insert Include tracks in your plot
      - Choose style of the track: X-axis

![](RackMultipart20210820-4-1t5hz5a_html_17fdc4281072b69a.png)

**Figure: pyGenomeTracks output**

As CTCF creates accessible regions, a region containing a peak with no corresponding CTCF peak or TSS could be a putative enhancer. In the pyGenomeTracks plot we see a region like this located in the intron of a gene and another one between genes.

**Conclusion**

ATAC-Seq is a method to investigate the chromatin accessibility and the genome is treated with a transposase (enzyme) called Tn5. It marks open chromatin regions by cutting and inserting adapters for sequencing. Low quality bases, adapter contamination, correct insert size and PCR duplicates (duplication level) were checked. Mapped the reads with  **Bowtie2** , filtered the reads for properly paired, good quality and reads that do not map to the mitochondrial genome. Open chromatin regions were found with  **MACS2** , a tool to find regions of genomic enrichment (peaks). The read coverage around TSS was investigated with the help of  **computeMatrix**  and  **plotHeatmap**. The peaks and other informative tracks, such as CTCF binding regions and hg38 genes were visualised with the help of  **pyGenomeTracks**. At the end, open chromatin regions that did not overlap with CTCF sites or TSS, which could be potential putative enhancer regions detected by the ATAC-Seq experiment.

**REFERENCES**

1. Lucille Delisle, Maria Doyle, Florian Heyl, 2021 **ATAC-Seq data analysis (Galaxy Training Materials)**. [https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/atac-seq/tutorial.html) Online; accessed Thu Aug 19 2021
2. Fu, Y., M. Sinha, C. L. Peterson, and Z. Weng, 2008  **The insulator binding protein CTCF positions 20 nucleosomes around its binding sites across the human genome**. PLoS genetics 4: e1000138. [10.1371/journal.pgen.1000138](https://doi.org/10.1371/journal.pgen.1000138)
3. Adey, A., H. G. Morrison, A. (no last name), X. Xun, J. O. Kitzman _et al._, 2010  **Rapid, low-input, low-bias construction of shotgun fragment libraries by high-density in vitro transposition**. Genome Biology 11: R119. [10.1186/gb-2010-11-12-r119](https://doi.org/10.1186/gb-2010-11-12-r119)
4. Green, B., C. Bouchier, C. Fairhead, N. L. Craig, and B. P. Cormack, 2012  **Insertion site preference of Mu, Tn5, and Tn7 transposons**. Mobile DNA 3: 3. [10.1186/1759-8753-3-3](https://doi.org/10.1186/1759-8753-3-3)
5. Buenrostro, J. D., P. G. Giresi, L. C. Zaba, H. Y. Chang, and W. J. Greenleaf, 2013  **Transposition of native chromatin for fast and sensitive epigenomic profiling of open chromatin, DNA-binding proteins and nucleosome position**. Nature Methods 10: 1213–1218. [10.1038/nmeth.2688](https://doi.org/10.1038/nmeth.2688)
6. Litzenburger, U. M., J. D. Buenrostro, B. Wu, Y. Shen, N. C. Sheffield _et al._, 2017  **Single-cell epigenomic variability reveals functional cancer heterogeneity**. Genome Biology 18: [10.1186/s13059-016-1133-7](https://doi.org/10.1186/s13059-016-1133-7)
7. Kia, A., C. Gloeckner, T. Osothprarop, N. Gormley, E. Bomati _et al._, 2017  **Improved genome sequencing using an engineered transposase**. BMC Biotechnology 17: [10.1186/s12896-016-0326-1](https://doi.org/10.1186/s12896-016-0326-1)
8. Corces, M. R., A. E. Trevino, E. G. Hamilton, P. G. Greenside, N. A. Sinnott-Armstrong _et al._, 2017  **An improved ATAC-seq protocol reduces background and enables interrogation of frozen tissues**. Nature Methods 14: 959–962. [10.1038/nmeth.4396](https://doi.org/10.1038/nmeth.4396)
</details>

