# metagenomics_gut_microbiome

The goal of this repository is to build a bioinformatics pipeline to assess a mock community of bacteria subsampled from the [Tara Oceans project](http://ocean-microbiome.embl.de/companion.html). The mock dataset can be found in the tutorial [here](https://www.hadriengourle.com/tutorials/meta_assembly/).

This pipeline will consist of the following steps.

- Check the raw data (fastq files) quality with FastQC


- Trim the dataset (if needed) with Trimmomatic, and check again with FastQC
- Assemble with MEGAHIT
- Binning with bowtie2, samtools, and metabat
- Plot!

## Quality control

Before doing anything, I checked the quality of the raw data using FastQC by running:

`fastqc data_raw/tara_reads_R1.fastq.gz`

`fastqc data_raw/tara_reads_R2.fastq.gz`

The quality of the raw data is already quite good: there aren't any over-represented sequences and there is no adapter content. The per base sequence quality is great, except for the bases near the end of the run. ![tara_reads_R1](tara_reads_R1.PNG) The quality score in the y axis is a log representation of the base call error probability e, according to: Q = -10log(e). This is sort of inevitable for biophysical/biochemical reasons (e.g. fluorophores bleach over time, nucleotides fail to incorporate, etc). To remedy these issues. I trimmed the datasets using trimmomatic, using the default parameters specified on their [page](http://www.usadellab.org/cms/?page=trimmomatic):

  `java -jar Trimmomatic-0.39/Trimmomatic-0.39/trimmomatic-0.39.jar  PE -threads 8 -trimlog trimmomatic_log/ data_raw/tara_reads_R1.fastq.gz data_raw/tara_reads_R2.fastq.gz -baseout data_trimmed/tara_reads_trimmed_R LEADING:3TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`
  
The result of this was a marginally better dataset, with none of the sequences having quality scores below 20 (0.01 error rate).
![tara_reads_R1_trimmed](tara_reads_R1_trimmed.PNG)

## Assembly

Now that I have cleaned up my reads, the next and most important question is: what do I actually want to learn from this data? Given that this is just a mock "mystery" dataset, I should probably figure out "Who is there?" To do this, I need to assemble all of my reads into contiguous sequences, or contigs. 

For this purpose, I used [megahit](https://academic.oup.com/bioinformatics/article/31/10/1674/177884?login=false), which is a NGS *de novo* assembler. I chose it because it is a competent but memory-efficient tool.
