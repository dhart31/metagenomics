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

```fastqc data_raw/tara_reads_R1.fastq.gz```

```fastqc data_raw/tara_reads_R2.fastq.gz```

The quality of the raw data is already quite good: there aren't any over-represented sequences and there is no adapter content. The per base sequence quality is great for most of the run, except for the last handful of bases ![shown here](tara_reads_R1.png)

