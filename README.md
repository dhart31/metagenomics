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

Now that I have cleaned up my reads, the next and most important question is: what do I actually want to learn from this data? Given that this is a mock dataset, I should probably figure out "Who is there?" To do this, I need to assemble all of my reads into contiguous sequences, or contigs. From these, I can infer some genomes

For this purpose, I used [megahit](https://academic.oup.com/bioinformatics/article/31/10/1674/177884?login=false), which is a NGS *de novo* assembler. I chose it because it is a memory-efficient tool that I can easily use on my local machine.

Despite its complicated algorithm under the hood, running megahit is very easy, and simply consists of writing this line:
`megahit -1 data_trimmed/tara_reads_trimmed_1P.gz -2 data_trimmed/tara_reads_trimmed_2P.gz -o assembly`

After 194 seconds, the code output the following line:

`2022-08-19 23:51:21 - 6046 contigs, total 22996494 bp, min 215 bp, max 2448145 bp, avg 3803 bp, N50 20909 bp`

The key measure of "assembly quality" given here is N50, which is equal to the contig size at the point of half of the mass of the contig distribution. In other words, the combined length of contigs with size N50 or larger is equal to the combined length of contigs smaller than N50. Since I do not know the genome sizes in this mock community, its difficult to say how good this number is. N50 has [several shortcomings](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5783553/), but alternatives like U50 rely on a reference genome, which I don't have.

## Binning

Now that I have assembled some contigs, I can try to reconstruct all of the genomes in my community. To do this, I will use [MetaBAT](https://bitbucket.org/berkeleylab/metabat/src/master/), which separates contigs into bins of taxonomic groups. This is achieved by comparing the tetranucleotide frequencies (TNR) and abundance (roughly, how many reads align to the contig) for a pair of contigs. If the TNR and abundance of a pair of contigs are similar, they will be paired together.

To use MetaBAT, I also need a sorted BAM file, which describes how all of my reads align to my contigs. To create the BAM file and then sort it, I used Bowtie2 and samtools in tandem, which are both popularly used for this task. 

 `bowtie2-build final.contigs.fa final.contigs`

```bowtie2 -x bam_files/final.contigs -1 data_trimmed_tara_reads_trimmed_1P.gz -2 data_trimmed/tara_reads_trimmed_2P.gz samtools view -bS -o bam_files/tara_needs_sorting.bam```

Doing these yielded the following output:

```
1499996 reads; of these:
  1499996 (100.00%) were paired; of these:
    1067699 (71.18%) aligned concordantly 0 times
    432069 (28.80%) aligned concordantly exactly 1 time
    228 (0.02%) aligned concordantly >1 times
    ----
    1067699 pairs aligned concordantly 0 times; of these:
      1047769 (98.13%) aligned discordantly 1 time
    ----
    19930 pairs aligned 0 times concordantly or discordantly; of these:
      39860 mates make up the pairs; of these:
        14917 (37.42%) aligned 0 times
        18023 (45.22%) aligned exactly 1 time
        6920 (17.36%) aligned >1 times
99.50% overall alignment rate
```
Given that I'm comparing my reads to sequences assembled from the reads, its not surprising that they align 99% percent of the time. The assembler worked. It is interesting to me that only 29% of my pairs aligned concordantly, meaning that their relative orientation and distance range was not what bowtie2 expected (illumina fragment sizes are supposed to be 200 to 500 bases). To me this suggests that there are structural variants of the same genome in this sample.

