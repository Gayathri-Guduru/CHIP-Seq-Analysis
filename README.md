# ChIP-Seq Analysis for Drosophila Spiked-In Datasets

## Overview:
In this workflow, we'll process Drosophila ChIP-Seq datasets that have been spiked-in with data from another species. The purpose of the spike-in is to account for variability in ChIP and sequencing efficiencies, allowing for more accurate normalization between samples.

## Table of Contents
• Quality Control \
• Alignment \
• Post-processing \
• Peak Calling \
• Differential Binding Analysis \
• Annotation & Functional Analysis 

## Commonly used file formats for ChIP-seq

• FASTA \
• FASTQ – Fasta with quality \ 
• SAM – Sequence Alignment or Map format \
• BAM – Binary Sequence Alignment/Map format \
• Bed – Basic genome interval \
• BedGraph \
• Wiggle (wig, bigwig) – tab-limited format to represent continuous value 

## ChIP-seq data have several characteristics:

• Histone modifications cover broader regions of DNA than TFs.
• Reads are trimmed to within a smaller number of bases.
• Fragments are quite large relative to binding sites of TFs.
• Measurements of histone modification often undulate following well-positioned nucleosomes.

## To extract meaningful data from the raw sequence reads, the ChIP-seq data analysis should:

• Identify genomic regions - ‘peaks’ - where TF binds or histones are modified.
• Quantify and compare levels of binding or histone modification between samples.
• Characterize the relationships among chromatin state and gene expression or splicing.

![image](https://github.com/Gayathri-Guduru/CHIP-Seq-Analysis/assets/98939664/41e438cb-a737-4d5a-9bae-a72bdf4f5097)


# Data pre-processing
## Quality control for Quality Check
.fastq files are quality-checked using FastQC
```
fastqc sample.fastq.gz -o output_directory/
```
 ## MultiQC for Aggregating Results
```
multiqc output_directory/ -o multiqc_report/
```

## Trimming (optional)
When needed, trimming is performed with ```trimmomatic```
```
Trimmomatic PE -threads 1 $input/file/1 $input/file/2 /$output/trimmed_paired/file/1 $output/trimmed_unpaired/file/1 $output/trimmed_paired/file/2 $output/trimmed_unpaired/file/2 ILLUMINACLIP:$Trimmomatic_Path/adapters/TruSeq3-PE-2.fa:2:40:12:8:true LEADING:10 SLIDINGWINDOW:4:15 MINLEN:50 2> read_processing.log
```

# Alignment to reference genome
## Building Combined Reference Genome

Clean reads are aligned to reference genome using ```Bowtie2```
```
# Combine Drosophila and spike-in (e.g., yeast) genomes
cat drosophila_genome.fasta yeast_genome.fasta > combined_genome.fasta
```

## bowtie2-build combined_genome.fasta combined_genome
```
bowtie2-build combined_genome.fasta combined_genome
```

## Aligning Reads to the Combined Genome
```
bowtie2 -x combined_genome -U sample.fastq.gz -S aligned_sample.sam
```
## Sorting
Output alignment files will be directly in .bam format
Aligned reads are sorted by genomic coordinates using ```samtools sort``` 

## Convert SAM to BAM
```
samtools view -bS aligned_sample.sam > aligned_sample.bam
```
## Sort BAM File
```
samtools sort aligned_sample.bam -o sorted_sample.bam
```

## Filtering
All unmapped reads, duplicates and multimappers are filtered out using ```sambamb```, ```PicardMarkDuplicates```
```
picard MarkDuplicates I=sorted_sample.bam O=marked_sample.bam M=marked_dup_metrics.txt
```

## Merging and indexing (optional)
Final filtered reads from each replicate are merged using ```samtools merge```, and indexed using ```samtools index```
```
samtools index marked_sample.bam
```

## Normalization using Spike-in:

Calculate the scaling factor based on the Drosophila spike-in reads.
Use the scaling factor to normalize the read counts in the target genome.

## Peak Calling: Use tools like MACS2 to identify regions of enrichment in the ChIP sample compared to the input control.
```
macs2 callpeak -t marked_sample.bam -c control_sample.bam -f BAM -g dm -n sample_name
```
## Visualization: Use tools like deepTools or the Integrative Genomics Viewer (IGV) to visualize the ChIP-Seq profiles and peaks.

## Differential Binding Analysis (if comparing different conditions): Use tools like DiffBind or DESeq2 to find regions with differential binding between conditions.
## In R
```library(DiffBind)
dba <- dba(sampleSheet="samples.csv")
dba <- dba.analyze(dba)
dba.report(dba)
```
Remember, the key difference when working with spike-ins is the normalization step. Proper normalization using the spike-in controls is crucial to obtaining accurate results.

## Annotation and Functional Analysis: Annotate the peaks using tools like ChIPseeker or HOMER. This will help identify the genes associated with the peaks and the potential functions of the protein being studied.
```
# Use bedtools to intersect peaks with gene annotations
bedtools intersect -a peaks.bed -b drosophila_annotations.bed -wa -wb > annotated_peaks.bed
```
## Functional Enrichment Analysis
```
# In R, using clusterProfiler
library(clusterProfiler)
genes <- read.delim("annotated_peaks.bed", header=FALSE)$V4
enrich <- enrichGO(genes, orgDb="org.Dm.eg.db", ont="BP")
dotplot(enrich)
```
