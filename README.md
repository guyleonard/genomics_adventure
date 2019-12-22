# genomics_workshop

## Set Up
### Programs

```
conda update -n base conda
conda create --name genomics_workshop
conda activate genomics_tutorial

conda install -c bioconda sra-tools #perl, zlib
conda install -c bioconda fastqc #openjdk
conda install -c bioconda ea-utils #blas
conda install -c bioconda bwa
conda install -c bioconda samtools #htslib
conda install -c bioconda qualimap
conda install -c bioconda igv
conda install -c bioconda igvtools
conda install -c bioconda bcftools
conda install -c bioconda vcftools
conda install -c bioconda bedtools
conda install -c bioconda spades #pip, python-3.8
conda install -c bioconda blast 
conda install -c bioconda emboss 
conda install -c bioconda pfam_scan #hmmer



```

## Deprecated / Outdated / No Conda?
 * samtools rmdup
 * ea-utils - fastq-mcf
 * bam2fastq
 * quast
 * filter_low_coverage_contigs.pl
 * reduce_fasta_10x.pl
 
