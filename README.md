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
conda install -c bioconda seqtk
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

### Data

#### Reference Data
```
mkdir reference_sequence && cd reference_sequence
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/196/095/GCF_000196095.1_ASM19609v1/GCF_000196095.1_ASM19609v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/196/095/GCF_000196095.1_ASM19609v1/GCF_000196095.1_ASM19609v1_genomic.gff.gz
gunzip *.gz

```

## Deprecated / Outdated / No Conda?
 * samtools rmdup
   * samtools fixmate -m | markdup -r
 * ea-utils - fastq-mcf
 * bam2fastq
   * samtools bam2fq or bedtools bamtofastq
 * quast
 * filter_low_coverage_contigs.pl
 * reduce_fasta_10x.pl
 
