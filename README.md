# A Genomics Adventure
Thanks to:
 * Konrad Paszkiewicz
 * Josie Paris
 * Sophie Shaw
 * Workshop on Genomics - evomics.org

## Getting Started
We will need to install some common software packages, and download some data to use for this tutorial. For the software we will use a program called 'conda', it will allow us to easily install lots of common bioinformatics software in a special 'environment' without the need for root/sudo access. For the data, we will use several methods - explained below.

### Software
You may copy and paste the commands below.
```
# Make sure we are up to date
conda update -n base conda
# Create our environment
conda create --name genomics_adventure
# Activate our environment
conda activate genomics_adventure
# Install the software
conda install -c bioconda bcftools bedtools blast bwa ea-utils emboss fastqc igv igvtools pfam_scan qualimap samtools seqtk spades sra-tools vcftools
```

### Data
#### Reference Data
We will be working with two different bacterial species for this adventure; Escherichia coli & Vibrio parahaemolyticus, they are two relatively small genomes, but the techniques you will learn here can be used with smaller & larger, and Eukaryotic genomes too.

```
# Create a directtory to store our data
mkdir reference_sequence && cd reference_sequence

# Download the E. coli reference genome in FASTA and GFF formats
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz

# Download the Vibrio reference genome in FASTA and GFF formats
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/196/095/GCF_000196095.1_ASM19609v1/GCF_000196095.1_ASM19609v1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/196/095/GCF_000196095.1_ASM19609v1/GCF_000196095.1_ASM19609v1_genomic.gff.gz

# Unzip the files
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
 
