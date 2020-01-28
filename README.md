# A Genomics Adventure
Thanks to:
 * Konrad Paszkiewicz
 * Josie Paris
 * Sophie Shaw
 * Workshop on Genomics - evomics.org

## Getting Started
We will need to install some common software packages, and download some data to use for this tutorial. For the software we will use a program called 'conda', it will allow us to easily install lots of common bioinformatics software in a special 'environment' without the need for root/sudo access. For the data, we will use several methods - explained below.

### Software
This section will create the 'environment' in which we will be running the tutorial, this allows us to keep all our software and data in on place for easy access and repeatability (e.g. you may wish to run different versions of software in other analyses). We won't explore which programs that are installed right now, but the tutorial will explain each as we get to them. 

You may copy and paste one-by-one the commands below.
```bash
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
We will be working with two different bacterial species for this adventure; *Escherichia coli* & *Vibrio parahaemolyticus*, as they are two relatively small genomes (which makes it easy for the timing of our tutorial), but the techniques you will learn here can be applied to any smaller or larger, and/or Eukaryotic genomes too.

We will access the data from the National Center for Biotechnology Information (NCBI), check out the links below:
 * *[E. coli](https://www.ncbi.nlm.nih.gov/genome/167?genome_assembly_id=161521)*
 * *[V. parahaemolyticus](https://www.ncbi.nlm.nih.gov/genome/691?genome_assembly_id=167995)*

There is a lot of information on these pages, but the main pieces of information we are interested in are; the genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format, and the gene annotations in [GFF](https://en.wikipedia.org/wiki/General_feature_format) format.

We will now download the data, as we are working with the command line we have already copied the links to the data below for you. Using the '[wget](https://www.gnu.org/software/wget/)' command we can download files directly from the web to our local dicrectory. The files are '*gzipped*', this means they are compressed to save space, it also allows us to make sure the data has not been corrupted. We will need to *unzip* them with the program '[gunzip](https://linux.die.net/man/1/gunzip)'.

```bash
# Create a directory to store our data
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

[Continue the adventure...]

# Deprecated / Outdated / No Conda?
 * samtools rmdup
   * samtools fixmate -m | markdup -r
 * ea-utils - fastq-mcf
 * bam2fastq
   * samtools bam2fq or bedtools bamtofastq
 * quast
 * filter_low_coverage_contigs.pl
 * reduce_fasta_10x.pl
 
