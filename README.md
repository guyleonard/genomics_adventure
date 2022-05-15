# A Genomics Adventure
This tutorial is based on a workshop that was written many years ago by Konrad Paszkiewicz, whilst he was the head of Exeter University's Sequencing Service. Over the years it has been extensivley updated and modified by several colleages (listed below), taught at Exeter Universtity, and also at the [evomics.org](https://evomics.org/) Workshop on Genomics since 2014. It has now been reproduced here, in a slightly new form - again with many updates.

Many thanks to:
 * [Konrad Paszkiewicz](https://scholar.google.com/citations?user=yrHDETIAAAAJ&hl=en), CTO Hummingbird Biosciences.
 * [Sophie Shaw](https://scholar.google.se/citations?user=_K3aFRYAAAAJ&hl=en), Bioinformatician, All Wales Medical Genomics Service. 
 * [Josie Paris](https://www.josieparis.com/), Postdoctoral Research Associate, University of L'Aquila, Italy.
 * Workshop on Genomics, [evomics.org](https://evomics.org/)

## A Few Notes on Styles
Throughout this adventure you will see various styles of text. Mostly they will be in the form of the main story - similar to what you are reading now - but some will be links to further and expanded reading (indicated with a magnifying glass - :mag: - this is mostly optional reading for the super curious :nerd_face:). However, there will also be commands for you to type. These will look something like:
```bash
# Comment line - do not type this
command -options
```

It will always be indicated when you need to run a command or just look at them. We also try to heavily discourage you from using 'copy and paste'. This is because we feel that you will learn much more by actually typing the commands yourself, and making mistakes. In fact, some of the commands are designed to intentionally fail if you copy and paste them. :stuck_out_tongue_closed_eyes: 

Often, you will see commands written in the style below, that is with a back-slash '\' at the end of lines, this is simply to break up the command for ease of reading. You would not normally type commands like this. You can see the actual command as you would type it below too.
```bash
# Large command that you can't fully see without scrolling
bwa mem -t 4 ~/genomics_adventure/ecoli/GCF_000005845.2_ASM584v2_genomic.fna ~/genomics_adventure/ecoli/XX_val1.fq.gz ~/genomics_adventure/ecoli/XX_val2.fq.gz > XXX.sam

# Much easier to read command
bwa mem -t 4 \
~/genomics_adventure/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
~/genomics_adventure/ecoli/XX_val1.fq.gz \
~/genomics_adventure/ecoli/XX_val2.fq.gz \
> XXX.sam
```

The tutorial (sorry adventure), like any good story, is designed to be long and you likely won't finish it in less than 3 hours or even in one sitting. Don't worry though, you can always come back to it anytime!

<p align="center">:dragon_face: Shall we begin? - Daenerys Targaryen :dragon_face:</p>

## Once Upon a Time...
Before our adventure begins we will need to install some common software, and download some data to help us on our way. For the software we will use a program called '[conda](https://docs.conda.io/en/latest/)':mag:, which will allow us to easily install lots of common bioinformatics software in a special 'environment' (think of it like a box :package:) without the need for admin access or other complications. For the data, we will download this from some publiic repositories.

NB - If you are attending the Workshop on Genomics then 'conda', the software and all data is already installed for you! So you can continue to the next section called ['Adventure Time'](https://github.com/guyleonard/genomics_adventure#adventure-time) :smiley:. Although you might like to read the next bit just for reference.

Firstly, we need to enter or create a directory called "workshop_materials" in your home directory and then clone this repository. All further commands will be run within the "genomics_adventure" directory.
```bash
cd workshop_materials
# or
mkdir workshop_materials && cd workshop_materials

# clone this repository
git clone https://github.com/guyleonard/genomics_adventure.git

# enter genomics_adventure
cd genomics_adventure
```

### Software
This section will create the 'environment' :package: in which we will be having our adventure, this allows us to keep all the software in one place for easy access and repeatability (e.g. you may wish to run different versions of software for other analyses, you can do that in other environments). We won't explore each of the programs that we will install right now, but the adventure will explain each as we get to them.

:squirrel: This time you may copy and paste, one-by-one, the commands below:
```bash
# Make sure we are up to date
conda update -n base conda

# Create our environment
conda create --name genomics_adventure python=3.7

# Activate our environment
conda activate genomics_adventure

# Install the software
conda install -c bioconda bcftools=1.12 bedtools blast bwa ea-utils emboss fastqc igv igvtools pfam_scan qualimap quast=5.0.2 samtools=1.12 seqtk spades sra-tools trim-galore vcftools
```

Make sure that the environment is manually activated everytime you come back to this adventure. You should see '(genomics_adventure)' before your normal terminal prompt. If it is not activated, use the 'activate' command from above.

### Data
We will need to retrieve several sets of data for our adventure, this is similar to how you may collate data for your own analyses.
 1) Sequence Data
  * Either directly from a Sequencing Service or from a public access database.
 2) Reference Data
  * If you are lucky to have a reference genome...
 3) Databases
  * PFam-A

We will be working with the bacterial species *Escherichia coli* as it is a relatively small genome (which makes it easy for the timings of this tutorial), but the techniques you will learn here can be applied to any smaller or larger, and/or Eukaryotic genomes too!

#### Sequencing Data
Back at your home institute you will likely retrieve your data from either the institute's sequencing service or a private outside provider. However, there is also a wealth :moneybag: of sequenced genomic data stored in publically accesible repositories such as NCBI's [SRA](https://www.ncbi.nlm.nih.gov/sra) or EMBL-EBI's [ENA](https://www.ebi.ac.uk/ena). These portals are also where you will be required to deposit your sequencing efforts during publication.

For this adventure we will be downloading and processing raw sequencing data. Please note that some sequencing services may provide trimmed or quality assessed reads as part of their standard service, however it is up to you whether you want to use that data directly or process the raw data yourself. Always ask: are their methods directly suited to your analysis needs?

The raw data that we will use for the *E. coli* genome is available from [NCBI](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=ERR2789854) or [EMBL-EBI](https://www.ebi.ac.uk/ena/data/view/ERR2789854) with the accession ERR2789854 (they are also archived at the [DDBJ-DRA](https://www.ddbj.nig.ac.jp/dra/index-e.html) too). This is the same data but it is mirrored between the sites, however each site has a different way of accessing the data. We will focus on NCBI and EMBL-EBI for now.

With NCBI you need to use a tool called '[fastq-dump](https://ncbi.github.io/sra-tools/fastq-dump.html)':mag:, which given an accession and several other options will download the 'fastq' data files - it can be notoriously tricky and difficult at times and has some issues with downloading PacBio data. You can give it a try below if you wish, however the EMBL-EBI downloads will be much faster for this tutorial, so we strongly suggest you start there.

At EMBL-EBI they provide direct links to the 'fastq' files that were submitted to the archive ("Submitted files (FTP)"), and so you can use tools such as 'wget' or 'curl' to retrieve them.

NB - These commands may take a little bit of time to complete depending on your connection (NCBI: ~15  minutes; EMBL-EBI: ~2 minutes), so you might want to skip ahead to the next chapter for some light reading about sequencing technologies and file formats whilst you wait... don't forget to come back soon!
```bash
# fastq-dump from NCBI - slow
fastq-dump --split-files --origfmt --gzip SRR857279

# or with wget from EMBL-EBI - faster
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/SRR857279/SRR857279_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR857/SRR857279/SRR857279_2.fastq.gz

# make a directory for them and move them there
mkdir -p sequencing_data/ecoli
mv *.gz sequencing_data/ecoli

#
chmod 444 sequencing_data/ecoli/*.gz

# Now do the same for Chapter 5's Pseudomonas data
cd .. && mkdir pseudomonas_gm41 && cd pseudomonas_gm41

# get the Illumina Data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR491/SRR491264/SRR491264_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR491/SRR491264/SRR491264_2.fastq.gz

# get the PacBio data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR104/006/SRR1042836/SRR1042836_subreads.fastq.gz

#
chmod 444 *.gz
```

#### Reference Data
We will access the reference data from the National Center for Biotechnology Information (NCBI), check out the links below by Ctrl or Cmd clicking the link to open in a new tab:

[*Escherichia coli* str. K-12 substr. MG1655](https://www.ncbi.nlm.nih.gov/genome/167?genome_assembly_id=161521)

There is a lot of information on this page, but the main pieces of information we are interested in are; the genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format, and the gene annotations in [GFF](https://en.wikipedia.org/wiki/General_feature_format) format. Can you see where these are? :eyes:

We will now download the data, as we are working with the command line we have already copied the links to the data below for you :slightly_smiling_face:. Using the '[wget](https://www.gnu.org/software/wget/)':mag: command we can download files directly from the web to our local dicrectory. The files are '*gzipped*', this means they are compressed to save space, it also allows us to make sure the data has not been corrupted during the transfer. We will also need to *unzip* them with the program '[gunzip](https://linux.die.net/man/1/gunzip)':mag:.

:squirrel: This time you may copy and paste, one-by-one, the commands below:
```bash
# Create a directory to store our data
mkdir reference_sequences && cd reference_sequences

# Download the E. coli reference genome in FASTA and GFF formats
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz

# Make an ecoli directory, and move the files there, and then unzip them
mkdir ecoli && mv *.gz ecoli && gunzip ecoli/*.gz

# Change write permissions, so that we can't edit them by accident
chmod -R 444 ecoli/*.fna
chmod -R 444 ecoli/*.gff
```

#### Databases
We will need to get the PFam-A database of Hidden Markov Models (HMMS) and an Active Site database from the [Pfam](https://pfam.xfam.org/) website. They are located in an ftp directory. Use the commands below. Make sure you are in the "genomics_adventure" directory.

```bash
# create a directory and a sub-directory and move there
mkdir -p db/pfam && cd db/pfam

# Download the HMMs and .dat files needed for Pfam-A
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz

# Uncompress the files
gunzip *.gz
```

When the downloads are finished you may continue on to the adventure by clicking the title below.

# [Adventure Time!](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_2/task_1.md)
