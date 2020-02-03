# Task 3 - Sequence Quality & Adaptor Trimming
From your terminal window, navigate to the 'sequencing_data' directory (you may be there already).
```bash
cd ~/genomics_adventure/sequencing_data
```
We are going to use the program '[Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)' :mag: to adaptor trim and QC our data. There are many other programs that can do this, see below for examples, but this one is our favourite. This is because it is actually a script that 'wraps' two programs together: '[cutadapt](https://cutadapt.readthedocs.io/en/stable/)' :mag: and '[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)' :mag:.

## Other QC Programs
A list (by no means exhaustive) of some of the other most common adaptor trimming and QC programs:

 * [fastp](https://github.com/OpenGene/fastp)
 * [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
 * [flexbar](https://github.com/seqan/flexbar)
 * [adapterremoval](https://github.com/MikkelSchubert/adapterremoval)
 * [fastq-mcf](https://expressionanalysis.github.io/ea-utils/)

Commands for installing these programs are shown below, however you do not need to complete them for this workshop. If you have extra time you may wish to journey back here at another time.

NB - Long read technologies may require a whole suite of different programs - especially Oxford Nanopore - due to the different technologies and errors involved in base calling.
 
```bash
conda install -c bioconda fastp
conda install -c bioconda trimmomatic
conda install -c bioconda flexbar
conda install -c bioconda adapterremoval
conda install -c bioconda ea-utils
```

## Running Trim Galore!
'Trim Galore!' can be run by typing the program name 'trim_galore' into the terminal. It won't display anything useful, so you will have to give it some instructions, e.g. '--help', to see what options are available for you to start.
```bash
trim_galore --help
```
You will see something similat to this:

[IMAGE]

There is always a lot of information to understand, but don't worry these walls of text will become your new best friends soon enough! :handshake:

### Input
First of all, let's ask the question - What do we know about our data? There's a few things:

 1. It is Illumina HiSeq
 2. It is in two fastq files
  * so we know it is paired end data. 
 3. It is using Illumina 1.9+ enconding.
 
This is probably the minimal set of information we should start with (mostly we knew this from FastQC), and as it happens this is  more or less the neccessary information that 'Trim Galore!' expects as input. So we could start with something like this:
```bash
trim_galore --illumina --paired --phred33 file_r1.fq.gz file_r2.fq.gz
```

However,  we may want to add some other options to control our output too. For example, we would like to run 'FastQC' on our trimmed reads, we want the output files to also be in a 'zipped' format, and we want to take advantage of the multi-threading capabilities of our computer. 

So we can actually type:
```bash
trim_galore --paired --fastqc --gzip --cores 4 file_r1.fq.gz file_r2.fq.gz
```

Notice that I have removed the '--illumina' and '--phred33' options, this is because 'Trim Galore!' is pretty smart and will now guess the encoding and adaptor type for you :crossed_fingers: (but if you are 100% sure, then you can leave them in). It will also output the results in a 'gzip' format and run 'FastQC' on them. But be careful, as other default parameters are set, such as the '-q' option which sets the quality trim option at a Phred score of 20. This is sufficient for our needs here, but may differ depending on your input libraries.




