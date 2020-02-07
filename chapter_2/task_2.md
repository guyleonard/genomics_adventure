# Task 2 -  Sequence Data Quality Control & Adaptor Trimming
In the next set of tasks we will filter the read data to ensure that any low quality reads are removed, and that any sequences containing adaptor sequences are either trimmed or removed altogether. Adaptors are not useful in an assembly or during mapping, and low quality reads can impair a genome assembler's ability to build contiguous sequences from reads, or suggest SNPs where there are none.

Note: Typically when submitting Illumina data to NCBI or EBI you would submit the raw unfiltered data, so don't delete your original FASTQ files!

We are going to use the program '[Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)' :mag: to adaptor trim and QC our data. There are many other programs that can do this, see below for examples, but this one is our current favourite - because it is actually a script that 'wraps' two programs together: '[cutadapt](https://cutadapt.readthedocs.io/en/stable/)' :mag: and '[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)' :mag:.

## Other QC Programs
A list (by no means exhaustive) of some of the other most common adaptor trimming and QC programs:

 * [fastp](https://github.com/OpenGene/fastp)
 * [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
 * [flexbar](https://github.com/seqan/flexbar)
 * [adapterremoval](https://github.com/MikkelSchubert/adapterremoval)
 * [fastq-mcf](https://expressionanalysis.github.io/ea-utils/)

Commands for installing these programs are shown below, however you do not need to complete them for this workshop. If you have extra time you may wish to journey back here at a later date to try them out.

PS - There are also programs for removing or demultiplexing special barcode data from your reads, but they are not covered in this tutorial.

PPS - Long reads require a whole suite of different programs - especially Oxford Nanopore - due to the different technologies and errors involved in base calling.
 
```bash
conda install -c bioconda fastp
conda install -c bioconda trimmomatic
conda install -c bioconda flexbar
conda install -c bioconda adapterremoval
conda install -c bioconda ea-utils
```

## Running Trim Galore!
'Trim Galore!' can be run by typing the program name 'trim_galore' into the terminal. However, it won't display anything useful to start with, so you will have to give it some instructions, e.g. '--help', this will show what options are available for you to start.
```bash
trim_galore --help
```
You will see something similat to this:

[IMAGE]

There is always a lot of information to take in when you look at the help pages of different programs, but don't worry these walls of text will become your new best friends soon enough! :handshake:

### Input Information
First of all, let's ask the question - What do we know about our data?

Here a few things to start:

 1. It is Illumina HiSeq
 2. It is in two FASTQ files - so we know it is probably paired-end data. 
 3. It is using Illumina 1.9+ enconding.
 
This is the minimal set of information we need to start with (in this case FastQC told us most of the information, but you would also know much more than this when you submit your samples to your sequencing facility). So, we could start with something like this (do not type this yet):
```bash
trim_galore --illumina --paired --phred33 read_1.fastq.gz read_2.fastq.gz
```

However, it is very useful to add some options to control our output too. For example, we should also run 'FastQC' on the trimmed reads, and we want the output files to be in a 'zipped' format, and we should also take advantage of the multi-threading capabilities of our computer (to speed up the processing). 

So we can actually type:
```bash
trim_galore --paired --fastqc --gzip --cores 4 --length 100 read_1.fastq.gz read_2.fastq.gz
```

Running the above command should take roughly 5 minutes. Read on below whilst you are waiting, or try and take in some of the output messages, what is the program telling you? :hourglass_flowing_sand:

Notice that I have removed the '--illumina' and '--phred33' options, this is because 'Trim Galore!' is pretty smart and will try to guess the encoding and adaptor types for you :crossed_fingers: (but if you are 100% sure, then you can leave them in). It will now also output the results in a 'gzip' format rather than plain text, and run 'FastQC' on those files automaically.

Most programs also have a bunch of default parameters which have been set for you, in this case the '-q' option is one that sets the quality trim option at a Phred score of 20. Luckily, this is sufficient for our needs here, but they may differ depending on your input libraries.

Did you notice the "--length 100" option was also enabled? The default is set at 20bp and that is very small - especially when we start with 150bp reads, if any reads were trimmed to this 20bp it is debatable how useful they would be in an assembly or mapping situation. Let's try and keep all our reads > 100bp in this case. It is likely that in the 'real world' :tm: you will end up trying several rounds of trimming until you are happy with the results. Remember, science is a process of trials and errors. :sweat:

### Output Files
Now let's see what the program has produced! Returning to your terminal, you can list the contents of the directory, and you should see something similar to this:

[IMAGE]

There are several sets of output files:
 1. two processed zipped 'fastq' files (read_1_val_1.fq.gz and read_2_val_2.fq.gz), containing the trimmed reads,
 2. two '.html' files which contain a visual report - much like the GUI,
 3. two '.txt' files which contain similar information as above, just without the graphs,
 4. and two '.zip' files which contain the same reports as the '.html' files.

You will notice that the 'read_2_val_2.fq.gz' file is a little bit smaller than 'read_1_val_1.fq.gz'. Why do you think this might be?

Now you should count the lines in all the files:
```bash
zcat read_1_val_1.fq.gz | wc -l

zcat read_2_val_2.fq.gz | wc -l
```

Although the reads have been trimmed differently - the number of reads in the 'read_1' and 'read_2' files are identical (2,951,293). This is required for all the tools we will use to analyse paired-end data.

Now you should check the quality scores and sequence distribution from the 'FastQC' outputs. Open the 'html' file in Firefox or your favourite browser. You should notice very little change (since comparatively few reads were filtered). However, you should notice a significant improvement in quality and the absence of adaptor sequences.

## [Task 3]()
