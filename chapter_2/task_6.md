# Task 6 - Mapping Reads to the Indexed Reference Sequence
Now we can begin to align read 1 and read 2 to the reference genome. First of all change back into raw reads directory where we made out trimmed and QC reads, and then create a subdirectory to contain our remapping results.
```bash
cd ~/genomics_adventure/sequencing_data/ecoli
mkdir mapping_to_reference
cd mapping_to_reference
```

## Mapping
Have a look at the output from typing 'bwa mem'. You should see something like this:

[IMAGE]

The basic format of the command is:

>Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]

From this we can see that we need to provide BWA with a set of FASTQ files containing the raw reads (denoted by
<in.fq> and [in2.fq]) to align to a reference file (unhelpfully this is listed as <idxbase>). There are also a number of options. The most important are the maximum number of differences in the seed ('-k' i.e. the first 32 bp of the sequence vs the reference) and the number of processors the program should use ('-t').
  
Our reference sequence is in:
>~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna

Our filtered reads are in:
>~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/XX_val1.fq.gz and
>~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/XX_val2.fq.gz

So, in order to align our paired reads using multi threading and output to a file XXX.sam:
```bash
bwa mem -t 4 ~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna ~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/XX_val1.fq.gz ~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/XX_val2.fq.gz > XXX.sam
```

This will take about 5 minutes to complete. There will be quite a lot of output but the end should look like:

[IMAGE]

## Viewing a SAM File- Basic
Once the alignment is complete, list the directory contents and check that the alignment file 'XXX.sam' is present. The raw alignment is stored in what is called 'SAM' format (Simple AlignMent format). It is in plain text and you can view it if you wish using the 'less', 'head', or 'tail' commands. Do not try to open the whole file in a text editor as you will likely run out of memory!

Each alignment line has 11 mandatory fields for essential alignment information such as mapping position, and a variable number of optional fields for flexible or aligner specific information. For further details as to what each field means see the PDF [here](http://samtools.sourceforge.net/SAM1.pdf) :mag:.


# [Task 7]()
