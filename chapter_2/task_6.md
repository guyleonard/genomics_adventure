# Task 6 - Mapping Reads to the Indexed Reference Sequence
Now we can begin to align 'read_1' and 'read_2' to the reference genome. First of all change back into raw reads directory where we made the trimmed and QC reads, and then create a subdirectory to contain our mapping results.
```bash
cd ~/workshop_materials/genomics_adventure/sequencing_data/ecoli

mkdir mapping_to_reference

cd mapping_to_reference
```

## Mapping
Have a look at the output from typing 'bwa mem'. You should see something like this:

[IMAGE]

The basic format of the command show as:

>Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]

From this we can see that we need to provide BWA with a set of FASTQ files containing the raw reads (denoted by
'\<in.fq>' meaning required and '[in2.fq]' as optional) and a reference file (unhelpfully this is listed as '\<idxbase>') and any other options we wish to change. The most important options that we should take notice of, are the maximum number of differences in the seed ('-k' i.e. the first 32 bp of the sequence vs the reference) and the number of threads/processes the program should use ('-t').
  
Our reference sequences are in the file:
>~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna

Our filtered reads are in the files:
>~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/reqd_1_val_1.fq.gz

>~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_2_val_2.fq.gz

So, in order to align our paired reads using multi threading and output to a file ecoli_mapped.sam:
```bash
bwa mem -t 4 \
~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_1_val_1.fq.gz \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_2_val_2.fq.gz \
-o ecoli_mapped.sam
```

NB - We are using the explicit '-o filename.sam' as our output but you could also use a simple 

This will take a couple of minutes to complete. You will be quite a lot of output on your terminal, but in the end it should look something similar to this:

[IMAGE]

Congratulations, you have just mapped reads to your first reference genome!! :trophy:

## SAM File Manipulation - Basic
Once the alignment is complete, list the directory contents and check that the alignment file 'ecoli_mapped.sam' is present. The raw alignment is stored in what is called 'SAM' format (Simple AlignMent format). It is a plain text file, and so you can view it if you wish using the 'less', 'head', or 'tail' commands. Do not try to open the whole file in a text editor as you will likely run out of memory, they can get very big!

Each alignment line has 11 mandatory fields for essential alignment information such as mapping position, and a variable number of optional fields for flexible or aligner specific information. For further details as to what each field means see the PDF [here](http://samtools.sourceforge.net/SAM1.pdf) :mag:.


# [Task 7]()
