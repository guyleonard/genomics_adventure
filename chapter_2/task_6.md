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
'\<in.fq>' meaning required and '[in2.fq]' as optional) and a reference file (unhelpfully this is listed as '\<idxbase>') and any other options we wish to change. The most important options that we should take notice of, are the output '-o' and the number of threads/processes the program should use ('-t').
  
Our reference sequences are in the file:
>~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna

Our filtered reads are in the files:
>~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/reqd_1_val_1.fq.gz

>~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_2_val_2.fq.gz

So, in order to align our paired reads using multi-threading, and output to a file called ecoli_mapped.sam, we can type:
```bash
bwa mem -t 4 \
~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_1_val_1.fq.gz \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_2_val_2.fq.gz \
-o ecoli_mapped.sam
```

NB - Here we are expllicity using the '-o filename.sam' option as our output, but you could also use a simple '> filename.sam' redirect instead (which might be more useful for inclusion in more advanced pipelines).

This process will take a couple of minutes to complete. You will see quite a lot of output on your terminal, but in the end it should look something similar to this:

[IMAGE]

Congratulations, you have just performed your first mapping of reads to a reference genome!! :trophy: But that only get's us so far... There's always more to do.

## SAM File Manipulation - Basic
Once the alignment is complete, list the directory contents and check that the alignment file 'ecoli_mapped.sam' is present. The raw alignment is stored in what is called 'SAM' format (Simple AlignMent format). It is a plain text file, and so you can view it using the 'less', 'head', or 'tail' commands, for example. It is best to not open the whole file in a text editor, as you will likely run out of memory, they can get very big! In our case it is about 2.2GB.

Each alignment line has 11 mandatory fields for essential alignment information including mapping position, and a variable number of optional fields for flexible or aligner specific information. For further details as to what each field means see the PDF [here](http://samtools.sourceforge.net/SAM1.pdf) :mag:. If you are on the Workshop on Genomics course then Mike Zody will have explained these files to you in gory detail! We will also explore this more in the next task.


# Go to [Task 7](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_2/task_7.md)
