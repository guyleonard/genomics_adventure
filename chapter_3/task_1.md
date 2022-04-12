# Chapter 3 - Assembly of Unmapped Reads
In this chapter of our adventure we will continue the analysis of a strain of ​*E. coli*.​ In the previous chapter we cleaned our data, checked QC metrics, mapped our data to a reference, obtained a list of variants and identyfied an overview of any missing regions.

Now, we will examine those reads which did not map to the reference genome. We want to know what these sequences represent. Are they novel genes, plasmids or just contamination? To do this we will extract the unmapped reads, evaluate their quality, prepare them for de-novo assembly, assemble them using SPAdes, generate assembly statistics and then produce some annotation via Pfam, BLAST and RAST. :cold_sweat:

## Task 1 - Extract the Unmapped Reads
Let's make sure we are in the right place to start out adventure and refresh our memories of what is in this directory.
```bash
cd ~/workhsop_materials/genomics_adventure/sequencing_data/ecoli

ls -lath
```

We are going to produce some new work, so remember we should keep things nice and tidy. Let's make a new directory and move to that folder.
```bash
mkdir unmapped_assembly

cd unmapped_assembly
```

Now we will use the [bam2fastq](https://gslweb.discoveryls.com/information/software/bam2fastq) :mag: program to extract from the BAM file just those reads which did NOT map to the reference genome. The
bam2fastq program has a number of options, most of which are self-explanatory. Whilst this tool has been discontinued, it still performs useful functions. However, you might like to learn the tool from the [Picard](http://picard.sourceforge.net/) :mag: package called SamToFastq at another time, which should perform a similar function.
```
bam2fastq --no-aligned -o unaligned#.fastq ../sequencing_data/ecoli/mapping_to_reference/ecoli_mapped_namesort_fixmate_sort_markdup.bam
```

