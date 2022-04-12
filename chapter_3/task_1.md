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

We want to extract all of the reads that do NOT map to the assembly. Luckily, in the SAM/BAM format there is a special 'bitwise flag' or code that identifies how the reads and their read-mates are aligned to a reference. They can be quite confusing, and there are many combinations. We will have a look at them now by viewing the first five lines our previously made BAM file.
```bash
samtools view ../sequencing_data/ecoli/mapping_to_reference/ecoli_mapped_namesort_fixmate_sort_markdup.bam | head -n 5
```

You should see a bunch a of text, numbers and sequence data on your screen. Don't panic. It is arranged in columns separated by a tab, and each row is one read. At this time we are only really interested in the second column (the flag), you can look up the rest [here](https://en.wikipedia.org/wiki/SAM_(file_format)#Format):mag:. You should see a number like "2147" on the first row. On its own this number doesn't tell us too much, but we can look up what it means [here](https://broadinstitute.github.io/picard/explain-flags.html). You can see that the tool tells us that this read is mapped, and that its read-mate is also mapped. It also tells us that it is mapped to the reverse strand and is the first read of the pair to be mapped.

Is this a read we are looking for? Using the "Decoding SAM flags" tool, can you figure out what flag number we want? Remember, We need reads that are unmapped and where their mates are also unmapped. Click below to reveal the answer.

<details>
  <summary>Did you guess correctly?</summary>
  The answer is "12".

  Some of you may have guessed 4 or 8. That's okay, but remember we wanted paired reads, not just one or the other of the pairs. The more astute of you will notice that 4 + 8 = 12, the flags are summative. However, they are not ontological - e.g. 12 will get unmapped pairs, but not the others associated with 4 (read unmapped) and 8 (mate unmapped).
    
</details>




Now we will use the [bamtofastq](https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html):mag: program from the [bedtools](https://bedtools.readthedocs.io/en/latest/index.html):mag: There are other tools that you can use too, for example in the [Picard](http://picard.sourceforge.net/):mag: program there is a tool called SamToFastq which provides a similar function. But we will not use this today. 
```
bedtools bamtofastq 
```


[IMAGE]


