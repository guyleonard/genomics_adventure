# Chapter 3 - Assembly of Unmapped Reads
In this chapter of our adventure we will continue the analysis of a strain of ​*E. coli*.​ In the previous chapter we cleaned our data, checked QC metrics, mapped our data to a reference, obtained a list of variants and identyfied an overview of any missing regions.

Now, we will examine those reads which did not map to the reference genome. We want to know what these sequences represent. Are they novel genes, plasmids or just contamination? To do this we will extract the unmapped reads, evaluate their quality, prepare them for de-novo assembly, assemble them using SPAdes, generate assembly statistics and then produce some annotation via Pfam, BLAST and RAST. :cold_sweat:

## Task 1 - Extract the Unmapped Reads
Let's make sure we are in the right place to start out adventure and refresh our memories of what is in this directory.
```bash
cd ~/workshop_materials/genomics_adventure/sequencing_data/ecoli

ls -lath
```

We are going to produce some new work, so remember we should keep things nice and tidy. Let's make a new directory and move to that folder.
```bash
mkdir unmapped_assembly

cd unmapped_assembly
```

We want to extract all of the reads that do NOT map to the assembly. Luckily, in the SAM/BAM format there is a special 'bitwise flag' or code that identifies how the reads and their read-mates are aligned to a reference. They can be quite confusing at first, and there are many combinations. We will have a look at them now by viewing the first five lines our previously made BAM file.
```bash
samtools view ../sequencing_data/ecoli/mapping_to_reference/ecoli_mapped_namesort_fixmate_sort_markdup.bam | head -n 5
```

You should see a bunch a of text, numbers and sequence data on your screen. Don't panic, this is just how the SAM format looks. It is arranged in columns separated by a tab, and each row is one read. At this time we are only really interested in the second column (the flag), you can look up the meaning for the rest [here](https://en.wikipedia.org/wiki/SAM_(file_format)#Format):mag:.

You should see a number like "2147" on the first row. On its own this number doesn't tell us too much, but we can use a tool to look up what it means [here](https://broadinstitute.github.io/picard/explain-flags.html). You can see that the tool tells us that this read is mapped, and that its read-mate is also mapped. It also tells us that it is mapped to the reverse strand, is the first read of the pair to be mapped, and that it is a supplementary alignment. So, ths is not a read we are looking for right now!

Now, using the "Decoding SAM flags" tool, can you figure out what flag number we need for reads that are unmapped and where their mates are also unmapped? Click below to reveal the answer.

<details>
  <summary>Did you guess correctly?</summary>
  The answer we were looking for is "12". :one::two:

  But some of you may have guessed 4 or 8 or even 13 or 15 or higher! :confused: So, why is it twelve?

  Let's talk about the "bit-flag" briefly. Brace yourselves! :grimacing: The number values we see are actually the summed positions of a binary code representing a set of outcomes for the reads and their pairs. Woah! Breathe. :nose: For example, we could have the binary code of "0000000100", which is equivalent to a decimal "4". Why? Well, each position from the right of the binary code can be represented in decimal as 1, 2, 4, 6, 8, 16...etc. So, a '1' in the third position from the right in binary is equivalent to a decimal "4". You can then see how this matches to each of the outcomes in the "Decoding SAM flags" tool, e.g. selecting the third box is equivalent to a value of 4! Easy huh!? :muscle:

  But why 12 and not 13 or some other combination? Remember we wanted "read unmapped" (4) AND "mate unmapped" (8), so selecting both gives us "12" (or 0000001100 in binary), that's all we need. Nonetheless, some of you may have also decided to include either "read paired" (1) or "read mapped in proper pair" (2) increasing the value. Well, the latter is not useful as we are looking for unmapped reads only. Secondly, even though "read paired" is what we are looking for it is not often a flag that is used on its own when reads are unmapped - you weren't to know. But, you can also think of what 13 represents as a subset of 12, and as we want to get all the reads we should use the lower number! :ok_woman:    
</details>

Now that we have identified the corect bit flag, we can go ahead with filtering of the BAM file from Chapter 2, here using samtools.
```bash
samtools view -b -f 12 ../sequencing_data/ecoli/mapping_to_reference/ecoli_mapped_namesort_fixmate_sort_markdup.bam -o unmapped.bam
```

This command outputs a BAM file "-b" and filters only those with a corresponding bit flag of "-f 12".

Have a look at some of the content of this new BAM file.
```bash
samtools view unmapped.bam | head -n 5
```

Oh no! Why do all the values say 77 and 141? :sob: I thought you said it was 12? :scream: Have a think why this might be. Remember, larger values are subsets of smaller values.

Okay, now we are happy again! :smiley: We can continue on our journey, remember we were hoping to extract the unmapped reads for assembly. However, we need our reads in FASTQ format, but right now they are trapped in BAM format.

To convert them we will use the [bamtofastq](https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html):mag: program from the [bedtools](https://bedtools.readthedocs.io/en/latest/index.html):mag: package. There are other tools that you can use too, for example in the [Picard](http://picard.sourceforge.net/):mag: package there is a tool called SamToFastq which provides a similar function. But we will not use this today. 
```
bedtools bamtofastq -i unmapped.bam -fq unmapped_r1.fastq -fq2 unmapped_r2.fastq
```

Nicely done! Now lets head over to Task 2.

# Go to [Task 2](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_3/task_2.md)
