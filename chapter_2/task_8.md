# Task 8 - Removing Suspected PCR Duplicates
'samtools' can do a reasonably good job of removing potential PCR duplicates, especially when using paired-end read data (refer back to the first part of this adventure if you are unsure what PCR duplicates are).

NB - In previous versions of 'samtools' it contained a command named 'rmdup' which removed suspected PCR duplicates, indeed you may see it still used in a lot of tutorials during your google adventures. However, the command has now been *deprecated* and so it may be removed in any future versions of the 'samtools' program. Therefore, it is best not to use it, in fact the authors of 'samtools' suggest that you do not.

There were a number of issues with the 'samtools rmdup' command, and it has now been replaced with the 'samtools markdup' command. This command marks duplicates, rather than removing them entirely - however, you will need to do two different rounds of sorting to make it work.

## samtools fixmate and samtools markdup
Firstly, we need to sort the BAM file by 'read name', (i.e., sort by the QNAME field in a SAM file), rather than by chromosomal coordinates (as we did in the previous task), simply because the command after requires this style of sorting.
```bash
samtools sort -n -o ecoli_mapped_namesort.bam ecoli_mapped.bam
```

Next we need 'samtools' to add 'ms' and 'MC' tags for the 'samtools markdup' program to use later on. This is because 'BWA' can sometimes leave unusual FLAG information in the 'SAM' file records. It is helpful when working with many downstream analysis tools to first clean up read pairing information and flags.
```bash
samtools fixmate -m ecoli_mapped_namesort.bam  ecoli_mapped_namesort_fixmate.bam
```

We are nearly ready to remove those duplicates, but first we need to re-sort our BAM file by chromosomal/position coordinates. Be very careful with the position of the output file, notice how that it is different to the previous command - bioinformatics tools can be frustrating like this at times...
```bash
samtools sort -o ecoli_mapped_namesort_fixmate_sort.bam ecoli_mapped_namesort_fixmate.bam
```

Finally we can mark the duplicates, and in this case actually remove them with the '-r' option. Again, watch that output file position!
```bash
samtools markdup -r ecoli_mapped_namesort_fixmate_sort.bam  ecoli_mapped_namesort_fixmate_sort_markdup.bam 
```

Each step should take about 2-3 minutes. There is a lot to take in here, so go back over the steps a few times in your head to make sure that you understand what is required. It's a bit of a laborious process, and some of the manipulations of the files may not seem as easy as one command, but at least it is a process that you can follow in your other analyses. Other programs can do the same step in different ways, for example the  '[MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)' :mag: in the 'Picard' suite of tools, but we won't go in to that here.

Now that we have a position-sorted, mate-fixed and removed-duplicates BAM file, we are very nearly ready for some actual analysis. Are you excited yet? :upside_down_face:

## Oh, One Last Thing...
Most programs used to 'view' BAM formatted data require an 'index file' in order to to locate the reads mapping to a particular location quickly. Remember that we did this for the reference genome file too. Again, you can think of this much like an index in a book. However, this time we'll use the 'samtools index' command to do this.
```bash
samtools index ecoli_mapped_namesort_fixmate_sort_markdup.bam
```

# [Task 9]()
