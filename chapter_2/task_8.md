# Task 8 - Remove Suspected PCR Duplicates
'samtools' can do a reasonably good job of removing potential PCR duplicates, especially when using paired-end reads (refer back to the first part of this adventure if you are unsure what PCR duplicates are).

In previous versions of 'samtools' it contained a command named 'rmdup' which removed suspected PCR duplicates, indeed you may see it still used in a lot of tutorials in your google adventures. However, the command has now been *deprecated* and so it may be removed in future versions of the 'samtools' program. Therefore, it is best not to use it, in fact the authors of 'samtools' suggest that you do not.

As there were a number of issues with the 'samtools rmdup', it has been replaced with the 'samtools markdup' command. This command marks duplicates, rather than removing them entirely - however, you will need to do some several rounds of sorting to make this work. It's more steps than 'rmdup', yes, but it is a better method of analysing your data and allows you to track what information is being removed.

## samtools fixmate and samtools markdup
Firstly, we need to sort the BAM file by 'read name', (i.e., the QNAME field in a SAM file), rather than by chromosomal coordinates (as we did in the previous task), because the next command requires this style of sorting.
```bash
samtools sort -n -o XXX_namesort.bam XXX.bam
```

Next we need 'samtools' to add 'ms' and 'MC' tags for the 'samtools markdup' program to use later. This is because BWA can sometimes leave unusual FLAG information on SAM records, it is helpful when working with many tools to first clean up read pairing information and flags.
```bash
samtools fixmate -m XXX_namesort.bam XXX_fixmate.bam
```

We are nearly ready to remove those duplicates, but first we need to re-sort our BAM file by chromosomal/position coordinates
```bash
samtools sort -o XXX_positionsort.bam XXX_fixmate.bam
```

Finally we can mark the duplicates, and in this case actually remove them with the '-r' option.
```bash
samtools markdup -r XXX_positionsort.bam XXX_markdup.bam
```

There's a lot to take in here, so go back over the steps a few times to make sure you understand what is required. It's a bit of a laborious process, and some of the manipulations of the files may not seem as easy as one command, but at least it is a process that you can follow in your other analyses. Different programs can do this in different ways, for example the  '[MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)' :mag: in the 'Picard' suite of tools, but we won't go in to that here.

Now that we have a position sorted, mate fixed and removed duplicates BAM file, we are nearly ready for some actual analysis. Are you excited yet? :upside_down_face:

It is generally good practice to keep your intermediary files whilst you are continuing your analysis, but in this case we don't really need the three extra files we have created during this process, we are only interested in the final file 'XX_markdup.bam'. Therefore, you can remove the others:
```bash
rm XXX_namesort.bam
rm XXX_fixmate.bam
rm XXX_positionsort.bam
```

## One Last Thing...
Most programs used to 'view' BAM formatted data require an 'index file' to locate the reads mapping to a particular location quickly. You can think of this as an index in a book, telling you where to go to find particular phrases or words. We'll use the 'samtools index' command to do this.
```bash
samtools index XXX_markdup.bam
```

# [Task 9]()
