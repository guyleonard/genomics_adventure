# Task 8 - Remove Suspected PCR Duplicates
Especially when using paired-end reads, samtools can do a reasonably good job of removing potential PCR duplicates (see the first part of this workshop if you are unsure what this means).

Previous versions of samtools contained a command 'rmdup' to automatically do this, you may see it still used in a lot of tutorials in your google searches. However, the command is now *deprecated* and may be removed at some point from future versions of the 'samtools' program. Therefore it is best not to use it, in fact the authors of samtools suggest you do not, unless you wish to compare a new analysis with an old analysis.

There were a number of issues with the 'samtools rmdup' and so it has been replaced with the 'samtools markdup' command. This marks duplicates, rather than removing them entirely - however, you will need to do some specific rounds of sorting to make this work. It's more steps, yes, but it is a better method of analysising your data that allows you to track what information is being removed and why.

## samtools fixmate and samtools markdup
First we need to sort the BAM file by 'read name', (i.e., the QNAME field in a SAM file), rather than by chromosomal coordinates, because the next command requires this style of sorting.
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

Now we have a position sorted, removed duplicates BAM file, nearly ready for actual ananlysis. Are you excited yet?

It is generally good practice to keep your intermediary files, but in this case we don't need the three extra files we have create here, we only need the final file 'XX_markdup.bam'. You can remove the others:
```bash
rm XXX_namesort.bam
rm XXX_fixmate.bam
rm XXX_positionsort.bam
```




