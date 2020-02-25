# Task 15 Locating Genes that are Missing Compared to the Reference
We can use a command from the [BEDTools](http://bedtools.readthedocs.org/en/latest/) package  to identify annotated genes which are not covered by reads across their full length.

For example, lets try the following:
```bash
bedtools coverage \
-a ~/workshop_materials/genomics_adventure/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.gff \
-b ecoli_mapped_namesort_fixmate_sort_markdup.bam > gene_coverage.txt
```

Oh no! :scream: What happenend? Did you get a message that says "Killed"? It's okay, Don't Panic! When you see a message similar to this it usually means that there isn't enough RAM available for this task - and so the Operating System has 'killed' the operation (sounds extreme I know) so that you can continue to use it. Limitations like these may well happen throughout your bioinformatics adventure! :weary: You've just got to learn to navigate around them, or go and do some mindfullness techniques for a while or both! :relaxed: You can do it! :muscle:

But we still want to see some results, right? So, let's try and subsample the BAM file to 50% of the data! Hopefully that should work! We can do that with the 'samtools view' command and the '-s' option. We need to give it a 'seed', in this case '42' :wink:, and a decimal representation of the percentage we want to subsample - 50% - so, '.5' and put them together into one as '42.5'.

```bash
samtools view -bs 42.5 \
ecoli_mapped_namesort_fixmate_sort_markdup.bam \
> ecoli_mapped_namesort_fixmate_sort_markdup_subsampled.bam
```

Now we can return to the 'bedtools coverage' command:
```bash
bedtools coverage \
-a ~/workshop_materials/genomics_adventure/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.gff \
-b ecoli_mapped_namesort_fixmate_sort_markdup_subsampled.bam > gene_coverage.txt
```

This should now only take a minute or so. Take a look and you will see that the output contains one row per annotated gene, whereby the 13th column contains the proportion of the gene that is covered by reads from our sequencing. 1.00 means the gene is 100% covered and 0.00 means there is no coverage. Of course, in our case, it could be anywhere up to 50% more coverage than is reported due to our subsampling, but that doesn't matter in order to find any missing genes!

Therefore, if we 'sort' our data by the 13th column we can see which genes are missing:
```bash
sort -t $'\t' -g -k 13 gene_coverage.txt | less -S
```
Can you find these missing genes in IGV?

# Congratulations! :tada:
That concludes the first part of the adveture. You have successfully, QC'd, filtered, mapped and analysed a whole bacterial genome! Well done! :trophy: :trophy: Time for a well deserved break! :coffee: :cookie:

In the next installment we will be looking at how to extract and assemble unmapped reads. This will enable us to look at material which may be present in the strain of interest but not in the reference sequence. See you soon!


# [Chapter 3]()
