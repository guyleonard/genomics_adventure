# Task 10 - QualiMap
['QualiMap'](http://qualimap.bioinfo.cipf.es/) is a program that summarises the mapped alignments in much more detail than the mapping stats file we produced previously. It’s a technical tool which allows you to assess the sequencing for any problems and biases in the sequencing and the alignment rather than a tool to deduce biological features.

There are a few options to the program, we want to run 'bamqc'
```bash
qualimap bamqc
```

To generate a 'QualiMap bamqc' report, you can run:
```bash
qualimap bamqc -outdir bamqc \
-bam ecoli_mapped_namesort_fixmate_sort_markdup.bam \
-gff ~/workshop_materials/genomics_adventure/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.gff
```

After a couple of minutes, this will create a subfolder called 'bamqc'. Change into this directory with 'cd' and view the '.html' file with your favourite browser, e.g.:
```bash
firefox qualimapReport.html
```

There is a lot of information in this report to explore, so we will give just a few highlights here, you are welcome to read the manual (RTM) for further information or discuss with your peers and TAs too.

## Coverage Across Reference
This shows the number of reads that 'cover' each section of the genome. The red line shows a rolling average around 175x - this means that on average every part of the genome was sequenced 175x. It is important to have sufficient depth of coverage in order to be confident that any features you find in your data are real and not a result of sequencing errors. What do you think the regions of low/zero coverage correspond to?

![Coverage Across Reference](https://github.com/guyleonard/genomics_adventure/blob/e219def38a300ab13dba1aa839c27b3fa7909c27/chapter_2/images/genome_coverage_across_reference.png)

## Insert Size Histogram
The Insert Size Histogram displays the range of sizes of the DNA fragments. It shows how well your DNA was size selected before sequencing. Note that the 'insert' refers to the DNA that was inserted between the sequencing adaptors, so equates to the size range of the DNA that was used. In this case we have 150bp paired end reads and our insert size varies around 120bp bases - so there should only be a small gap between the reads that were not sequenced.

![Insert Size Histogram](https://github.com/guyleonard/genomics_adventure/blob/e219def38a300ab13dba1aa839c27b3fa7909c27/chapter_2/images/genome_insert_size_histogram.png)

You can have a look at some of the other graphs produced if you like and refer to the 'genome_results.txt' for similar statistics, but let's move on to a more interesting way of looking at our alignments...

# Go to [Task 11](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_2/task_11.md)
