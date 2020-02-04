# Task 10 - QualiMap
['QualiMap'](http://qualimap.bioinfo.cipf.es/) is a program that summarises the mapped alignments in much more detail than the mapping stats file we produced previously. It’s a technical tool which allows you to assess the sequencing for any problems and biases in the sequencing and the alignment rather than a tool to deduce biological features.

There are a few options to the program, we want to run 'bamqc'
```bash
qualimap bamqc
```

That will show you some helpful information on this command.

To generate a 'QualiMap bamqc' report, you can run:
```bash
qualimap bamqc -outdir bamqc -bam XXX_rmdup.bam -gff ~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.gff
```

this creates a subfolder called 'bamqc'. Change into this directory and view the 'html' file with your favourite browser, e.g.:
```bash
firefox qualimapReport.html
```

There is a lot in this report, so we will give just a few highlights here, you are welcome to read the manual for further information.

## Coverage Across Reference
This shows the number of reads that 'cover' each section of the genome. The red line shows a rolling average around 50x - this means that on average every part of the genome was sequenced 50X. It is important to have sufficient depth of coverage in order to be confident that any features you find in your data are real and not a result of sequencing errors. What do you think the regions of low/zero coverage correspond to?

[IMAGE]

## Insert Size Histogram
The Insert Size Histogram displays the range of sizes of the DNA fragments. It shows how well your DNA was size selected before sequencing. Note that the 'insert' refers to the DNA that was inserted between the sequencing adaptors, so equates to the size range of the DNA that was used. In this case we have XXX paired end reads and our insert size varies around XXX bases - so there should only be a small gap between the reads that was not sequenced.

[IMAGE]

You can have a look at some of the other graphs produced if you like and refer to the qualimap.txt for similar statistics, but let's move on to a more interesting way of looking at our alignments...

# [Task 11]()
