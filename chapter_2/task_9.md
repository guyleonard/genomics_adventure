# Task 9 - Mapping Statistics
Finally, we can generate some summary statistics!
```bash
samtools flagstat ecoli_mapped_namesort_fixmate_sort_markdup.bam > mappingstats.txt
```
This should only take a few seconds. Once complete you can view the 'mappingstats.txt' file using your favourite text editor or command (e.g. nano, vi, more or cat).

[IMAGE]

>5645403 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
4516 + 0 supplementary
0 + 0 duplicates
5531549 + 0 mapped (97.98% : N/A)
5640887 + 0 paired in sequencing
2820442 + 0 read1
2820445 + 0 read2
5503328 + 0 properly paired (97.56% : N/A)
5526970 + 0 with itself and mate mapped
63 + 0 singletons (0.00% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)


Here we can see we have XXX reads in total, none of which failed QC. XX % of reads mapped to the reference genome and XX % mapped with the expected XXX-XXX bp distance between them. XXXX reads could not have their read-pair mapped.

XX reads have mapped to a different chromosome than their pair (X has a mapping quality > X – this is a Phred scaled quality score much as we say in the FASTQ files). If there were any such reads they would likely be due to repetitive sequences (e.g phage insertion sites) or an insertion of plasmid or phage DNA into the main chromosome.

## Cleaning Up
It is generally good practice to keep your intermediary files whilst you are continuing your analysis, but in this case we don't really need some of the extra files we have created during this process, we are only really interested in the final file 'XX_markdup.bam'. Therefore, you can remove the others:
```bash
rm XXX.sam
rm XXX_sorted.bam
rm XXX_namesort.bam
rm XXX_fixmate.bam
rm XXX_positionsort.bam
```

In case you get asked if you are sure to remove the files just type in “yes” and hit enter. You should now be left with just the processed alignment file, the index file and the mapping stats. Well done!

# [Task 10]()
