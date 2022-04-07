# Task 9 - Mapping Statistics
Finally, we can generate some summary statistics! :boom:
```bash
samtools flagstat ecoli_mapped_namesort_fixmate_sort_markdup.bam > mappingstats.txt
```
This should only take a few seconds. Once complete you can view the 'mappingstats.txt' file using your favourite text editor or command (e.g. nano, vi, more or cat).

[IMAGE]

On the first line we can see that we have 5,645,403 reads in total, and none of which failed QC. Next, 97.9% of the reads mapped to the reference genome and 97.56% mapped with the expected distance between them, as paired reads. And only 63 reads could not have their other read-pair mapped.

0 reads have mapped to a different chromosome than their pair (0 have a mapping quality > 5 – this is a Phred scaled quality score much as we say in the FASTQ files). If there were any such reads they would likely be due to repetitive sequences (e.g phage insertion sites) or an insertion of plasmid or phage DNA into the main chromosome.

## Cleaning Up
It is generally good practice to keep your intermediary files whilst you are continuing your analysis, but in this case we really don't need some of the extra files we have created during this process. Indeed, we are only really interested in the final file 'ecoli_mapped_namesort_fixmate_sort_markdup.bam'. Therefore, you can remove the others:
```bash
rm ecoli_mapped.bam
rm ecoli_mapped.sam
rm ecoli_mapped_namesort.bam
rm ecoli_mapped_namesort_fixmate.bam
rm ecoli_mapped_namesort_fixmate_sort.bam
rm ecoli_mapped_sorted.bam
```

In case you get asked if you are sure to remove the files, just say 'yes'. You should now be left with just the processed alignment file, the index file and the mapping stats. Well done! :1st_place_medal:

# Go to [Task 10](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_2/task_10.md)
