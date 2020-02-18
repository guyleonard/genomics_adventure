# Task 15 Locating Genes that are Missing Compared to the Reference
We can use a command from the [BEDTools](http://bedtools.readthedocs.org/en/latest/) package  to identify annotated genes which are not covered by reads across their full length.

For example, lets try the following:
```bash
bedtools coverage \
-a ~/workshop_materials/genomics_adventure/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.gff \
-b ecoli_mapped_namesort_fixmate_sort_markdup.bam > gene_coverage.txt
```
This should only take a minute or so. The output contains one row per annotated gene - the 13th column contains the proportion of the gene that is covered by reads from our sequencing. 1.00 means the gene is 100% covered and 0.00 means no coverage.

