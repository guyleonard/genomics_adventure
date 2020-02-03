# Task 2 - Evaluating the Quality of Illumina Data Continues...
Do the same for read 2 as we have for read 1. Open fastqc and analyse the read 2 file. Look at the
various plots and metrics which are generated. How similar are they?

Note that the number of reads reported in both files is identical. Overall, both read 1 and read 2 can be regarded as 'good' data-sets.

## Quality control – filtering of Illumina data
In the next set of tasks we will be filtering the data to ensure that any low quality reads are removed, and that any
sequences containing adaptor sequences are either trimmed or removed altogether. To do this we will
use the 'Trime Galore!' program which is a program that integrates both 'fastqc' and 'cuadapt'. This package is remarkably fast and ensures that after filtering both read 1 and read 2 files are in the correct order, it is also nice to view your trimmed/cleaned sequencing data in fastqc.

Note: Typically when submitting raw Illumina data to NCBI or EBI you would submit unfiltered data, so
don't delete your original fastq files!

## Contaminant Checking
A number of tools are available which also enable to you to quickly search through your reads and assign them to particular taxa or taxonomic group. These can serve as a quick check to make sure your samples or libraries are not contaminated with DNA from other sources. If you are performing a de-novo assembly, for example, and have DNA sequences present from multiple
organisms, you will risk poor results and chimeric contigs.

Some ‘contaminants’ may turn out to be inevitable by-products of sampling and DNA extraction, and this is often the case with algae, and/or other symbionts but some groups have made amazing discoveries such as the discovery of a third symbiont (which turned out to be a yeast) in lichen, [here](http://science.sciencemag.org/content/353/6298/488.full).

Some tools you can use to check the taxonomic classification of reads include:
 * [Kraken](https://ccb.jhu.edu/software/kraken2/)
 * [Centrifuge](https://ccb.jhu.edu/software/centrifuge/)
 * [Blobology](https://blobtoolkit.genomehubs.org/)
 * Blast (in conjunction with sub-sampling your reads) and Krona to plot results
 * and many more!

## [Task 3]()
