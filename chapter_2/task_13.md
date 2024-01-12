## Task 13 - Identify SNPs and Indels Manually

Now, can you find any SNPs? Which genes (if any) are they in? How reliable do they look? (Hint – look at the number of reads mapping, their orientation - i.e. which strand they are on and how 'bright' the base-calls are).

Zoom into the maximum magnification at the site of the SNP. Can you determine whether a SNP results in a synonymous (i.e. silent) or non-synonymous change in the amino acid? Can you use the [Protein Data Bank](https://www.rcsb.org/search/)'s search, or other resources to determine whether or not this occurs in a catalytic site, or another functionally crucial region? (NB - this may not always be possible depending on how well the gene you select has been studied).

What effect do you think this would have on the cell?

### Example: Identifying Variants Manually

Here are some regions where there are differences in the organism sequenced and the reference: Can you interpret what has happened to the genome of our strain? Try to work out what is going on yourself before looking at the comments.

Paste each of the genomic locations in this box to IGV and click go:

```
NC_000913.3:2,108,392-2,133,153
NC_000913.3:3,662,049-3,663,291
NC_000913.3:4,296,249-4,296,510
NC_000913.3:565,965-566,489
```

#### Region NC_000913.3:2,108,392-2,133,153

This area corresponds to a drop in coverage identified by QualiMap. It looks like a fairly large region of about 17 kbases which was present in the reference and is missing from our sequenced genome. It looks like about 12 genes from the reference strain are not present in our strain.

Q: Is this real or an artefact?

A: We have coverage of about 60x either side of the deletion and nothing at all within. There are nice clean edges to the start and end of the deletion. We also have paired reads which span the deletion. This is exactly what you would expect if the two regions of coverage were actually joined together. Therefore, it is most likely a deletion event.

#### Region NC_000913.3:3,662,049-3,663,291

Zoom right in until you can see the reference sequence and protein sequence at the bottom of the display.

The first thing to note is that only discrepancies with respect to the reference are shown. If a read is entirely the same as its reference, it will appear compltely grey. Blue and red blocks indicate the presence of an 'abnormal' distance between paired-end reads. Note, unless this is consistent across most of the reads at a given position, it is not likely to be significant.

Here we can observe a C -> T SNP. This changes the codon from CAG->TAG (remember to check what strand the gene is on - this one is on the forward strand, if it was on the reverse strand you would have to take the reverse complement of the codon to interpret the amino acid it codes for) and results in a Gln -> Stop mutation in the final protein product which is very likely to change the effect of the protein product. 

Hover over the gene to get some more information from the annotation... Since it is a drug resistance protein it could be very significant!!


One additional check we can perform is whether the SNPs occur when reading the forward strand. We can check this by looking at the direction of the grey reads, or by hovering over the coverage graph. We can see that approximately half of the bases reporting the C -> T mutation occur in read 1 (forward arrow), and half in read 2 (reverse arrow). This adds confidence to the base-call as it reduces the likelihood of this SNP being the result of a PCR duplication error.

Note, sequencing errors in Illumina data can be quite common (We rely on depth of sequencing to average out these errors. This is why people often mention that a minimum median coverage of 20-30x across the genome is required for accurate SNP-calling with Illumina data). This is not always true for 'simple' organisms such as prokaryotes, but for diploid and polyploid organisms it becomes much more important because each position may have one, two or multiple alleles changed.

#### Regions U00096.3:4,296,332-4,296,428

Much of the same guidelines apply for indels as they do for SNPs. Here we have an insertion of two bases 'CG' in our sample compared to the reference. Again, we can see how much confidence we have that the insertion is real by checking that the indel is found on both read 1 and read 2 and on both strands.  

The insertion is signified by the presence of a purple bar. Hover your mouse over the reference sequence to get details of the gene. We can see that it occurs in a repeat region and is unlikely to have very significant effects.

One can research the effect that a SNP or Indel may have by finding the relevant gene at http://www.uniprot.org/ (or google 'mdtF uniprot' in the previous case).

### Conclusions

It should be clear from this quick exercise that trying to work out where SNPs and Indels are manually is a fairly tedious process especially where there are many mutations. As such, the next section will look at how to obtain a 'spreadsheet friendly' summary of these.

## Recap - Tips for SNP/Indel Identification
 * Only changes from the reference sequence are displayed in IGV
 * Genuine SNPs/Indels should be present on both read 1 and read 2
 * Genuine SNPs/Indels should be present on both strands
 * Genuine SNPs/Indels should be supported by a good (i.e. 20-30x) depth of coverage
 * Very important mutations (i.e. ones relied upon in a paper) should be confirmed via PCR/Sanger sequencing.

# Go to [Task 14](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_2/task_14.md)
