# Task 14 - Automated Analyses
Viewing alignments is useful when convincing yourself or others that a particular mutation is real rather than an artefact and for getting a feel for short read sequencing datasets. However, if we want to quickly and easily find variants we need to be able to generate lists of variants, in which gene they occur (if any) and what effect they have. We also need to know which (if any) genes are missing (i.e. have zero coverage).

## Automated Variant Calling
To call (predict) variants we can use a number of packages (e.g. VarScan, GATK). However here, we will show you how to use the 'bcftools' package. First we need to generate a 'pileup' file which contains only the locations with the variants.

## Identify SNPs and Indels using Automated Variant Callers
Make sure you are in the directory.
> ~/workshop_materials/genomics_adventure/sequencing_data/ecoli/mapping_to_reference

Then type the following
```bash
bcftools mpileup
```

You should see a screen similar to the following:

[IMAGE]

If you are run 'bcftools' on large numbers of datasets with limited coverage where recombination is a factor, you can obtain increased sensitivity by passing all the BAM files to the variant caller simultaneously (hence the multiple BAM file options in bcftools).

Now, lets type the following:
```bash
bcftools mpileup -O v -P Illumina --threads 4 \
-f ~/workshop_materials/genomics_adventure/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
ecoli_mapped_namesort_fixmate_sort_markdup.bam > var.raw.vcf
```

This process may take 20 minutes or so, and will generate a '[VCF](https://en.wikipedia.org/wiki/Variant_Call_Format)' file containing the raw unfiltered variant calls for each position in the genome. Note that we are asking 'bcftools mpileup' to generate an uncompressed VCF output with the '-O v' option. '-P' tells 'bcftools' that it is dealing with Illumina data, so that it can apply to the correct model to help account for mis-calls or indels. 

This output by itself is not super useful on its own, as it contains information on each position in the genome. So let’s use 'bcftools' again to call what it thinks are the variant sites.

```bash
bcftools call -c -v --ploidy 1 -O v -o var.called.vcf var.raw.vcf
```

Note that we are asking 'bcftools' to call assuming a ploidy of 1, and to output only the variant sites in the 'VCF' format. Using the tool 'gre'p we can count how many sites were identified as being variant sites (i.e. sites with a potential mutation). We can ask grep not to count lines beginning with a comment (#) also.

```bash
grep -v -c  "^#" var.called.vcf
```

You should find around 174 sites. Don't worry if your number isn't exactly the same. Now we should filter this further, and with a tool made specifically to work with vcf files, in order to ensure we only retain regions where we have >90% allele frequency - we can do this with a tool called 'vcftools'.

```bash
vcftools --minDP 10 --min-alleles 2 --max-alleles 2 \
--non-ref-af 0.9 --vcf var.called.vcf --recode --recode-INFO-all \
--out var.called.filt 
```

You are safe to ignore any warnings. This command creates a file called 'var.called.filt.recode.vcf'. Once complete, you can view the file using the 'more' command (or your favourite editor). You should see something similar to below, (lines beginning with a '#' are just comment lines explaining the output):

[IMAGE]

You can see the chromosome, position, reference and alternate allele as well as a quality score for the SNP. This is a VCF file (Variant Call File),  a standard developed for the 1000 genomes project. The full specification is given [here](http://samtools.github.io/hts-specs/VCFv4.2.pdf) :PDF: :mag:.

The lines starting with 'DP' and 'INDEL' contain various details concerning the variants. For haploid organisms, most of these details are not necessary. This forms our definitive list of variants for this sample.

You can load the VCF file to IGV:

[IMAGE]

## Compare the Variants Found using this Method to Those You Found in the Manual Section
Can you see any variants which may have been missed? Often variants within a few bp of indels are filtered out as they could be spurious SNPs thrown up by a poor alignment. This is especially the case if you use non-gapped aligners such as Bowtie.

# [Task 15]()
