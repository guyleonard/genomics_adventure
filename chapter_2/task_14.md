# Task 14 - Automated Analyses
Viewing alignments is useful when convincing yourself, and your peers, that a particular mutation is real rather than an artefact, and for getting a feel for short read sequencing datasets. However, if we want to quickly and easily find variants we need to be able to generate lists of variants, and in which gene they occur (if any) and what effect that they might have. We also need to know which (if any) genes are missing (i.e. have zero coverage).

## Automated Variant Calling
To 'call' (predict) variants we can use a number of packages (e.g. VarScan, GATK, Picard Tools). However here, we will show you how to use the 'bcftools' package. First we need to generate a 'pileup' file which contains only the locations with the variants.

## Identify SNPs and Indels using Automated Variant Callers
Make sure you are in the directory.
> ~/workshop_materials/genomics_adventure/sequencing_data/ecoli/mapping_to_reference

Then type the following
```bash
bcftools mpileup
```

You should see a screen similar to the following:

![bcf tools](https://github.com/guyleonard/genomics_adventure/blob/f6986943cc8b37de06003550777a771068e7fbee/chapter_2/images/chapter_2_task_14_image_1.png)

If you run 'bcftools' on large numbers of datasets with limited coverage where recombination is a factor, you can obtain increased sensitivity by passing all the BAM files to the variant caller simultaneously (hence the multiple BAM file options in bcftools).

Now, lets type the following:
```bash
bcftools mpileup -O v -P Illumina --threads 4 \
-f ~/workshop_materials/genomics_adventure/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
ecoli_mapped_namesort_fixmate_sort_markdup.bam > var.raw.vcf
```

This process may take 20 minutes or so, and will generate a '[VCF](https://en.wikipedia.org/wiki/Variant_Call_Format)' :mag: (read up a little and read on whilst you wait) file containing the raw unfiltered variant calls for each position in the genome. Note that we are asking 'bcftools mpileup' to generate an uncompressed VCF output with the '-O v' option. The '-P Illumina' tells 'bcftools' that it is dealing with Illumina data, so that it can apply to the correct model to help account for mis-calls and/or indels. 

This output by itself is not super useful on its own, as it contains information on each position in the genome. So let’s use 'bcftools' again to 'call' what it thinks are the variant sites.

```bash
bcftools call -c -v --ploidy 1 -O v -o var.called.vcf var.raw.vcf
```

Here we have asked 'bcftools' to 'call' variants assuming a ploidy of 1, and to output only the variant sites in the 'VCF' format. Now, using the tool 'grep' we can count how many sites were identified as being variant sites (i.e. sites with a potential mutation). We can ask grep not to count lines beginning with a comment (#) also. Please ask a TA (or read grep's manual if you want to find out what the command is doing).

```bash
grep -v -c  "^#" var.called.vcf
```

You should find around 174 sites. Don't worry if your number isn't exactly the same for this adventure (different versions of tools may contain bugs or updates that can change output between them, this is why it is important to document your work including version numbers of the programs you have used). Now we can filter this further, and with a tool made specifically to work with VCF files, in order to ensure we only retain regions where we have >90% allele frequency - we can do this with a tool called 'vcftools'. (are these tools pronounced the same in Spanish? :stuck_out_tongue_winking_eye: answers on a postcard :love_letter:)

```bash
vcftools --minDP 10 --min-alleles 2 --max-alleles 2 \
--non-ref-af 0.9 --vcf var.called.vcf --recode --recode-INFO-all \
--out var.called.filt 
```

You are safe to ignore any warnings that you see. This command creates a file called 'var.called.filt.recode.vcf'. Once complete, you can try viewing the file using the 'more' command (or your favourite text editor). You should see something similar to below, (lines beginning with a '#' are just comment lines explaining the output):

![vcf output](https://github.com/guyleonard/genomics_adventure/blob/f6986943cc8b37de06003550777a771068e7fbee/chapter_2/images/chapter_2_task_14_image_2.png)

You can see the chromosome, position, reference and alternate allele as well as a quality score for the SNP. This is a VCF file (Variant Call File), a standard developed for the [1000 Genomes Project](https://en.wikipedia.org/wiki/1000_Genomes_Project). The full specification is given [here](http://samtools.github.io/hts-specs/VCFv4.2.pdf) :PDF: :mag:, but you don't need to fully understand or read the whole document - just use it for reference.

The lines starting with 'DP' and 'INDEL' contain various details concerning the variants. For haploid organisms, most of these details are not necessary. This forms our definitive list of variants for this sample.

You can now try loading the VCF file into IGV:

![igv with vcf](https://github.com/guyleonard/genomics_adventure/blob/f6986943cc8b37de06003550777a771068e7fbee/chapter_2/images/chapter_2_task_14_image_3.png)

## Compare the Variants Found using this Method to Those You Found in the Manual Section
Can you see any variants which may have been missed? Often variants within a few bp of indels are filtered out as they could be spurious SNPs thrown up by a poor alignment. This is especially the case if you use non-gapped aligners such as Bowtie.

# Go to [Task 15](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_2/task_15.md)
