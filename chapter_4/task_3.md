# Map Reads Back to Assembly
Here we will use BWA again to index the contigs.fasta file and remap the reads. This is almost identical to the procedure we followed during the alignment section, the only difference is that instead of aligning to the reference genome, we are aligning to our newly created reference.

Make sure that you are in the directory "~/workshop_materials/genomics_adventure/denovo_assembly/". Let's create a new directory to keep our work seprate and organised, we will also create a link to our contig data. You may prefer to copy your data, but this way we can see where the contigs.fasta file has come from when we come back to our analyses at a later point.

```bash
mkdir mapping_to_assembly

cd mapping_to_assembly

ln -s ../assembly/contigs.fasta .
```

We now need to index and map our files, this is similar to the steps from Chapter 2... Indeed we will use the QC reads we created in Task 2! We will go through the steps one-by-one, although don't forget you can combine them all in to one step if you are feeling brave!
```bash
# Index the contigs
bwa index contig.fasta

# align QC reads to contigs and output SAM file
bwa mem -t 2 contigs.fasta \
../../sequencing_data/ecoli/read_1_val_1.fq.gz \
../../sequencing_data/ecoli/read_2_val_2.fq.gz \
> contigs_mapped.sam
```

Once complete we can convert the SAM file to a BAM file:
```bash
samtools view -bS contigs_mapped.sam > contigs_mapped.bam
```

And then we can sort the BAM file:
```bash
samtools sort -o contigs_mapped_sorted.bam contigs_mapped.bam
```

Once completed, we can index the BAM file:
```
samtools index contigs_mapped_sorted.bam
```

<details>
  <summary>Advanced: All in one command</summary>
  ```bash
  bwa index contigs.fasta && \
  bwa mem -t 2 contigs.fasta \
  ../../sequencing_data/ecoli/read_1_val_1.fq.gz \
  ../../sequencing_data/ecoli/read_2_val_2.fq.gz \
  | samtools sort -O bam -o contigs_mapped_sorted.bam && \
  bwa index contigs_mapped_sorted.bam
  ```
</details>

We can then (at last!) obtain some basic summary statistics using the samtools flagstat command:
```
samtools flagstat contigs_mapped_sorted.bam

5904122 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
1536 + 0 supplementary
0 + 0 duplicates
5903994 + 0 mapped (100.00% : N/A)
5902586 + 0 paired in sequencing
2951293 + 0 read1
2951293 + 0 read2
5846610 + 0 properly paired (99.05% : N/A)
5902444 + 0 with itself and mate mapped
14 + 0 singletons (0.00% : N/A)
36952 + 0 with mate mapped to a different chr
36263 + 0 with mate mapped to a different chr (mapQ>=5)
```

We can see that very few of the reads do not map back to the conigs. Importantly 99% of the reads are properly paired, which gives us some indication that there are not too many mis-assemblies.

We can run 'qualimap' to get some more detailed information (and some images too), it'll take a couple of minutes:
```bash
qualimap bamqc -outdir bamqc -bam contigs_mapped_sorted.bam

firefox bamqc/qualimapReport.html
# or
nano bamqc/genome_results.txt
```

Go to either the "Chromosome stats" section if you opened the '.html' or the "Coverage per contig" section from the text file. We can see that the larger of our contigs have a mean coverage of around 160 - which is what we would expect from our original alignment.

If you notice very carefully :wink:, there is one contig which has a size of 46899 - this is very very close to the size (46850) of the main contig we found in the unmapped reads assembly - another good indication that it is a separate sequence (remember we suspected it was a plasmid) and not integrated into a chromosome. We can double check this with a quick blast search...

```bash
blastn -subject contigs.fasta \
-query ../../unmapped_assembly/spades_assembly/contigs.fasta \
-outfmt 6 -out check_plasmid.blastn
```

Opening the 'check_plasmid.blastn' we can see the top hit as:
```
NODE_1_length_46850_cov_69.441152       NODE_33_length_46899_cov_70.549934      100.000 46850   0       0       1       46850   50      46899   0.0     86516
```

This shows us that this contig exactly almost matches that in the unmapped assembly, strongly supporting that this is a plasmid sequence and not integrated into the chromosomes.

# Go to [Task 4](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_4/task_3.md)
