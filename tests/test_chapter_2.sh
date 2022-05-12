#!/bin/bash

# task 1
cd ~/workshop_materials/genomics_adventure/sequencing_data/ecoli
ls -lath
ln -s SRR857279_1.fastq.gz read_1.fastq.gz
ln -s SRR857279_2.fastq.gz read_2.fastq.gz
zcat read_1.fastq.gz | tail | grep @SRR
zcat read_2.fastq.gz | tail | grep @SRR
zcat read_1.fastq.gz | grep @SRR | wc -l
zcat read_2.fastq.gz | grep @SRR | wc -l

fastqc -o fastqc read_1.fastq.gz read_2.fastq.gz

# task 2
trim_galore --help
trim_galore --illumina --paired --phred33 read_1.fastq.gz read_2.fastq.gz
zcat read_1_val_1.fq.gz | wc -l
zcat read_2_val_2.fq.gz | wc -l

# task 3
seqtk
seqtk sample
seqtk sample read_1_val_1.fq.gz 0.1 > read_1_val_1_subsample_one.fq
seqtk sample read_1_val_1.fq.gz 0.1 > read_1_val_1_subsample_two.fq
head read_1_val_1_subsample*
seqtk sample -s 1234 read_1_val_1.fq.gz 0.1 > read_1_val_1_subsample_three.fq
seqtk sample -s 5678 read_1_val_1.fq.gz 0.1 > read_1_val_1_subsample_four.fq
head read_1_val_1_subsample_three.fq
head read_1_val_1_subsample_four.fq
rm *.fq
seqtk sample -s 628 read_1_val_1.fq.gz 0.5 > read_1_val_1_subsampled.fq
seqtk sample -s 628 read_2_val_2.fq.gz 0.5 > read_2_val_2_subsampled.fq

# task 4
# task 5
cd ~/workshop_materials/genomics_adventure/reference_sequences/ecoli
ls -lath
bwa index GCF_000005845.2_ASM584v2_genomic.fna

# task 6
cd ~/workshop_materials/genomics_adventure/sequencing_data/ecoli
mkdir mapping_to_reference
cd mapping_to_reference
bwa mem -t 4 \
~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_1_val_1.fq.gz \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_2_val_2.fq.gz \
-o ecoli_mapped.sam

# task 7
samtools view
samtools view -b -S \
-T ~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
ecoli_mapped.sam -o ecoli_mapped.bam
samtools sort ecoli_mapped.bam -o ecoli_mapped_sorted.bam

bwa mem -t 4 \
~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_1_val_1.fq.gz \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_2_val_2.fq.gz \
-o ecoli_mapped.sam
samtools view -b -S \
-T ~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
ecoli_mapped.sam -o ecoli_mapped.bam
samtools sort ecoli_mapped.bam -o ecoli_mapped_sorted.bam

bwa mem -t 4 \
~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_1_val_1.fq.gz \
~/workshop_materials/genomics_tutorial/sequencing_data/ecoli/read_2_val_2.fq.gz \
| samtools sort -O bam -o ecoli_mapped_sorted_onecommand.bam

# task 8
samtools sort -n -o ecoli_mapped_namesort.bam ecoli_mapped.bam
samtools fixmate -m ecoli_mapped_namesort.bam  ecoli_mapped_namesort_fixmate.bam
samtools sort -o ecoli_mapped_namesort_fixmate_sort.bam ecoli_mapped_namesort_fixmate.bam
samtools markdup -r ecoli_mapped_namesort_fixmate_sort.bam  ecoli_mapped_namesort_fixmate_sort_markdup.bam 
samtools index ecoli_mapped_namesort_fixmate_sort_markdup.bam

# task 9
samtools flagstat ecoli_mapped_namesort_fixmate_sort_markdup.bam > mappingstats.txt
rm ecoli_mapped.bam
rm ecoli_mapped.sam
rm ecoli_mapped_namesort.bam
rm ecoli_mapped_namesort_fixmate.bam
rm ecoli_mapped_namesort_fixmate_sort.bam
rm ecoli_mapped_sorted.bam

# task 10
qualimap bamqc
qualimap bamqc -outdir bamqc \
-bam ecoli_mapped_namesort_fixmate_sort_markdup.bam \
-gff ~/workshop_materials/genomics_tutorial/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.gff

# task 14
bcftools mpileup
bcftools mpileup -O v -P Illumina --threads 4 \
-f ~/workshop_materials/genomics_adventure/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
ecoli_mapped_namesort_fixmate_sort_markdup.bam > var.raw.vcf
bcftools call -c -v --ploidy 1 -O v -o var.called.vcf var.raw.vcf
grep -v -c  "^#" var.called.vcf
vcftools --minDP 10 --min-alleles 2 --max-alleles 2 \
--non-ref-af 0.9 --vcf var.called.vcf --recode --recode-INFO-all \
--out var.called.filt

# task 15
bedtools coverage \
-a ~/workshop_materials/genomics_adventure/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.gff \
-b ecoli_mapped_namesort_fixmate_sort_markdup.bam > gene_coverage.txt

samtools view -bs 42.5 \
ecoli_mapped_namesort_fixmate_sort_markdup.bam \
> ecoli_mapped_namesort_fixmate_sort_markdup_subsampled.bam
bedtools coverage \
-a ~/workshop_materials/genomics_adventure/reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.gff \
-b ecoli_mapped_namesort_fixmate_sort_markdup_subsampled.bam > gene_coverage.txt
sort -t $'\t' -g -k 13 gene_coverage.txt | less -S
sort -t $'\t' -g -k 13 gene_coverage.txt | cut -f1-8,10-13 | less -S


