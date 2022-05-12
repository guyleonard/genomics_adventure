# task 1
cd ~/workhsop_materials/genomics_adventure/sequencing_data/ecoli
ls -lath
mkdir unmapped_assembly && cd unmapped_assembly
samtools view ../sequencing_data/ecoli/mapping_to_reference/ecoli_mapped_namesort_fixmate_sort_markdup.bam | head -n 5
samtools view -b -f 12 ../sequencing_data/ecoli/mapping_to_reference/ecoli_mapped_namesort_fixmate_sort_markdup.bam -o unmapped.bam
samtools view unmapped.bam | head -n 5
bedtools bamtofastq -i unmapped.bam -fq unmapped_r1.fastq -fq2 unmapped_r2.fastq

# task 2
grep -c "^@SRR" unmapped_r1.fastq unmapped_r2.fastq
tail -n 4 unmapped_r1.fastq unmapped_r2.fastq

# task 3
spades.py --careful -o spades_assembly -1 unmapped_r1.fastq -2 unmapped_r2.fastq

# task 4
quast.py --output-dir quast contigs.fasta
cat quast/report.txt

# task 5
blastn -db nt \
-query contigs.fasta \
-evalue 1e-06 \
-num_threads 4 \
-num_alignments 10 \
-num_descriptions 10 \
-out contigs.fasta.blastn
nano ../../precomputed/unmapped_assembly/spades_assembly/contigs.fasta.blastn
blastn -db /databases/ncbi/nt/2020-10-15/nt \
-query contigs.fasta \
-evalue 1e-06 \
-num_threads 4 \
-num_alignments 10 \
-num_descriptions 10 \
-outfmt '6 std stitle' \
-out contigs.fasta.blastn6

# task 6
getorf -table 11 -circular N -minsize 300 -sequence contigs.fasta -outseq contigs.orf.fasta
blastp -db nr \
-query contigs.orf.fasta \
-evalue 1e-06 \
-num_threads 4 \
-num_alignments 10 \
-num_descriptions 10 \
-outfmt '6 std stitle' \
-out contigs.orf.blastp6
blastn -subject ../../reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
-query contigs.fasta \
-outfmt 6

# task 7
hmmpress ../../db/pfam/Pfam-A.hmm
pfam_scan.pl -fasta contigs.orf.fasta \
-dir ~/workshop_materials/genomics_tutorial/db/pfam/ \
-outfile contigs.orf.pfam \
-cpu 4 \
-as


