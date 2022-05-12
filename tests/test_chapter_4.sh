# task 1
cd ~/workshop_materials/genomics_adventure
mkdir denovo_assembly && cd denovo_assembly
ln -s ../../precomputed/denovo_assembly/assembly .
spades.py -o assembly -1 ../sequencing_data/ecoli/read_1_val_1.fq.gz -2 ../sequencing_data/ecoli/read_2_val_2.fq.gz

# task 2
cd assembly
quast.py --output-dir quast contigs.fasta

# task 3
mkdir mapping_to_assembly
cd mapping_to_assembly
ln -s ../assembly/contigs.fasta .
bwa index contig.fasta
bwa mem -t 2 contigs.fasta \
../../sequencing_data/ecoli/read_1_val_1.fq.gz \
../../sequencing_data/ecoli/read_2_val_2.fq.gz \
> contigs_mapped.sam
samtools view -bS contigs_mapped.sam > contigs_mapped.bam
samtools sort -o contigs_mapped_sorted.bam contigs_mapped.bam
samtools index contigs_mapped_sorted.bam
bash bwa index contigs.fasta && \ bwa mem -t 2 contigs.fasta \ ../../sequencing_data/ecoli/read_1_val_1.fq.gz \ ../../sequencing_data/ecoli/read_2_val_2.fq.gz \ | samtools sort -O bam -o contigs_mapped_sorted.bam && \ bwa index contigs_mapped_sorted.bam
samtools flagstat contigs_mapped_sorted.bam
qualimap bamqc -outdir bamqc -bam contigs_mapped_sorted.bam
blastn -subject contigs.fasta \
-query ../../unmapped_assembly/spades_assembly/contigs.fasta \
-outfmt 6 -out check_plasmid.blastn

# task 4

