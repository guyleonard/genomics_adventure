# Task 5 - Annotation of *de novo* Assembled Contigs
We will now annotate the contigs using BLAST and Pfam as with the unmapped contigs.

As before, we’ll 'call' open reading frames within the *de novo* assembly. We will use codon table 11 which defines the bacterial [codon usage table](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)🔍 and state that the sequences we are dealing with are not circular. We will also restrict the ORFs to just those sequences longer than 300 nucleotides (i.e. 100 amino acids). We will store the results in file contigs.orf.fasta.

```bash
getorf -table 11 -circular N -minsize 300 -sequence contigs.fasta -outseq contigs.orf.fasta
```

As with the unmapped reads we will search the open reading frames against the Pfam HMM database of protein families. Later on we will be able to use these results to identify Pfam domains which are unique to a particular strain.

This will take around 5 hours on an average laptop (or the workshop VMs) so it is recommended that you don’t run this command right now. This is mainly for future reference.

⚠️ do not run this ⚠️
```bash
pfam_scan.pl -fasta contigs.orf.fasta \
-dir  ~/workshop_meterials/genomics_adventure/db/pfam/ -outfile contigs.orf.pfam \
-cpu 2 -as
```

Precomputed results are available here:

> ~/workshop_materials/genomics_adventure/denovo_assembly/mapping_to_assembly/contigs.orf.pfam

# [Chapter 5](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_5/task_1.md) awaits..!
