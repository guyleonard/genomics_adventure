# Task 6 - Obtain Open Reading Frames
The first task is to call open reading frames within the contigs. These are designated by canonical start and stop codons and are usually identified by searching for regions free of stop codons. We will use the EMBOSS package program 'getorf' to call these.

We will use codon table 11 which defines the bacterial codon usage [table](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) :mag: and state that the sequences we are dealing with are not circular (they are nowhere near long enough!). We will also restrict the ORFs to just those sequences longer than 300 nucleotides (i.e. 100 amino acids). We will store the results in file contigs.orf.fa.

PS - You can go back to typing the commands now! :smile:
```bash
getorf -table 11 -circular N -minsize 300 -sequence contigs.fasta -outseq contigs.orf.fasta
```

If we look at the output file we can see that it is a FASTA formatted file containing the name of the contig on which the ORF occurs, followed by an underscore and a number (e.g. _1) to indicate the number of the ORF on that contig. The numbers in square brackets indicate the start and end position of the ORF on the contig (i.e. in nucleotide space).

Taking a look through the output, you can see that the first ORF occurs on NODE_1 and is between position 345 and 644 and the third ORF occurs between positions 983 and 1900. Scrolling down through the file, we can see that the 81st ORF can be found between base positions of 46505 and 45828 but on the reverse strand (labelled (REVERSE SENSE)).

Also note that many ORFs do not start with a Methionine. This is because by default the getorf program calls ORFs between stop codons rather than start and stop codons. Primarily this is to avoid spurious ORFs due to Met residues within a protein sequence and to ensure untranslated regions are captured.

## Search Open Reading Frames against NCBI non-redundant Database
The first thing we can do with these open reading frames is to search them against the NCBI non-redundant database of protein sequences to see what they may match.

Here we will step through the process to perform a BLAST search using the non-redundant (nr) database, using the 'blastp' program and store the results in contigs.orf.blastp. We'll apply an [e-value](http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html) :mag: cutoff of 1e-06 to limit ourselves to statistically significant hits (i.e. in this case 1 in 1 million likelihood of a hit to a database of this size by a sequence of this length). The '-num_alignments' and '-num_descriptions' options tell blastp to only display the top 10 results for each hit, the "-num_threads" tells blastp to use 4 CPU cores.

:warning::no_entry: PLEASE DO NOT RUN THIS! :no_entry_sign::warning: It would take about 18+ hours to run on 4 CPUs. We did it for you though! :smile:
```bash
blastp -db nr \
-query contigs.orf.fasta \
-evalue 1e-06 \
-num_threads 4 \
-num_alignments 10 \
-num_descriptions 10 \
-outfmt '6 std stitle' \
-out contigs.orf.blastp6
```

Open the results file with your favourite text editor (e.g. nano) and search for plasmid in the text. You should find a number of hits to plasmid related proteins, can you find any others? The evidence is not conclusive, but combined with the coverage, it is starting to look like this contig is most likely a plasmid.

## Additional Checks
Check that the contigs do not appear in the reference sequence. In theory, the unmapped reads used to generate the contigs should not assemble into something which will map against the genome. However, it is always possible (especially with more complex genomes), that this might happen.

To do this, we can use the BLAST+ package in a different mode to compare two sequences against each other. Unlike the previous examples where we have searched against a database of sequences, here we are doing a simple search of the contigs against the reference genome. You can run this one, it will be very fast. We won't save the results to a file, instead they will show on your terminal.
```
blastn -subject ../../reference_sequences/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
-query contigs.fasta \
-outfmt 6
```

You should see something similar to below, it shows that some of the contigs have hit the reference sequence.
```bash
NODE_1_length_46850_cov_69.441152	NC_000913.3	100.000	28	0	0	1	28	392940	392967	2.72e-05	52.8
NODE_1_length_46850_cov_69.441152	NC_000913.3	96.774	31	1	0	1	31	566803	566773	2.72e-05	52.8
NODE_1_length_46850_cov_69.441152	NC_000913.3	100.000	28	0	0	1	28	1095476	1095503	2.72e-05	52.8
NODE_1_length_46850_cov_69.441152	NC_000913.3	100.000	28	0	0	1	28	2170197	2170170	2.72e-05	52.8
NODE_5_length_1363_cov_71.788491	NC_000913.3	100.000	28	0	0	1	28	4499486	4499513	7.79e-07	52.8
```

In the records above the e-values are very significant; however, since our subject sequence is now very small and we know that it is related to our strain, it is not surprising that there are some hits. We should only be concerned about whole contigs that map closely to the reference genome.

# Go to [Task 7](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_3/task_7.md)

