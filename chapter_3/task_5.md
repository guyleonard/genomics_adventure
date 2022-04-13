# Task 5 - Analysing the *de novo* Assembled Reads
Now that we have assembled the reads and have a feel for how much (or in this case, how little) data we have, we can set about analysing it. By analysing, we mean identifying which genes are present, which organism they are from and whether they form part of the main chromosome or are an independent unit (e.g. plasmid).

We are going to take a 3-prong approach:

1. Firstly, we will search the nucleotide sequences of the contigs against the NCBI non-redundant database. This will enable us to identify the species to which a given contig matches best (or most closely).
2. Secondly, we will call open reading frames within the contigs and search those against the Swissprot database of manually curated (i.e. high quality) annotated protein sequences.
3. And finally, we will search the open reading frames against the Pfam database of protein families [Pfam](​http://pfam.xfam.org).

Why not only search the NCBI blast database? Well, remember nearly all of our biological knowledge is based on homology – if two proteins are similar they probably share an evolutionary history and may thus share functional characteristics. Metrics to define whether two sequences are homologous are notoriously difficult to define accurately. If two sequences share 90% sequence identity over their length, you can be pretty sure they are homologous. If they share 2% they probably aren't. But what if they share 30%? This is the notorious twilight zone of 20-30% sequence identity where it is very difficult to judge whether two proteins are homologous based on sequence alone.

To help overcome this searching more subtle signatures may help – this is where Pfam comes in. Pfam is a database which contains protein families identified by particular signatures or patterns in their protein sequence. These signatures are modeled by Hidden Markov Models ([HMM](https://en.wikipedia.org/wiki/Hidden_Markov_model)s :mag:) and used to search query sequences. These can provide a high level annotation where BLAST might otherwise fail. It also has the advantage of being much faster than BLAST.

## Search Contigs against NCBI nucleotide Database
To do this we will use the command line version of the Basic Local Alignment Search Tool ([BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)), which finds regions of similarity between biological sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance.

We will use the database commongly denoted as 'nt' - the nucleotide database is a collection of sequences from several sources, including GenBank, RefSeq, TPA and PDB including genome, gene and transcript sequences. 

The following command executes a nucleotide BLAST search (using 'blastn' for nucleotides) of the sequences in the 'contigs.fasta' file against the 'nt' database. :warning::no_entry: PLEASE DO NOT RUN THIS! :no_entry_sign::warning:. Remember the notation with "\" is just to make it easier to read here, you don't have to type them. Not that we will this time anyway.
```
blastn -db nt \
-query contigs.fasta \
-evalue 1e-06 \
-num_threads 4 \
-num_alignments 10 \
-num_descriptions 10 \
-out contigs.fasta.blastn
```

The reason I stress to not run this is that blast can be somewhat slow, and would require too much downtime in our adventure. Nevertheless, the command is there so you can understand what you would do to run one of your own. We will provide a pre-computed output below. But first let's explore the command a little bit.

* '-db'​ - the prepared blast database to search
* '-​evalue' - apply an e-value (expectation value) cutoff (read more [here](http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html​):mag:), here a cutoff of 1e-06 to limit ourselves to statistically significant hits (i.e. in this case 1 in 1 million likelihood of a hit to a database of this size by a sequence of this length).
* '-num_alignments' and '-num_descriptions' - options tell blastn to only display the top 10 results for each hit
* '-num_threads'​ - use 4 CPU cores
* '-out' - the file in which to place the output

There is a lot of further information on running blast from the command line in NCBI's handbook, [here](http://www.ncbi.nlm.nih.gov/books/NBK1763/) :mag:.

Let's take a look at the results of the blast in your favourite text editor. As mentioned before we have pre-computed this data for you to save some time. It should be in the directory "~/workshop_materials/genomics_adventure/precomputed/unmapped_assembly/spades_assembly". We can access it like below. You may like to copy it to your local directory also.
```
nano ../../precomputed/unmapped_assembly/spades_assembly/contigs.fasta.blastn
```

You should see something similar to this: 
```
Query= NODE_1_length_46850_cov_69.441152

Length=46850
                                                                      Score     E
Sequences producing significant alignments:                          (Bits)  Value

CP048440.1 Escherichia coli strain NBRC 3301 plasmid putative_pEc...  86516   0.0
MF370216.1 Escherichia coli strain J53 plasmid pOX38-Gen, complet...  86516   0.0
CP053606.1 Escherichia coli strain NEB_Turbo plasmid F', complete...  86505   0.0
CP014273.1 Escherichia coli K-12 strain C3026 plasmid F128-(C3026...  86505   0.0
CP014271.1 Escherichia coli K-12 strain DHB4 plasmid F128-(DHB4),...  86505   0.0
CP008801.1 Escherichia coli KLY, complete genome                      86505   0.0
CP053608.1 Escherichia coli strain NEB5-alpha_F'Iq plasmid F'Iq, ...  86500   0.0
CP040643.1 Escherichia coli strain BE104 chromosome, complete genome  86500   0.0
CP018801.1 Escherichia coli strain tolC-, complete genome             86500   0.0
LR134083.1 Escherichia coli strain NCTC12655 genome assembly, chr...  86492   0.0


>CP048440.1 Escherichia coli strain NBRC 3301 plasmid putative_pEcol1, complete
sequence
Length=98786

 Score = 86516 bits (46850),  Expect = 0.0
 Identities = 46850/46850 (100%), Gaps = 0/46850 (0%)
 Strand=Plus/Plus

Query  1      CGTGTCCACTATTGCTGGGTAAGATCAGCTTTCTGAACGGGTGATGGCCTTCAGCGCCCT  60
              ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  37222  CGTGTCCACTATTGCTGGGTAAGATCAGCTTTCTGAACGGGTGATGGCCTTCAGCGCCCT  37281
```

This is the standard blast output format. It shows you the hits for each query sequence, and the alignments of the sequences. It is somewhat similar to the webserver version but in text format. The first query should be "Node_1". There are a number of good hits, do they give you a clue to what it might be?

The output is quite 'human-friendly' but it isn't the most useful for parsing with bioinformatics tools, or even opening in a spreadsheet. Luckily blast has some different output options that we can use, for example '6', which is a tabulated format. :warning::no_entry: PLEASE DO NOT RUN THIS! :no_entry_sign::warning:
```
blastn -db /databases/ncbi/nt/2020-10-15/nt \
-query contigs.fasta \
-evalue 1e-06 \
-num_threads 4 \
-num_alignments 10 \
-num_descriptions 10 \
-outfmt '6 std stitle' \
-out contigs.fasta.blastn6
```

Here you can see we have added the option "-outfmt '6 std stitle'" which changes the output to a table and includes the standard information of e-value and query hit along with the title of the blast hit.

Take a look at the new output, much nicer to read and parse!
```
nano ../../precomputed/unmapped_assembly/spades_assembly/contigs.fasta.blastn6
```

# Go to [Task 6](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_3/task_6.md)
