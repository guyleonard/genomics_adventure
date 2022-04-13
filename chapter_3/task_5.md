# Task 5 - Analysing the *de novo* Assembled Reads
Now that we have assembled the reads and have a feel for how much (or in this case, how little) data we have, we can set about analysing it. By analysing, we mean identifying which genes are present, which organism they are from and whether they form part of the main chromosome or are an independent unit (e.g. plasmid).

We are going to take a 3-prong approach:

1. Firstly we will search the nucleotide sequences of the contigs against the NCBI non-redundant database. This will enable us to identify the species to which a given contig matches best (or most closely).
2. Secondly will call open reading frames within the contigs and search those against the Swissprot database of manually curated (i.e. high quality) annotated protein sequences.
3. And finally, we will search the open reading frames against the Pfam database of protein families [Pfam](​http://pfam.xfam.org).

Why not just search the NCBI blast database? Well, remember nearly all of our biological knowledge is based on homology – if two proteins are similar they probably share an evolutionary history and may thus share functional characteristics. Metrics to define whether two sequences are homologous are notoriously difficult to define accurately. If two sequences share 90% sequence identity over their length, you can be pretty sure they are homologous. If they share 2% they probably aren't. But what if they share 30%? This is the notorious twilight zone of 20-30% sequence identity where it is very difficult to judge whether two proteins are homologous based on sequence alone.

To help overcome this searching more subtle signatures may help – this is where Pfam comes in. Pfam is a database which contains protein families identified by particular signatures or patterns in their protein sequence. These signatures are modeled by Hidden Markov Models ([HMM](https://en.wikipedia.org/wiki/Hidden_Markov_model)s :mag:) and used to search query sequences. These can provide a high level annotation where BLAST might otherwise fail. It also has the advantage of being much faster than BLAST.

## Search Contigs against NCBI non-redundant Database
To do this we will use the command line version of the Basic Local Alignment Search Tool ([BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)), which finds regions of similarity between biological sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance.

We will use the database commongly denoted as 'nt' - the nucleotide database is a collection of sequences from several sources, including GenBank, RefSeq, TPA and PDB including genome, gene and transcript sequences. 

The following command executes a nucleotide BLAST search (using 'blastn' for nucleotides) of the sequences in the 'contigs.fasta' file against the 'nt' database. :warning::no_entry: PLEASE DO NOT RUN THIS! :no_entry_sign::warning:
```
blastn -db nt -query contigs.fasta -evalue 1e-06 -num_threads 4 -num_alignments 10 -num_descriptions 10 -out contigs.fasta.blastn
```

The reason I stress to not run this is that blast can be somewhat slow, and would require too much downtime in our adventure. Nevertheless the command is there so you can understand what you would do to run one of your own. We will provide a pre-computed output below. But first let's explore the command a little bit.

* '-db'​ - the prepared blast database to search
* '-​evalue' - apply an e-value (expectation value) cutoff (read more [here](http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html​):mag:), here a cutoff of 1e-06 to limit ourselves to statistically significant hits (i.e. in this case 1 in 1 million likelihood of a hit to a database of this size by a sequence of this length).
* '-num_alignments' and '-num_descriptions' - options tell blastn to only display the top 10 results for each hit
* '-num_threads'​ - use 4 CPU cores
* '-out' - the file in which to place the output

There is lots of information on running blast from the command line in NCBI's handbook, [here](http://www.ncbi.nlm.nih.gov/books/NBK1763/).

Take a look at the results of the blast in your favourite text editor, e.g.
```
nano contigs.fasta.blastn
```


