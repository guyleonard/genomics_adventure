# Task 7 Run Open Reading Frames Through pfam_scan
Pfam is a database of protein families. They are grouped together using a number of criteria based on their function. For more information you can read more [here](http://en.wikipedia.org/wiki/Pfam). Pfam is grouped into several databases depending on the level of curation. Pfam-A is high-quality manual curation and consists of around 20,000 families. Pfam-B is full of automated predictions which may be informative but should not be relied upon without additional evidence. Pfam will also search for signatures of active-sites if you specify the correct option.

Before we can use Pfam-A we will need to make sure the database is ready to go, if you are on one of the workshops then the database will be downloaded already, and if you are following at home you will have downloaded it during the setup. This step won't take too long, and it's a bit like when we index BAM files, except we use the program 'hmmpress'.

```bash
hmmpress ../../db/pfam/Pfam-A.hmm
```

Now we want to search the Pfam database of Hidden Markov Models to see which protein families are contained within this contig. You'll notice that this runs considerably faster than BLAST. We will search using
the contigs.orf.fasta file against the Pfam databases and output the results to contigs.orf.pfam. We'll use 4 CPU cores for the search and state that we want to search active site residues. This should only take a few minutes.

```bash
pfam_scan.pl -fasta contigs.orf.fasta \
-dir ~/workshop_materials/genomics_tutorial/db/pfam/ \
-outfile contigs.orf.pfam \
-cpu 4 \
-as
```

Now take a look at the output, it should look similar to the data shown below (note that there are about 25 lines of information/comments at the top of the file, so a simple 'head' command won't be usefule here):
```
# <seq id> <alignment start> <alignment end> <envelope start> <envelope end> <hmm acc> <hmm name> <type> <hmm start> <hmm end> <hmm length> <bit score> <E-value> <significance> <clan> <predicted_active_site_res>

NODE_1_length_46850_cov_69.441152_12      18    138     17    139 PF01850.24  PIN               Domain     2   121   122     60.6   2.2e-16   1 CL0280
NODE_1_length_46850_cov_69.441152_54      76    185     75    193 PF01464.23  SLT               Domain     2   108   117     76.2   1.7e-21   1 CL0037   predicted_active_site
NODE_1_length_46850_cov_69.441152_81      89    224     88    226 PF06924.14  DUF1281           Family     2   140   179     90.4   1.4e-25   1 No_clan
NODE_1_length_46850_cov_69.441152_82      53    125     53    140 PF18406.4   DUF1281_C         Domain     1    72    87     62.4   3.2e-17   1 No_clan
NODE_1_length_46850_cov_69.441152_86      53    232     53    235 PF01555.21  N6_N4_Mtase       Family     1   228   231    110.6   9.7e-32   1 CL0063
NODE_1_length_46850_cov_69.441152_87      52    187     52    187 PF07128.15  DUF1380           Family     1   137   137    170.6   1.9e-50   1 No_clan
NODE_1_length_46850_cov_69.441152_91      51    147     50    147 PF03230.16  Antirestrict      Family     2    92    92     94.2   4.2e-27   1 No_clan
NODE_1_length_46850_cov_69.441152_94      41    145     41    147 PF00436.28  SSB               Domain     1   102   104    125.0   1.3e-36   1 CL0021
NODE_1_length_46850_cov_69.441152_96      61    154     59    155 PF02195.21  ParBc             Family     3    89    90     55.8   4.4e-15   1 CL0248
NODE_1_length_46850_cov_69.441152_97      22    158     21    158 PF06290.14  PsiB              Family     2   138   138    201.7   3.4e-60   1 No_clan
```

The 8th column shows the type of entry that was hit in the pfam database. Let's take a look at Pfam domain "SLT" (accession number PF01464.23). Go to [http://pfam.xfam.org](http://pfam.xfam.org​)​ and enter the accession number for this Pfam domain in the search box.

[IMAGE]

There are a lot of hits to phage domains and domains that manipulate DNA. You might expect this as these sequences have presumably been incorporated into our strain since it diverged from the reference. Also look at Family (the most specific type of hit) from our large contig NODE_2_... is there any evidence for it being a plasmid?

Examine one or two more domains from your results file - is there anything interesting?

# Go to [Task 8](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_3/task_8.md)
