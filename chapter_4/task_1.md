# Chapter 4 ​-  *De novo* Assembly Using Short Reads

## Introduction
In this section of our adventure we will continue the analysis of a strain of ​* E.coli*. In the previous chapter we extracted the reads which did not map to the reference genome and assembled them. However, it is often necessary to be able to perform a *de novo* assembly of a genome. In this case, rather than doing any mapping, we will start with the filtered reads we obtained in Chapter 3 of the adventure.

To do this we will use the program SPAdes to try to get the best possible assembly for a given genome. We will then generate assembly statistics and produce some annotation via Pfam and BLAST. Easy stuff! You should be a natural at this by now...

## Stop. Assembly Time!

Well, the assembly will take a long time - so we have kindly precomputed it for you - nonetheless we shall step through the process. Here we go!

Let's get our files and folders in order, return to the main directory and create a new one for this work, and make a link to the pre-computed assembly.
```bash
cd ~/workshop_materials/genomics_adventure

mkdir denovo_assembly && cd denovo_assembly

ln -s ../../precomputed/denovo_assembly/assembly .
```

The command you would use to run SPAdes is this (do not run, it would take XX hours):
```bash
spades.py -o assembly -1 ../sequencing_data/ecoli/read_1_val_1.fq.gz -2 ../sequencing_data/ecoli/read_2_val_2.fq.gz
```

### Recap on Assemby Theory
We will be using an assembler called [SPAdes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342519/). It generally performs pretty well with a variety of genomes. It can also incorporate longer reads produced from PacBio sequencers that we will use later in our adventure.

One big advantage is that it is not just a pure assembler - it is a suite of programs that prepare the reads you have, assembles them and then refines the assembly.

SPAdes runs the modules that are required for a particular dataset and it produces the assembly with a minimum of preparation and parameter selection - making it very straightforward to produce a decent assembly. As with everything in bioinformatics you should try to assess the results critically and understand the implications for further analysis.

A pretty old overview of how SPAdes differs from 'velvet', a very old assembly program, can be found [here](http://thegenomefactory.blogspot.co.uk/2013/08/how-spades-differs-from-velvet.html). Nonetheless, it outlines the overall process quite nicely:

1. Read error correction based on k-mer frequencies using ​BayesHammer.
2. De Bruijn graph assembly at ​multiple ​k-mer sizes, not just a single fixed one.
3. Merging of different k-mer assemblies (good for varying coverage).
4. Scaffolding of contigs from paired end/mate pair reads.
5. Repeat resolution from paired end/mate pair data using rectangle graphs.
6. Contig error correction based on aligning the original reads with ​BWA​ back to contigs.

#### K-mers
Rather than store all reads individually which would be unfeasible for Illumina type datasets, de Bruijn assemblers convert each read to a series of k-mers and stores each k-mer once, along with information about how often it occurs and which other k-mers it links to. A short k-mer length (e.g. 21) reduces the chance that data will be missed from an assembly (e.g. due to reads being shorter than the k-mer length or sequencing errors in the k-mer), but can result in shorter contigs as repeat regions cannot be resolved.

For a genomic assembly you want to try to obtain the lowest number of contigs, with the longest length, with the fewest errors. However, although numbers of contigs and longest lengths are easy to evaluate, it is extremely difficult to know what is or isn't an error when sequencing a genome for the first time.

SPAdes allows you to choose more than one k-mer length - it then performs an assembly for each k-mer and merges the result - trying to get the best of both worlds. It actually has some pre-calculated k-mer settings based on the length of reads you have, so you don't even have to choose them.

Let's look at the assembly process in more detail, let's say you have a single read "AACTAACGACGCGCATCAAAA". The set of k-mers, with length 6 (i.e. 6-mers), obtained from this read, would be created by taking the first six bases, then moving the window along one base, taking the next 6 bases and so-on until the end of the read. For example, "AACTAA", followed by "ACTAAC", then "CTAACG", "TAACGA", "TAACGAC" and so on... You may well ask, “So what? How does that help”? For a single read, it really doesn't help. However, let's say that you have another read which is identical except for a single base. Rather than represent both reads separately, we need only store the k-mers which differ and the number of times they occur. Note the 'bubble' like structure which occurs when a single base-change occurs. This kind of representation of reads is called a 'k-mer graph' (sometimes inaccurately referred to as a de Bruijn graph).

[IMAGE]

Now let's see what happens when we add in a third read. This is identical to the first read except for a change at another location. This results in an extra dead-end being added to the path.

[IMAGE]

The job of any k-mer based assembler is to find a path through the k-mer graph which correctly represents the genome sequence.

#### Coverage Cut-off
In the figure above, you can see that the coverage of various k-mers varies between 1x and 3x. The question is which parts of the graph can be trimmed or removed so that we avoid any errors. As the graph stands, we could output three different contigs as there are three possible paths through the graph. However, we might wish to apply a coverage cutoff and remove the top right part of the graph because it has only 1x coverage and is more likely to be an error than a genuine variant.

In a real graph you would have millions of k-mers and thousands of possible paths to deal with. The best way to estimate the coverage cutoff in such cases is to look at the frequency plot of contig (node) coverage, weighted by length. In the example below you can see that contigs with a coverage below 7x or 8x occur very infrequently. As such it is probably a good idea to exclude those contigs which have coverage less than this – they are likely to be errors.

[IMAGE]

In the example below you can see a stretch of DNA with many reads mapping to it. There are two repetitive regions A1 and A2 which have identical sequence. If we try to assemble the reads without any knowledge of the true DNA sequence, we will end up with an assembly that is split into two or more contigs rather than one.

One contig will contain all the reads which did not fall into A1 and A2. The other will contain reads from both A1 and A2. As such the coverage of the repetitive contig will be twice as high as that of the non-repetitive contig.

If we had 5 repeats we would expect 5x more coverage relative to the non-repetitive contig. As such, provided we know what level of coverage we expect for a given set of data, we can use this information to try and resolve the number of repeats we expect.

# Now go to [Task 2]()
