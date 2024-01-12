# Task 3 - Random Subsampling and Digital Normalisation
Sometimes you may find that you have sequenced too much data! :open_mouth: However, this is not as bad a place to be in as it first seems, e.g it may put limits on your infrastructure; for example, your computer/server/HPC may not be able to keep all the data in memory to assemble it. Or you may want to do some quick analyses to sanity check your data without doing a full analysis.

There are several methods of reducing your data, here we will discuss two of them:

### 1. Random Subsampling
As easy as the name suggests, we take a random subsample of the original dataset, e.g. 10% of the data and then we can use that data to perform an assembly, although we should really use all the original data to do any mapping. Why do you think this might be the case?

For this task we will use the program ['seqtk'](https://github.com/lh3/seqtk) :mag: which is an excellent little toolkit to do FASTA/Q processing.

Check it out:
```bash
seqtk

seqtk sample
```

Now let's try randomly subsampling our *E. coli* 'read_1' dataset to 10% of the data. We need to express 10% as a fraction for 'seqtk', so 0.1, and supply that as an option to the 'seqtk sample' program. Let's do it twice, to see if we truly do get a random subsampling of the data...
```bash
seqtk sample read_1_val_1.fq.gz 0.1 > read_1_val_1_subsample_one.fq

seqtk sample read_1_val_1.fq.gz 0.1 > read_1_val_1_subsample_two.fq

head read_1_val_1_subsample*
```

What do the headers look like, do they look like a random selection to you? Ha ha! :stuck_out_tongue_closed_eyes: I played another trick on you! Don't worry though, this is a common pitfall when you start out learning to use bioinformatic software.

So, what is going on? Well 'seqtk sample', by default, sets a special 'seed' value of '11' - it is sneakily there in the help section - and so in this case both files happen to be the same random subsampled selection! :dizzy_face:

This may seem somewhat counter-intertuitive at first - how can something be random if you can repeat it exactly - however it is immensely useful, as you will come to see. The program uses a special trick (called [Reservoir sampling](https://en.wikipedia.org/wiki/Reservoir_sampling) :mag: - you don't need to know this, just to be aware) that will take the same 'random' sample when given a starting 'seed' - in this case '11' - allowing it to start the random selection process at the same place each time. Change the 'seed', and you will get a different subsampled selection! Neat. :sunglasses:

Using the '-s' option we can set the 'seed' starting position of the random subsampling. So, if we change the value to '1234' and '5678' we should have two different random samples. For this next bit you can choose any numbers you like! :metal:
```bash
seqtk sample -s 1234 read_1_val_1.fq.gz 0.1 > read_1_val_1_subsample_three.fq

seqtk sample -s 5678 read_1_val_1.fq.gz 0.1 > read_1_val_1_subsample_four.fq

head read_1_val_1_subsample_three.fq

head read_1_val_1_subsample_four.fq
```

How do they look now? Much more random right!? :thumbsup: Of course, if you are doing this in groups (or on multiple machines) the random samples should be the same give the same seed.

Let's tidy up our testing files now, and then repeat the process one last time - we will need this next set of files later on in our adventure. Please use the 'seed' value that is indicated, or your analyses will look very different to the examples shown here.
```bash
rm *.fq

seqtk sample -s 628 read_1_val_1.fq.gz 0.5 > read_1_val_1_subsampled.fq

seqtk sample -s 628 read_2_val_2.fq.gz 0.5 > read_2_val_2_subsampled.fq
```

There is just one final thing we need to do, although 'seqtk' can read 'gzipped' files, it does not produce them. So we should 'gzip' them now. You can either use the program 'gzip' or it's speedier cousin 'pigz'.
```bash
pigz *.fq
```

One thing to be careful of when randomly subsampling data like this is that we may, by accident, subsample regions in our genome that are more, or less, sequenced than others - especially if our coverage is not even across the sequencing libraries. This can lead to issues with downstream analysis, especially if for example, the libraries come from Multiple Displacement Amplification (e.g. Single Cells). So we should always be careful with this approach. One method to minimise this potential source of bias in our data is called Digital Normalisation.

### 2. Digital Normalisation

In all seriousness, I forgot to write this section. Mea culpa. So as it's 2024 and we are all using ChatGPT anyway, I gave it the prompt "write me a few short paragraphs about digital normalisation of fastq illumina reads".

Let's have a poquito look at what it came up with.

>Digital normalization is an intriguing concept in the realm of bioinformatics, particularly when dealing with high-throughput sequencing data, such as that produced by Illumina platforms. The primary goal of digital normalization is to reduce the redundancy in sequencing data, thereby improving both computational efficiency and downstream analysis quality.

In the context of fastq files from Illumina reads, digital normalization operates by selectively down-sampling the dataset. It's akin to strategically thinning out a dense forest so that each tree (read) can be examined more effectively. The process involves analyzing the k-mer frequency distribution within the reads. K-mers are subsequences of length 'k' derived from the longer DNA sequences. By assessing the abundance of these k-mers across the dataset, digital normalization identifies and discards reads that are over-represented, presumably because they originate from highly abundant sequences.

This technique is particularly beneficial when dealing with large datasets, a common scenario given the high output of Illumina platforms. By reducing the number of reads, digital normalization lessens the computational resources needed for tasks like assembly and alignment, without significantly compromising the quality or coverage of the data. It's a bit like efficiently packing a suitcase for a trip; you want to include everything necessary but avoid unnecessary bulk.

A noteworthy aspect of digital normalization is its impact on error reduction. Redundant reads in high-throughput sequencing often include errors. By reducing redundancy, digital normalization can indirectly decrease the error rate, leading to cleaner, more accurate datasets. This is crucial for bioinformatic analyses, where precision is paramount.

However, it's important to approach digital normalization with a bit of caution. While it effectively reduces data size and computational load, it might inadvertently eliminate low-abundance reads that could be biologically significant. Thus, the parameters for k-mer thresholding and read retention must be chosen carefully to strike a balance between efficiency and data integrity.

In summary, digital normalization is a smart and efficient way to handle the deluge of data from Illumina sequencing, ensuring that the focus remains on the most informative and relevant parts of the genomic jigsaw puzzle.


Not too bad! We won't use this technique today but a good tool is [BBNorm](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbnorm-guide/) 🔍 from the BBTools package.

## Contaminant Checking
A number of tools are available which also enable to you to quickly search through your reads and assign them to particular taxa or taxonomic groups. These can serve as a quick check to make sure your samples or libraries are not contaminated with DNA from other sources. If you are performing a de-novo assembly, for example, and have DNA sequences present from multiple organisms, you will risk poor results and chimeric contigs.

Some ‘contaminants’ may turn out to be inevitable by-products of sampling and DNA extraction, and this is often the case with algae, and/or other symbionts but some groups have made amazing discoveries such as the discovery of a third symbiont (which turned out to be a yeast) in lichen, see [here](http://science.sciencemag.org/content/353/6298/488.full) :mag:.

We won't cover this topic in our current adventure, but here is a list of some tools you can use to check the taxonomic classification of reads:
 * [Kraken](https://ccb.jhu.edu/software/kraken2/)
 * [Centrifuge](https://ccb.jhu.edu/software/centrifuge/)
 * [Blobology](https://blobtoolkit.genomehubs.org/)
 * Blast (in conjunction with sub-sampling your reads) and Krona to plot results
 * and many more!

# Go to [Task 4](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_2/task_4.md)
