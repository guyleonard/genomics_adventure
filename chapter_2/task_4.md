# Task 4 - Aligning Illumina Data to a Reference Sequence
Now that we have checked the quality of our raw data, we can begin the process of aligning the trimmed reads against a reference sequence. In this way we can compare how the reference sequence, and the strain we are interested in compare.

To do this we will be using a program called 'BWA' (Program: [here](https://github.com/lh3/bwa) & Citation: [Li. H, & Durbin. R. (2009)](https://www.ncbi.nlm.nih.gov/pubmed/19451168) :mag:). This uses an algorithm called 'Burrows-Wheeler' to rapidly map reads to a reference. 'BWA' allows for a certain number of mismatches to account for variants which may be present in the sequences vs the reference. Unlike other alignment packages such as Bowtie (version 1), BWA allows for insertions or deletions as well.

By mapping reads against a reference, what we mean is that we want to go from a FASTQ file of lots of reads and a FASTA file with a set of contiguous sequences, to another type of file (described later) that lists the reads AND where/if those reads map (align) against the reference. Bioinformatics is as much about converting one data/file type to another as wet-lab biology is moving colourless liquids between different plastic tubes. :sweat_smile:

The figure below illustrates what we are trying to achieve. Along the top in grey is the reference sequence. The coloured sequences below indicate individual sequences (reads) and how they are positioned (mapped) to the reference. If there is a real variant in a bacterial genome we would expect that (nearly) all the reads would contain the variant at the relevant position rather than the same base as the reference genome. Remember, the error rates for any single read on "next generation" platforms tend to be around 0.5-1%. Therefore a 300bp read, is on average, likely to contain around 2-3 errors. Not bad!

[IMAGE]

Now let's examine at two potential sources of artefacts:
## 1. Sequencing Error
The region highlighted in green on the image above shows that most reads agree with the reference sequence (i.e. the C-base). However, two reads near the bottom show an A-base. In this situation we can safely assume that the A-bases are due to a sequencing error rather than a genuine variant since the ‘variant’ has only one read supporting it. If this occurred at a higher frequency however, we would struggle to determine whether it was a genuine variant or an error.

## 2. PCR Duplication
The highlighted region red in the image above shows what appears to be a variant. A C-base is present in the reference and half the reads, whilst an A-base is present in a set of reads which all start at the same position.
 
So we must ask, is this a genuine difference or a sequencing or sample preparation error? What do you think? If this was a real sample, would you expect all the reads containing an A-base to start at the same location?

The answer is probably not. This 'SNP' is in fact probably an artefact of PCR duplication - i.e. the same fragment of DNA has been replicated many times more than the average and happens to contain an error at the first position. We can filter out such reads after mapping to the reference (we will do this shortly).

Note that the entire region in the image above seems to contain lots of PCR duplicates with reads starting at the same location. In the case of the region highlighted in red, this will likely cause a false SNP to be called. The area in green also contains PCR duplicates – the A-bases at these positions are probably either sequencing errors or errors introduced during PCR.

It's always important to think critically about any findings - don't assume that whatever bioinformatic tools you are using are perfect. Or that you have used them perfectly.

# Go to [Task 5](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_2/task_5.md)
