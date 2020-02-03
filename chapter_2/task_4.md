# Task 4 - Aligning Illumina Data to a Reference Sequence

Now that we have checked the quality of our raw data, we can begin to align the reads against a reference sequence. In this way we can compare how the reference sequence and the strain we have sequenced compare.

To do this we will be using a program called 'BWA' (Program: [here](https://github.com/lh3/bwa) & Citation: [Li. H, & Durbin. R. (2009)](https://www.ncbi.nlm.nih.gov/pubmed/19451168)). This uses an algorithm called (unsurprisingly) 'Burrows-Wheeler' to rapidly map reads to a reference genome. 'BWA' also allows for a certain number of mismatches to account for variants which may be present in the sequenced strain vs the reference genome. Unlike other alignment packages such as Bowtie (version 1), BWA allows for insertions or deletions as well.

By mapping reads against a reference, what we mean is that we want to go from a FASTQ file of lots of reads, to another type of file (which we'll describe later) which lists the reads AND where/if that read maps (aligns) against the reference genome. The figure below illustrates what we are trying to achieve here. Along the top in grey is the reference sequence. The coloured sequences below indicate individual sequences and how they map to the reference. If there is a real variant in a bacterial genome we would expect that (nearly) all the reads would contain the variant at the relevant position rather than the same base as the reference genome. Remember that error rates for any single read on second generation platforms tend to be around 0.5-1%. Therefore a 300bp read is on average likely to contain at 2-3 errors.

Now let's look at two potential sources of artefacts.

## 1. Sequencing Error
The region highlighted in green on the image below shows that most reads agree with the reference sequence (i.e. C-base). However, 2 reads near the bottom show an A-base. In this situation we can safely assume that the A-bases are due to a sequencing error rather than a genuine variant since the ‘variant’ has only one read supporting it. If this occurred at a higher frequency however, we would struggle to determine whether it was a genuine variant or an error.

## 2. PCR Duplication
The highlighted region red in the image below shows what appears to be a variant. A C-base is present in the reference and half the reads, whilst an A-base is present in a set of reads which all start at the same position.
 
[IMAGE]

Is this a genuine difference or a sequencing or sample prep error? What do you think? If this was a real sample, would you expect all the reads containing an A to start at the same location?

The answer is probably not. This 'SNP' is in fact probably an artefact of PCR duplication. I.e. the same fragment of DNA has been replicated many times more than the average and happens to contain an error at the first position. We can filter out such reads during after alignment to the reference (see later).

Note that the entire region above seems to contain lots of PCR duplicates with reads starting at the same location. In the case of the region highlighted in red, this will likely cause a false SNP call. The area in green also contains PCR duplicates – the As at these positions are probably either sequencing errors or errors introduced during PCR.

It's always important to think critically about any finding - don't assume that whatever bioinformatic tools you are using are perfect. Or that you have used them perfectly.

## Indexing a Reference Genome
Before we can start aligning reads to a reference genome, the genome sequence needs to be 'indexed' (similarly like a book). This means sorting the genome into easily searched chunks.

# [Task 5]()
