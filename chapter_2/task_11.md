# Task 11 - IGV - Graphical View of Your Alignments
The Integrative Genome Viewer (IGV) is a tool developed by the Broad Institute for browsing interactively the alignment data you produced. It has a huge wealth of features :dizzy_face: and we will only cover some of the basics to get you started. You can read more more information [here](http://www.broadinstitute.org/igv/).

To run IGV:
```bash
igv
```

The IGV viewer should appear, as below, and you will notice that by default a human genome has been selected in the first drop down box.

[IMAGE]

## Importing the *E. coli* Reference Sequence
By default IGV does not contain our reference genome. We'll need to import it. Therefore, click on 'Genomes -> Load Genome from File...'

Navigate to the location of your sequence file, as in the image below, select it and click "Open". You should see the drop down menu change to the name of your file.

[IMAGE]

## Importing the *E. coli* Reference Annotation
We can also load the annotation (.gff). Click on 'File -> Load from File...'. Navigate to the location of your sequence file, as in the image below, select it and click "Open". You should see a new 'track' of blue boxes with white arrows appear.

[IMAGE]

## Load the BAM Alignment
Load the alignment file (.bam). Note that IGV requires the .bai index file to also be in the same directory (we generated that earlier with samtools index). Select 'File...' and 'Load From File...'. Select your BAM file and click 'Open'.

[IMAGE]

Once loaded your screen should look similar to the image below (you may need to use the zoom '+' tool to see the features you just loaded). Note that you can load more BAM files if you wish to compare different samples, sequencing technologies or the results of different mapping programs.

[IMAGE]

Use the +/- keys to zoom in or use the zoom bar at the top right of the screen to zoom into about 1-2kbases as above.

Right click on the main area and select view as pairs.

[IMAGE]

The gray graph at the top of the figure indicates the coverage of the genome.

[IMAGE]

The more reads mapping to a certain location, the higher the peak on the graph. You'll see a coloured line of blue, green or red in this coverage plot if there are any SNPs (single-nucleotide polymorphisms) present. If there are any regions in the genome which are not covered by the reads, you will see these as gaps in the coverage graph. Sometimes these gaps are caused by repetitive regions; others are caused by genuine insertions/deletions in your new strain with respect to the reference.

Below the coverage graph is a representation of each read pair as it is mapped to the genome. One pair is highlighted.

[IMAGE]

This pair consists of 2 reads with a gap (there may be no gap if the reads overlap) Any areas of mismatch either due to inconsistent distances between paired-end reads or due to differences between the reference and the read and are highlighted by a colour. The brighter the colour, the higher the base-calling quality is estimated to be. Differences in a single read are likely to be sequencing errors. Differences consistent in all reads are likely to be mutations.

Hover over a read to get detailed information about the reads' alignment.

[IMAGE]

You don't need to understand every value, but compare this to the SAM format to get an idea of what is there.

# Go to [Task 12](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_2/task_12.md)
