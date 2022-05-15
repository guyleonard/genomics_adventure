# Task 4 - View Assembly in IGV

You are by now experts at this bit, load up IGV from the console or desktop adn click "Genomes --> Load Genome from File..."

We are going to import the contigs we have assembled as the reference. Unlike the reference genome though, we have no annotation available. Make sure you select the contigs.fasta file for
the complete de novo assembly (not the unmapped reads assembly).

Once loaded, click on "File --> Load From File..." and select the contigs_mapped.sorted.bam file. Again mak sure you load the file in the "remapping_to_assembly" directory.

Once loaded, explore some of the contigs in IGV. See if you can find anything unusual or interesting  in any of the contigs.

Here's one to get you started... Zoom in on "NODE_1..." until you can see the reads. Then go to the start of the contig and have a look at the reads.

[IMAGE]

What do he colours of the reads mean? Can you guess before you read the [manual](https://software.broadinstitute.org/software/igv/interpreting_insert_size) 🔍?

So, at the start of the contig we have a lot of reads that map to one contig but their read-pair maps to another contig. This is suggestive of repetetive regions, and if an assembler cannot resolve these repetitive regions with paired-end reads or coverage information, it will generally be unable to assemble any further sequence for that contig. Therefore it is quite common to see contigs which start and end in sequence which is repeated elsewhere.

This part is open-ended and depends on how much you want to explore your genome. There's another task waiting though!

# Go to [Task 5](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_4/task_5.md)
