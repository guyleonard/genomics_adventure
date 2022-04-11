# Task 12 - SNPs and Indels

*The following 3 tasks are open-ended. Please take your time with these but don't forget there are more adventures to be had in the other chapters following this!*

## Alignment Display Format
Visit  the BROAD's website [here](http://www.broadinstitute.org/software/igv/AlignmentData) :mag: and have a brief read of the format.

## Manually Identify a Region Without any Reads Mapping
These regions can be quite difficult to find even with a very small genome. Zoom out as far as you can, but where you are still able to see the read mappings. You could try and use the coverage plot from QualiMap [Task 10](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_2/task_10.md) to try to find one... Are there any genes associated with this region?

Due to the way IGV handles BAM files, it will not display coverage information if you zoom out too far. To be able to see this coverage information across the entire genome, regardless of how far you are zoomed out, you’ll need to create a 'TDF file' which contains coverage information across windows of 'X' number of bases on the genome.

You can do this within IGV:

> Select Tools->Run igvtools

[IMAGE]

> Load the BAM alignment file in the Input field and click Run

[IMAGE]

Once completed, close the igvtools window and then you can load this TDF file:

> Select 'File -> Load...' from the 'File' menu.

[IMAGE]

You should now see the extra coverage track which remains visible even after you zoom out! What a great trick! :ok_hand: 

## Manually Identify a Region Containing Repetitive Sequences.
You could also try using the QualiMap reports from before [Task 10](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_2/task_10.md) to give you an idea where to look. What could these regions represent? Are there any genes close-by? (Think about repetitive sequences, what does BWA do if a region in the genome has been duplicated?).

# Go to [Task 13](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_2/task_13.md)
