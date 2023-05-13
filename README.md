# A Genomics Adventure
This tutorial is based on a workshop that was written many years ago by Konrad Paszkiewicz, whilst he was the head of Exeter University's Sequencing Service. Over the years it has been extensivley updated and modified by several colleages (listed below), taught at Exeter Universtity, and also at the [evomics.org](https://evomics.org/) Workshop on Genomics since 2014. It has now been reproduced here, in a slightly new form - again with many updates.

Many thanks to:
 * [Konrad Paszkiewicz](https://scholar.google.com/citations?user=yrHDETIAAAAJ&hl=en), CTO Hummingbird Biosciences.
 * [Sophie Shaw](https://scholar.google.se/citations?user=_K3aFRYAAAAJ&hl=en), Bioinformatician, All Wales Medical Genomics Service. 
 * [Josie Paris](https://www.josieparis.com/), MSCA Postdoctoral Fellow, Universita' Politecnica delle Marche, Italy.
 * Workshop on Genomics, [evomics.org](https://evomics.org/)

## A Few Notes on Styles
Throughout this adventure you will see various styles of text. Mostly they will be in the form of the main story - similar to what you are reading now - but some will be links to further and expanded reading (indicated with a magnifying glass - :mag: - this is mostly optional reading for the super curious :nerd_face:). However, there will also be commands for you to type. These will look something like:
```bash
# Comment line - do not type this
command -options
```

It will always be indicated when you need to run a command or just look at them. We also try to heavily discourage you from using 'copy and paste'. This is because we feel that you will learn much more by actually typing the commands yourself, and making mistakes. In fact, some of the commands are designed to intentionally fail if you copy and paste them. :stuck_out_tongue_closed_eyes: 

Often, you will see commands written in the style below, that is with a back-slash '\' at the end of lines, this is simply to break up the command for ease of reading. You would not normally type commands like this. You can see the actual command as you would type it below too.
```bash
# Large command that you can't fully see without scrolling
bwa mem -t 4 ~/genomics_adventure/ecoli/GCF_000005845.2_ASM584v2_genomic.fna ~/genomics_adventure/ecoli/XX_val1.fq.gz ~/genomics_adventure/ecoli/XX_val2.fq.gz > XXX.sam

# Much easier to read command
bwa mem -t 4 \
~/genomics_adventure/ecoli/GCF_000005845.2_ASM584v2_genomic.fna \
~/genomics_adventure/ecoli/XX_val1.fq.gz \
~/genomics_adventure/ecoli/XX_val2.fq.gz \
> XXX.sam
```

The tutorial (sorry adventure), like any good story, is designed to be long and you likely won't finish it in less than 3 hours or even in one sitting. Don't worry though, you can always come back to it anytime!

<p align="center">:dragon_face: Shall we begin? - Daenerys Targaryen :dragon_face:</p>

## Chapter 1: Once Upon a Time...
Before our adventure begins we will need to choose our starting point:

1. [Workshop on Genomics](workshop_on_genomics.md)
2. [Home](home.md)
