# Task 3 - Random Subsampling and Digital Normalisation
Sometimes you may find that you have sequenced too much data! :open_mouth: However, this is not as bad a place to be in as it first seems, but it may put limits on your infrastructure - for example, your computer/server may not be able to keep all the data in memory to assemble it.

There are several methods of reducing your data, here we will discuss two:

### 1. Random Subsampling
As easy as the naming suggests, we take a random subsampling of the original dataset, e.g. 10% of the data and then we use that data to perform an assembly, but importantly, the original data to do mapping.

For this we will use the program ['seqtk'](https://github.com/lh3/seqtk) :mag: which is an excellent little toolkit of all kinds of FASTA/Q processing ability.

Check it out..
```bash
seqtk

# and
seqtk sample
```

Now let's try randomly subsampling our *E. coli* 'read_1' dataset to 10% of the data. We need to express 10% as a fraction for 'seqtk', so 0.1, and supply that as an option to the 'seqtk sample' program. Let's do it twice, to see if we truly do get a random subsampling of the data...
```bash
seqtk sample read_1.fq.gz 0.1 > read_1_subsample_one.fq

seqtk sample read_1.fq.gz 0.1 > read_1_subsample_two.fq

head read_1_subsample*
```

What do the headers look like? Do they look random? It doesn't look like it huh? So what is going on? 'seqtk sample', by default, sets a 'seed' value of '11' - so in this case both files are the same random subsampled selection! dizzy_face:

This may seem somewhat counter-intertuitive at first - how can something be random if you can repeat it exactly? Well, the program uses a special trick (called [Reservoir sampling](https://en.wikipedia.org/wiki/Reservoir_sampling) :mag:) that can take the same 'random' sample when given a starting 'seed' - in this case '11' - which allows it to start the random selection at the same place each time. Change the 'seed' and you will get a different subsampled selection! Neat. :sunglasses:

Using the '-s' option allows us to set the 'seed' starting position of the random subsampling. So, if we change the value to '1234' and '5678' we should have two different random samples. But you can choose any numbers you like.
```bash
seqtk sample -s 1234 read_1.fq.gz 0.1 > read_1_subsample_three.fq

seqtk sample -s 5678 read_1.fq.gz 0.1 > read_1_subsample_four.fq

head read_1_subsample*
```

How do they look now? Much better right!? :thumbsup:

There's just one final thing we need to do, although 'seqtk' can read 'gzipped' files, it does not produce them. So we should 'gzip' them now. You can use the program 'gzip' or it's speedier cousin 'pigz'.
```bash
pigz XX.fq
```
### 2. Digital Normalisation


## Contaminant Checking
A number of tools are available which also enable to you to quickly search through your reads and assign them to particular taxa or taxonomic group. These can serve as a quick check to make sure your samples or libraries are not contaminated with DNA from other sources. If you are performing a de-novo assembly, for example, and have DNA sequences present from multiple organisms, you will risk poor results and chimeric contigs.

Some ‘contaminants’ may turn out to be inevitable by-products of sampling and DNA extraction, and this is often the case with algae, and/or other symbionts but some groups have made amazing discoveries such as the discovery of a third symbiont (which turned out to be a yeast) in lichen, [here](http://science.sciencemag.org/content/353/6298/488.full).

Some tools you can use to check the taxonomic classification of reads include:
 * [Kraken](https://ccb.jhu.edu/software/kraken2/)
 * [Centrifuge](https://ccb.jhu.edu/software/centrifuge/)
 * [Blobology](https://blobtoolkit.genomehubs.org/)
 * Blast (in conjunction with sub-sampling your reads) and Krona to plot results
 * and many more!

# [Task 4]()
