# Chapter Two
## Task 1 - Evaluating the Quality of Illumina Data
From your terminal (command line), navigate to the 'sequencing_data/ecoli' directory (you may be there already) and list the contents of the directory.
```bash
cd genomics_adventure/sequencing_data/ecoli

ls -lath
```
[IMAGE]

These are the paired-end read data that were previously downloaded in the introduction section. The filenames are a bit long and boring to continue to type out - so let's make them a little easier by making a symbolic link to them (this is a type of special file that acts as a pointer to the real file).
```bash
ln -s SRR857279_1.fastq.gz read_1.fastq.gz
ln -s SRR857279_1.fastq.gz read_2.fastq.gz
```

Now we have nice and easy file names to work with :smile:. You should always record or note down what your original filenames are, so that you can refer to the correct data in the future! Just like keeping notes in a lab book. :open_book:

Most programs that work with sequence data require that the 'read 1' and 'read 2' files have the reads in the same order. You can identify reads from the same pair because they will have the same header followed by either a "1" or a "2". We will now look at the raw reads to make sure they look 'correct'. To view the first few headers we can use the 'zcat' command (this is similar to 'cat' but works with zipped files), 'head' to see the top three lines, and then 'grep' to catch the header lines (in this case we know they start with "@HWI" - you may need to check your own files for a suitable header).
```bash
zcat read_1.fastq.gz | head | grep @SRR
zcat read_2.fastq.gz | head | grep @SRR
```

[IMAGE]

The only difference in the headers for the two reads is the read number. Of course, this is no guarantee that all the headers in the file are consistent. To get some more confidence, lets repeat the above commands using 'tail' instead of 'head' to compare reads at the end of the files. This will take a little longer, as tail has to read all the way through to the end of the file!
```bash
zcat read_1.fastq.gz | tail | grep @SRR
zcat read_2.fastq.gz | tail | grep @SRR
```

You can also check that there is an identical number of reads in each file using 'zcat', 'grep' and 'wc -l'. This should take about 4 minutes for each.
```bash
zcat read_1.fastq.gz | grep @SRR | wc –l
zcat read_2.fastq.gz | grep @SRR | wc –l
```

Oops! :trollface: Did you just copy and paste that and receive and error saying "wc: –l: No such file or directory" :stuck_out_tongue_winking_eye:. That's okay! Remember, typing the commands so that you get used to using them, and so that you understand what the options do is a much better way of learning! Don't forget you can use 'Tab complete' to automatically complete filenames. Anyway, the answer you should have received is '4273258. Try again, with the correct command, and see what you get! :hugs: Hint: the "–" was wrong in the above command...

Now, let's run the 'fastqc' program on the raw reads. Runing 'fastqc' as below will bring up the Graphical User Interface version.
```bash
fastqc &
```

Load the **read_1.fastq.gz** file from the ***~/workshop_materials/genomics_adventure/sequencing_data/ecoli*** directory. This should take approximately 6 minutes, so please continue reading below. Many steps in your bioinformatics analysis will take time, you will become a master of multi-tasking, or drinking too much coffee :coffee:.

The 'fastqc' program performs a number of tests which determines whether a green tick (pass) :heavy_check_mark:, exclamation mark (warning) :exclamation:, or a red cross (fail) :x: is displayed. However, it is important to realise that fastqc has no knowledge of what your specific library is or should look like. All of its tests are based on a completely random library with 50% GC content. Therefore, if you have a sample which does not match these assumptions, it may 'fail' the library.

For example, if you have a high AT or high GC organism it may fail the 'per sequence GC content' test. If you have any barcodes or low complexity libraries (e.g. small RNA libraries, RAD-Seq, Amplicons) they may also fail some of the sequence complexity tests. The bottom line is that you need to be aware of what your library is and whether what 'fastqc' is reporting makes sense for that type of library. Know your data from study design to sequencing and beyond!

[IMAGE]

In this case we have a number of errors and warnings which at first sight suggest that there has been a problem - but don't worry too much yet. Let's go through them in turn.

### Quality Scores
This is one of the most important metrics. If the quality scores are poor, either the wrong FASTQ encoding has been guessed by fastqc (see the title of the chart), or the data itself is poor quality. This view shows an overview of the range of quality values across all bases at each position in the FASTQ file.  Generally anything with a median quality score greater than Q20 is regarded as acceptable; anything above Q30 is regarded as 'good'. For more details, see the help documentation in fastqc.

[IMAGE]

In this case this check is red - and it is true that the quality drops off at the end of the reads. It is normal for read quality to get worse towards the end of the read. You can see that at 250 bases the quality is still very good.

### Per tile Sequence Quality
This is a purely technical view on the sequencing run, it is more important for the team running the sequencer. The sequencing flow cell is divided up into areas called cells. The colour of the tiles indicate the read quality and you can see that the quality drops off in some cells faster than others. This maybe because of the way the sample flowed over the flow cell or a mark or smear on the lens of the optics.

[IMAGE]

### Per base Sequence Content
For a completely randomly generated library with a GC content of 50% one expects that at any given position within a read there will be a 25% chance of finding an A,C,T or G base. Here we can see that our library satisfies these criteria, although there appears to be some minor bias at the beginning of the read. This may be due to PCR duplicates during amplification or during library preparation. It is unlikely that one will ever see a perfectly uniform distribution. See [here](http://sequencing.exeter.ac.uk/guide-to-your-data/quality-control/) for examples of good vs bad runs as well as the fastqc help for more details :mag:.

[IMAGE]

### Sequence Duplication Levels
In a library that covers a whole genome uniformly most sequences will occur only once in the final set. A low level of duplication may indicate a very high level of coverage of the target sequence, but a high level of duplication is more likely to indicate some kind of enrichment bias (e.g. PCR over-amplification).
This module counts the degree of duplication for every sequence in the set and creates a plot showing the relative number of sequences with different degrees of duplication. 

[IMAGE]

### Overrepresented Sequences
This checks for sequences that occur more frequently than expected in your data. It also checks any sequences it finds against a small database of known sequences. In this case it has found that a small number of reads 4000 out of 600000 appear to contain a sequence used in the preparation for the library. A typical cause is that the original DNA was shorter than the length of the read - so the sequencing overruns the actual DNA and runs into the adaptors used to bind it to the flow cell.

[IMAGE]

### Other Reports
Have a look at them and at what the author of FastQC has to say [here](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/):mag:. Or check out their youtube tutorial video [here](https://www.youtube.com/watch?v=bz93ReOv87Y).

Remember the error and warning flags are his (albeit experienced) judgement of what typical data should look like. It is up to you to use some initiative and understand whether what you are seeing is typical for your dataset and how that might affect any analysis you are performing.

## Task 1 - Evaluating the Quality of Illumina Data Continued...
Do the same for 'read 2' as you did for 'read 1', have a look at the various plots and metrics which are generated. How similar are they?

Note that the number of reads reported in both files is identical. Overall, both 'read 1' and 'read 2' can be regarded as 'good' data-sets.

## Got to [Task 2]()
