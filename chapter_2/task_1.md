# Chapter Two
## Task 1 - Evaluating the Quality of Illumina Data
From your terminal window, navigate to the 'sequencing_data' directory (you may be there already) and list the contents of the directory.
```bash
cd genomics_adventure/sequencing_data

ls -lath
```
[IMAGE]

These are paired-end read data that we downloaded in the introduction. Many programs require that the read 1 and read 2 files have the reads in the same order. You can identify reads from the same pair because they will have the same header followed by either a "1" or a "2".  We will now look at the raw reads. To view the first few headers we can use the 'zcat', 'head' and 'grep' commands:
```bash
zcat ERR2789854_1.fastq.gz | head | grep @HWI
zcat ERR2789854_1.fastq.gz | head | grep @HWI
```
[IMAGE]

The only difference in the headers for the two reads is the read number. Of course this is no guarantee that all the headers in the file are consistent. To get some more confidence, repeat the above commands using 'tail' instead of 'head' to compare reads at the end of the files. 
```bash
zcat ERR2789854_1.fastq.gz | tail | grep @HWI
zcat ERR2789854_1.fastq.gz | tail | grep @HWI
```
You can also check that there is an identical number of reads in each file using 'cat', 'grep' and 'wc -l':
```bash
zcat ERR2789854_1.fastq.gz | grep @MISEQ | wc –l
zcat ERR2789854_2.fastq.gz | grep @MISEQ | wc –l
```
Ha ha! :trollface: Did you just copy and paste that?! Did it not work? :stuck_out_tongue_winking_eye: Remember, typing the commands so that you get used to using them, and so that you understand what they do is a much better way of learning! Anyway, the answer is XXX. Try again and see what you get! :hugs:

Now, let's run the 'fastqc' program on the data. Unlike the QC lab, we will open up a Graphical User Interface (GUI) and load the data this way. To do this, run:
```bash
fastqc &
```
Load the **ERR2789854_1.fastq.gz** file from the ***~/workshop_materials/genomics_adventure/sequencing_data directory***.
The fastqc program performs a number of tests which determines whether a green tick (pass) :heavy_check_mark:, exclamation mark (warning) :grey_exclamation:, or red cross (fail) :x: is displayed. However, it is important to realise that fastqc has no knowledge of what your library is or should look like. All of its tests are based on a completely random library with 50% GC content. Therefore if you have a sample which does not match these assumptions, it may 'fail' the library. For example, if you have a high AT or high GC organism it may fail the 'per sequence GC content' test. If you have any barcodes or low complexity libraries (e.g. small RNA libraries, RAD-Seq, Amplicons) they may also fail some of the sequence complexity tests.
The bottom line is that you need to be aware of what your library is and whether what fastqc is reporting makes sense for that type of library. Know your data from study design to sequencing!
[IMAGE]

In this case we have a number of errors and warnings which at first sight suggest that there has been a problem - but don't worry too much yet. Let's go through them in turn.

### Quality Scores
This is one of the most important metrics. If the quality scores are poor, either the wrong FASTQ encoding has been guessed by fastqc (see the title of the chart), or the data itself is poor quality. This view shows an overview of the range of quality values across all bases at each position in the FASTQ file.  Generally anything with a median quality score greater than Q20 is regarded as acceptable; anything above Q30 is regarded as 'good'. For more details, see the help documentation in fastqc.

