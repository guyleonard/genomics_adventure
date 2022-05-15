# Task 2 - QC the Data!

It is always important to check and understand the quality of the data you are working with:

Change to the directory and run fastqc:

```bash
cd ~/workshop_meterials/genomics_adventure/sequencing_data/pseudomonas_gm41

fastqc &
```

Open the files SRR1042836_subreads.fastq, SRR491287_1.fastq & SRR491287_2.fastq and then look at the reports generated

![FastQC](https://github.com/guyleonard/genomics_adventure/blob/ed2dd272b2241fc00fcae6c224fb3552360230eb/chapter_5/images/chapter_5_task_2_image_1.png)

Note that the quality of the PacBio reads (SRR1042836_subreads.fastq) is much lower than the Illumina reads with a greater than 1 chance in 10 of there being a mistake for most reads.

![FastQC](https://github.com/guyleonard/genomics_adventure/blob/ed2dd272b2241fc00fcae6c224fb3552360230eb/chapter_5/images/chapter_5_task_2_image_2.png)

However, importantly, the length of the PacBio reads is much longer. Yay! :smile:.

# Go to [Task 3](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_5/task_3.md)
