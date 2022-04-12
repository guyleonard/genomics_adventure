# Task 2

## Task 2a - Check Reads
Check that the number of entries in both fastq files is the same. Also check that the last few entries in the read 1 and read 2 files have the same header (i.e. that they have been correctly paired).
<details>
  <summary>Try before you check here</summary>
  There are multiple ways to do this. Here's how we did it:
  
  ```
  grep -c "^@SRR" unmapped_r1.fastq unmapped_r2.fastq 
  unmapped_r1.fastq:56710
  unmapped_r2.fastq:56710

  tail -n 4 unmapped_r1.fastq unmapped_r2.fastq 
  ==> unmapped_r1.fastq <==
  @SRR857279.4273239/1
  GTATAAATCTTGCCGTCATTCTGATCAGTTTGTAACATTCTGTAATGATCACCATTGGCTGGCGATTTTTCTGTTCAGTAATGTAATTAACCTTATCTGATGCGCTGGCCACTATTCCATCAGCTGTACTGATGGCAGGCTCCCTGTTG
  +
  ?????B/B?BBBBBBBCEFFFFHHHHFFFHHHFHHHHHHHHFGHHHHHHHHHHHHHDGHHHHHHHHHHFHHHHHHHHHDGEGGGGHHIHHHHHHHHHHGFHHHHHHBCEHHHHFFFHFHHHFHFFF,@DFFD,BDFFD=@DDEEEEEE@

  ==> unmapped_r2.fastq <==
  @SRR857279.4273239/2
  CGTCATTGCCGCCCCTGCCAGGGACTATATCGACCTGACCCTTGATCAGTTTCCAGGCTATCATAACCGGATTGTGGCAGAGCCCGTTGAGTCCGGCGGACAGCTCGCGGCAGACCTCAACAGGGAGCCTGCCATCAGTACAGCTGATGG
  +
  ?A???BBBDDBDDDDDGGGGGGHHHHIIHFHFHHHHFHHHHIIIIIIIIHHIIIFGHIHHHIIIIIHIHHHHHIHHHIHHHHGHHHHHHHHHHHHHGGGEGGGGGGED2;DGGGGGEGGGCC?CEGCCEEGGEEEGGECCC??C:CEGEC
  ```

  How did your methods compare? 
</details>

## Task 2b - Evaluate QC
Back to some familiar territory! Use the fastqc program to look at the statistics and QC for the ​unaligned_r1.fastq and unaligned_r2.fastq​ files.

Do these look reasonably good? Remember, some reads will fail to map to the reference because they are poor quality, so the average scores will be lower than the initial fastqc report we did in the mapping task. The aim here is to see if it looks as though there are reads of reasonable quality which did not map.

Assuming these reads look ok to you, we will proceed with preparing them for de novo assembly!

# Go to [Task 3](https://github.com/guyleonard/genomics_adventure/blob/release/chapter_3/task_3.md)
