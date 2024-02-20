# MyGenome
Analyses for ABT480/CS485G genome assembly

## 1. Analysis of sequence quality
The F1 and R1 sequence datasets were analyzed using FASTQC:
```bash
ssh -Y linkblueID@linkblueID.cs.uky.edu
cd MyGenome
fastqc &
```

Load F1 and R1 datasets into GUI interface.
Take screenshots of output files:

![F1screenshot.png](/data/F1screenshot.png)

## 2. Ran trimmomatic
```bash
java -jar ../sequences/trimmomatic-0.38.jar PE -threads 2 -phred33 -trimlog UFVPY218_errorlog.txt UFVPY218_1.fq.gz UFVPY218_2.fq.gz UFVPY218_1_paired.fastq UFVPY218_1_unpaired.fastq UFVPY218_2_paired.fastq UFVPY218_2_unpaired.fastq SLIDINGWINDOW:20:20 MINLEN:120 ILLUMINACLIP:adaptors.fasta:2:30:10
```

## 3. Count number of forward reads remaining
```bash
grep '@A00261' -c UFVPY218_1_paired.fastq
```
