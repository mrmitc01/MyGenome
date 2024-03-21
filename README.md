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

![F1screenshot.png](/data/UFVP218_1_screenshot.PNG)
![R1screenshot.png](/data/UFVP218_2_screenshot.PNG)

## 2. Ran trimmomatic
```bash
java -jar ../sequences/trimmomatic-0.38.jar PE -threads 2 -phred33 -trimlog UFVPY218_errorlog.txt UFVPY218_1.fq.gz UFVPY218_2.fq.gz UFVPY218_1_paired.fastq UFVPY218_1_unpaired.fastq UFVPY218_2_paired.fastq UFVPY218_2_unpaired.fastq SLIDINGWINDOW:20:20 MINLEN:120 ILLUMINACLIP:adaptors.fasta:2:30:10
```

## 3. Count number of forward reads remaining
```bash
grep '@A00261' -c UFVPY218_1_paired.fastq
```

## 4. Assemble myGenome
```bash
ssh linkblueID@mcc.uky.edu
cd /project/farman_s24cs485g/linkblueID
cp ../SLURM_SCRIPTS/velvetoptimiser_noclean.sh .
nano velvetoptimiser_noclean.sh
```
Use nano to add your email address to the mail-user line of the slurm script

#SBATCH --mail-user farman@uky.edu,linkBlueID@uky.edu
```bash
sbatch velvetoptimiser_noclean.sh UFVPY218 61 131 10
cat slurm-21274784.out

sbatch velvetoptimiser_noclean.sh UFVPY218 91 111 2
cat slurm-21274802.out
```
