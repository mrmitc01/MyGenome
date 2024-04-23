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

## 2. Run trimmomatic
```bash
$env:DISPLAY = "localhost:0"
ssh -Y linkblueID@linkblueID.cs.uky.edu
cd MyGenome
java -jar ../sequences/trimmomatic-0.38.jar PE -threads 2 -phred33 -trimlog UFVPY218_errorlog.txt UFVPY218_1.fq.gz UFVPY218_2.fq.gz UFVPY218_1_paired.fastq UFVPY218_1_unpaired.fastq UFVPY218_2_paired.fastq UFVPY218_2_unpaired.fastq SLIDINGWINDOW:20:20 MINLEN:120 ILLUMINACLIP:adaptors.fasta:2:30:10
fastqc &
```

Load four trimmed output files (UFVPY218_1_paired.fastq, UFVPY218_1_unpaired.fastq, UFVPY218_2_paired.fastq, and UFVPY218_2_unpaired.fastq) into GUI interface. Take screenshots of output files:

![Trim1PairedScreenshot.png](/data/UFVPY218_1_paired_fastq_screenshot.PNG)
![Trim1UnpairedScreenshot.png](/data/UFVPY218_1_unpaired_fastq_screenshot.PNG)
![Trim2PairedScreenshot.png](/data/UFVPY218_2_paired_fastq_screenshot.PNG)
![Trim2UnpairedScreenshot.png](/data/UFVPY218_2_unpaired_fastq_screenshot.PNG)

## 3. Count number of forward reads remaining
```bash
grep '@A00261' -c UFVPY218_1_paired.fastq
```
Result: 6955655

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

cd /project/farman_s24cs485g/linblueID/UFVPY218/velvet_UFVPY218_91_111_2_noclean
mv contigs.fa UFVPY218.fasta
perl ../../../SCRIPTs/SimpleFastaHeaders.pl UFVPY218.fasta UFVPY218
```
Output files:
[slurm-21274784.out](/data/slurm-21274784.out), [slurm-21274802.out](/data/slurm-21274802.out), and [19-03-2024-15-17-15_Logfile.txt](/data/19-03-2024-15-17-15_Logfile.txt)

## 5. Check genome completeness using BUSCO
```bash
ssh linkblueID@mcc.uky.edu
cd /project/farman_s24cs485g/linblueID/
cp ../SLURM_SCRIPTS/BuscoSingularity.sh .
nano BuscoSingularity.sh
```
Use nano to add your email address to the mail-user line of the slurm script

#SBATCH --mail-user farman@uky.edu,linkBlueID@uky.edu
```bash
cd /UFVPY218/velvet_UFVPY218_91_111_2_noclean
sbatch ../../BuscoSingularity.sh UFVPY218.fasta
cat slurm-21290995.out
cat busco_916945.log
scp UFVPY218_nh.fasta mrmi248@mrmi248.cs.uky.edu:~/blast/UFVPY218_nh.fasta
```

## 6. BLAST myGenome
```bash
ssh linkblueID@linkblueID.cs.uky.edu
cd blast
blastn -subject UFVPY218_nh.fasta -query MoRepeats.fasta -out MoRepeats.UFVPY218_genomeBLASTn6 -evalue 1e-20 -outfmt 6

scp linkblueID@mcc.uky.edu:/project/farman_s24cs485g/SCRIPTs/CullShortContigs.pl ./CullShortContigs.pl
perl CullShortContigs.pl UFVPY218_nh.fasta
scp linkblueID@mcc.uky.edu:/project/farman_s24cs485g/SCRIPTs/SequenceLengths.pl ./SequenceLengths.pl
perl SequenceLengths.pl UFVPY218_final.fasta
blastn -query MoMitochondrion.fasta -subject UFVPY218_nh.fasta -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid slen length qstart qend sstart send btop' -out MoMitochondrion.UFVPY218.BLAST
awk '$4/$3 > 0.9 {print $2 ",mitochondrion"}' MoMitochondrion.UFVPY218.BLAST > UFVPY218_mitochondrion.csv

ssh linkblueID@mcc.uky.edu
cp /project/farman_s24cs485g/B71v2sh_masked.fasta /project/farman_s24cs485g/linkblueID/B71v2sh_masked.fasta
cd /project/farman_s24cs485g/linkblueID
scp linkblueID@linkblueID.cs.uky.edu:~/blast/UFVPY218_final.fasta ./UFVPY218_final.fasta
scp ./B71v2sh_masked.fasta linkblueID@linkblueID.cs.uky.edu:~/blast/B71v2sh_masked.fasta
mkdir UFVPY218_BLAST
logout
blastn -query B71v2sh_masked.fasta -subject UFVPY218_final.fasta -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' -out B71v2sh.UFVPY218.BLAST

mkdir UFVPY218_BLAST
mv B71v2sh.UFVPY218.BLAST ./UFVPY218_BLAST
scp ./UFVPY218_BLAST/B71v2sh.UFVPY218.BLAST linkblueID@mcc.uky.edu:/project/farman_s24cs485g/linkblueID/UFVPY218_BLAST/B71v2sh.UFVPY218.BLAST
ssh linkblueID@mcc.uky.edu
cd /project/farman_s24cs485g/linkblueID
cp ../SLURM_SCRIPTS/CallVariants.sh ./CallVariants.sh

nano CallVariants.sh
```
Use nano to add your email address to the mail-user line of the slurm script

#SBATCH --mail-user linkblueID@uky.edu
```bash
sbatch CallVariants.sh ./UFVPY218_BLAST
```

## 7. Gene Prediction
```bash
ssh linkblueID@linkblueID.cs.uky.edu
cd genes/snap
echo '##FASTA' | cat B71Ref2_a0.3.gff3 - B71ref2.fasta > B71Ref2.gff3
fathom genome.ann genome.dna -gene-stats
fathom genome.ann genome.dna -categorize 1000
fathom uni.ann uni.dna -gene-stats
fathom uni.ann uni.dna -export 1000 -plus
fathom export.ann export.dna -gene-stats
forge export.ann export.dna
hmm-assembler.pl Moryzae . > Moryzae.hmm
cp ../../blast/UFVPY218_final.fasta UFVPY218_final.fasta
snap-hmm Moryzae.hmm UFVPY218_final.fasta  > UFVPY218_final-snap.zff

fathom  UFVPY218_final-snap.zff  UFVPY218_final.fasta -gene-stats
snap-hmm Moryzae.hmm UFVPY218_final.fasta -gff > UFVPY218_final-snap.gff2
cd ../augustus
screen
augustus --species=magnaporthe_grisea --gff3=on --singlestrand=true --progress=true ../snap/UFVPY218_final.fasta > UFVPY218_final-augustus.gff3

cd ../maker
maker -CTL
nano maker_opts.ctl
```
Use nano to change the following fields:

genome=/home/yourusername/genes/snap/MyGenome.fasta

model_org= #must be set to blank

repeat_protein= #must be set to blank

snaphmm=/home/yourusername/genes/snap/Moryzae.hmm

augustus_species=magnaporthe_grisea

keep_preds=1

protein=/home/yourusername/genes/maker/genbank/ncbi-protein-Magnaporthe_organism.fasta
```bash
ls -l ~ | tee listing.txt
screen
maker 2>&1 | tee maker.log

gff3_merge -d UFVPY218_final.maker.output/UFVPY218_final_master_datastore_index.log -o UFVPY218_final-annotations.gff
fasta-merge -d UFVPY218_final.maker.output/UFVPY218_final_master_datastore_index.log -o UFVPY218_final-genes.fasta
```
