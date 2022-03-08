#!/bin/bash
# Quality control of raw data with fastqc and removal of adapters.
for f1 in *.fastq
do
     f2=${f1%%.fastq}"_trim.fastq"
   
     echo $f1
     echo $f2
     
     java -jar /home/soft/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 6 -phred33 -trimlog trim.log data/$f1 $f2 ILLUMINACLIP:adapters:2:30:7 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:16
     mv $f2 fastq_trim 
     fastqc $f1
done

# Reads are mapped against the reference genome.
for f1 in *_trim.fastq
do	
	f2=${f1%%_trim.fastq}".fa"
	f3=${f1%%_trim.fastq}".arf"
	f4=${f1%%_trim.fastq}".log"

	mapper.pl fastq_trim/$f1 -e -h -j -m -p ovis_aries_rambouillet_v1 -q -s reads/$f2 -t reads_against_genome/$f3 &> mapper_logs/$f4
done

# Merge the outputs of the mapper.pl
cat reads/*.fa > all_samples_reads.fa
cat reads_against_genome/*.arf > all_samples_reads_against_genome.arf

#The inputs for this command line are the merged outputs of the mapper.pl (the merged fasta of reads of all samples and the merged file of reads against the genome of all samples). Run de core command of mirdeep2.
miRDeep2.pl all_samples_reads.fa Ovis_aries_rambouillet.Oar_rambouillet_v1.0.dna.toplevel.fa all_samples_reads_against_genome.arf sheep_mature_mirna.fa mammals_mature_mirna.fa sheep_premirna.fa -P 2> mirdeep2.log 	

