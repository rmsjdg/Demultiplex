#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=rj.dxNone

DB=/projects/bgmp/shared/2017_sequencing

/usr/bin/time -v /projects/bgmp/roj/bioinfo/Bi622/Demultiplex/Assignment-the-third/demultiplex.py \
 -a ${DB}/1294_S1_L008_R1_001.fastq.gz \
 -b ${DB}/1294_S1_L008_R2_001.fastq.gz \
 -c ${DB}/1294_S1_L008_R3_001.fastq.gz \
 -d ${DB}/1294_S1_L008_R4_001.fastq.gz \
 -i ${DB}/indexes.txt \
 -o /scratch/bgmp/roj/demux/qNone/