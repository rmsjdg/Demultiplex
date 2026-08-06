#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=rj.dx30

/usr/bin/time -v /projects/bgmp/roj/bioinfo/Bi622/Demultiplex/Assignment-the-third/demultiplex.py \
 -a /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
 -b /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
 -c /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
 -d /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
 -i /projects/bgmp/shared/2017_sequencing/indexes.txt \
 -o /scratch/bgmp/roj/demux/q30/ \
 -q 30