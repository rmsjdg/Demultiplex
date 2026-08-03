#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --job-name=qscores-R4


/usr/bin/time -v pixi run ./qscores.py -a /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -l 101 -r R4