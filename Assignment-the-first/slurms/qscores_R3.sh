#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8
#SBATCH --job-name=qscores-R3


/usr/bin/time -v pixi run ./qscores.py -a /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -l 8 -r R3