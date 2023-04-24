#!/bin/bash

#SBATCH -p ipuq # partition (queue)
#SBATCH -t 0-24:00 # time limit (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

sh scripts/run_demo.sh 24192
