#!/bin/sh
#SBATCH -J tmapp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joelfiddes@gmail.com
#SBATCH --workdir=/home/caduff/src/tmapp2/
#SBATCH --ntasks=39       # tasks requested
#SBATCH --mem-per-cpu=2000
#SBATCH -o outfile  # send stdout to outfile
#SBATCH -e errfile  # send stderr to errfile
#SBATCH -t 2:00:00  # time requested in hour:minute:second


python slurm.py 39 /home/caduff/sim/tscale "1979-09-01" "2018-09-01"



