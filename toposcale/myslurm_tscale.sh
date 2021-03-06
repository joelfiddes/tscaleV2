#!/bin/sh
#SBATCH -J tscale
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joelfiddes@gmail.com
#SBATCH --workdir=/home/caduff/src/tscaleV2/toposcale
#SBATCH --ntasks=39       # tasks requested
#SBATCH --mem-per-cpu=2000
#SBATCH -o outfile  # send stdout to outfile
#SBATCH -e errfile  # send stderr to errfile
#SBATCH -t 2:00:00  # time requested in hour:minute:second

for i in {1..10}
do
echo "computing member" $i
python slurm_tscale.py 39 /home/caduff/sim/tscale 1979-09-01 2018-09-01 EDA $i
done
#Args
# 39 - number of cores/jobs
# wdir
# start date
# end date
# optional ensemble memeber number


# run as: sbatch -J "tscale"  myslurm_tscale.sh

