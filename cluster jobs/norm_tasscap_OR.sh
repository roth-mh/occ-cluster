#!/bin/csh
#SBATCH -J norm_and_reproj                       # name of job

#SBATCH -A eecs                        # name of my sponsored account, e.g. class or research group

#SBATCH -p gpu                                # name of partition or queue

#SBATCH -o norm_and_reproj.out                    # name of output file for this submission script

#SBATCH -e err_norm_and_reproj.err                    # name of error file for this submission script

#SBATCH --mail-type=BEGIN,END,FAIL                # send email when job begins, ends or aborts

#SBATCH --mail-user=rothmark@oregonstate.edu        # send email to this address

# load any software environment module required for app
module load gdal/2.4.4

# run my job
Rscript norm_and_reproj.R
