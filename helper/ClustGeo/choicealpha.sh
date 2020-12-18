#!/bin/csh

# use current working directory for input and output
# default is to use the users home directory
#$ -cwd

# name this job
#$ -N choicealpha_cluster

# send stdout and stderror to this file
#$ -o choicealpha_cluster.out
#$ -j y

# select queue - if needed 
#$ -q eecs3,share,share2,share3,share4 

# set mail notificatoin - set who to send mail to
#$ -M rothmark@engr.orst.edu

# next, set when to send email, currently beginning and end of job
#$ -m be

# see where the job is being run
hostname

# print date and time
date

# Sleep for 20 seconds
sleep 20

# run job
#Rscript sample_smote.R
Rscript choicealpha_cluster.R

# print date and time again
date
