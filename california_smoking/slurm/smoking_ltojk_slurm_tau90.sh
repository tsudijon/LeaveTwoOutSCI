#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 00-5:00 # Runtime in D-HH:MM
#SBATCH -n 20 # number of cores used 
#SBATCH --mem=15G
#SBATCH -o smokinglto_tau90_sim.out # File to which STDOUT will be written
#SBATCH -e smokinglto_tau90_sim.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tsudijon@stanford.edu # Email to which notifications will be sent
#SBATCH -J smokinglto_tau90_ # name of the job

ml R
ml gcc

# size.resampled.dataset, alpha, tau, mc_samples
Rscript smoking_ltojk_poweranalysis_slurm.R 30 0.05 -90 20
