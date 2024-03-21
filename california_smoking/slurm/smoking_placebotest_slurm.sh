#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 00-15:00 # Runtime in D-HH:MM
#SBATCH -n 20 # number of cores used 
#SBATCH --mem=15G
#SBATCH -o smoking_placebo_poweranalysis.out # File to which STDOUT will be written
#SBATCH -e smoking_placebo_poweranalysis.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tsudijon@stanford.edu # Email to which notifications will be sent
#SBATCH -J smoking_placebotest_slurm # name of the job

ml R/3.5.1
ml gcc
# size.resampled.dataset, alpha, tau, mc_samples

Rscript smoking_placebotest_poweranalysis_slurm.R 30 0.02 0 20

Rscript smoking_placebotest_poweranalysis_slurm.R 30 0.02 -30 20

Rscript smoking_placebotest_poweranalysis_slurm.R 30 0.02 -60 20

Rscript smoking_placebotest_poweranalysis_slurm.R 30 0.02 -90 20