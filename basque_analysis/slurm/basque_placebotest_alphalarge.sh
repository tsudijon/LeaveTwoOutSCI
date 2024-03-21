#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 00-15:00 # Runtime in D-HH:MM
#SBATCH -n 20 # number of cores used 
#SBATCH --mem=15G
#SBATCH -o basque_placebo_poweranalysis.out # File to which STDOUT will be written
#SBATCH -e basque_placebo_poweranalysis.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tsudijon@stanford.edu # Email to which notifications will be sent
#SBATCH -J basque_simulation_placebo_tau0 # name of the job

ml R/3.5.1
ml gcc
# size.resampled.dataset, alpha, tau, mc_samples

Rscript basque_placebotest_poweranalysis_slurm.R 15 0.1 0 20

Rscript basque_placebotest_poweranalysis_slurm.R 15 0.1 -1.5 20

Rscript basque_placebotest_poweranalysis_slurm.R 15 0.1 -3.0 20

Rscript basque_placebotest_poweranalysis_slurm.R 15 0.1 -4.5 20