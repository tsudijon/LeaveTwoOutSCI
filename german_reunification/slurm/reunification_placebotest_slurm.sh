#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 00-15:00 # Runtime in D-HH:MM
#SBATCH -n 20 # number of cores used 
#SBATCH --mem=15G
#SBATCH -o reunification_placebo_sim.out # File to which STDOUT will be written
#SBATCH -e reunification_placebo_sim.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tsudijon@stanford.edu # Email to which notifications will be sent
#SBATCH -J reunification_simulation_placebotest # name of the job

ml R
ml gcc
# size.resampled.dataset, alpha, tau, mc_samples

Rscript reunification_placebotest_poweranalysis_slurm.R 14 0.05 0 20

Rscript reunification_placebotest_poweranalysis_slurm.R 14 0.05 -1800 20

Rscript reunification_placebotest_poweranalysis_slurm.R 14 0.05 -3600 20

Rscript reunification_placebotest_poweranalysis_slurm.R 14 0.05 -7200 20

Rscript reunification_placebotest_poweranalysis_slurm.R 14 0.05 -10800 20