#!/bin/bash
#
#SBATCH -N 1
#SBATCH -t 00-15:00 # Runtime in D-HH:MM
#SBATCH -n 20 # number of cores used 
#SBATCH --mem=15G
#SBATCH -o reunification_loo_sim.out # File to which STDOUT will be written
#SBATCH -e reunification_loo_sim.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=tsudijon@stanford.edu # Email to which notifications will be sent
#SBATCH -J reunification_simulation_loo # name of the job

ml R/3.5.1
ml gcc
# size.resampled.dataset, alpha, tau, mc_samples

Rscript reunification_loo_poweranalysis_slurm.R 14 0.1 0 20

Rscript reunification_loo_poweranalysis_slurm.R 14 0.1 -3600 20

Rscript reunification_loo_poweranalysis_slurm.R 14 0.1 -7200 20

Rscript reunification_loo_poweranalysis_slurm.R 14 0.1 -10800 20