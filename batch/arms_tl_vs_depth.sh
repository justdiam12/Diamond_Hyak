#!/bin/bash

#SBATCH --job-name=arms_tl_vs_depth
#SBATCH --account=stf
#SBATCH --partition=cpu-g2

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SBATCH --mem=64G
#SBATCH --time=24:00:00

#SBATCH --chdir=/mmfs1/gscratch/stf/jdiam12/Diamond_Hyak
#SBATCH --output=run_output/arms_tl_vs_depth.txt
#SBATCH --error=run_output/arms_tl_vs_depth_error.txt

module load conda
module load matlab/r2023a

source $(conda info --base)/etc/profile.d/conda.sh
conda activate acoustics

matlab -batch "cd('/mmfs1/gscratch/stf/jdiam12/Diamond_Hyak'); run('runs/arms_tl_vs_depth.m')"