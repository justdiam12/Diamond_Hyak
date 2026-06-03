#!/bin/bash

#SBATCH --job-name=rampe_21
#SBATCH --account=stf
#SBATCH --partition=cpu-g2

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64

#SBATCH --mem=64G
#SBATCH --time=24:00:00

#SBATCH --chdir=/mmfs1/gscratch/stf/jdiam12/Diamond_Hyak
#SBATCH --output=run_output/rampe_21.txt
#SBATCH --error=run_output/rampe_21_error.txt

module load conda
module load gcc/12.3.0

source $(conda info --base)/etc/profile.d/conda.sh
conda activate acoustics

python runs/rampe_21.py