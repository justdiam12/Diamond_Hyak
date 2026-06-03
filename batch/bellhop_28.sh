#!/bin/bash

#SBATCH --job-name=bellhop_28
#SBATCH --account=stf
#SBATCH --partition=cpu-g2

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64

#SBATCH --mem=64G
#SBATCH --time=24:00:00

#SBATCH --chdir=/mmfs1/gscratch/stf/jdiam12/Diamond_Hyak
#SBATCH --output=run_output/bellhop_28.txt
#SBATCH --error=run_output/bellhop_28_error.txt

module load conda
module load gcc/12.3.0

source $(conda info --base)/etc/profile.d/conda.sh
conda activate acoustics

python runs/bellhop_28.py