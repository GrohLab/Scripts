#!/bin/bash
#SBATCH --partition=cpu-single
#SBATCH --ntasks=32
#SBATCH --time=24:00:00
#SBATCH --mem=128gb

module load math/matlab/R2023a

matlab -nodisplay -r poolEphBeh_regression > results.out 2>&1