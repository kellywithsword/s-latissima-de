#!/bin/bash
#SBATCH -p cegs
#SBATCH --mem=1g
#SBATCH --time=1-0

fruitcake=$2

source activate trimgalore
conda activate trimgalore

