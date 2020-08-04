#!/bin/bash

#SBATCH --partition=shas
#SBATCH --qos=normal
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --job-name=ENDs
#SBATCH --output=OUT-%j.out
#SBATCH --mail-type=all

module load python/2.7.11

python /projects/alge9397/code/python/rebound_runs_alro/end_aleksey_config_b.py --config OUT/config
