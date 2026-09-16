#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node 24
#SBATCH -t 4:00:00
#SBATCH --mail-user=haltomj@chop.edu
#SBATCH --mail-type=ALL


module load R/4.6.1
Rscript scripts/DGE.r
