#!/bin/bash
#SBATCH --job-name=SparseGRM
#SBATCH --output=logs/SparseGRM.out
#SBATCH --error=logs/SparseGRM.err
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu=2000
#SBATCH --partition=tier3q

module load gcc/11.3.0
module load openblas/0.3.20
module load python/3.10.5
module load plink/2.0
module load R

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile_UKBB
snakemake --cores all --rerun-incomplete -s snakefile_UKBB    
#snakemake -j1 --rerun-incomplete 
