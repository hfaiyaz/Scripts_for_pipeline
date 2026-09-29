#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --time=04:00:00
#SBATCH --mem=12G
#SBATCH --ntasks=4
#SBATCH --account=amc-general
#SBATCH --job-name=cdd_hmm_maker
#SBATCH --output=./data/cdd_hmm_maker_%J.out
#SBATCH --array=1-200
#SBATCH --mail-user=holly.vose@cuanschutz.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load hmmer
module load perl

perl cdd_hmm_maker.pl cdd_fasta/all_fasta_list.txt $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT

