#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --time=04:00:00
#SBATCH --mem=48G
#SBATCH --ntasks=8
#SBATCH --account=amc-general
#SBATCH --job-name=pipeline_1
#SBATCH --output=./data/pipeline_1_%J.out
#SBATCH --array=1-4
#SBATCH --mail-user=faiyaz.hasan@cuanschutz.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load sra-toolkit
module load biobloom
module load megahit

python3 Pipeline_1.py C_albicans_SRR_test.txt $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT $SLURM_NTASKS /projects/fhasan1@xsede.org/bloom_references/C_albicans_SC5314_fkh.bf