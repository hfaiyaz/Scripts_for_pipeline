#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --time=2:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --account=amc-general
#SBATCH --job-name=bloomfiltertest
#SBATCH --output=/scratch/alpine/fhasan1@xsede.org/data/bloomfiltertest_%J.out
#SBATCH --mail-user=faiyaz.hasan@cuanschutz.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load biobloom

python3 biobloom1.py /scratch/alpine/fhasan1@xsede.org/C_auris_B11244.fna $SLURM_NTASKS /projects/fhasan1@xsede.org/bloom_references/C_auris_B8441.bf