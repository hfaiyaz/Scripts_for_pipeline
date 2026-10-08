#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --time=03:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --account=amc-general
#SBATCH --job-name=getorf_C_albicans
#SBATCH --output=/scratch/alpine/fhasan1@xsede.org/data/getorf_C_albicans_%J.out
#SBATCH --mail-user=faiyaz.hasan@cuanschutz.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load emboss

python3 getorf.py /pl/active/Viralogue/final_contigs/C_albicans/C_albicans_negative_contigs.fna /pl/active/Viralogue/C_albicans_trasnlated_333aa_contigs.fa