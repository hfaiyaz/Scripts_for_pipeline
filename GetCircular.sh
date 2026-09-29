#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --time=2:00:00
#SBATCH --mem=24G
#SBATCH --ntasks=2
#SBATCH --account=amc-general
#SBATCH --job-name=GetCicrular_C_albicans_negative
#SBATCH --output=./data/GetCicrular_C_albicans_negative_%J.out
#SBATCH --mail-user=faiyaz.hasan@cuanschutz.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load perl

perl GetCircular.pl /pl/active/Viralogue/final_contigs/C_albicans/C_albicans_negative_contigs.fna

