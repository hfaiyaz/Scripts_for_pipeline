#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --time=06:00:00
#SBATCH --mem=6G
#SBATCH --ntasks=4
#SBATCH --account=amc-general
#SBATCH --job-name=single_hmm_scanner_C_albicans
#SBATCH --output=/scratch/alpine/fhasan1@xsede.org/data/single_hmm_scanner_C_albicans_%J.out
#SBATCH --mail-user=faiyaz.hasan@cuanschutz.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load hmmer

python3 single_hmm_scanner.py /projects/fhasan1@xsede.org/viral_hmms/ALL_VIRAL_HMMs.hmm /pl/active/Viralogue/C_albicans_translated_333aa_contigs.fa $SLURM_NTASKS /scratch/alpine/fhasan1@xsede.org/C_albicans_hmm.txt