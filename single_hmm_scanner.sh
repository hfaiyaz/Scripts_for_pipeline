#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --time=03:00:00
#SBATCH --mem=2G
#SBATCH --ntasks=4
#SBATCH --account=amc-general
#SBATCH --job-name=single_hmm_scanner_C_auris_RNA
#SBATCH --output=/scratch/alpine/fhasan1@xsede.org/data/single_hmm_scanner_C_auris_RNA_%J.out
#SBATCH --mail-user=faiyaz.hasan@cuanschutz.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

module load hmmer

python3 single_hmm_scanner.py /projects/fhasan1@xsede.org/viral_hmms/ALL_VIRAL_HMMs.hmm /pl/active/Viralogue/final_contigs/C_auris_RNA/C_auris_negative_contigs.fna $SLURM_NTASKS /scratch/alpine/fhasan1@xsede.org/C_auris_RNA_hmm.txt