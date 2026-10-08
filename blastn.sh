#!/bin/sh

#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=amc-general
#SBATCH --job-name=blastn_C_auris
#SBATCH --output=./data/blastn_C_auris_%J.out
#SBATCH --mail-user=faiyaz.hasan@cuanschutz.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END

set -euo pipefail

QUERY="/scratch/alpine/fhasan1@xsede.org/C_auris_RNA_negative_CircularContigs.fa"
DB="/pl/active/Viralogue/blastdb/circulardb/C_auris_RNA/C_auris_RNA_circular"
OUT="C_auris_RNA_circular_self_blastn.csv"

# Write CSV header
echo "qseqid,sseqid,pident,length,qlen,qcovs,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore" > "$OUT"

# Run BLASTN
~/software/ncbi-blast-2.17.0+/bin/blastn \
    -query "$QUERY" \
    -db "$DB" \
    -outfmt '10 qseqid sseqid pident length qlen qcovs mismatch gapopen qstart qend sstart send evalue bitscore' \
    -num_threads "$SLURM_CPUS_PER_TASK" \
    >> "$OUT"

echo "BLAST finished."
echo "Results written to: $OUT"