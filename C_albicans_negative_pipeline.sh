#!/bin/bash

#SBATCH --nodes=1
#SBATCH --qos=cpu-normal
#SBATCH --partition=acpu
#SBATCH --time=04:00:00
#SBATCH --mem=128G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=amc-general
#SBATCH --job-name=C_albicans_negative
#SBATCH --output=C_albicans_negative_%j.out
#SBATCH --mail-user=faiyaz.hasan@cuanschutz.edu
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# -----------------------------
# Paths
# -----------------------------

INPUT_DIR="/pl/active/Viralogue/final_contigs/C_albicans"

DB="/pl/active/Viralogue/blastdb/candida_albicans/Candida"

BLASTN="$HOME/software/ncbi-blast-2.17.0+/bin/blastn"
BLASTDBCMD="$HOME/software/ncbi-blast-2.17.0+/bin/blastdbcmd"

WORKDIR="/scratch/alpine/fhasan1@xsede.org/C_albicans_negative_pipeline"

COMBINED="$WORKDIR/C_albicans_all_contigs_unique.fna"
BLAST_OUT="$WORKDIR/C_albicans_blastn_results_unique.csv"
NEGATIVE_FASTA="$WORKDIR/C_albicans_negative_contigs.fna"

mkdir -p "$WORKDIR"

# Temporary working files
TMPDIR=$(mktemp -d "$WORKDIR/tmp.XXXXXX")
trap 'rm -rf "$TMPDIR"' EXIT

ALL_IDS="$TMPDIR/all_ids.txt"
HIT_IDS="$TMPDIR/hit_ids.txt"
HIGH_IDS="$TMPDIR/high80_ids.txt"
LOW_IDS="$TMPDIR/low80_ids.txt"
LOW_NO_HIGH_IDS="$TMPDIR/low80_no_high80_ids.txt"
NO_HIT_IDS="$TMPDIR/no_hit_ids.txt"
NEGATIVE_IDS="$TMPDIR/negative_ids.txt"

echo "========================================"
echo "Starting C. albicans negative-contig pipeline"
echo "========================================"


# ============================================================
# 1. COMBINE FASTA FILES AND ADD ACCESSION TO EVERY HEADER
# ============================================================

echo
echo "Combining FASTA files..."

: > "$COMBINED"

for file in "$INPUT_DIR"/*_final.contigs.fa; do

    [ -e "$file" ] || {
        echo "ERROR: No *_final.contigs.fa files found in $INPUT_DIR"
        exit 1
    }

    acc=$(basename "$file" _final.contigs.fa)

    awk -v acc="$acc" '
        /^>/ {
            sub(/^>/, ">" acc "_")
        }
        {print}
    ' "$file" >> "$COMBINED"

done

echo "Combined FASTA created:"
echo "$COMBINED"


# ============================================================
# 2. VERIFY THAT ALL CONTIG IDS ARE UNIQUE
# ============================================================

echo
echo "Checking contig ID uniqueness..."

TOTAL_CONTIGS=$(grep -c '^>' "$COMBINED")

grep '^>' "$COMBINED" \
    | sed 's/^>//' \
    | awk '{print $1}' \
    | LC_ALL=C sort -u \
    > "$ALL_IDS"

UNIQUE_CONTIGS=$(wc -l < "$ALL_IDS")

echo "Total contigs:  $TOTAL_CONTIGS"
echo "Unique IDs:     $UNIQUE_CONTIGS"

if [ "$TOTAL_CONTIGS" -ne "$UNIQUE_CONTIGS" ]; then

    echo
    echo "ERROR: CONTIG IDS ARE NOT UNIQUE."
    echo "Pipeline stopped before BLAST."
    echo

    exit 1
fi

echo "PASS: All contig IDs are unique."


# ============================================================
# 3. VERIFY BLAST DATABASE
# ============================================================

echo
echo "Checking BLAST database..."

"$BLASTDBCMD" -db "$DB" -info > /dev/null

echo "BLAST database found:"
echo "$DB"


# ============================================================
# 4. RUN BLASTN
# ============================================================

echo
echo "Running BLASTN..."

echo "qseqid,sseqid,pident,length,qlen,qcovs,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore" \
    > "$BLAST_OUT"

"$BLASTN" \
    -query "$COMBINED" \
    -db "$DB" \
    -outfmt '10 qseqid sseqid pident length qlen qcovs mismatch gapopen qstart qend sstart send evalue bitscore' \
    -num_threads "$SLURM_CPUS_PER_TASK" \
    >> "$BLAST_OUT"

echo "BLAST completed."


# ============================================================
# 5. FIND ALL CONTIGS THAT HAD ANY BLAST HIT
# ============================================================

tail -n +2 "$BLAST_OUT" \
    | cut -d',' -f1 \
    | LC_ALL=C sort -u \
    > "$HIT_IDS"

NUM_HITS=$(wc -l < "$HIT_IDS")

echo
echo "Contigs with at least one BLAST hit: $NUM_HITS"


# ============================================================
# 6. FIND CONTIGS WITH NO BLAST HIT
# ============================================================

comm -23 "$ALL_IDS" "$HIT_IDS" > "$NO_HIT_IDS"

NUM_NO_HITS=$(wc -l < "$NO_HIT_IDS")

echo "Contigs with no BLAST hit: $NUM_NO_HITS"


# ============================================================
# 7. FIND CONTIGS WITH >80 IDENTITY AND >80 COVERAGE
# ============================================================

awk -F',' '
    NR > 1 && $3 > 80 && $6 > 80 {
        print $1
    }
' "$BLAST_OUT" \
    | LC_ALL=C sort -u \
    > "$HIGH_IDS"

NUM_HIGH=$(wc -l < "$HIGH_IDS")

echo "Contigs with at least one >80/>80 hit: $NUM_HIGH"


# ============================================================
# 8. FIND CONTIGS WITH <80 IDENTITY AND <80 COVERAGE
# ============================================================

awk -F',' '
    NR > 1 && $3 < 80 && $6 < 80 {
        print $1
    }
' "$BLAST_OUT" \
    | LC_ALL=C sort -u \
    > "$LOW_IDS"

NUM_LOW=$(wc -l < "$LOW_IDS")

echo "Contigs with at least one <80/<80 hit: $NUM_LOW"


# ============================================================
# 9. REMOVE ANY LOW-80 CONTIG THAT ALSO HAS A >80/>80 HIT
# ============================================================

comm -23 "$LOW_IDS" "$HIGH_IDS" > "$LOW_NO_HIGH_IDS"

NUM_LOW_NO_HIGH=$(wc -l < "$LOW_NO_HIGH_IDS")

echo "Contigs with <80/<80 hit AND no >80/>80 hit: $NUM_LOW_NO_HIGH"


# ============================================================
# 10. COMBINE THE TWO NEGATIVE GROUPS
#
#     Negative =
#       no BLAST hit
#       OR
#       <80/<80 hit with no >80/>80 hit
# ============================================================

cat "$NO_HIT_IDS" "$LOW_NO_HIGH_IDS" \
    | LC_ALL=C sort -u \
    > "$NEGATIVE_IDS"

NUM_NEGATIVE=$(wc -l < "$NEGATIVE_IDS")

echo
echo "Total negative contigs: $NUM_NEGATIVE"


# ============================================================
# 11. EXTRACT ACTUAL NEGATIVE CONTIG SEQUENCES
# ============================================================

echo
echo "Extracting negative contig sequences..."

awk '
    NR == FNR {
        ids[$1] = 1
        next
    }

    /^>/ {
        header = substr($0, 2)
        split(header, a, /[[:space:]]+/)
        id = a[1]
        keep = (id in ids)
    }

    keep {
        print
    }
' "$NEGATIVE_IDS" "$COMBINED" > "$NEGATIVE_FASTA"


# ============================================================
# 12. VERIFY FINAL NEGATIVE FASTA
# ============================================================

FINAL_COUNT=$(grep -c '^>' "$NEGATIVE_FASTA")

echo
echo "Negative IDs expected:  $NUM_NEGATIVE"
echo "Negative FASTA records: $FINAL_COUNT"

if [ "$NUM_NEGATIVE" -ne "$FINAL_COUNT" ]; then

    echo
    echo "ERROR: Negative FASTA count does not match negative ID count."
    echo "Pipeline stopped."
    exit 1

fi


# ============================================================
# FINISHED
# ============================================================

echo
echo "========================================"
echo "PIPELINE COMPLETED SUCCESSFULLY"
echo "========================================"
echo
echo "Total input contigs:             $TOTAL_CONTIGS"
echo "Contigs with any BLAST hit:      $NUM_HITS"
echo "Contigs with no BLAST hit:       $NUM_NO_HITS"
echo "Contigs with >80/>80 hit:        $NUM_HIGH"
echo "Contigs with <80/<80 hit:        $NUM_LOW"
echo "<80/<80 with no >80/>80 hit:    $NUM_LOW_NO_HIGH"
echo "Final negative contigs:          $FINAL_COUNT"
echo
echo "Combined FASTA:"
echo "$COMBINED"
echo
echo "BLAST results:"
echo "$BLAST_OUT"
echo
echo "Final negative FASTA:"
echo "$NEGATIVE_FASTA"
echo