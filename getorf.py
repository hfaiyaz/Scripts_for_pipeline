#!/usr/bin/env python3

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


MIN_AA_LENGTH = 333


def filter_fasta_by_length(input_fasta: Path, output_fasta: Path, min_length: int) -> None:
    kept = 0
    removed = 0

    with input_fasta.open("r") as infile, output_fasta.open("w") as outfile:
        header = None
        sequence_parts = []

        def write_if_long_enough(header, sequence_parts):
            nonlocal kept, removed

            if header is None:
                return

            sequence = "".join(sequence_parts).strip()

            if len(sequence) >= min_length:
                outfile.write(header + "\n")

                # Write sequence in 60-aa lines
                for i in range(0, len(sequence), 60):
                    outfile.write(sequence[i:i + 60] + "\n")

                kept += 1
            else:
                removed += 1

        for line in infile:
            line = line.strip()

            if line.startswith(">"):
                write_if_long_enough(header, sequence_parts)
                header = line
                sequence_parts = []
            else:
                sequence_parts.append(line)

        # Process final sequence
        write_if_long_enough(header, sequence_parts)

    print(f"Filter complete: {kept} ORFs kept, {removed} ORFs removed.")
    print(f"Minimum ORF length: {min_length} amino acids")


def run_getorf(input_file: Path, output_file: Path) -> None:
    getorf_path = shutil.which("getorf")

    if getorf_path is None:
        sys.exit(
            "Error: getorf was not found. Activate the EMBOSS environment "
            "before running this script."
        )

    if not input_file.is_file():
        sys.exit(f"Error: input FASTA file does not exist: {input_file}")

    # Temporary file for unfiltered getorf output
    with tempfile.NamedTemporaryFile(
        suffix=".fa",
        delete=False
    ) as temp:
        temp_output = Path(temp.name)

    command = [
        getorf_path,
        "-sequence", str(input_file),
        "-outseq", str(temp_output),
        "-auto",
    ]

    print(f"Input:  {input_file}")
    print("Running EMBOSS getorf...")

    try:
        subprocess.run(
            command,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )

    except subprocess.CalledProcessError as error:
        temp_output.unlink(missing_ok=True)

        print(f"getorf failed for {input_file}", file=sys.stderr)
        print(error.stderr, file=sys.stderr)
        sys.exit(1)

    print("Translation finished.")
    print(f"Filtering ORFs shorter than {MIN_AA_LENGTH} amino acids...")

    filter_fasta_by_length(
        temp_output,
        output_file,
        MIN_AA_LENGTH,
    )

    # Delete unfiltered temporary FASTA
    temp_output.unlink(missing_ok=True)

    print(f"Final output: {output_file}")


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(
            f"Usage: {sys.argv[0]} <input_fasta> <output_fasta>"
        )

    input_file = Path(sys.argv[1]).resolve()
    output_file = Path(sys.argv[2]).resolve()

    run_getorf(input_file, output_file)


if __name__ == "__main__":
    main()