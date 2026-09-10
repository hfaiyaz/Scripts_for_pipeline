#!/usr/bin/env python3

from pathlib import Path
import argparse
import shutil
import subprocess
import sys
import time
import re


def run_getorf(search_root: Path, overwrite: bool = False) -> None:
    getorf_path = shutil.which("getorf")

    if getorf_path is None:
        sys.exit(
            "Error: getorf was not found. Activate the EMBOSS environment "
            "before running this script."
        )

    out_directories = sorted(
        path
        for path in search_root.rglob("*")
        if path.is_dir() and "_out" in path.name
    )

    if not out_directories:
        sys.exit(f"No directories containing '_out' found under {search_root}")

    processed = 0
    skipped = 0
    failed = 0

    for out_directory in out_directories:
        input_file = out_directory / "final.contigs.fa"
        output_file = out_directory / "output.fa"

        if not input_file.is_file():
            print(f"Skipping: {input_file} does not exist")
            skipped += 1
            continue

        if output_file.exists() and not overwrite:
            print(f"Skipping: {output_file} already exists")
            skipped += 1
            continue

        command = [
            getorf_path,
            "-sequence", str(input_file),
            "-outseq", str(output_file),
            "-auto",
        ]

        print(f"Processing: {input_file}")

        try:
            subprocess.run(
                command,
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
            )
            processed += 1

        except subprocess.CalledProcessError as error:
            print(f"getorf failed for {input_file}", file=sys.stderr)
            print(error.stderr, file=sys.stderr)
            failed += 1

    print(
        f"\nFinished: {processed} processed, "
        f"{skipped} skipped, {failed} failed."
    )

    if failed:
        sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Find directories containing '_out' and run EMBOSS getorf "
            "on each final.contigs.fa file."
        )
    )
    parser.add_argument(
        "search_root",
        nargs="?",
        default=".",
        type=Path,
        help="Directory to search recursively (default: current directory)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing output.fa files",
    )

    args = parser.parse_args()
    search_root = args.search_root.resolve()

    if not search_root.is_dir():
        parser.error(f"Search root is not a directory: {search_root}")

    run_getorf(search_root, args.overwrite)


if __name__ == "__main__":
    main()