 #!/usr/bin/env python3

import os
import sys
import subprocess
import time
import re
import shutil
from pathlib import Path
import glob

Dir = '/scratch/alpine/fhasan1@xsede.org/C_auris'


BLASTNucleotide = sys.argv[1]
Catch = int(sys.argv[2]) - 1
NumJobs = int(sys.argv[3])

with open(BLASTNucleotide, 'r') as BLASTN:
    j = 0

    for line in BLASTN:
    
        bloom_too_big = False

        if j % NumJobs == Catch:
            print(f"Numer={j}. Numjobs={NumJobs}, Catch={Catch}")
            print(f"{j%NumJobs}")

            line = line.rstrip('\n')
            Array = line.split('\t')

            accession = Array[0]

            print(f"{accession} ", end='')

            # ----------------------------------------------------------
            # Check whether this accession has already been completed
            # ----------------------------------------------------------

            megahit_out = Path(Dir) / accession
            final_contigs = megahit_out / "final.contigs.fa"

            if final_contigs.is_file():
                print("final.contigs.fa already exists")
                rm_pattern_1 = glob.glob("/scratch/alpine/fhasan1@xsede.org/C_auris/" + accession + "*.*")
                my_og_delete_cmd = ["rm"] + rm_pattern_1
                try:
                    subprocess.run(my_og_delete_cmd, check = True)
                except subprocess.CalledProcessError as bloom_e:
                    print(bloom_e, "Warning: Did not find files to delete", accession)
                
                j += 1
                continue
            
            else:
                rm_pattern_1 = glob.glob("/scratch/alpine/fhasan1@xsede.org/C_auris/" + accession + "*")
                my_og_delete_cmd = ["rm", "-r"] + rm_pattern_1
                try:
                    subprocess.run(my_og_delete_cmd, check = True)
                except subprocess.CalledProcessError as bloom_e:
                    print(bloom_e, "Warning: Did not find files and/or directories to delete", accession)
                


            # ----------------------------------------------------------
            # FASTQ locations
            # ----------------------------------------------------------

            Output = Dir

            pe_item = Output + '/' + accession + '_1.fastq'
            se_item = Output + '/' + accession + '.fastq'


            # ----------------------------------------------------------
            # Download from SRA if FASTQ is not already present
            # ----------------------------------------------------------

            if os.path.exists(pe_item) or os.path.exists(se_item):

                print("FASTQ already exists")

            else:

                try:

                    subprocess.run(
                        [
                            'fasterq-dump',
                            accession,
                            '--split-files', '-e', sys.argv[4],
                            '-O',
                            Output
                        ],
                        check=True
                    )

                except subprocess.CalledProcessError as e:

                    print(e, "Failed at SRADump", accession)
                    j += 1
                    continue


            print("SRADump done")


            # ==========================================================
            # Bloom filtering
            # ==========================================================

            start_time_bf = time.time()


            # ----------------------------------------------------------
            # Paired-end
            # ----------------------------------------------------------

            if Path(pe_item).is_file():

                item_out = pe_item + "_out"

                item2 = re.sub("_1.fastq$", "_2.fastq", pe_item)

                my_bf_cmd = [
                    "biobloomcategorizer",
                    "-d",
                    "-n",
                    "-t",
                    sys.argv[4],
                    "-e",
                    "-p",
                    item_out,
                    "-f",
                    sys.argv[5],
                    pe_item,
                    item2
                ]


                # IMPORTANT:
                # Explicitly define both filtered paired-end files.
                #
                # These are the exact filenames created by awk.

                item_bloom = item_out + "noMatch_1.fq"
                item_bloom_2 = item_out + "noMatch_2.fq"


                awk_cmd = [
                    "awk",
                    "-v",
                    f"out1={item_bloom}",
                    "-v",
                    f"out2={item_bloom_2}",
                    '{print > (int((NR-1)/4)%2==0 ? out1 : out2)}'
                ]


                try:

                    with subprocess.Popen(
                        my_bf_cmd,
                        stdout=subprocess.PIPE
                    ) as bf_out:

                        subprocess.run(
                            awk_cmd,
                            stdin=bf_out.stdout,
                            check=True
                        )

                        bf_out.stdout.close()

                        return_code = bf_out.wait()


                        if return_code != 0:

                            print(
                                f"{accession} failed with exit code "
                                f"{return_code}"
                            )

                            j += 1
                            continue


                except subprocess.CalledProcessError as e:

                    print(
                        e,
                        "Failed at paired-end Bloomfiltering",
                        accession
                    )

                    j += 1
                    continue


            # ----------------------------------------------------------
            # Single-end
            # ----------------------------------------------------------

            else:

                item_out = se_item + "_out"

                my_bf_cmd = [
                    "biobloomcategorizer",
                    "-d",
                    "-n",
                    "-t",
                    sys.argv[4],
                    "-p",
                    item_out,
                    "-f",
                    sys.argv[5],
                    se_item
                ]


                item_bloom = item_out + "noMatch.fq"


                with open(item_bloom, "w") as f:

                    try:

                        subprocess.run(
                            my_bf_cmd,
                            check=True,
                            stdout=f
                        )

                    except subprocess.CalledProcessError as e:

                        print(
                            e,
                            "Failed at se Bloomfiltering",
                            accession
                        )

                        j += 1
                        continue


            end_time_bf = time.time()

            runtime_bf = end_time_bf - start_time_bf

            print(runtime_bf, "Bloomfilter done")
            
            # ==========================================================
            # Remove bloom filter output that is larger than a 20% hit
            # ==========================================================
            
            if Path(pe_item).is_file():
                bloom_summary = "/scratch/alpine/fhasan1@xsede.org/C_auris/" + accession + "_1.fastq_out_summary.tsv"
                with open(bloom_summary, "r") as bloomf:
                    for line  in bloomf:
                        if "noMatch" in line:
                            line = line.strip().split("\t")
                            print(accession, "noMatch rate:", line[4])
                            if float(line[4]) > 0.2 or float(line[4]) == 0:
                                bloom_too_big = True
                                rm_pattern_2 = glob.glob("/scratch/alpine/fhasan1@xsede.org/C_auris/" + accession + "*.*")
                                print(rm_pattern_2)
                                my_bloom_delete_cmd =  ["rm"] + rm_pattern_2
                                my_sra_delete_cmd = ["rm",pe_item, item2]
                                try:
                                    subprocess.run(my_sra_delete_cmd, check = True)
                                    subprocess.run(my_bloom_delete_cmd, check = True)
                                    print("Bloom files larger than 20% matches and SRA deleted")
                                except subprocess.CalledProcessError as bloom_e:
                                    print(bloom_e, "Warning: Bloom files largesr than 20% matches could not be deleted", accession)
                                    
                                    
                                    
                                
            
            else:
                bloom_summary = "/scratch/alpine/fhasan1@xsede.org/C_auris/" + accession + ".fastq_out_summary.tsv"
                with open(bloom_summary, "r") as bloomf:
                    for line  in bloomf:
                        if "noMatch" in line:
                            line = line.strip().split("\t")
                            print(accession, "noMatch rate:", line[4])
                            if float(line[4]) > 0.2 or float(line[4]) == 0:
                                bloom_too_big = True
                                rm_pattern_2 = glob.glob("/scratch/alpine/fhasan1@xsede.org/C_auris/" + accession + "*.*")
                                print(rm_pattern_2)
                                my_bloom_delete_cmd = ["rm"] + rm_pattern_2
                                my_sra_delete_cmd = ["rm", se_item]
                                try:
                                    subprocess.run(my_sra_delete_cmd, check = True)
                                    subprocess.run(my_bloom_delete_cmd, check = True)
                                    print("Bloom files larger than 20% matches and SRA deleted")
                                except subprocess.CalledProcessError as bloom_e:
                                    print(bloom_e, "Warning: Bloom files largesr than 20% matches could not be deleted", accession)
            if bloom_too_big:
                j += 1
                continue
                                
                            


            # ==========================================================
            # Delete original SRA FASTQ files
            # ==========================================================

            if Path(pe_item).is_file():

                my_sra_delete_cmd = [
                    "rm",
                    pe_item,
                    item2
                ]

            else:

                my_sra_delete_cmd = [
                    "rm",
                    se_item
                ]


            try:

                subprocess.run(
                    my_sra_delete_cmd,
                    check=True
                )

            except subprocess.CalledProcessError as e:

                print(
                    e,
                    "Warning: could not delete original FASTQ",
                    accession
                )


            print("SRA deleted", j)
            
            
            # ==========================================================
            # Remove incomplete MEGAHIT directory from previous attempt
            # ==========================================================

            if megahit_out.exists() and not final_contigs.is_file():

                print(
                    "Removing incomplete MEGAHIT directory:",
                    megahit_out
                )

                shutil.rmtree(megahit_out)


            # ==========================================================
            # MEGAHIT
            # ==========================================================

            if Path(pe_item).exists():

                my_mh_cmd = [
                    "megahit",
                    "-1",
                    item_bloom,
                    "-2",
                    item_bloom_2,
                    "-o",
                    str(megahit_out)
                ]

            else:

                my_mh_cmd = [
                    "megahit",
                    "-r",
                    item_bloom,
                    "-o",
                    str(megahit_out)
                ]


            start_time_mh = time.time()


            try:

                result_mh = subprocess.run(
                    my_mh_cmd,
                    capture_output=True,
                    text=True,
                    check=True
                )


            except subprocess.CalledProcessError as e:

                print(
                    "Failed at Megahit",
                    accession
                )

                print("Command:")
                print(" ".join(my_mh_cmd))

                print("MEGAHIT STDOUT:")
                print(e.stdout)

                print("MEGAHIT STDERR:")
                print(e.stderr)

                j += 1
                continue


            end_time_mh = time.time()

            runtime_mh = end_time_mh - start_time_mh


            print(
                runtime_mh,
                result_mh.stdout,
                result_mh.stderr,
                "Megahit done"
            )


            # ==========================================================
            # Delete MEGAHIT intermediate files
            # Keep only final.contigs.fa
            # ==========================================================

            my_mh_delete_cmd = [
                "find",
                str(megahit_out) + "/.",
                "-mindepth",
                "1",
                "-not",
                "-name",
                "final.contigs.fa",
                "-delete"
            ]


            try:

                subprocess.run(
                    my_mh_delete_cmd,
                    check=True
                )

            except subprocess.CalledProcessError as e:

                print(
                    e,
                    "Warning: MEGAHIT cleanup failed",
                    accession
                )
                continue


             # ==========================================================
            # Delete SRA cache
            # ==========================================================
            cache1 = Path("/scratch/alpine/fhasan1@xsede.org/sra/sra/" + accession + ".sra")
            cache2 = Path("/scratch/alpine/fhasan1@xsede.org/sra/sra/" + accession + ".sra.cache")
            
            if not cache1.is_file() and not cache2.is_file():
                pass
            else:
                if cache1.is_file():
                    my_sra_cache_delete_cmd = ["rm", cache1]
                elif cache2.is_file():
                    my_sra_cache_delete_cmd = ["rm", cache2]
                else:
                    print("SRA cache not found")
                try:
                    subprocess.run(my_sra_cache_delete_cmd, check=True)
                    print("SRA cache deleted")
                except subprocess.CalledProcessError as e:
                    print(e,"Warning: SRA cache delete failed", accession)
                    continue

                    
            
            
                        # ==========================================================
            # BLASTN final.contigs.fa and create negative contigs FASTA
            # ==========================================================

            blast_db = sys.argv[6]

            negative_output_dir = Path(sys.argv[7])
            negative_output_dir.mkdir(parents=True, exist_ok=True)

            blast_out = megahit_out / "final.contigs_blastn.csv"

            negative_contigs = (
                negative_output_dir /
                f"{accession}_negative_final_contigs.fa"
            )


            # ----------------------------------------------------------
            # Run BLASTN
            # ----------------------------------------------------------

            print("Running BLASTN:", accession)

            try:

                with open(blast_out, "w") as blastf:

                    blastf.write(
                        "qseqid,sseqid,pident,length,qlen,qcovs,"
                        "mismatch,gapopen,qstart,qend,sstart,send,"
                        "evalue,bitscore\n"
                    )

                    my_blast_cmd = [
                        os.path.expanduser(
                            "~/software/ncbi-blast-2.17.0+/bin/blastn"
                        ),
                        "-query",
                        str(final_contigs),
                        "-db",
                        blast_db,
                        "-outfmt",
                        "10 qseqid sseqid pident length qlen qcovs "
                        "mismatch gapopen qstart qend sstart send "
                        "evalue bitscore",
                        "-num_threads",
                        sys.argv[4]
                    ]

                    subprocess.run(
                        my_blast_cmd,
                        stdout=blastf,
                        check=True
                    )

            except subprocess.CalledProcessError as e:

                print(
                    e,
                    "Failed at BLASTN",
                    accession
                )

                j += 1
                continue


            print("BLASTN done:", accession)


            # ==========================================================
            # Read BLAST results
            # ==========================================================

            hit_ids = set()
            low80_ids = set()
            high80_ids = set()


            with open(blast_out, "r") as blastf:

                next(blastf)

                for blast_line in blastf:

                    blast_line = blast_line.rstrip("\n")

                    if not blast_line:
                        continue

                    blast_array = blast_line.split(",")

                    qseqid = blast_array[0]
                    pident = float(blast_array[2])
                    qcovs = float(blast_array[5])


                    # Contig had at least one BLAST alignment
                    hit_ids.add(qseqid)


                    # Contig had a hit with:
                    # pident < 80 AND qcovs < 80
                    if pident < 80 and qcovs < 80:

                        low80_ids.add(qseqid)


                    # Contig had a hit with:
                    # pident > 80 AND qcovs > 80
                    if pident > 80 and qcovs > 80:

                        high80_ids.add(qseqid)


            # ==========================================================
            # Get all contig IDs from final.contigs.fa
            # ==========================================================

            all_contig_ids = set()


            with open(final_contigs, "r") as contigf:

                for contig_line in contigf:

                    if contig_line.startswith(">"):

                        contig_id = contig_line[1:].split()[0]

                        all_contig_ids.add(contig_id)


            # ==========================================================
            # Define negative contigs
            #
            # Negative contigs are:
            #
            # 1. Contigs with NO BLAST alignment
            #
            # OR
            #
            # 2. Contigs with at least one <80/<80 hit
            #    AND no >80/>80 hit
            # ==========================================================

            no_hit_ids = all_contig_ids - hit_ids

            low80_no_high80_ids = low80_ids - high80_ids

            negative_ids = no_hit_ids | low80_no_high80_ids


            print(
                accession,
                "total contigs:",
                len(all_contig_ids)
            )

            print(
                accession,
                "contigs with BLAST hits:",
                len(hit_ids)
            )

            print(
                accession,
                "contigs with no BLAST hits:",
                len(no_hit_ids)
            )

            print(
                accession,
                "contigs with <80/<80 and no >80/>80 hit:",
                len(low80_no_high80_ids)
            )

            print(
                accession,
                "total negative contigs:",
                len(negative_ids)
            )


            # ==========================================================
            # Write negative contigs FASTA to sys.argv[7]
            # ==========================================================

            negative_written = 0
            keep_contig = False


            with open(final_contigs, "r") as infile, \
                 open(negative_contigs, "w") as outfile:

                for contig_line in infile:

                    if contig_line.startswith(">"):

                        contig_id = contig_line[1:].split()[0]

                        keep_contig = contig_id in negative_ids

                        if keep_contig:

                            negative_written += 1


                    if keep_contig:

                        outfile.write(contig_line)


            # ==========================================================
            # Verify negative FASTA before deleting anything
            # ==========================================================

            if negative_written != len(negative_ids):

                print(
                    "ERROR:",
                    accession,
                    "negative FASTA count does not match "
                    "negative ID count"
                )

                print(
                    "Expected:",
                    len(negative_ids),
                    "Written:",
                    negative_written
                )

                j += 1
                continue


            print(
                accession,
                "negative contigs written:",
                negative_written
            )

            print(
                accession,
                "negative FASTA:",
                negative_contigs
            )


            # ==========================================================
            # Final cleanup
            #
            # The negative FASTA has been successfully written and
            # verified in sys.argv[7].
            #
            # Delete the entire MEGAHIT accession directory.
            # This removes:
            #
            #   final.contigs.fa
            #   final.contigs_blastn.csv
            #   megahit_out itself
            #
            # No intermediate ID files are created.
            # ==========================================================

            try:

                shutil.rmtree(megahit_out)

                print(
                    accession,
                    "MEGAHIT directory deleted:",
                    megahit_out
                )

            except OSError as e:

                print(
                    e,
                    "Warning: Could not delete MEGAHIT directory",
                    accession
                )

                j += 1
                continue


            print(
                accession,
                "BLAST filtering and final cleanup done"
            )
            print("Pipeline part 1 done")


        j += 1