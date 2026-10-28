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
            