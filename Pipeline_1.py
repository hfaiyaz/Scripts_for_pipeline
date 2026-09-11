#!/usr/bin/env python3

import os
import sys
import subprocess
import time
import re
from pathlib import Path

Dir = '/scratch/alpine/fhasan1@xsede.org/sra_python'


BLASTNucleotide = sys.argv[1]
Catch = int(sys.argv[2]) - 1
NumJobs = int(sys.argv[3])

with open(BLASTNucleotide, 'r') as BLASTN:
    j = 0
    for line in BLASTN:
        if j % NumJobs == Catch:
            line = line.rstrip('\n')
            Array = line.split('\t')
            print(f"{Array[0]} ", end='')
            
            if Path(Dir + '/' + Array[0] + '_1.fastq_out_noMatch_1.fq_out/final.contigs.fa').is_file():
                print("folder exists")
                continue
            else:
                pass

            Output = Dir
            pe_item = Output + '/' + Array[0] + '_1.fastq'
            se_item = Output + '/' + Array[0] + '.fastq'

            if (os.path.exists(pe_item)) or (os.path.exists(se_item)):
                j += 1
            else:
                try:
                    subprocess.run(['fasterq-dump', Array[0], '--split-files', '-O', Output])
                except subprocess.CalledProcessError as e:
                        print(e, "Failed at SRADump", Array[0])
                        continue
            print("SRADump done")
            
            ##Bloomfiltering
            start_time_bf = time.time()
            if Path(pe_item).is_file():
                item_out = pe_item + "_out"
                item2 = re.sub("_1", "_2", pe_item)
                my_bf_cmd = ["biobloomcategorizer", "-d", "-n", "-t", sys.argv[4], "-e", "-p", item_out, "-f", sys.argv[5], pe_item, item2]
                item_bloom = item_out + "noMatch_1.fq"
                awk_cmd = ["awk", '{print > (int((NR-1)/4)%2==0 ? item_bloom : item_out + "noMatch_2.fq")}']
                with subprocess.Popen(my_bf_cmd, stdout = subprocess.PIPE) as bf_out:
                    subprocess.run(awk_cmd, stdin = bf_out.stdout, check = True)
                    bf_out.stdout.close()
                    return_code = bf_out.wait()
                    if return_code != 0:
                        print(f"{Array[0]} failed with exit code {return_code}")
                        continue
            else:
                item_out = se_item + "_out"
                my_bf_cmd = ["biobloomcategorizer", "-d", "-n", "-t", sys.argv[4], "-p", item_out, "-f", sys.argv[5], se_item]
                item_bloom = item_out + "noMatch.fq"
                with open(item_bloom, "w") as f:
                    try:
                        subprocess.run(my_bf_cmd, check = True, stdout = f)
                    except subprocess.CalledProcessError as e:
                        print(e, "Failed at se Bloomfiltering", Array[0])
                        continue
            end_time_bf = time.time()
            runtime_bf = end_time_bf - start_time_bf
            print(runtime_bf, "Bloomfilter done")
            
            if Path(pe_item).is_file():
                my_sra_delete_cmd = ["rm", pe_item, item2]
            else:
                 my_sra_delete_cmd = ["rm", se_item]
            subprocess.run(my_sra_delete_cmd, check = True)
            print("SRA deleted")

    
            if "noMatch_1" in item_bloom:
                item_bloom_2 = re.sub("_1", "_2", item_bloom)
                my_mh_cmd = ["megahit", "-1", item_bloom, "-2", item_bloom_2, "-o", Array[0]] 
            else:    
                my_mh_cmd = ["megahit", "-r", item_bloom, "-o", Array[0]] 
            start_time_mh = time.time()
            try:
                result_mh = subprocess.run(my_mh_cmd, capture_output = True, text = True, check = True)
            except subprocess.CalledProcessError as e:
                        print(e, "Failed at Megahit", Array[0])
                        continue
            end_time_mh = time.time()
            runtime_mh = end_time_mh - start_time_mh
            print(runtime_mh, result_mh.stdout, result_mh.stderr, "Megahit done")
            
            
            my_mh_delete_cmd = ["find", Array[0] + "/.", "-mindepth", "1", "-not", "-name", "final.contigs.fa", "-delete"]
            subprocess.run(my_mh_delete_cmd)
            print("Pipeline part 1 done")
        j += 1
        
            