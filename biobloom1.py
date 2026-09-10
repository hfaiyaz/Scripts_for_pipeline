#!/usr/bin/env python3
import time
import subprocess
import sys
import re

with open(sys.argv[1], "r") as f1:
    for index, item in enumerate(f1):
        item = item.strip()
        if index%int(sys.argv[2]) == int(sys.argv[3]):
            item_out = item + "_out"
            if "_" in item:
                item2 = re.sub("_1", "_2", item)
                my_cmd = ["biobloomcategorizer", "--fq", "-t", sys.argv[4], "-e", "-p", item_out, "-f", "/projects/fhasan1@xsede.org/bloom_references/mm39_input.bf", item, item2] 
            else:    
                my_cmd = ["biobloomcategorizer", "--fq", "-t", sys.argv[4], "-p", item_out, "-f", "/projects/fhasan1@xsede.org/bloom_references/mm39_input.bf", item]
            start_time = time.time()
            result = subprocess.run(my_cmd, capture_output = True, text = True, check = True)
            end_time = time.time()
            runtime = end_time - start_time
            print(runtime, result.stdout, result.stderr)