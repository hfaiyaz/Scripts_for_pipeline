#!/usr/bin/env python3
import time
import subprocess
import sys
import re

with open(sys.argv[1], "r") as f1:
    for index, item in enumerate(f1):
        item = item.strip()
        item = "/scratch/alpine/fhasan1@xsede.org/sra/" + item
        if index%int(sys.argv[2]) == int(sys.argv[3]):
            item_out = item + "_out"
            if "noMatch_1" in item:
                item2 = re.sub("_1", "_2", item)
                my_cmd = ["megahit", "-1", item, "-2", item2, "-o", item_out] 
            else:    
                my_cmd = ["megahit", "-r", item, "-o", item_out] 
            start_time = time.time()
            result = subprocess.run(my_cmd, capture_output = True, text = True, check = True)
            end_time = time.time()
            runtime = end_time - start_time
            print(runtime, result.stdout, result.stderr)