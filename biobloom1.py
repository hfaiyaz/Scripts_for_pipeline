#!/usr/bin/env python3
import time
import subprocess
import sys
import re

#with open(sys.argv[1], "r") as f1:
    #for index, item in enumerate(f1):
        #item = item.strip()
        #if index%int(sys.argv[2]) == int(sys.argv[3]):
item = sys.argv[1]
item_out = item + "_out"
my_cmd = ["biobloomcategorizer", "--fq", "-t", sys.argv[2], "-p", item_out, "-f", sys.argv[3], item]
start_time = time.time()
result = subprocess.run(my_cmd, capture_output = True, text = True, check = True)
end_time = time.time()
runtime = end_time - start_time
print(runtime, result.stdout, result.stderr)