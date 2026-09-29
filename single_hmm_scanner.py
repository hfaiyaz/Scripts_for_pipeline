#!/usr/bin/env python3

import subprocess
import sys

hmm_reference = sys.argv[1]
hmm_file = sys.argv[2]
num_cpus = sys.argv[3]
out_file = sys.argv[4]

hmm_cmd = ["hmmscan","--incE", "0.00001", "-o", out_file, "--cpu",str(num_cpus), hmm_reference, hmm_file]

try:
    subprocess.run(hmm_cmd,check=True)
except subprocess.CalledProcessError as e:
    print("Could not finish file, see specific error for more details:",e)
    sys.exit()