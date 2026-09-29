#!/usr/bin/env python3

import os
import sys
import subprocess

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

            Output = Dir
            OutFile = Output + '/' + Array[0] + '_1.fastq'

            if os.path.exists(OutFile):
                j += 1
                continue

            subprocess.run(
                ['fasterq-dump', Array[0], '--split-files', '-O', Output]
            )

        j += 1
