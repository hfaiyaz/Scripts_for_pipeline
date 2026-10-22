#!/usr/bin/env python3

import my_utils

"""Grabs SRA files from a provided list."""

def _sra_dump(self,stop_here:bool=False):
        sra_cmd = ['fasterq-dump',self.accession,'--split-files','-e',str(self.cores),'-O',self.DIR]

        sra_run = my_utils.run_command(
            sra_cmd,
            step="SRADump",
            check_state=True,
            extra_error=f"SRADump could not be completed, see specific error for more details.")
        
        self._check_file_exists(prefix=self._accession,
                                ending=".fastq",
                                list_to_add=self.fq)
        if stop_here:
            return False
        return sra_run

