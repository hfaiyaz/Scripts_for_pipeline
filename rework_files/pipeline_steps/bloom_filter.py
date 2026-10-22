#!/usr/bin/env python3

"""Filters given reads from a fastq file through a Bloom filter."""

import glob
import subprocess
import my_utils
from pathlib import Path

def _bloom_filter(self):
        """Filters samples based on provided bloom filter reference."""
        if self.paired:
            fq_1, fq_2 = self.get_files(self.fq)
            item_bloom_1 = Path(str(self.DIR) + '/' + self._accession + "_outnoMatch_1.fq")
            item_bloom_2 = Path(str(self.DIR) + '/' + self._accession + "_outnoMatch_2.fq")
            self._add_files(files=[item_bloom_1,item_bloom_2],
                           list_to_add=self.bloom_out)
            bloom_cmd = [
                "biobloomcategorizer",
                "-d",
                "-n",
                "-t",
                str(self.cores),
                "-e",
                "-p",
                self.bloom_prefix,
                "-f",
                self.bloom_filter,
                fq_1,
                fq_2
                ]
                
            print(f"passing command: {bloom_cmd}")    

            #command to split single stdout stream into two noMatch files
            awk_cmd = [
                "awk",
                "-v",
                f"out1={item_bloom_1}",
                "-v",
                f"out2={item_bloom_2}",
                '{print > (int((NR-1)/4)%2==0 ? out1 : out2)}'
            ]

            try:
                with subprocess.Popen(
                    bloom_cmd,
                    stdout=subprocess.PIPE
                ) as bf_out:

                    bloom_run = my_utils.run_command(awk_cmd,
                                check_state=True,
                                step="paired-end Bloom filtering",
                                stdin_redirect=bf_out.stdout,
                                extra_error=self._accession
                                )

                    bf_out.stdout.close()

                    return_code = bf_out.wait()

                    if (return_code != 0) or (not bloom_run):

                        print(f"{self.accession} failed with exit code {return_code}")
                        return False
                    return True
            except subprocess.CalledProcessError as e:
                print(f"Failed at separation of paired-end Bloom filters.", e, self._accession)
                return False

            
        else:
            fq_se = self.get_files(self.fq)
            item_bloom_se = Path(str(self.DIR) + '/' + self._accession + "_outnoMatch.fq")
            self._add_files(item_bloom_se,self.bloom_out)
            print(self.bloom_out)
            bloom_cmd = [
                    "biobloomcategorizer",
                    "-d",
                    "-n",
                    "-t",
                    str(self.cores),
                    "-p",
                    self.bloom_prefix,
                    "-f",
                    self.bloom_filter,
                    fq_se
                ]
            print(f"passing command {bloom_cmd}")    
                
            with open(self.bloom_out[0], "w") as f:
                se_bloom_run = my_utils.run_command(bloom_cmd,
                step="single end Bloom filtering",
                check_state=True,
                stdout_redirect=f)
                return se_bloom_run
            
def _bloom_threshold_check(self,thresh:float=0.20):
    tsv_file = Path(str(self.DIR) + "/" + self._accession + "_out_summary.tsv")
    with open(tsv_file,"r",encoding="utf-8") as summary:
        for line in summary:
            if "noMatch" in line:
                nm_line = line.strip().split("\t")
                print(f"{self._accession} noMatch rate: {nm_line[4]}")
                if float(nm_line[4]) > thresh:
                    #deletes ALL files associated with the accession, will not touch folders in current state
                    rm_pattern = glob.glob(self.DIR + self._accession + "*.*")
                    my_bloom_delete_cmd = ["rm"] + rm_pattern
                    print(f"Files to be deleted: {rm_pattern}")
                    bloom_thresh_delete = my_utils.run_command(my_bloom_delete_cmd,
                                step="Large Bloom Filter Deletion",
                                check_state=True,
                                extra_error=f"WARNING! Large noMatch files (over {thresh}%) could not be deleted.")
                    return bloom_thresh_delete
                else:
                    return True