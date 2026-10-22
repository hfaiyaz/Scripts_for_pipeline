#!/usr/bin/env python3

import glob
import io
import os
import sys
import subprocess
import time
from pathlib import Path

### GLOBAL FUNCTIONS
### NOW IN: utils.py
def run_command(
    cmd:list,
    step:str,
    cap_out:bool=False,
    check_state:bool=False,
    stdin_redirect:None|io.BufferedWriter=None,
    stdout_redirect:io.TextIOWrapper|None=None,
    extra_error:str|None=None
) -> bool|str:
    """Runs a subprocess command. Returns True if the run was successful, and false otherwise."""
    try:
        if cap_out:
            cmd_run = subprocess.run(cmd,
                            stdin=stdin_redirect,
                            stdout=stdout_redirect,
                            check=check_state,
                            capture_output=cap_out,
                            text=True)
        else:
            cmd_run = subprocess.run(cmd,
                            stdin=stdin_redirect,
                            stdout=stdout_redirect,
                            check=check_state,
                            capture_output=cap_out,
                            text=True)
        if cap_out:
            return cmd_run.stdout
        return True
    except subprocess.CalledProcessError as run_error:
        print(run_error,f"failed at {step}")
        if extra_error is not None:
            print(extra_error)
        return False

### BIG CLASS

class SampleAccession:
    """Class to contain sample accession paths and methods."""

    #TODO: get way to set automatically regardless of system
    BASE = Path(os.path.expandvars('/scratch/alpine/$USER'))
    DIR = Path(os.path.expandvars('/scratch/alpine/$USER/c_albicans_benchmark')) 

    def __init__(self, acc, n_cores, b_filter, n_filter, overwrite:bool=False):
        self._accession: str = acc
        self._paired = False
        self.fq = []
        self.bloom_filter = b_filter
        self.blastn_filter = n_filter
        self.bloom_prefix = str(self.DIR) + '/' + self._accession + "_out"
        self.bloom_out = []
        self.final_contigs = Path(self.DIR / self._accession / "final.contigs.fa")
        self.cores = n_cores
        self.overwrite = overwrite

    @property
    def accession(self):
        return self._accession

    @property
    def paired(self):
        return self._paired
    
    @paired.setter
    def paired(self,value:bool):
        self._paired = value

    def init_file_check(self):
    # check file existence 
        if self.final_contigs.is_file():
            if self.overwrite:
                print("final contigs file found, overwriting files...")
                self.delete_files(delete_step="pipeline_overwrite",
                                    glob_accession=True,
                                    delete_dirs=True)
            else:
                print(f"final contigs file found, skipping {self._accession}")
                return True
        return False
    
    def _add_files(self,
                 files: list|str,
                 list_to_add: list) -> None:
        """Adds a file to a provided list. Should only be used on self list
        objects such as fq or bloom_out."""
        if type(files) is list:
            for added_file in files:
                list_to_add.append(added_file)
        else:
            list_to_add.append(files)
        return None

    def get_files(self,
                fq_list:list|None=None):
        "Returns the contents of self-attribute file lists. Automatically"
        "assess for paired-end status."
        if self.paired:
            return fq_list[0],fq_list[1]
        else:
            return fq_list[0]

    def _check_file_exists(self,
        prefix:str,
        ending:str="",
        list_to_add:list|None=None) -> bool:
        """Checks if the given accession is single or paired end, and checks if the files exist.
        Requires built-in continue logic if neither file exists but sra_run is True.
        Will exit unless an iterator is provided."""


        # Define the directory and the prefix you want to match
        files = [f for f in self.DIR.glob(f"{prefix}*{ending}") if f.is_file()]

        if len(files) > 2:
            print(f"Deleting all files associated with sample {self._accession}")
            self._delete_files(delete_step="Removing Failed Run Contents",
            glob_accession=True,
            delete_dirs=True)
            return False
            
        if len(files) == 2:
            self.paired = True
            self._add_files(files,list_to_add)
        
        elif len(files) == 1:
            self._add_files(files,list_to_add)
            
        else:
            print(f"No files found on disk.")
            return False

        return True

    def _delete_files(self,
                    delete_step:str,
                    files_to_delete:list|str=[],
                    glob_accession:bool=False,
                    delete_dirs:bool=False,
                    cmd_override:list|None=None) -> bool|None: 
        """Deletes files provided in list or str format. Will exit loop unless iterator
        is provided."""
        if cmd_override is not None:
            print("WARNING! Using command override. Assuming defaults for run_command.")
            rm_run = run_command(cmd_override,
                                step=delete_step,
                                check_state=True,
                                extra_error=f"WARNING: Failed at {delete_step}. \
                                Accumulation of files may cause unwanted IO/server errors.")
            
        if glob_accession:
            #TODO: include prefix option so _accession isn't default
            files_to_delete = [f for f in self.DIR.glob(f"{self._accession}*") if f.is_file()]
        if delete_dirs:
            rm_base = ["rm","-r"]
        else:
            rm_base = ["rm"]
        if type(files_to_delete) is list:
            rm_cmd = rm_base + files_to_delete
        else:
            rm_cmd = rm_base.append(files_to_delete)

        rm_run = run_command(rm_cmd,
                    step=delete_step,
                    check_state=True,
                    extra_error=f"WARNING: Failed at {delete_step}. \
                    Accumulation of files may cause unwanted IO/server errors.")

        return rm_run

    ### PIPELINE STEPS ###

    def _sra_dump(self,stop_here:bool=False):
        sra_cmd = ['fasterq-dump',self.accession,'--split-files','-e',str(self.cores),'-O',self.DIR]

        sra_run = run_command(
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

                    bloom_run = run_command(awk_cmd,
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
                se_bloom_run = run_command(bloom_cmd,
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
                        bloom_thresh_delete = run_command(my_bloom_delete_cmd,
                                    step="Large Bloom Filter Deletion",
                                    check_state=True,
                                    extra_error=f"WARNING! Large noMatch files (over {thresh}%) could not be deleted.")
                        return bloom_thresh_delete
                    else:
                        return True

    def _megahit_assembly(self):
        """Runs the MEGAHIT assembly."""
        megahit_out = Path(str(self.DIR) + "/" + self._accession)
        if self._paired:
            mh_1, mh_2 = self.get_files(self.bloom_out)
            my_mh_cmd = [
                "megahit",
                "-1",mh_1,
                "-2",mh_2,
                "-t",str(self.cores),
                "-o",megahit_out
            ]
        else:
            mh_fq = self.get_files(self.bloom_out)
            my_mh_cmd = [
                "megahit",
                "-r",mh_fq,
                "-t",str(self.cores),
                "-o",megahit_out
            ]

        result_mh = run_command(
            my_mh_cmd,
            step="MEGAHIT assembly",
            check_state=True
            )

        return result_mh

    def get_contig_ids(self) -> set:
        """Looks in the final.contigs.fa file and returns all IDs present."""
        all_contig_ids = set()
        with open(self.final_contigs, "r") as contigf:
            for contig_line in contigf:
                if contig_line.startswith(">"):
                    contig_id = contig_line[1:].split()[0]
                    all_contig_ids.add(contig_id)
        return all_contig_ids


    def _delete_sra_cache(self):
        """Checks for valid SRA cache files and deletes the corresponding files to help save space.
        Should only be run after the entire pipeline is completed to avoid file I/O errors."""
        cache1 = Path(str(self.BASE) + "/sra/sra/" + self._accession + ".sra")
        cache2 = Path(str(self.BASE) + "/sra/sra/" + self._accession + ".sra.cache")
                    
        if not cache1.is_file() and not cache2.is_file():
            print("No SRA cache to delete.")
            pass
        else:
            if cache1.is_file():
                my_sra_cache_delete_cmd = ["rm", cache1]
            elif cache2.is_file():
                my_sra_cache_delete_cmd = ["rm", cache2]
            else:
                print("SRA cache not found")
            sra_cache_delete = run_command(my_sra_cache_delete_cmd,
                        step="SRA Cache Deletion",
                        check_state=True)

            return sra_cache_delete

    def _delete_folder_subcontents(self,dir_path:str|None=None,exclude:str="final.contigs.fa"):
        """Deletes subfolder contents generated by MEGAHIT."""
        if dir_path is None:
            dir_path = Path(str(self.DIR) + "/" + self._accession + "/.")
        my_mh_delete_cmd = [
            "find",
            dir_path,
            "-mindepth",
            "1",
            "-not",
            "-name",
            exclude,
            "-delete"
        ]
        mh_delete = run_command(my_mh_delete_cmd,
                                step="Megahit delete",
                                check_state=True)
        return mh_delete

    # def self_blastn(self):
    #     """Performs a self blastn against a specified database."""
    #     blast_output = Path(self.DIR / self._accession / "final.contigs_blastn.csv")
        
    #     my_blast_cmd = [os.path.expanduser("~/software/ncbi-blast-2.17.0+/bin/blastn"),
    #                     "-query",
    #                     self.final_contigs,
    #                     "-db",
    #                     self.blastn_filter,
    #                     "-outfmt",
    #                     "10 qseqid sseqid pident length qlen qcovs "
    #                     "mismatch gapopen qstart qend sstart send "
    #                     "evalue bitscore",
    #                     "-num_threads",
    #                     self.cores]
    #     self_blastn = run_command(my_blast_cmd,step="self blastn",cap_out=True,check_state=True)
    #     if type(self_blastn) == bool:
    #         return False
    #     with open(blast_output,"w",encoding="utf-8") as blastf:
    #                 blastf.write(
    #                             "qseqid,sseqid,pident,length,qlen,qcovs,"
    #                             "mismatch,gapopen,qstart,qend,sstart,send,"
    #                             "evalue,bitscore\n"
    #                         )
    #                 blastf.write(self_blastn)
    #     print("BLASTN done:", self.accession)

    #     return True

    # def header_id(header_line: str) -> str:
    #     """First whitespace-delimited token of a FASTA header (without '>')."""
    #     parts = header_line[1:].split()
    #     return parts[0] if parts else ""

    # def combine_fastas() -> list[str]:
    #     """Combines fasta files and adds accessions to every header"""
    #     print("\nCombining FASTA files...")
    
    #     files = sorted(INPUT_DIR.glob(f"*{FILE_SUFFIX}"))
    #     if not files:
    #         raise(f"No *{FILE_SUFFIX} files found in {INPUT_DIR}")
    
    #     ids: list[str] = []
    #     with open(COMBINED, "w") as out:
    #         for path in files:
    #             acc = path.name[: -len(FILE_SUFFIX)]
    #             with open(path) as fh:
    #                 for line in fh:
    #                     if line.startswith(">"):
    #                         line = ">" + acc + "_" + line[1:]
    #                         ids.append(header_id(line))
    #                     # make sure the last line of a file ends in a newline
    #                     # so the next file's header never gets glued onto it
    #                     if not line.endswith("\n"):
    #                         line += "\n"
    #                     out.write(line)
    
    #     print("Combined FASTA created:")
    #     print(COMBINED)
    #     return ids
    
    # def check_unique(ids: list[str]) -> set[str]:
    #     """Checks the uniqueness of given IDs"""
    #     print("\nChecking contig ID uniqueness...")
    
    #     total = len(ids)
    #     all_ids = set(ids)
    #     print(f"Total contigs:  {total}")
    #     print(f"Unique IDs:     {len(all_ids)}")
    
    #     if total != len(all_ids):
    #         raise("CONTIG IDS ARE NOT UNIQUE.\nPipeline stopped before BLAST.")
    
    #     print("PASS: All contig IDs are unique.")
    #     return all_ids
    
    # def check_blast_db() -> None:
    #     """Verifies the BLAST database"""
    #     print("\nChecking BLAST database...")
    
    #     try:
    #         subprocess.run(
    #             [str(BLASTDBCMD), "-db", DB, "-info"],
    #             check=True,
    #             stdout=subprocess.DEVNULL,
    #         )
    #     except (subprocess.CalledProcessError, FileNotFoundError) as exc:
    #         raise(f"BLAST database check failed: {exc}")
    
    #     print("BLAST database found:")
    #     return DB
    
    # def run_blastn() -> None:
    #     """Runs BLASTN on a provided query file"""
    #     print("\nRunning BLASTN...")
    
    #     # Falls back to 1 thread if not running under SLURM
    #     threads = os.environ.get(num_cores, "1")
    
    #     with open(BLAST_OUT, "w") as fh:
    #         fh.write(",".join(BLAST_FIELDS) + "\n")
    #         fh.flush()
    #         subprocess.run(
    #             [
    #                 str(BLASTN),
    #                 "-query", str(COMBINED),
    #                 "-db", DB,
    #                 "-outfmt", "10 " + " ".join(BLAST_FIELDS),
    #                 "-num_threads", str(threads),
    #             ],
    #             stdout=fh,
    #             check=True,
    #         )
    
    #     print("BLAST completed.")
    
#     def classify_hits() -> tuple[set[str], set[str], set[str]]:
#         """Return (contigs_with_any_hit, high_80_contigs, low_80_contigs)."""
#         hit_ids: set[str] = set()
#         high_ids: set[str] = set()
#         low_ids: set[str] = set()
    
#         # Split from the right so a comma inside sseqid cannot shift the columns:
#         # everything after sseqid is a fixed number of numeric fields.
#         n_trailing = len(BLAST_FIELDS) - 2
    
#         with open(BLAST_OUT) as fh:
#             next(fh, None)  # skip header
#             for line in fh:
#                 line = line.rstrip("\n")
#                 if not line:
#                     continue
    
#                 parts = line.rsplit(",", n_trailing)
#                 qid = parts[0].split(",", 1)[0]
#                 pident = float(parts[1])  # field 3 (pident)
#                 qcovs = float(parts[4])   # field 6 (qcovs)
    
#                 hit_ids.add(qid)
    
#                 if pident > THRESHOLD and qcovs > THRESHOLD:
#                     high_ids.add(qid)
#                 if pident < THRESHOLD and qcovs < THRESHOLD:
#                     low_ids.add(qid)
    
#         return hit_ids, high_ids, low_ids
 
#     def extract_negatives(negative_ids: set[str]) -> int:
#         """Extracts negative contigs from the combined file post-BLAST"""
#         print("\nExtracting negative contig sequences...")
    
#         count = 0
#         keep = False
#         with open(COMBINED) as fin, open(NEGATIVE_FASTA, "w") as fout:
#             for line in fin:
#                 if line.startswith(">"):
#                     keep = header_id(line) in negative_ids
#                     if keep:
#                         count += 1
#                 if keep:
#                     fout.write(line)
#         return count
 

#     def get_negative_contigs(self,blastdb,output_folder):
#         """Returns contigs that only pass a certain threshold out of the BLAST call."""
#         negative_contigs = (negative_output_dir / f"{self.accession}_negative_final_contigs.fa")

#         INPUT_DIR = Path("/pl/active/Viralogue/final_contigs/C_albicans") #where to source the contigs from
#         DB = "/pl/active/Viralogue/blastdb/candida_albicans/Candida" #blast_db
        
#         BLAST_BIN = Path.home() / "software" / "ncbi-blast-2.17.0+" / "bin" #find BLAST software
#         BLASTN = BLAST_BIN / "blastn"
#         BLASTDBCMD = BLAST_BIN / "blastdbcmd"
        
#         COMBINED = output_folder / "C_albicans_all_contigs_unique.fna"
#         BLAST_OUT = output_folder / "C_albicans_blastn_results_unique.csv"
#         NEGATIVE_FASTA = output_folder / "C_albicans_negative_contigs.fna"
        
#         FILE_SUFFIX = "_final.contigs.fa"
#         THRESHOLD = 80.0  # percent identity / query coverage cutoff
        
#         BLAST_FIELDS = [
#             "qseqid", "sseqid", "pident", "length", "qlen", "qcovs",
#             "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
#             "evalue", "bitscore",
# ]

#         print("Running BLASTN:", self.accession)

#         ids = combine_fastas()
#         all_ids = check_unique(ids)
#         total_contigs = len(ids)
        
#         check_blast_db()
#         run_blastn()
        
#         hit_ids, high_ids, low_ids = classify_hits()
#         num_hits = len(hit_ids)
#         print(f"\nContigs with at least one BLAST hit: {num_hits}")
        
#         no_hit_ids = all_ids - hit_ids
#         num_no_hits = len(no_hit_ids)
#         print(f"Contigs with no BLAST hit: {num_no_hits}")
        
#         num_high = len(high_ids)
#         print(f"Contigs with at least one >80/>80 hit: {num_high}")
        
#         num_low = len(low_ids)
#         print(f"Contigs with at least one <80/<80 hit: {num_low}")
        
#         # Remove any low-80 contig that also has a >80/>80 hit
#         low_no_high_ids = low_ids - high_ids
#         num_low_no_high = len(low_no_high_ids)
#         print(f"Contigs with <80/<80 hit AND no >80/>80 hit: {num_low_no_high}")
        
#         # Negative = no BLAST hit OR (<80/<80 hit with no >80/>80 hit)
#         negative_ids = no_hit_ids | low_no_high_ids
#         num_negative = len(negative_ids)
#         print(f"\nTotal negative contigs: {num_negative}")
        
#         final_count = extract_negatives(negative_ids)
        
#         print(f"\nNegative IDs expected:  {num_negative}")
#         print(f"Negative FASTA records: {final_count}")
        
#         if num_negative != final_count:
#             raise("Negative FASTA count does not match negative ID count.\nPipeline stopped.")
        
#         print(f"\nTotal input contigs:             {total_contigs}")
#         print(f"Contigs with any BLAST hit:      {num_hits}")
#         print(f"Contigs with no BLAST hit:       {num_no_hits}")
#         print(f"Contigs with >80/>80 hit:        {num_high}")
#         print(f"Contigs with <80/<80 hit:        {num_low}")
#         print(f"<80/<80 with no >80/>80 hit:    {num_low_no_high}")
#         print(f"Final negative contigs:          {final_count}")
#         print(f"\nCombined FASTA:\n{COMBINED}")
#         print(f"\nBLAST results:\n{BLAST_OUT}")
#         print(f"\nFinal negative FASTA:\n{NEGATIVE_FASTA}\n")

#         pass

    def hmmer_scan(self):
        pass
    

    def pipeline_run_phase_1(self,iter:int,blast_database:Path) -> int|Path:
        """Runs the pipeline steps and associated cleanup: \
            1. SRA acquisition \
            2. Bloom filtering \
            3. MEGAHIT assembly"""
        contigs_exist = self.init_file_check()
        if contigs_exist:
            return iter + 1
            
        SRA_exists = self._check_file_exists(prefix=self._accession, ending=".fastq",list_to_add=self.fq)
        bloom_exists = self._check_file_exists(prefix=self._accession,ending="noMatch",list_to_add=self.bloom_out)
        
        
        print(f"STARTING RUN:")
        if not SRA_exists:
            start_sra = time.time()
            sra = self._sra_dump(stop_here=True)
            if not sra:
                print(f"SRA failed in {time.time() - start_sra} seconds.")
                return iter + 1
            print(f"SRA completed in {time.time() - start_sra} seconds.")
        else:
            print("Found SRA file, continuing...")

        if not bloom_exists:
            start_bloom = time.time()
            bloom = self._bloom_filter()
            if not bloom:
                print(f"Bloom filtering failed in {time.time() - start_bloom} seconds.")
                return iter + 1
            
            bloom_thresh = self._bloom_threshold_check()
            print(f"bloom thresh: {bloom_thresh}")
            if not bloom_thresh:
                print(f"Bloom thresholding failed in {time.time() - start_bloom} seconds.")
                return iter + 1
            else:
                print(f"did we make it here???")
            sra_delete = self._delete_files("SRA FQ Delete",self.fq)
            if not sra_delete:
                return iter + 1
                
            print(f"Complete Bloom filtering and thresholding completed in {time.time() - start_bloom} seconds")
        else:
            print("Found bloom filtering files, continuing...")
        
        print("SRA Files deleted, continuing with MEGAHIT")

        start_megahit = time.time()
        megahit = self._megahit_assembly()
        if not megahit:
            print(f"MEGAHIT failed in {time.time() - start_megahit} seconds.")
            return iter + 1

        print(f"MEGAHIT completed in {time.time() - start_megahit} seconds.")
        cache_delete = self._delete_sra_cache()
        if not cache_delete:
            return iter + 1
        delete_mh_folder = self._delete_folder_subcontents()
        if not delete_mh_folder:
            return iter + 1
        delete_extraneous_files = self._delete_files("Bloom and summary delete",glob_accession=True)
        #concatenate??? incorporate later, not necessary right now
        return self.final_contigs

        start_blast = time.time()
        blast_run = self.get_negative_contigs(self,blast_database,)

    def pipeline_run_hmmer_phase(self):
        pass

    def pipeline_run_circular_phase(self):
        pass


### MAKE THESE INTO ARGPARSE ARGUMENTS
sample_list = sys.argv[1]
catch = int(sys.argv[2]) - 1
num_jobs = int(sys.argv[3])
num_cores = int(sys.argv[4])
bloom_filter_reference = sys.argv[5]
blastn_reference =  sys.argv[6]
negative_output_dir = Path(sys.argv[7])
negative_output_dir.mkdir(parents=True, exist_ok=True)

with open(sample_list, 'r') as acc_list:
    j = 0
    for line in acc_list:
        if j % num_jobs == catch:
            ### MAKE PRINT STATEMENTS INTO VERBOSE OUTPUT AS PART OF ARGPARSE
            print(f"Number={j}. NumJobs={num_jobs}, Catch={catch}")
            print(f"{j%num_jobs}")

            line = line.rstrip('\n')
            array = line.split('\t')

            sample = SampleAccession(acc=array[0],
                                     n_cores=num_cores,
                                     b_filter=bloom_filter_reference,
                                     n_filter=blastn_reference)

            print(f"Now processing sample {sample.accession}\n")

            pipeline = sample.pipeline_run_phase_1(j)  
            if type(pipeline) == int:
                j = pipeline
                continue  
            else:
                print(f"Pipeline part 1 complete: final.contigs.fa file stored at {pipeline}") 



        j += 1
        
####################################################################
#Fix underneath here
####################################################################
 
 
