#!/usr/bin/env python3

"""All-module script utilities"""

import io
import subprocess

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
