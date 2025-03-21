import subprocess
from pathlib import Path
import yaml
import pandas as pd


class Rrvgo_submit():
    def __init__(self, r_exe_path, rrvgo_script_path, inp_file, op_dir, ontology, taxon_id, fdr_thresh=0.25, verbose=False):
        self.r_exe_path = r_exe_path  # Removed trailing comma
        self.r_script_path = rrvgo_script_path  # Removed trailing comma
        self.inp_file = str(Path(inp_file))
        self.op_dir = str(Path(op_dir))
        self.taxon_id = str(taxon_id)  # Ensure this is a string
        self.ontology = str(ontology)  # Ensure this is a string
        self.fdr_thresh = str(fdr_thresh)
        self.verbose = verbose
        self.reduced_full_df = None

    def build_Rrvgo_command(self):
        # Construct command as a list

        command = [
            self.r_exe_path,
            self.r_script_path,
            self.inp_file,
            self.op_dir,
            self.ontology,
            self.taxon_id,
            self.fdr_thresh
        ]

        # Ensure all arguments are strings
        command = [str(arg) for arg in command]
        return command

    def run_Rrvgo(self):
        """
        Run the Rrvgo command using subprocess with a timeout and logging.
        """
        command_1 = self.build_Rrvgo_command()
        if self.verbose:
            print("Running command:", " ".join(command_1))  # Print the full command for debugging

        # Redirect output and errors to log files
        log_file = Path(self.op_dir) / "rrvgo_output.log"
        error_file = Path(self.op_dir) / "rrvgo_error.log"

        try:
            with open(log_file, "w") as log, open(error_file, "w") as err:
                # Run the command with a timeout
                result = subprocess.run(
                    command_1,
                    stdout=log,
                    stderr=err,
                    text=True,
                    timeout=1200  # Timeout in seconds (e.g., 20 minutes)
                )

            # Check return code
            if result.returncode != 0:
                print(f"R script returned an error (see {error_file} for details).")
            else:
                print(f"R script executed successfully (see {log_file} for details).")

        except subprocess.TimeoutExpired:
            print(f"R script timed out after 20 minutes. Check {log_file} and {error_file} for details.")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

    
    def run(self):
        self.run_Rrvgo()