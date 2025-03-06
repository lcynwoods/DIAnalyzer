import subprocess
from pathlib import Path
import yaml
import pandas as pd


class Rrvgo_submit():
    def __init__(self, inp_file, op_dir, ontology, taxon_id, fdr_thresh=0.25):

        self.inp_file = str(Path(inp_file))
        self.op_dir = str(Path(op_dir))
        self.taxon_id = str(taxon_id)  # Ensure this is a string
        self.ontology = str(ontology)  # Ensure this is a string
        self.fdr_thresh = str(fdr_thresh)

        config_file = "./config.yaml"
        config_path = str(Path(config_file).resolve())
        print(f"Looking for config file at: {config_path}")
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)

        self.r_exe_path = str(Path(config['r_exe_path']))
        self.r_script_path = str(Path(config['r_rrvgo_path']))  

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
        command_1 = self.build_Rrvgo_command()
        print("Running command:", " ".join(command_1))  # Print the full command for debugging

        # Run the command
        result = subprocess.run(command_1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
        
        # Print the output and errors, if any
        print("Output:\n", result.stdout)
        print("Errors:\n", result.stderr)

        # Check return code
        if result.returncode != 0:
            print("R script returned an error:", result.returncode)

        else:
            print("R script executed successfully.")

    
    def run(self):
        self.run_Rrvgo()