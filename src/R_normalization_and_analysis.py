import subprocess
import datetime

class R_normalization_and_analysis():
    def __init__(self,
        r_exe_path,
        r_script_path,
        op_file,
        raw_quant_file,
        metadata_file,
        comparisons = "All",
        MS_data_type="DIA",
        biomolecule_level='protein',
        DAPAR_norm_type = "LOESS",
        DAPAR_normalization_within ='overall',
        th_filter = ">=0.5",
        th_push = ">=0.5",
        batch_col = "No"
        ):
        
        self.data_file = raw_quant_file
        self.metadata_file = metadata_file
        self.op_file = op_file
        self.r_exe_path = r_exe_path
        self.r_script_path = r_script_path
        self.comparisons = comparisons
        self.MS_data_type = MS_data_type
        self.biomolecule_level = biomolecule_level
        self.DAPAR_norm_type = DAPAR_norm_type
        self.DAPAR_normalization_within = DAPAR_normalization_within
        self.th_filter = th_filter
        self.th_push = th_push
        self.batch_col = batch_col

    def build_DAPAR_command(self):
        command = [
            self.r_exe_path,
            self.r_script_path,
            self.op_file,
            self.data_file,
            self.metadata_file,
            self.comparisons,
            self.MS_data_type,
            self.biomolecule_level,
            self.DAPAR_norm_type,
            self.DAPAR_normalization_within,
            self.th_filter,
            self.th_push,
            self.batch_col
                   ]
        command = [str(arg) for arg in command]
        return command

    def run_DAPAR(self):
        command_1 = self.build_DAPAR_command()
        print(command_1)
        result = subprocess.run(command_1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        # Print the output and errors, if any
        print("Output:\n", result.stdout)
        print("Errors:\n", result.stderr)
        # Check return code
        if result.returncode != 0:
            print("R script returned an error:", result.returncode)
        print(command_1)