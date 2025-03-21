# -*- coding: utf-8 -*-
"""
Created on Wed Aug 30 15:34:12 2023

...
"""

import sys
import os
import getpass
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()
sys.path.append(str(Path(__file__).parent / "src"))

from src.common_classes import *
from src.analysis_loop import AnalysisLoop
from src.setup_utils import check_inputs, create_folders_and_logger
from src.data_prep_utils import perform_data_prep, sort_data_and_assign_colours

class DIAnalyzer:
    """
    A class to analyze Spectronaut or DIANN reports, perform data processing, stats, and push results to GitHub.
    """

    def __init__(self, spec_report, conditions_file, analyses_file, op_folder,
                 data_source='DIANN', remove_outliers=False, POC=False,
                 q1=0.05, q3=0.95, use_DIANN_maxLFQ=False,
                 reattribute_peptides=False, peptides_to_remove=None,
                 min_pepts=2, biomolecule_level='Protein', species_id=9606,
                 reference_condition=None, contaminants_fasta=None,
                 mark_contaminants=False, trispecies_test=False,
                 gene_sets =['GO_Biological_Process_2023', 'MSigDB_Hallmark_2020'],
                 log_file=None):
        """
        Initialize the DIAnalyzer class.
        
        Args:
            spec_report (str): Path to the Spectronaut or DIANN report file.
            conditions_file (str): Path to the conditions file.
            analyses_file (str): Path to the analyses file.
            op_folder (str): Output folder path.
            data_source (str, optional): Source format ('Spectronaut' or 'DIANN'). Defaults to 'DIANN'.
            remove_outliers (bool, optional): Whether to remove outliers. Defaults to False.
            POC (bool, optional): Proof of Concept flag. Defaults to False. If True, skips limma testing and plots MA and Rank plots instead.
            q1 (float, optional): Lower quantile threshold for outlier removal. Defaults to 0.05.
            q3 (float, optional): Upper quantile threshold for outlier removal. Defaults to 0.95.
            use_DIANN_maxLFQ (bool, optional): Whether to use DIANN maxLFQ. Defaults to False. DIANN maxLFQ supercedes all other normalizations.
            reattribute_peptides (bool, optional): Whether to reassign peptides to proteins of interest (POIs) where peptides are shared. Defaults to False.
            peptides_to_remove (str, optional): List separated by commas. Defaults to None. 
            min_pepts (int, optional): Minimum peptides required per protein. Defaults to 2. 
            biomolecule_level (str, optional): Biomolecule level ('Protein' or 'Peptide'). Defaults to 'Protein'.
            species_id (int, optional): Species ID. Should be radio buttons for Human (ID:9606) or Mouse (ID:10090). Defaults to 9606. Options depend on GO term availability (TODO)
            reference_condition (str, optional): Condition against which all other conditions are compared in relative POI abundance. Defaults (None) to first available.
            contaminants_fasta (str, optional): Path to contaminants FASTA file. Defaults to None.
            mark_contaminants (bool, optional): Whether to mark contaminants. Defaults to False. If True but no contaminants file provided, DIAnalyzer will identify contaminants as "CON__"
            trispecies_test (bool, optional): Whether the pipeline is being run to assess quantification performance with trispecies ratio test. Defaults to False. If True, no comparative analyses performed.
            log_file (str, optional): Path to log file. Defaults to None.
        """

        config_file = "./config/config.yaml"
        config_path = str(Path(config_file).resolve())
        print(f"Looking for config file at: {config_path}")
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)

        self.spec_report = spec_report
        self.conditions_file = conditions_file
        self.analyses_file = analyses_file
        self.op_folder = op_folder
        self.q1 = q1
        self.q3 = q3
        self.use_DIANN_maxLFQ = use_DIANN_maxLFQ
        self.reattribute_peptides = reattribute_peptides
        self.peptides_to_remove = peptides_to_remove
        self.remove_outliers = remove_outliers
        self.min_pepts = min_pepts
        self.biomolecule_level = biomolecule_level
        self.species_id = species_id
        self.reference_condition = reference_condition
        self.POC = POC
        self.contaminants_fasta = contaminants_fasta
        self.mark_contaminants = mark_contaminants
        self.trispecies_test = trispecies_test
        self.gene_sets = gene_sets
        self.log_file = log_file

        # Retrieve from config
        self.html_template = config['html_template']
        self.colour_palette = config['colour_palette']
        self.r_exe_path = str(Path(config['r_exe_path']))
        self.r_DAPAR_path = str(Path(config['r_DAPAR_path']))
        self.r_rrvgo_path = str(Path(config['r_rrvgo_path']))
        self.format = data_source  # 'Spectronaut' or 'DIANN'
        self.column_mapping = config[self.format]  # Get the mapping for the selected format
        self.gh_repo = config['gh_repo']
        self.local_repo_parent = config['local_repo_parent']
        self.gh_token = os.getenv('GH_TOKEN')
        self.username = getpass.getuser()
        self.colour_palette = config['colour_palette']  

        # Initialize some vars
        self.report_data = {
            'data_filtering': {},
            'summary_data': {},
            'tables': [],
            'plots': []
        }

        self.spec_full_df = None
        self.conditions_df = None
        self.sorted_conditions = None
        self.colour_mapper = None
        self.POIs=None
        self.POI_colour_map=None
        #self.intensity_per_prot_df = None
        #self.pivoted=None

    def setup(self):
        check_inputs(
            self.spec_report, self.conditions_file, self.analyses_file, self.op_folder,
            self.r_exe_path, self.r_DAPAR_path, self.r_rrvgo_path, self.species_id, self.gh_token
        )
        create_folders_and_logger(self.op_folder, self.log_file)

    def perform_data_prep(self):
        self.spec_full_df, self.conditions_df = perform_data_prep(
            self.spec_report, self.conditions_file, self.column_mapping, self.reattribute_peptides,
            self.peptides_to_remove, self.remove_outliers, self.min_pepts, self.log_file, self.report_data
        )

    def sort_data_and_assign_colours(self):
        self.sorted_conditions, self.reference_condition, self.colour_mapper, self.POIs, self.POI_colour_map = sort_data_and_assign_colours(
            self.conditions_df, self.reference_condition, self.colour_palette
        )
    
    def perform_quantification_data_analysis(self):
        print(self.use_DIANN_maxLFQ)
        job_name = Path(self.op_folder).stem
        plot_folder = Path(self.op_folder) / f'{job_name}_data' / 'plots'
        intermediate_folder = Path(self.op_folder) / f'{job_name}_data' / 'intermediate_files'
        summary_data_folder = plot_folder / 'summary_data'
        summary_intermediate_op_folder = intermediate_folder / 'summary_data'
        for folder in [summary_data_folder, summary_intermediate_op_folder]:
            Path(folder).mkdir(parents=True, exist_ok=True)
        quantification_data = PrepareQuantificationData(
            filtered_df=self.spec_full_df,
            conditions_df=self.conditions_df,
            sorted_conditions=self.sorted_conditions,
            colour_map=self.colour_mapper,
            log_file=self.log_file,
            op_folder=summary_data_folder,
            intermediate_op_folder=summary_intermediate_op_folder,
            use_DIANN_maxLFQ=self.use_DIANN_maxLFQ,
            biomolecule_level=self.biomolecule_level,
            POIs=self.POIs,
            POI_colour_map=self.POI_colour_map,
            #is_IP=self.is_IP,
            #is_POC=self.is_POC,
            report_data=self.report_data
        )
        self.spec_full_df = quantification_data.run()

    def perform_overview_data_analysis(self):
        print(self.use_DIANN_maxLFQ)
        job_name = Path(self.op_folder).stem
        plot_folder = Path(self.op_folder) / f'{job_name}_data' / 'plots'
        plot_folder.mkdir(exist_ok=True, parents=True)
        summary_data_folder = plot_folder / 'summary_data'
        overview_data = OverviewAnalysis(
            filtered_df=self.spec_full_df,
            conditions_df=self.conditions_df,
            sorted_conditions=self.sorted_conditions,
            column_mapping=self.column_mapping,     
            colour_map=self.colour_mapper,
            op_folder=summary_data_folder,
            use_DIANN_maxLFQ=self.use_DIANN_maxLFQ,
            log_file=self.log_file,
            report_data=self.report_data
        )
        overview_data.run()

    def loop_through_analyses(self):
        job_name = Path(self.op_folder).stem
        intermediate_folder = Path(self.op_folder) / f'{job_name}_data' / 'intermediate_files' / 'summary_data'

        my_analyses = AnalysisLoop(
            self.analyses_file,
            intermediate_folder / f'{self.biomolecule_level}_quantification_for_R.xlsx',
            self.conditions_df,
            self.op_folder,
            self.r_exe_path,
            self.r_DAPAR_path,
            self.r_rrvgo_path,
            self.colour_mapper,
            sorted_conditions=self.sorted_conditions,
            reference_condition=self.reference_condition,
            POIs=self.POIs,
            POI_colour_map=self.POI_colour_map,
            use_DIANN_maxLFQ=self.use_DIANN_maxLFQ,
            biomolecule_level=self.biomolecule_level,
            species_id=self.species_id,
            POC=self.POC,
            contaminants_fasta=self.contaminants_fasta,
            mark_contaminants=self.mark_contaminants,
            gene_sets=self.gene_sets
        )
        my_analyses.run()
    
    def push_to_github(self):
        repo_url = self.gh_repo
        folder_path = self.op_folder
        local_repo_path = Path(self.local_repo_parent) / self.username / self.gh_repo.split('/')[-1]
        token = self.gh_token
        upload_html_to_github(repo_url, folder_path, local_repo_path, token)

    def run(self):
        self.setup()
        self.perform_data_prep()
        self.sort_data_and_assign_colours()
        self.perform_quantification_data_analysis()
        self.perform_overview_data_analysis()
        if not self.trispecies_test:
            self.loop_through_analyses()
        if not self.trispecies_test:
            self.push_to_github()
        #self.generate_html_report()

if __name__ == "__main__":
    SPEC_REPORT = r"C:\Users\lwoods\Documents\LW_Projects_folder\01_MO\G12_Eva_Gonzalez\R03_Mario_Rodriguez\P02\myE01_RANKL1_IP_triplicates\20250320122458_51676a6c-4e27-4917-b4e6-5a6dfc4a47c6_subset\diannsummary\diann_report.parquet"
    CONDITIONS_FILE = r"\\bespin\Biotechnology\Proteomics\01_MO\G12_Eva_Gonzalez\R03_Mario Rodriguez\P02_RANKL1_IP\E02_RANKL1_Nt_Turbo_Triplicates\conditions.xlsx"
    ANALYSES_FILE = r"\\bespin\Biotechnology\Proteomics\01_MO\G12_Eva_Gonzalez\R03_Mario Rodriguez\P02_RANKL1_IP\E02_RANKL1_Nt_Turbo_Triplicates\analyses.xlsx"
    OP_FOLDER =  r"C:\Users\lwoods\Documents\LW_Projects_folder\01_MO\G12_Eva_Gonzalez\R03_Mario_Rodriguez\P02\myE01_RANKL1_IP_triplicates"
    REMOVE_OUTLIERS = False
    POC = False
    USE_DIANN_MAXLFQ = False
    SPECIES_ID = 9606
    # integer, no range
    MIN_PEPTS = 2
    REFERENCE_CONDITION = "EV"
    MARK_CONTAMINANTS = True
    CONTAMINANTS_FASTA = r"C:\Users\lwoods\Documents\LW_Projects_folder\general\20210914_525_CONTAMINANTS_UNIPROT_FORMAT.fasta"
    TRISPECIES_TEST = False

    # Create an instance of the SpectronautAnalyzer class
    analyzer = DIAnalyzer(
        spec_report=SPEC_REPORT,
        conditions_file=CONDITIONS_FILE,
        analyses_file=ANALYSES_FILE,
        op_folder=OP_FOLDER,
        remove_outliers=REMOVE_OUTLIERS,
        use_DIANN_maxLFQ=USE_DIANN_MAXLFQ,
        species_id=SPECIES_ID,
        min_pepts=MIN_PEPTS,
        reference_condition=REFERENCE_CONDITION,
        POC=POC,
        contaminants_fasta=CONTAMINANTS_FASTA,
        mark_contaminants=MARK_CONTAMINANTS,
        trispecies_test=TRISPECIES_TEST
    )
    analyzer.run()
    logging.shutdown()