import gseapy as gp
import pandas as pd
import re
import os
import multiprocessing
import matplotlib.pyplot as plt
from pathlib import Path

from Rrvgo_submit import Rrvgo_submit
from GO_dot_plot import GO_term_dot_plot

class GSEAPy_Prerank():
    """
    A structured pipeline for running GSEAPrerank with error handling, RRVGO support, and dot plot visualization.
    """

    def __init__(self, 
                 logFC_data_file, 
                 analysis_group, 
                 op_dir,
                 intermediate_op_dir=None,
                 table_op_dir=None,
                 taxon_id=9606,
                 gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'MSigDB_Hallmark_2020'],
                 r_exe_path="C:\\Program Files\\R\\R-4.3.1\\bin\\Rscript.exe", 
                 rrvgo_script_path="C:\\Users\\lwoods\\Documents\\LW_scripts\\R Scripts\\rrvgo_analyser_v2.R"):

        self.op_dir = Path(op_dir)
        self.intermediate_op_dir = Path(intermediate_op_dir) if intermediate_op_dir else self.op_dir
        self.table_op_dir = Path(table_op_dir) if table_op_dir else self.op_dir
        self.taxon_id = taxon_id
        if self.taxon_id == 9606:
            self.species = 'human'
        if self.taxon_id == 10090:
            self.species = 'mouse'
        self.gene_sets = gene_sets

        self.r_exe_path = r_exe_path
        self.rrvgo_script_path = str(Path(rrvgo_script_path))

        # GSEA Parameters
        self.logFC_data_file = logFC_data_file
        self.analysis_group = analysis_group
        self.df = None
        self.results = None

    def prepare_files(self):
        """Loads logFC data from multiple formats."""
        try:
            if isinstance(self.logFC_data_file, str):
                ext = Path(self.logFC_data_file).suffix.lower()
                if ext in ['.txt', '.tsv']:
                    self.df = pd.read_csv(self.logFC_data_file, sep='\t')
                elif ext == '.csv':
                    self.df = pd.read_csv(self.logFC_data_file)
                elif ext == '.xlsx':
                    self.df = pd.read_excel(self.logFC_data_file)
            elif isinstance(self.logFC_data_file, pd.DataFrame):
                self.df = self.logFC_data_file

            if self.df is None or self.df.empty:
                raise ValueError("Input data file is empty or invalid.")
        except Exception as e:
            print(f"Error loading logFC file: {e}")
            self.df = None

    def gsea_for_groups(self):
        """Runs GSEAPrerank for groups with error handling."""
        try:

            df = self.df[['Genes', f'{self.analysis_group}_logFC']].dropna()
            df = df.rename(columns={'Genes': 'Gene', f'{self.analysis_group}_logFC': 'log2FC'})

            # Auto-detect available CPUs
            if "SLURM_CPUS_PER_TASK" in os.environ:
                allocated_cpus = int(os.environ["SLURM_CPUS_PER_TASK"])
                num_threads = max(8, min(allocated_cpus, 12))  # Limit to 8-12 threads
            else:
                total_cpus = multiprocessing.cpu_count()
                num_threads = max(1, int(total_cpus * 0.75))  # Use 75% of available CPUs

            print(f"Using {num_threads} threads for GSEAPrerank ({self.analysis_group})")

            self.results = gp.prerank(
                rnk=df,
                gene_sets=self.gene_sets,
                permutation_num=1000,
                outdir=str(self.table_op_dir),
                min_size=10,
                max_size=500,
                threads=num_threads
            )
            gsea_result_file = Path(self.table_op_dir) / 'gseapy.gene_set.prerank.report.csv'

            if self.results.res2d.empty:
                print(f"No enrichment results found for {self.analysis_group}.")
                return

            self.write_rank_plot(self.results, self.analysis_group)
            df = pd.read_csv(gsea_result_file)

            # Extract the gene set category (first part before "__")
            df["Category"] = df["Term"].apply(lambda x: x.split("__")[0] if "__" in x else "Unknown")


            #{category}.{self.species}.enrichr.reports.txt'
            # Save separate files for each gene set category
            for category, sub_df in df.groupby("Category"):
                output_file = self.table_op_dir / f'{category}.prerank.report.csv'
                sub_df.to_csv(output_file, index=False)

        except Exception as e:
            print(f"Error running GSEAPrerank for {self.analysis_group}: {e}")

    def write_rank_plot(self, gs_res, group):
        """Generates enrichment ranking plots."""
        try:
            terms = gs_res.res2d.Term[1:5]
            hits = [gs_res.results[t]['hits'] for t in terms]
            runes = [gs_res.results[t]['RES'] for t in terms]

            figs = gp.gseaplot2(terms=terms, hits=hits, RESs=runes,
                                rank_metric=gs_res.ranking,
                                legend_kws={'loc': (1.2, 0)},
                                figsize=(5, 5))

            fig = figs[0] if isinstance(figs, list) else figs
            rank_fig_file = Path(self.op_dir) / f"{group}_gsea_rankplot.png"
            fig.figure.savefig(rank_fig_file, dpi=300, bbox_inches='tight')
        except Exception as e:
            print(f"Error generating rank plot for {group}: {e}")

    def run_rrvgo(self):
        """Run RRVGO to reduce GO terms. Continue plotting even if RRVGO isn't available or applicable."""
        if self.results is None:
            print("No results available for term reduction.")
            return
        
        #Name,Term,ES,NES,NOM p-val,FDR q-val,FWER p-val,Tag %,Gene %,Lead_genes
        #gseapy.gene_set.prerank.report.csv

        for category in self.gene_sets:
            try:
                #inp_file = str(Path(self.table_op_dir) / 'gseapy.gene_set.prerank.report.csv')
                inp_file = Path(self.table_op_dir) / f'{category}.prerank.report.csv'
                if not inp_file.is_file():
                    print(f"Input file for {category} not found, skipping.")
                    continue

                inp_df = pd.read_csv(inp_file, sep=",", usecols=['Term', 'FDR q-val', 'NES', 'ES','Gene %', 'Lead_genes'])
                inp_df = inp_df.rename(columns={'Term': 'term', 'FDR q-val': 'fdr', 'Lead_genes':'Genes'})
                inp_df['description'] = inp_df['term']
                inp_df['term'] = inp_df['term'].apply(lambda x: re.sub(r".*(GO:\d+).*", r"\1", x))

                # Check if RRVGO is applicable (for GO categories)
                if 'GO_' in category:
                    rrvgo_inp = str(Path(self.intermediate_op_dir) / f'{category}.gsea_prerank_for_rrvgo.csv')
                    inp_df.to_csv(rrvgo_inp, index=False)
                    ontology = None

                    if "GO_Biological_Process" in category:
                        ontology = "BP"
                    elif "GO_Molecular_Function" in category:
                        ontology = "MF"
                    elif "GO_Cellular_Component" in category:
                        ontology = "CC"

                    try:
                        rrvgo = Rrvgo_submit(
                            inp_file=rrvgo_inp,
                            op_dir=self.intermediate_op_dir,
                            ontology=ontology,
                            taxon_id=self.taxon_id,
                            fdr_thresh=1        # No threshold
                        )
                        rrvgo.run()

                        rrvgo_result_file = Path(self.intermediate_op_dir) / f'rrvgo_{ontology}_reduced_terms.csv'
                        if rrvgo_result_file.is_file():
                            rrvgo_result_df = pd.read_csv(rrvgo_result_file, usecols=['parent', 'parentTerm'])
                            rrvgo_result_df.drop_duplicates(inplace=True)
                            reduced_full_df = rrvgo_result_df.merge(inp_df, how='left', left_on='parent', right_on='term')
                            reduced_term_file = Path(self.table_op_dir) / f'{ontology}_rrvgo.csv'
                            reduced_full_df.to_csv(reduced_term_file)
                        else:
                            print(f"No RRVGO results for {category}, skipping.")
                            continue

                        go_plot = GO_term_dot_plot(
                            data_file=str(reduced_term_file),
                            op_dir=self.op_dir,
                            test_type="GSEA",
                            category=category,
                            analysis_group=self.analysis_group,
                            sig_thresh=0.25,
                            reduced = True
                        )
                        try:
                            go_plot.run()
                        except FileNotFoundError:
                            print(f"Plotting failed for {category}.")
                    except Exception as e:
                        print(f"Error during RRVGO processing for {category}: {e}")
                        # Continue to plot without RRVGO reduction
                        go_plot = GO_term_dot_plot(
                            data_file=str(rrvgo_inp),
                            op_dir=self.op_dir,
                            test_type="GSEA",
                            category=category,
                            analysis_group=self.analysis_group,
                            sig_thresh=0.25
                        )
                        go_plot.run()
                else:
                    # If RRVGO is not applicable, still plot
                    print(f"Rrvgo not available for {category}, proceeding with plot.")
                    new_inp_file = str(Path(self.intermediate_op_dir) / f'{category}.csv')
                    inp_df.to_csv(new_inp_file, index=False)
                    go_plot = GO_term_dot_plot(
                        data_file=str(new_inp_file),
                        op_dir=self.op_dir,
                        test_type="GSEA",
                        category=category,
                        analysis_group=self.analysis_group,
                        sig_thresh=0.25,
                    )
                    go_plot.run()
            except Exception as e:
                print(f"Error while processing {category}: {e}")
                continue  # Continue with the next category

    def run(self):
        """Runs the full pipeline."""
        self.prepare_files()
        if self.df is None:
            print("Skipping GSEA due to missing input data.")
            return
        
        self.gsea_for_groups()
        self.run_rrvgo()

if __name__ == '__main__':
    LOGFC_DATA_FILE = r"C:\Users\lwoods\Documents\LW_Projects_folder\01_MO\G07_Mendez\R07_SergioMuñoz\P03\myE02\myE02_data\tables\DA_data\consolidated_results.xlsx"
    ANALYSIS_GROUP = "Re-replicating_vs_Control"
    OP_DIR = r"C:\Users\lwoods\Documents\LW_Projects_folder\01_MO\G07_Mendez\R07_SergioMuñoz\P03\myE02\myE02_data"

    gene_sets = [
            'GO_Biological_Process_2023','MSigDB_Hallmark_2020'
        ]
    GSEA = GSEAPy_Prerank(logFC_data_file=LOGFC_DATA_FILE,
        analysis_group=ANALYSIS_GROUP,
        op_dir=OP_DIR)
    GSEA.run()