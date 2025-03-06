from common_classes import *

class EnrichrORA():

    """Note: Enrichr uses a list of Entrez gene symbols as input. You should convert all gene names to uppercase."""

    def __init__(self,
                 gene_list,
                 op_dir,
                 background,
                 analysis_group,
                 intermediate_op_dir=None,
                 table_op_dir=None,
                 taxon_id=9606,
                 gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'MSigDB_Hallmark_2020'],
                 r_exe_path="C:\\Program\\\ Files\\R\\R-4.3.1\\bin\\Rscript.exe", 
                 rrvgo_script_path="C:\\Users\\lwoods\\Documents\\LW_scripts\\R Scripts\\rrvgo_analyser_v2.R"):

        self.gene_list = [gene.upper() for gene in gene_list]
        self.op_dir = op_dir
        # Put outputs together if optional folder(s) not specified
        if not table_op_dir:
            self.table_op_dir = self.op_dir
        else:
            self.table_op_dir = table_op_dir
        if not intermediate_op_dir:
            self.intermediate_op_dir = self.op_dir
        else:
            self.intermediate_op_dir = intermediate_op_dir

        self.intermediate_op_dir = intermediate_op_dir
        self.background = background
        self.analysis_group = analysis_group
        self.gene_sets = gene_sets
        self.taxon_id = taxon_id
        self.results = None  # to store the enrichment results
        if self.taxon_id == 9606:
            self.species = 'human'
        if self.taxon_id == 10090:
            self.species = 'mouse'
        self.r_exe_path = str(Path(r_exe_path))
        self.rrvgo_script_path = str(Path(rrvgo_script_path))

    def enrichr_api_call(self):
        try:
            enr_bg = gp.enrichr(gene_list=self.gene_list,
                                gene_sets=self.gene_sets,
                                background=self.background,
                                organism=self.species,
                                cutoff=1,
                                outdir=self.table_op_dir)  # write to disk
            self.results = enr_bg.results  # store the results
            if self.results.empty:
                print("No enrichment results found.")
                return None
            return self.results.head()  # return the top results
        except ValueError as e:
            print(f"Error in Enrichr API call: {e}")
            return None

    def plot_enrichr(self):
        if self.results is None or self.results.empty:
            print("No enrichment results found--skipping plot!")
            return

        try:
            op_file = str(Path(self.intermediate_op_dir) / 'ORA_barplot.png')
            barplot(self.results,
                column="Adjusted P-value",
                group='Gene_set',
                ofname=op_file,
                size=10,
                top_term=5,
                figsize=(3, 5),
                cutoff=0.25,
                color='red')
        except Exception as e:
            print(f"Error while plotting: {e}")

    def run_rrvgo(self):
        """Run RRVGO to reduce GO terms. Continue plotting even if RRVGO isn't available or applicable."""
        if self.results is None or self.results.empty:
            print("No results available for term reduction.")
            return

        for category in self.gene_sets:
            try:
                inp_file = str(Path(self.table_op_dir) / f'{category}.{self.species}.enrichr.reports.txt')
                if not Path(inp_file).is_file():
                    print(f"Input file for {category} not found, skipping.")
                    continue

                inp_df = pd.read_csv(inp_file, sep="\t", usecols=['Term', 'Adjusted P-value', 'Odds Ratio', 'Combined Score', 'Genes'])
                inp_df = inp_df.rename(columns={'Term': 'term', 'Adjusted P-value': 'fdr'})
                inp_df['description'] = inp_df['term']
                inp_df['term'] = inp_df['term'].apply(lambda x: re.sub(r".*(GO:\d+).*", r"\1", x))

                # Check if RRVGO is applicable (for GO categories)
                if 'GO_' in category:
                    rrvgo_inp = str(Path(self.intermediate_op_dir) / f'{category}.{self.species}.enrichr_for_rrvgo.csv')
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
                            fdr_thresh=1
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
                            test_type="ORA",
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
                            test_type="ORA",
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
                        test_type="ORA",
                        category=category,
                        analysis_group=self.analysis_group,
                        sig_thresh=0.25,
                    )
                    go_plot.run()
            except Exception as e:
                print(f"Error while processing {category}: {e}")
                continue  # Continue with the next category

    def run(self):
        self.enrichr_api_call()
        self.plot_enrichr()
        self.run_rrvgo()

if __name__ == "__main__":

    GENE_LIST = ['TP53', 'EGFR', 'BRCA1', 'BRCA2', 'PTEN', 'AKT1', 'PIK3CA', 'BRAF', 'KRAS', 'CTNNB1']

    overlaps_file = r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\LOESS\comparison_overlaps_for_upset.csv"
    BACKGROUND=pd.read_csv(overlaps_file)['Genes'].values.tolist()

    OP_DIR = r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\LOESS\UpSet_ORA\Decreased_368_413"
    TAXON_ID = 9606
    ANALYSIS_GROUP = "Upset Group: Decreased in 368 and 417"

    # Initialize the class
    enrichr = EnrichrORA(
        gene_list=GENE_LIST,
        op_dir=r"C:\Users\lwoods\Downloads\TEEST_plot",
        table_op_dir=r"C:\Users\lwoods\Downloads\TEEST_table",
        intermediate_op_dir=r"C:\Users\lwoods\Downloads\TEEST_intermediate",
        background=BACKGROUND,
        analysis_group=ANALYSIS_GROUP,
        taxon_id=TAXON_ID,
        #gene_sets=['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'MSigDB_Hallmark_2020'],
        #colour_map={'GO_Biological_Process_2023': 'salmon', 'GO_Cellular_Component_2023': 'navy', 'GO_Molecular_Function_2023': 'green', 'MSigDB_Hallmark_2020': 'magenta'}
        #gene_sets=['MSigDB_Hallmark_2020'],
        #colour_map={'MSigDB_Hallmark_2020': 'magenta'}
    )

    enrichr.run()
