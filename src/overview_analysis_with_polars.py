
import logging
from pathlib import Path
import pandas as pd
import plotly.express as px
import polars as pl

class OverviewAnalysis():
    def __init__(self,
        filtered_df,
        conditions_df,
        column_mapping,
        sorted_conditions,
        colour_map,
        op_folder,
        use_DIANN_maxLFQ=False,
        log_file=None,
        report_data = None
    ):
                
        self.spec_full_df = filtered_df
        self.conditions_df = conditions_df
        self.column_mapping = column_mapping
        self.sorted_conditions = sorted_conditions
        self.colour_map = colour_map
        self.op_folder = op_folder
        self.use_DIANN_maxLFQ = use_DIANN_maxLFQ
        self.log_file = log_file
        self.report_data = report_data

        self.logging_setup()

        # Generated vars
        self.file_name_col = None
        self.protein_group_col = None
        self.protein_ids_col = None
        self.genes_col = None
        self.precursor_id_col = None
        self.precursor_quantity_col = None
        self.modified_sequence_col = None
        self.protein_names_col = None
        self.run_col = 'Run'
        self.species_col = 'Species'
        self.prot_by_run_df = None
        self.pept_by_run_df = None
        self.avg_pept_per_prot_df = None
        self.runs = None
    
    def logging_setup(self):
        logging.basicConfig(filename=self.log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', force=True)
    
    def enforce_mapped_columns(self):
        self.file_name_col = self.column_mapping['file_name_col']
        self.protein_group_col = 'Reattributed.Protein'
        #self.protein_ids_col = self.column_mapping['protein_ids_col']
        self.protein_ids_col = 'Reordered.Protein.Ids'
        #self.genes_col = self.column_mapping['genes_col']
        self.genes_col = 'Reattributed.Genes'
        self.precursor_id_col = self.column_mapping['precursor_id_col']
        self.precursor_quantity_col = self.column_mapping['precursor_quantity_col']
        self.modified_sequence_col = self.column_mapping['modified_sequence_col']
        #self.protein_names_col = self.column_mapping['protein_names_col']
        self.protein_names_col = 'Reattributed.Names'

    def get_plotting_data(self):
        if isinstance(self.spec_full_df, pd.DataFrame):
            self.spec_full_df = pl.from_pandas(self.spec_full_df)

        # Ensure columns exist in the DataFrame
        required_cols = ['CombinedName', self.protein_group_col, 'Condition', 'Bio.Rep', 'Tech.Rep']
        missing_cols = [col for col in required_cols if col not in self.spec_full_df.columns]
        if missing_cols:
            raise ValueError(f"Missing columns in DataFrame: {missing_cols}")
        
        print(1)
        print(self.spec_full_df.columns)

        # Get the total number of unique samples
        total_samples = self.spec_full_df.select('CombinedName').n_unique()

        # Count the number of unique samples each protein appears in
        protein_counts = (
            self.spec_full_df
            .group_by(self.protein_group_col)
            .agg(pl.col('CombinedName').n_unique().alias('sample_count'))
        )

        # Filter proteins that appear in at least 80% of the samples
        required_count = int(total_samples * 0.8)
        common_proteins = protein_counts.filter(pl.col('sample_count') >= required_count)[self.protein_group_col]

        # Filter the original DataFrame to keep only valid proteins
        filtered_spec_df = self.spec_full_df.filter(pl.col(self.protein_group_col).is_in(common_proteins))

        # Find number of proteins by run
        prot_by_run_df = (
            self.spec_full_df
            .group_by('CombinedName')
            .agg(pl.col(self.protein_group_col).n_unique().alias('Number of protein groups'))
        )

        common_prot_by_run_df = (
            filtered_spec_df
            .group_by('CombinedName')
            .agg(pl.col(self.protein_group_col).n_unique().alias('Number of protein groups'))
        )

        prot_by_run_df = prot_by_run_df.join(
            self.spec_full_df.select(['CombinedName', 'Condition', 'Bio.Rep', 'Tech.Rep']).unique(),
            on='CombinedName',
            how='left'
        )

        logging.info("PROTEIN SUMMARY")
        # Find number of proteins by run
        
        self.prots_sum = prot_by_run_df['Number of protein groups'].sum()
        logging.info(f'Number of proteins across all experiments: {self.prots_sum}')
        # Number of different proteins (at least one run) in experiment is
        self.num_prots_uniq = self.spec_full_df.select(pl.col(self.protein_group_col)).n_unique()
        logging.info(f'Number of different proteins in at least one experiment: {self.num_prots_uniq}')
        self.num_common_prots = filtered_spec_df.select(pl.col(self.protein_group_col)).n_unique()
        logging.info(f'Number of different proteins in at least one experiment: {self.num_prots_uniq}')
        # Number of proteins by run
        logging.info('Number of proteins by run:\n')
        logging.info(prot_by_run_df)

        logging.info("PEPTIDE SUMMARY")
        logging.info("Peptides here refer distinct modified AA sequences")
        pept_by_run_df = (
            filtered_spec_df
            .group_by('CombinedName')
            .agg(pl.col(self.modified_sequence_col).n_unique().alias('Number of peptides'))
        )
        pept_by_run_df = pept_by_run_df.join(
            self.spec_full_df.select(['CombinedName', 'Condition', 'Bio.Rep', 'Tech.Rep']).unique(), 
            on='CombinedName', 
            how='left'
        )
        self.pept_sum = pept_by_run_df['Number of peptides'].sum()
        self.num_pepts_uniq = self.spec_full_df.select(pl.col(self.modified_sequence_col)).n_unique()

        logging.info(f'Number of peptides overall: {self.pept_sum}')
        logging.info(f'Number of different peptides in at least one experiment: {self.num_pepts_uniq}')
        logging.info('Number of peptides by run:')
        logging.info(pept_by_run_df)

        logging.info("PROTEIN-LEVEL INFORMATION")
        logging.info("Average number of peptides per protein")
        avg_pept_per_prot_df = (
            self.spec_full_df
            .group_by(['CombinedName', 'Reattributed.Protein'])
            .agg(pl.col(self.modified_sequence_col).count().alias('peptide_count'))
            .group_by('CombinedName')
            .agg(pl.col('peptide_count').mean().alias('Average number of peptides per protein'))
        )
        avg_pept_per_prot_df = avg_pept_per_prot_df.join(self.conditions_df, on='CombinedName', how='left')

        logging.info(avg_pept_per_prot_df)

        # Storage for report vars
        if self.report_data:
            self.report_data['summary_data']['total_proteins'] = self.prots_sum
            self.report_data['summary_data']['unique_proteins'] = self.num_prots_uniq
            self.report_data['summary_data']['total_peptides'] = self.pept_sum
            self.report_data['summary_data']['unique_peptides'] = self.num_pepts_uniq

            self.report_data['tables'].append({
                'name': 'Proteins per run',
                'data': prot_by_run_df.to_pandas()
            })
            self.report_data['tables'].append({
                'name': 'Peptides per run',
                'data': pept_by_run_df.to_pandas()
            })
            self.report_data['tables'].append({
                'name': 'Peptides per protein',
                'data': avg_pept_per_prot_df.to_pandas()
            })
        
        self.prot_by_run_df = prot_by_run_df.to_pandas()
        self.common_prot_by_run_df = common_prot_by_run_df.to_pandas()
        self.pept_by_run_df = pept_by_run_df.to_pandas()
        self.avg_pept_per_prot_df = avg_pept_per_prot_df.to_pandas()
        self.runs = self.prot_by_run_df['CombinedName'].tolist()
        return self.spec_full_df

    
    def save_plot(self, fig, title, y_title, file_name):
        # HTML Configuration
        config = {
            f'toImageButtonOptions': {
                'format': 'png',  # one of png, svg, jpeg, webp
                'filename': f'{file_name}.png'},
                'scale': 3  # Multiply title/legend/axis/canvas sizes by this factor
        }

        fig.update_layout(
            title=title,
            template='plotly_white',
            title_font_color='black',
            xaxis=dict(
                title="Sample",
                titlefont_size=20,
                color='black'
            ),
            yaxis=dict(
                title=y_title,
                titlefont_size=16,
                color='black'
            ),
            legend=dict(
                font=dict(
                    color='black'
                )
            )
        )
        plot_folder = Path(self.op_folder)
        plot_folder.mkdir(exist_ok=True)
        output_path = plot_folder / f'{file_name}.html'
        fig.write_html(output_path, config=config)

        if self.report_data:
            self.report_data['plots'].append({
                'name': f'{title}',
                'fig': fig.to_html()
            })
 

    def write_plots(self):
        print(self.prot_by_run_df[self.prot_by_run_df['Condition'].isin(self.sorted_conditions)])

        fig1 = px.bar(self.prot_by_run_df[self.prot_by_run_df['Condition'].isin(self.sorted_conditions)],
                    y='Number of protein groups', x='CombinedName',
                    color='Condition', hover_data=['Condition', 'Number of protein groups'],
                    color_discrete_map=self.colour_map,
                    category_orders={'Condition': self.sorted_conditions})  # Apply sorted order
        title = f'Identified proteins by run<br><sup>Num. proteins across samples: <b>{self.num_prots_uniq}</b>; Num. proteins in ≥80% of samples: <b>{self.num_common_prots}</b>.'
        self.save_plot(fig1, title=title, y_title='Identified proteins', file_name="Proteins_by_run")

        fig2 = px.bar(self.pept_by_run_df[self.pept_by_run_df['Condition'].isin(self.sorted_conditions)],
                    y='Number of peptides', x='CombinedName',
                    color='Condition', hover_data=['Condition', 'Number of peptides'],
                    color_discrete_map=self.colour_map,
                    category_orders={'Condition': self.sorted_conditions})
        title = f'Identified peptides by run<br><sup> Num. peptides across samples: <b>{self.num_pepts_uniq}</b>.'
        self.save_plot(fig2, title=title, y_title='Identified peptides', file_name="Peptides_by_run")

        fig3 = px.bar(self.avg_pept_per_prot_df[self.avg_pept_per_prot_df['Condition'].isin(self.sorted_conditions)],
                    x='CombinedName', y='Average number of peptides per protein',
                    color='Condition', hover_data=['Condition', 'Average number of peptides per protein'],
                    color_discrete_map=self.colour_map,
                    category_orders={'Condition': self.sorted_conditions})
        title = f'Average peptides per protein by run'
        self.save_plot(fig3, title=title, y_title='Average peptides per protein', file_name="Average_peptides_per_prot")
    

    def run(self):
        logging.info("Organising columns.")
        self.enforce_mapped_columns()
        logging.info("Organising overview tables")
        self.get_plotting_data()
        logging.info("Writing overview plots")
        self.write_plots()
