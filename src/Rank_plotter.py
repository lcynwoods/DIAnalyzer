import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path

class PlotRank:
    def __init__(self, bkgd_file, metadata_file, op_folder, comparison_columns, condition_colour_map, mark_contams=True, contaminants_fasta=None, poi_dict=None):
        self.bkgd_file = bkgd_file
        self.metadata_file = metadata_file
        self.op_folder = op_folder
        self.comparison_columns = comparison_columns
        self.condition_colour_map = condition_colour_map
        self.mark_contams = mark_contams
        self.contaminants_fasta = contaminants_fasta
        self.poi_dict = poi_dict if poi_dict else {}

        self.ordered_dataframe = None
        self.contaminant_ids = None
        self.intensity_columns = None

        or_rd_r_colors = [
            "#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84",
            "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000"
        ]

        self.light_colour= or_rd_r_colors[3]
        self.medium_colour= or_rd_r_colors[5]
        self.dark_colour= or_rd_r_colors[8]

    def load_file_inputs(self):
        # Load input files
        if '.tsv' in str(self.bkgd_file):
            self.ordered_dataframe = pd.read_csv(self.bkgd_file, sep="\t")
        elif '.xlsx' in str(self.bkgd_file):
            self.ordered_dataframe = pd.read_excel(self.bkgd_file)
        self.metadata = pd.read_excel(self.metadata_file, sheet_name="Sheet1")
        
        if self.contaminants_fasta:
            with open(self.contaminants_fasta, 'r') as f:
                lines = f.readlines()

            # Extract protein IDs (assuming uniprot ids)   
            self.contaminant_ids = [line.strip().lstrip('>').split("|")[1] for line in lines if line.startswith('>')]

    def mark_contaminants_dataframe(self):
        def mark_contaminants(row, contaminants_list):
            # Split the Protein_Group by ";"
            proteins = row['Protein_Group'].split(";")
            row['is_contaminant'] = "No"
            for protein in proteins:
                if contaminants_list and protein in contaminants_list:
                    row['Protein_Group'] = "CON__" + row['Protein_Group']
                    row['is_contaminant'] = "Yes"
                    break
                elif "CON__" in protein:
                    row['is_contaminant'] = "Yes"
                    break
            return row

        if self.mark_contams:
            self.ordered_dataframe = self.ordered_dataframe.apply(mark_contaminants, args=(self.contaminant_ids,), axis=1)

    def calculate_log2FC(self):
        # Calculate log2 fold changes (M) and average log2 intensities (A) for each comparison
        for comparison in self.comparison_columns:
            cond1, cond2 = comparison.split('_vs_')

            # Extract relevant intensity columns
            cond1_columns = [f'Intensity_{sample}' for sample in self.metadata[self.metadata['Condition'] == cond1]['Sample.Name']]
            cond2_columns = [f'Intensity_{sample}' for sample in self.metadata[self.metadata['Condition'] == cond2]['Sample.Name']]

            # Calculate average intensity per condition and log2 fold change, skipping NA values
            self.ordered_dataframe[f'{comparison}_log2_mean_{cond1}'] = np.log2(self.ordered_dataframe[cond1_columns].mean(axis=1, skipna=True) + 1e-10)
            self.ordered_dataframe[f'{comparison}_log2_mean_{cond2}'] = np.log2(self.ordered_dataframe[cond2_columns].mean(axis=1, skipna=True) + 1e-10)

            # Log2 fold change: M = log2(cond1/cond2)
            self.ordered_dataframe[f'{comparison}_log2FC'] = self.ordered_dataframe[f'{comparison}_log2_mean_{cond1}'] - self.ordered_dataframe[f'{comparison}_log2_mean_{cond2}']

            # Average log2 intensity: A = (log2(cond1) + log2(cond2)) / 2
            self.ordered_dataframe[f'{comparison}_log2_mean_intensity'] = (self.ordered_dataframe[f'{comparison}_log2_mean_{cond1}'] + self.ordered_dataframe[f'{comparison}_log2_mean_{cond2}']) / 2

    def create_rank_plot(self):
        for comparison in self.comparison_columns:
            cond1, _ = comparison.split('_vs_')
            log2_mean_col = f'{comparison}_log2_mean_{cond1}'

            # Sort by average log2 intensity for condition 1 and get the rank
            self.ordered_dataframe = self.ordered_dataframe.sort_values(by=log2_mean_col, ascending=False).reset_index(drop=True)
            self.ordered_dataframe['rank'] = self.ordered_dataframe.index + 1

            # Define color and size for each point based on POIs and percentile
            log2_mean_95th = self.ordered_dataframe[log2_mean_col].quantile(0.95)
            log2_mean_75th = self.ordered_dataframe[log2_mean_col].quantile(0.75)

            def color_mapper(row):
                if row['Genes'] in self.poi_dict:
                    return self.poi_dict[row['Genes']]
                elif row[log2_mean_col] > log2_mean_95th:
                    return self.dark_colour
                elif row[log2_mean_col] > log2_mean_75th:
                    return self.medium_colour
                else:
                    return self.light_colour

            def size_mapper(row):
                return 20 if row['Genes'] in self.poi_dict else 6

            colors = self.ordered_dataframe.apply(color_mapper, axis=1)
            sizes = self.ordered_dataframe.apply(size_mapper, axis=1)

            # Define hover text excluding log2FC
            hover_text = self.ordered_dataframe.apply(
                lambda row: f"<b>Gene:</b> {row['Genes']}<br>"
                            f"<b>Protein ID:</b> <a href='https://www.uniprot.org/uniprot/{row['Protein_Group'].split(';')[0]}' target='_blank' style='background-color: white;'>{row['Protein_Group']}</a><br>"
                            f"<b>Contaminant:</b> {row['is_contaminant']}<br>"
                            f"<b>Average log2 Intensity:</b> {row[log2_mean_col]:.2f}<br>"
                            f"<b>Rank:</b> {row['rank']}",
                axis=1
            )

            # Create the rank vs average log2 intensity line plot
            fig = go.Figure()

             # Plot line that changes color based on percentile cutoffs
            segments = []
            current_color = 'lightblue'
            start_index = 0

            for index, row in self.ordered_dataframe.iterrows():
                if row[log2_mean_col] > log2_mean_95th:
                    new_color = self.dark_colour
                elif row[log2_mean_col] > log2_mean_75th:
                    new_color = self.medium_colour
                else:
                    new_color = self.light_colour

                if new_color != current_color:
                    # To smooth out transitions, add a small blend between segments
                    if start_index != index:
                        blend_index = min(index + 2, len(self.ordered_dataframe))  # Adding slight overlap for blending
                        segments.append((start_index, blend_index, current_color))
                    start_index = index
                    current_color = new_color

            segments.append((start_index, len(self.ordered_dataframe), current_color))

            for start, end, color in segments:
                fig.add_trace(go.Scatter(
                    x=self.ordered_dataframe['rank'][start:end],
                    y=self.ordered_dataframe[log2_mean_col][start:end],
                    mode='lines',
                    line=dict(color=color, width=5),
                    hoverinfo='skip',
                    showlegend=False
                ))

            # Add scatter points for POIs
            poi_df = self.ordered_dataframe[self.ordered_dataframe['Genes'].isin(self.poi_dict.keys())]
            hover_text = poi_df.apply(
                lambda row: f"<b>Gene:</b> {row['Genes']}<br>"
                            f"<b>Protein ID:</b> {row['Protein_Group']}<br>"
                            f"<b>Contaminant:</b> {row['is_contaminant']}<br>"
                            f"<b>Average log2 Intensity:</b> {row[log2_mean_col]:.2f}<br>"
                            f"<b>Rank:</b> {row['rank']}",
                axis=1
            )

            fig.add_trace(go.Scatter(
                x=poi_df['rank'],
                y=poi_df[log2_mean_col],
                mode='markers',
                marker=dict(
                    color=[self.poi_dict[gene] for gene in poi_df['Genes']],
                    size=20
                ),
                text=hover_text,
                hoverinfo='text',
                showlegend=False
            ))


            # Include the number of proteins in the title
            num_proteins = len(self.ordered_dataframe[log2_mean_col].dropna())
            fig.update_layout(
                title=f"Rank vs Average log2 Intensity Plot: {cond1} (N={num_proteins})",
                xaxis_title="Rank",
                yaxis_title="Average log2 Intensity",
            )

            # Manually add legend items for different percentiles and POIs
            fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                marker=dict(size=8, color=self.light_colour),
                name='0-75% Percentile'
            ))

            fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                marker=dict(size=8, color=self.medium_colour),
                name='75-95% Percentile'
            ))

            fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                marker=dict(size=8, color=self.dark_colour),
                name='95%+ Percentile'
            ))

            # Add legend items for POIs
            for poi, color in self.poi_dict.items():
                fig.add_trace(go.Scatter(
                    x=[None],
                    y=[None],
                    mode='markers',
                    marker=dict(size=20, color=color),
                    name=f'POI: {poi}'
                ))

            # Save plot as HTML
            plot_path = f"{self.op_folder}/Rank_Intensity_plot_{cond1}.html"
            fig.update_layout(template='plotly_white',font=dict(color='black'))
            fig.write_html(plot_path)

    def run_analysis(self):
        self.load_file_inputs()
        self.mark_contaminants_dataframe()
        self.calculate_log2FC()
        self.create_rank_plot()


# Example Usage:

if __name__ == "__main__":
    condition_colour_map = {
        "DMSO": "red",
        "413": "blue"
    }

    poi_dict = {
        "GeneA": "green",
        "GeneB": "purple"
    }

    plot_rank = PlotRank(
        bkgd_file=r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\LOESS\DAPARv1.30.6_DA_Sheet1_output.tsv",
        metadata_file=r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\LOESS\metadata_for_DAPAR.xlsx",
        op_folder=r"C:\Users\lwoods\Documents\LW_Projects_folder\TESTS\DIA_summarizer",
        comparison_columns=["413_vs_DMSO"],
        condition_colour_map=condition_colour_map,
        mark_contams=True,
        contaminants_fasta=r"C:\Users\lwoods\Documents\LW_Projects_folder\general\20210914_525_CONTAMINANTS_UNIPROT_FORMAT.fasta",
        poi_dict=poi_dict
    )

    plot_rank.run_analysis()