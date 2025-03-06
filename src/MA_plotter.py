import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path

class PlotMA:
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
                    #row['Protein_Group'] = "CON__" + row['Protein_Group']
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

    def create_ma_plot(self):

        or_rd_r_colors = [
            "#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84",
            "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000"
        ]

        for comparison in self.comparison_columns:
            log2FC_col = f'{comparison}_log2FC'
            log2_mean_col = f'{comparison}_log2_mean_intensity'

            cond1, cond2 = comparison.split('_vs_')

            # Define color map for each point based on logFC values with different percentiles or POI status
            log2FC_95th = self.ordered_dataframe[log2FC_col].abs().quantile(0.95)
            log2FC_75th = self.ordered_dataframe[log2FC_col].abs().quantile(0.75)

            def color_mapper(row):
                if row['Genes'] in self.poi_dict:
                    return self.poi_dict[row['Genes']]
                if self.mark_contams:
                    if row['is_contaminant']=='Yes':
                        return '#8D7467'
                if abs(row[log2FC_col]) > log2FC_95th:
                    return or_rd_r_colors[8]
                elif abs(row[log2FC_col]) > log2FC_75th:
                    return or_rd_r_colors[5]
                else:
                    return or_rd_r_colors[3]

            colors = self.ordered_dataframe.apply(color_mapper, axis=1)

            # Define hover text including average log2 intensity
            hover_text = self.ordered_dataframe.apply(
                lambda row: f"<b>Gene:</b> {row['Genes']}<br>"
                            f"<b>Protein ID:</b> {row['Protein_Group'].split(';')[0]}<br>"
                            f"<b>Contaminant?:</b> {row['is_contaminant']}<br>"
                            f"<b>log2FC:</b> {row[log2FC_col]:.2f}<br>"
                            f"<b>Average log2(Intensity):</b> {row[log2_mean_col]:.2f}<br><br>"
                            f"Click to link to Uniprot for this protein.",
                axis=1
            )

            # Create the MA plot
            fig = go.Figure(
                data=[go.Scatter(
                    x=self.ordered_dataframe[log2_mean_col],
                    y=self.ordered_dataframe[log2FC_col],
                    mode='markers',
                    marker=dict(
                        color=colors,
                        size=6
                    ),
                    line=dict(color='white'),
                    opacity=0.7,
                    text=hover_text,  # Hover text with gene, protein ID, contaminant, log2FC, and average log2 intensity
                    hoverinfo='text',
                    customdata=self.ordered_dataframe['Protein_Group'].apply(lambda x: x.split(';')[0]),
                    showlegend=False
                )]
            )

            # Include the number of proteins in the title
            num_proteins = len(self.ordered_dataframe[log2FC_col].dropna())
            fig.update_layout(
                title=f"MA Plot: {cond1} vs {cond2} (N={num_proteins})",
                xaxis_title="Average log2(Intensity)",
                yaxis_title="log2 Fold Change",
                template='plotly_white',
                font=dict(color='black'),
                legend=dict(
                    itemsizing='constant'
                )
            )

            # Manually add legend items for different percentiles and POIs
            fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                marker=dict(size=10, color=or_rd_r_colors[3]),
                name='0-75% Percentile'
            ))

            fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                marker=dict(size=10, color=or_rd_r_colors[5]),
                name='75-95% Percentile'
            ))

            fig.add_trace(go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                marker=dict(size=10, color=or_rd_r_colors[8]),
                name='95%+ Percentile'
            ))

            # Add legend items for POIs
            for poi, color in self.poi_dict.items():
                fig.add_trace(go.Scatter(
                    x=[None],
                    y=[None],
                    mode='markers',
                    marker=dict(size=10, color=color),
                    name=f'POI: {poi}'
                ))
            
            if self.mark_contams:
                fig.add_trace(go.Scatter(
                    x=[None],
                    y=[None],
                    mode='markers',
                    marker=dict(size=10, color='#8D7467'),
                    name="Contaminants"
                ))
            
            # Save plot as HTML with a specific div ID
            plot_path = f"{self.op_folder}/MA_plot_{cond1}_vs_{cond2}.html"
            fig.write_html(plot_path, include_plotlyjs='cdn', full_html=True, div_id='ma-plot-div')

            # Append JavaScript to handle click events to open the corresponding URL
            with open(plot_path, "a") as html_file:
                html_file.write("""
<script>
  document.addEventListener('DOMContentLoaded', function() {
    var plot = document.getElementById('ma-plot-div');

    if (plot) {
      console.log('Plot element found:', plot);

      plot.on('plotly_click', function(data) {
        console.log('Click event triggered:', data);
        var point = data.points[0];
        var proteinId = point.customdata;
        if (proteinId) {
          var url = 'https://www.uniprot.org/uniprot/' + proteinId;
          console.log('Opening URL:', url);
          window.open(url, '_blank');
        }
      });
    } else {
      console.error('Plot element not found. Check if the ID is correct.');
    }
  });
</script>
                """)






    def run_analysis(self):
        self.load_file_inputs()
        self.mark_contaminants_dataframe()
        self.calculate_log2FC()
        self.create_ma_plot()


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

    plot_ma = PlotMA(
        bkgd_file=r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\LOESS\DAPARv1.30.6_DA_Sheet1_output.tsv",
        metadata_file=r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\LOESS\metadata_for_DAPAR.xlsx",
        op_folder=r"C:\Users\lwoods\Documents\LW_Projects_folder\TESTS\DIA_summarizer",
        comparison_columns=["413_vs_DMSO"],
        condition_colour_map=condition_colour_map,
        mark_contams=True,
        contaminants_fasta=r"C:\Users\lwoods\Documents\LW_Projects_folder\general\20210914_525_CONTAMINANTS_UNIPROT_FORMAT.fasta",
        poi_dict=poi_dict
    )

    plot_ma.run_analysis()