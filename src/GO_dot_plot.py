import pandas as pd
from pathlib import Path
import plotly.express as px
import numpy as np
import textwrap  # For wrapping the text

class GO_term_dot_plot():

    def __init__(self, data_file, op_dir, test_type, category, analysis_group, n_terms=20, sig_thresh=0.10, reduced=False):
        self.data_file = data_file
        self.op_dir = op_dir
        self.test_type = test_type
        self.category = category
        self.analysis_group = analysis_group
        self.n_terms = n_terms
        self.sig_thresh = sig_thresh
        self.reduced = reduced

        self.id_col="Genes"
        self.sig_col='fdr'
        if self.test_type == 'ORA':
            print("It's an ORA")
            self.score = 'Odds Ratio'
        elif self.test_type == 'GSEA':
            print("It's a GSEA")
            self.score = 'NES'
        
        self.data_df = None
    
    def load_data(self):
        # Load data and conditions
        if isinstance(self.data_file, str):
            if '.xlsx' in self.data_file:
                self.data_df = pd.read_excel(self.data_file)
            elif '.csv' in self.data_file:
                print("CSV")
                self.data_df = pd.read_csv(self.data_file)
            elif ('.tsv' in self.data_file or '.txt' in self.data_file):
                self.data_df = pd.read_csv(self.data_file, sep='\t')
        else:
            return

    def wrap_text(self, text, width=75, truncate=True):
        """Wrap the text into lines of max `width` characters and truncate with '...' if too long."""
        wrapped_text = "<br>".join(textwrap.wrap(text, width=width))
        if truncate and len(wrapped_text) > 150:  # Truncate only if the total length exceeds 150 characters
            return wrapped_text[:147] + '...'  # Truncate to 147 characters, and add '...'
        else:
            return wrapped_text

    def plot_GO(self):
        
        def percent_to_float(percent_str):
            return float(percent_str.rstrip('%')) / 100
        
        def count_genes(gene_list):
            return len(gene_list.split(";"))

        # Make a copy of the DataFrame to work with
        df_plot = self.data_df.copy()

        # Handle infinite values first (for both the score and Odds Ratio columns)
        df_plot.replace([np.inf, -np.inf], np.nan, inplace=True)  # Replace inf with NaN first

        # Convert necessary columns to numeric
        df_plot[self.score] = pd.to_numeric(df_plot[self.score])
        df_plot[self.sig_col] = pd.to_numeric(df_plot[self.sig_col])

        max_thresh = 0.5  # Maximum allowable FDR threshold
        min_required_results = 3  # Minimum number of results you want

        # Keep loosening the threshold if there are fewer than `min_required_results`
        while True:
            filtered_df = df_plot[df_plot[self.sig_col] < self.sig_thresh]  # Filter by FDR q-value
            if self.test_type == "ORA":
                filtered_df = filtered_df[filtered_df[self.score] > 0]  # Ensure score is positive
            
            # If enough rows meet the criteria or we've hit the max threshold, stop
            if len(filtered_df) >= min_required_results or self.sig_thresh >= max_thresh:
                break
            
            # Loosen the FDR threshold incrementally
            self.sig_thresh += 0.05

        # Now apply the final filtering to `df_plot`
        df_plot = filtered_df

        # Replace NaN in score with the maximum finite value
        max_score = df_plot[self.score].dropna().max()
        df_plot[self.score].fillna(max_score, inplace=True)

        # Wrap the 'description' and 'Genes' text for better readability
        df_plot['description'] = df_plot['description'].apply(lambda x: self.wrap_text(x, truncate=True))
        df_plot['Genes'] = df_plot['Genes'].apply(lambda x: self.wrap_text(x, truncate=False))  # Wrap 'Genes' text without truncation

        title_size_label = None
        # Only enriched
        if self.test_type == "GSEA":
            df_plot['size'] = df_plot['Gene %'].apply(percent_to_float)
            hover_size_label = "Gene %"
            # Get the actual min and max sizes used in the plot
            min_size = df_plot['size'].min()
            max_size = df_plot['size'].max()
            title_size_label = f'% in gene set (min: {min_size*100:.2f}; max: {max_size*100:.2f})'
            size_description = '% Genes in Set'
            score_description = "Normalized Enrichment Score"
        
        elif self.test_type == 'ORA':
            df_plot['size'] = df_plot['Genes'].apply(count_genes)
            df_plot = df_plot[df_plot['size'] > 1]
            df_plot["Number of genes"] = df_plot['size']
            # Get the actual min and max sizes used in the plot
            min_size = df_plot['size'].min()
            max_size = df_plot['size'].max()
            hover_size_label = "Number of genes"
            title_size_label = f'number of genes (min: {min_size}; max: {max_size})'
            size_description = 'Number of genes'
            score_description = "Odds Ratio*<br><sub>*As calculated by EnrichR: https://maayanlab.cloud/Enrichr/help#background"

        plot_title = f"Annotation terms from gene set \"{self.category}\" significant in {self.analysis_group}<sub><br>FDR < {self.sig_thresh:.2f}; marker size corresponds to {title_size_label}"

        if self.reduced:
            plot_title += f"</sub><sup><br>GO terms reduced with rrvgo (threshold=0.5)"

        # Create the scatter plot
        fig = px.scatter(df_plot,
                        x=self.score,
                        y="description",
                        size="size",
                        color="fdr",
                        hover_data={self.score: True, 'fdr': True, hover_size_label: True, 'size': False, 'Genes': True},
                        title=plot_title,
                        labels={"description": "Description", "fdr": "FDR q-value", hover_size_label: size_description},
                        color_continuous_scale="OrRd_r",
                        )

        # Add marker borders
        fig.update_traces(marker=dict(line=dict(width=0.5, color='black')))

        # Adjust colorbar size and position
        fig.update_layout(
            xaxis_title=score_description,
            yaxis_title="Gene Set Term",
            font=dict(size=12),
            title_font_size=13,  # Make the main title a bit smaller
            coloraxis_colorbar=dict(
                title="FDR q-value",  # Title for the color bar
                thickness=15,  # Adjust thickness of the color bar
                len=0.5,  # Adjust length of the color bar (in proportion to the plot height)
                x=1,  # Shift the color bar horizontally (to move it away from the marker legend)
                y=0.5,  # Adjust vertical positioning (centered in the plot)
                yanchor="middle"  # Set y position anchor in the middle of the plot
            ),
            legend_title_text=size_description,  # Optional: customize the legend title
        )

        # Increase the size of the y-axis tick labels and wrap the descriptions
        fig.update_layout(
            template='plotly_white',
            title_font_color='black',
            xaxis=dict(
                titlefont_size=12,
                color='black'
            ),
            yaxis=dict(
                titlefont_size=12,
                color='black',
                tickfont=dict(size=11),  # Increase y-axis tick label size
            ),
            legend=dict(
                font=dict(
                    color='black'
                )
            ),
        )

        # Show the plot
        #fig.show()

        # Save the plot to an HTML file
        dot_fig_filestem = Path(self.data_file).stem
        dot_fig_file = Path(self.op_dir) / f'{dot_fig_filestem}.html'
        fig.write_html(dot_fig_file)

    def run(self):
        self.load_data()
        self.plot_GO()