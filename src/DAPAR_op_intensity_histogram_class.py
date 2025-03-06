import pandas as pd
import plotly.figure_factory as ff
import plotly.graph_objs as go
import os
import numpy as np

class PlotNormCurve:
    def __init__(self, op_folder, prot_file, metadata_file, condition_colour_map, normalization_method="None", require_log_transform=True):
        self.op_folder = op_folder
        self.prot_file = prot_file
        self.metadata_file = metadata_file
        self.condition_colour_map = condition_colour_map
        self.normalization_method = normalization_method
        self.require_log_transform = require_log_transform

        self.data = None
        self.metadata = None
        self.intensity_columns = None
        self.intensity_prefix = None
        self.intensity_suffix = None
        self.fig1 = None
        
    def load_data(self):
        # Load the data into a Pandas DataFrame
        if '.tsv' in self.prot_file:
            self.data = pd.read_csv(self.prot_file, sep='\t', low_memory=False)
        elif '.txt' in self.prot_file:
            self.data = pd.read_csv(self.prot_file, sep='\t', low_memory=False)
        elif '.csv' in self.prot_file:
            self.data = pd.read_csv(self.prot_file, low_memory=False)
        elif '.xlsx' in self.prot_file:
             self.data = pd.read_excel(self.prot_file)
        self.metadata = pd.read_excel(self.metadata_file, sheet_name="Sheet1")
    
    def extract_data(self):
        # Set intensity column prefixes and suffixes based on normalization method
        if self.normalization_method == "LOESS":
            intensity_prefix = 'Intensity_'
            intensity_suffix = '_normed'
        elif self.normalization_method == "None":
            intensity_prefix = 'Intensity_'
            intensity_suffix = ''
        else:
            intensity_prefix = 'Intensity_'
            intensity_suffix = '_normed'

        self.intensity_columns = [f'{intensity_prefix}{sample}{intensity_suffix}' for sample in self.metadata['Sample.Name']]
        self.intensity_prefix = intensity_prefix
        self.intensity_suffix = intensity_suffix

        # Apply log transformation if needed
        if self.require_log_transform:
            for col in self.intensity_columns:
                self.data[col] = self.data[col].apply(lambda x: np.log2(x, out=np.zeros_like(x), where=(x != 0)))

    def plot_intensity(self):
        # Title string for the plot
        if self.normalization_method == 'DIANN_maxLFQ':
            title_string = "Intensity distribution by run<br><sup>Values showing are normalized with maxLFQ within diaNN.</sup>"
        elif self.normalization_method == "None":
            title_string = "Intensity distribution by run"
        elif self.normalization_method == 'LOESS':
            title_string = "Intensity distribution by run<br><sup>Following LOESS normalization.</sup>"
        elif self.normalization_method:
            title_string = f"Intensity distribution by run<br><sup>Following {self.normalization_method} normalization.</sup>"
        else:
            print("No normalization selected")

        plot_data = []
        sample_labels = [col.replace(self.intensity_prefix, '').replace(self.intensity_suffix, '') for col in self.intensity_columns]

        # Prepare the data for the plot
        for col in self.intensity_columns:
            sample_data = self.data[[col, 'Protein_Group','Genes']].dropna(subset=[col])  # Assuming 'Protein' column has the protein names
            plot_data.append(sample_data[col].values)

        sample_condition_map = pd.Series(self.metadata['Condition'].values, index=self.metadata['Sample.Name']).to_dict()

        # Assign colors to each sample based on their condition
        colors = [self.condition_colour_map[sample_condition_map.get(label)] for label in sample_labels]

        # Create the distplot (disable hover by setting hoverinfo to 'skip')
        self.fig1 = ff.create_distplot(plot_data, sample_labels, show_hist=False, show_rug=False, colors=colors)

        # Disable hover for the markers and lines
        for trace in self.fig1['data']:
            trace['hoverinfo'] = 'skip'  # Disable hover
            trace['showlegend'] = False  # Hide legend for sample data

        # Create a dummy legend for the conditions
        condition_labels = list(self.condition_colour_map.keys())
        condition_colors = [self.condition_colour_map[condition] for condition in condition_labels]

        # Create traces for the dummy legend
        legend_traces = []
        for condition, color in zip(condition_labels, condition_colors):
            legend_traces.append(go.Scatter(
                x=[None], y=[None],
                mode='markers',
                marker=dict(size=10, color=color),
                name=condition
            ))

        # Add dummy legend traces
        for trace in legend_traces:
            self.fig1.add_trace(trace)

        # Update layout for better visual presentation
        self.fig1.update_layout(
            xaxis_title="log2(Intensity)",
            yaxis_title="Density",
            title=title_string,
            template='plotly_white',
            font=dict(color='black'),
            legend_title_text="Condition"
        )

    def save_plot(self, fig):
        # Save as HTML and PNG
        config = {
            'toImageButtonOptions': {
                'format': 'png',  # one of png, svg, jpeg, webp
                'filename': 'Normed_intensity_plot',
                'scale': 2  # Multiply title/legend/axis/canvas sizes by this factor
            }
        }
        os.makedirs(self.op_folder, exist_ok=True)
        if self.normalization_method:
            fig.write_html(os.path.join(self.op_folder, f'log2FC_intensity_dist_with_{self.normalization_method}_normalization.html'), config=config)  # Save as HTML
        else:
            fig.write_html(os.path.join(self.op_folder, f'log2FC_intensity_dist_with_no_normalization.html'), config=config)

    def run(self):
        self.load_data()
        self.extract_data()
        self.plot_intensity()
        self.save_plot(self.fig1)


# Example usage:
op_folder = r"C:\Users\lwoods\Documents\LW_Projects_folder\test_DIA_summarizer"
prot_file = r"C:\Users\lwoods\Documents\LW_Projects_folder\test_DIA_summarizer\DAPAR_Protein_None_diff-abundance_DAPARv1.30.6.tsv"
metadata_file = r"C:\Users\lwoods\Documents\LW_Projects_folder\test_DIA_summarizer\metadata_for_DAPAR.xlsx"

if __name__ == "__main__":
    condition_colour_map = {
        "A": "red",
        "B": "blue"
    }

    intensity_plotter = PlotNormCurve(op_folder, prot_file, metadata_file, condition_colour_map, require_log_transform=True)
    intensity_plotter.run()