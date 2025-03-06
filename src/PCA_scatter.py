import pandas as pd
import seaborn as sns
import plotly.express as px
from sklearn.decomposition import PCA

class PCAComparison:
    def __init__(self, data_file, op_file, col_suffix="log_transformed", conditions_file=None, colour_map=None, title="PCA Projection", id_col='Genes'):
        self.data_file = data_file
        self.op_file = op_file
        self.col_suffix = col_suffix  
        self.conditions_file = conditions_file  
        self.colour_map = colour_map  
        self.id_col = id_col
        self.title = title

        self.data_df = None
        self.filtered_data = None
        self.conditions_df = None
        self.condition_color_map = None
        self.samples = None
        self.sample_conditions = None
        self.pca = None  # Store PCA results

    def load_data(self):
        if isinstance(self.data_file, str):
            if '.xlsx' in self.data_file:
                self.data_df = pd.read_excel(self.data_file).set_index(self.id_col)
            elif '.csv' in self.data_file:
                self.data_df = pd.read_csv(self.data_file).set_index(self.id_col)
            elif '.tsv' in self.data_file or '.txt' in self.data_file:
                self.data_df = pd.read_csv(self.data_file, sep='\t').set_index(self.id_col)
        elif isinstance(self.data_file, pd.DataFrame):
            self.data_df = self.data_file.set_index(self.id_col)

    def filter_data(self):
        relevant_columns = [col for col in self.data_df.columns if (col.startswith("Intensity_") and col.endswith(self.col_suffix))]
        self.filtered_data = self.data_df[relevant_columns].dropna()
        
        if self.conditions_file:
            if '.xlsx' in self.conditions_file:
                self.conditions_df = pd.read_excel(self.conditions_file)
            elif '.csv' in self.conditions_file:
                self.conditions_df = pd.read_csv(self.conditions_file)
            condition_map = dict(zip(self.conditions_df['Sample.Name'], self.conditions_df['Condition']))

            self.samples = [col.replace("Intensity_", "").replace(f"_{self.col_suffix}", "") for col in self.filtered_data.columns]
            self.sample_conditions = [condition_map.get(sample, "Unknown") for sample in self.samples]

            if self.colour_map:
                self.condition_color_map = self.colour_map
            else:
                unique_conditions = self.conditions_df['Condition'].unique()
                palette = sns.color_palette("Set2", len(unique_conditions))
                self.condition_color_map = dict(zip(unique_conditions, [f'rgb({r*255},{g*255},{b*255})' for r, g, b in palette]))

    def perform_pca(self):
        data_transposed = self.filtered_data.T
        pca_model = PCA(n_components=2)
        self.pca = pca_model.fit_transform(data_transposed)

    def plot_pca(self):
        self.perform_pca()
        
        pca_df = pd.DataFrame(self.pca, columns=['PC1', 'PC2'])
        pca_df['Condition'] = self.sample_conditions
        pca_df['Sample'] = self.samples

        color_discrete_map = {condition: self.condition_color_map[condition] for condition in pca_df['Condition'].unique()}

        fig = px.scatter(pca_df, x='PC1', y='PC2', color='Condition', title=self.title,
                         labels={'PC1': 'Principal Component 1', 'PC2': 'Principal Component 2'},
                         hover_data={'Sample': True, 'Condition': True},
                         template='plotly_white', color_discrete_map=color_discrete_map)

        fig.update_traces(marker={'size': 12})
        fig.update_layout(font=dict(color='black'))
        fig.write_html(self.op_file)

    def run(self):
        self.load_data()
        self.filter_data()
        self.plot_pca()

