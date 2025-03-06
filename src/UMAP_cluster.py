import pandas as pd
import seaborn as sns
from umap import UMAP
import plotly.express as px
from sklearn.metrics import silhouette_score
from sklearn.model_selection import ParameterGrid

class UMAPComparison:
    def __init__(self, data_file, op_file, col_suffix="log_transformed", conditions_file=None, colour_map=None, title=None, id_col='Genes', param_grid=None):
        self.data_file = data_file
        self.op_file = op_file
        self.col_suffix = col_suffix  # Suffix to select relevant columns (e.g., log_transformed, with_imputed)
        self.conditions_file = conditions_file  # Path to conditions file
        self.colour_map = colour_map  # Custom color map (optional)
        self.id_col = id_col
        self.title = title

        self.data_df = None
        self.filtered_data = None
        self.imputed_data_df = None
        self.conditions_df = None
        self.condition_color_map = None
        self.samples = None
        self.sample_conditions = None
        self.umap = None

        if not param_grid:
            self.param_grid = {
                'n_neighbors': [5, 15, 30],
                'min_dist': [0.01, 0.1, 0.3],
                'metric': ['euclidean', 'cosine'],
                'n_components': [2]
            }
        else:
            self.param_grid = param_grid

        self.best_params = None
        self.best_score = -1

    def load_data(self):
        # Load data and conditions
        incoming_data = None
        print(self.data_file)
        if isinstance(self.data_file, str):
            print("STRING")
            if '.xlsx' in self.data_file:
                incoming_data = pd.read_excel(self.data_file).set_index(self.id_col)
            elif '.csv' in self.data_file:
                print("CSV")
                incoming_data = pd.read_csv(self.data_file).set_index(self.id_col)
            elif ('.tsv' in self.data_file or '.txt' in self.data_file):
                incoming_data = pd.read_csv(self.data_file, sep='\t').set_index(self.id_col) 
        elif isinstance(self.data_file, pd.DataFrame):
            incoming_data = self.data_file.set_index(self.id_col)
        else:
            return
        print("HELLO")
        self.data_df = incoming_data

    def filter_data(self):

        # Select columns based on the suffix provided

        relevant_columns = [col for col in self.data_df.columns if (col.startswith("Intensity_") and col.endswith(self.col_suffix))]
        self.filtered_data = self.data_df[relevant_columns]
        self.filtered_data.dropna(inplace=True)
        print(f"Selected columns with suffix '{self.col_suffix}': {self.filtered_data.shape}")

        # Load conditions file and map sample names to conditions
        if self.conditions_file:
            if '.xlsx' in self.conditions_file:
                self.conditions_df = pd.read_excel(self.conditions_file)
            elif '.csv' in self.conditions_file:
                self.conditions_df = pd.read_csv(self.conditions_file)
            condition_map = dict(zip(self.conditions_df['Sample.Name'], self.conditions_df['Condition']))

            # Map the columns (samples) to their respective conditions
            self.samples = [col.replace(f"Intensity_", '').replace(f'_{self.col_suffix}', '') for col in self.filtered_data.columns]
            self.sample_conditions = [condition_map[col.replace(f"Intensity_", '').replace(f'_{self.col_suffix}', '')] for col in self.filtered_data.columns]

            # Create or use a custom color map
            if self.colour_map:
                self.condition_color_map = self.colour_map
            else:
                # Default color map if none is provided
                unique_conditions = self.conditions_df['Condition'].unique()
                palette = sns.color_palette("Set2", len(unique_conditions))
                self.condition_color_map = dict(zip(unique_conditions, [f'rgb({r*255},{g*255},{b*255})' for r, g, b in palette]))
        else:
            raise ValueError("No conditions file provided!")

    def perform_umap(self, data, n_neighbors, min_dist, metric, n_components, random_state=42):
        data_transposed = data.T
        umap_model = UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric, n_components=n_components, random_state=random_state)
        embedding = umap_model.fit_transform(data_transposed)
        return embedding

    def evaluate_clusters(self, embedding):
        score = silhouette_score(embedding, self.conditions_df['Condition'])
        return score

    def grid_search(self):
        grid = ParameterGrid(self.param_grid)
        for params in grid:
            print(f"Testing with params: {params}")
            embedding = self.perform_umap(data=self.filtered_data, **params)
            score = self.evaluate_clusters(embedding)
            print(f"Silhouette Score: {score}")
            
            if score > self.best_score:
                self.best_score = score
                self.best_params = params

        print(f"Best params: {self.best_params}, Best Silhouette Score: {self.best_score}")

    def plot_umap(self):
        self.umap = self.perform_umap(self.filtered_data, n_neighbors=self.best_params['n_neighbors'],
                                      min_dist=self.best_params['min_dist'],
                                      metric=self.best_params['metric'],
                                      n_components=self.best_params['n_components'])

        print(self.umap)
        # Create a DataFrame for Plotly
        umap_df = pd.DataFrame(self.umap, columns=['UMAP1', 'UMAP2'])
        umap_df['Condition'] = self.sample_conditions
        umap_df['Sample'] = self.samples

        # Use condition_color_map to assign colors
        color_discrete_map = {condition: self.condition_color_map[condition] for condition in umap_df['Condition'].unique()}

        if not self.title:
            self.title = "UMAP Projection of Samples"
        
        self.title = (
            f"{self.title}<sup><br>Parameters: N neighbours={self.best_params['n_neighbors']}; "
            f"Min distance={self.best_params['min_dist']}; "
            f"Metric={self.best_params['metric']}; "
            f"N components={self.best_params['n_components']}"
        )

        # Plot using Plotly
        fig = px.scatter(umap_df, x='UMAP1', y='UMAP2', color='Condition', title=self.title,
                         labels={'UMAP1': 'UMAP Dimension 1', 'UMAP2': 'UMAP Dimension 2'},
                         hover_data={'Sample':True,'Condition':True,'UMAP1':True, 'UMAP2':True},
                         template='plotly_white', color_discrete_map=color_discrete_map,)

        # Save and show the plot
        fig.update_traces(marker={'size': 15})

        fig.update_layout(font=dict(color='black'))

        fig.write_html(self.op_file)
        #fig.show()

    def run(self):
        self.load_data()
        self.filter_data()
        self.grid_search()
        self.plot_umap()

if __name__ == "__main__":

    import pandas as pd
    import numpy as np
    from sklearn.datasets import make_blobs
    import random

    n_samples = 10
    n_features = 1000
    centers = 2  # Two clusters: Control and Treatment

    # Create blobs (clusters)
    X, y = make_blobs(n_samples=n_samples, n_features=n_features, centers=centers, cluster_std=5, random_state=42)

    # Convert to a DataFrame where rows are samples and columns are proteins
    samples = [f"Intensity_sample_{i}" for i in range(1, n_samples + 1)]
    proteins = [f"Gene_{i}" for i in range(1, n_features + 1)]
    clustered_data = pd.DataFrame(X.T, columns=samples, index=proteins)  # Transpose to make samples as columns

    # Create log-transformed data
    log_transformed_data = np.log1p(np.abs(clustered_data))  # Log transformation

    # Introduce some missing values to simulate imputed data
    with_imputed_data = clustered_data.copy()
    for col in with_imputed_data.columns:
        missing_indices = random.sample(range(n_features), k=int(n_features * 0.1))  # 10% missing values
        with_imputed_data.iloc[missing_indices, clustered_data.columns.get_loc(col)] = np.nan

    # Save original, log-transformed, and imputed data as a single CSV
    combined_data = pd.concat([
        clustered_data.add_suffix('_log_transformed'),
        with_imputed_data.add_suffix('_with_imputed')
    ], axis=1)
    combined_data['Genes'] = proteins
    combined_data.set_index('Genes', inplace=True)
    combined_data.to_csv('clustered_demo_data.csv')
    print("Clustered data saved as 'clustered_demo_data.csv'")

    # Conditions file (Control and Treatment alternating for simplicity)
    conditions_df = pd.DataFrame({
        'Sample.Name': [f"sample_{i}" for i in range(1, n_samples + 1)],
        'Condition': ['Control' if i % 2 == 0 else 'Treatment' for i in range(1, n_samples + 1)]
    })
    conditions_df.to_csv('demo_conditions.csv', index=False)
    print("Conditions file saved as 'demo_conditions.csv'")

    # Create a corresponding conditions file
    samples = [f"sample_{i}" for i in range(1, 11)]  # Sample names without 'Intensity_' prefix
    conditions = ['Control', 'Treatment'] * 5  # Alternating Control and Treatment conditions

    # Create a DataFrame for the conditions
    conditions_df = pd.DataFrame({
        'Sample.Name': samples,
        'Condition': conditions
    })

    # Save the conditions file
    conditions_df.to_csv('demo_conditions.csv', index=False)
    print("Conditions data saved as 'demo_conditions.csv'")

    DATA_FILE = r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\LOESS\DAPARv1.30.6_DA_Sheet1_output.tsv"  # Path to the data file with various suffixes
    OUTPUT_FILE = 'umap.html'  # Path to save the plot
    COL_SUFFIX = 'with_imputed'  # Suffix for the columns you want to use (log_transformed, with_imputed, etc.)
    CONDITIONS_FILE = r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\LOESS\metadata_for_DAPAR.xlsx"  # Path to your conditions file
    TITLE = "UMAP Projection of Samples Intensities"

    umap_comp = UMAPComparison(data_file=DATA_FILE, op_file=OUTPUT_FILE, col_suffix=COL_SUFFIX, conditions_file=CONDITIONS_FILE, title=TITLE)
    umap_comp.run()

