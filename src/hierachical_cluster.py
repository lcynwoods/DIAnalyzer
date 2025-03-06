

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.cluster.hierarchy import linkage
import scipy.spatial.distance as ssd
import matplotlib.patches as mpatches

class HierarchicalCluster():
    def __init__(self, data_df, op_file, col_suffix="log_transformed", conditions_file=None, colour_map=None, prot_ids=None):
        self.data_df = data_df
        self.op_file = op_file
        self.col_suffix = col_suffix
        self.conditions_file = conditions_file
        self.colour_map = colour_map
        self.prot_ids = prot_ids


        self.id_col = 'Genes'
        print("HELLO")

        print(self.data_df)
    
    def load_and_prepare_data(self):
        # Load data and conditions
        if isinstance(self.data_df, str):
            if '.xlsx' in self.data_df:
                incoming_data = pd.read_excel(self.data_df).set_index(self.id_col)
            elif '.csv' in self.data_df:
                incoming_data = pd.read_csv(self.data_df).set_index(self.id_col)
            elif ('.tsv' in self.data_df or '.txt' in self.data_df):
                incoming_data = pd.read_csv(self.data_df, sep='\t').set_index(self.id_col) 
        elif isinstance(self.data_df, pd.DataFrame):
            incoming_data = self.data_df.set_index(self.id_col)
        else:
            return
        # Filter to proteins of interest, if required
        if self.prot_ids:
            incoming_data = incoming_data[incoming_data.index.isin(self.prot_ids)]

        # Select relevant columns based on the column suffix
        self.cols_of_interest = [col for col in incoming_data.columns if self.col_suffix in col]
        self.data_df = incoming_data[self.cols_of_interest]
        new_cols = [col.replace(f'_{self.col_suffix}', '').replace('Intensity_', '') for col in self.data_df.columns]
        self.data_df.columns = new_cols
        print("Initial Data Shape:", self.data_df.shape)

        self.data_df = self.data_df.apply(pd.to_numeric, errors='coerce')        
        self.data_df = self.data_df.apply(lambda x: (x - np.mean(x)) / np.std(x), axis=1)
        
        # Replace non-finite values with zero (or another value if appropriate)
        self.data_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        self.data_df.fillna(0, inplace=True)

        # Load and map conditions from the conditions file
        self.conditions_df = pd.read_excel(self.conditions_file)

        # Create a dictionary mapping Sample.Name to Condition
        condition_map = dict(zip(self.conditions_df['Sample.Name'], self.conditions_df['Condition']))

        # Define a color palette and map conditions to colors
        condition_palette = sns.color_palette("Set2", len(self.conditions_df['Condition'].unique()))
        condition_color_map = dict(zip(self.conditions_df['Condition'].unique(), condition_palette))

        # Map the column conditions to their respective colors
        self.col_colors = self.data_df.columns.map(lambda x: condition_color_map[condition_map[x]])

    def perform_cluster(self):

        self.data_df = self.data_df.apply(pd.to_numeric, errors='coerce')        
        self.data_df = self.data_df.apply(lambda x: (x - np.mean(x)) / np.std(x), axis=1)
        # Replace non-finite values with zero (or another value if appropriate)
        self.data_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        self.data_df.dropna(inplace=True)

        untransformed_data = self.data_df  # Revert log2 transformation?
        
        # Ensure numeric data types
        untransformed_data = untransformed_data.apply(pd.to_numeric, errors='coerce')
        
        # Optionally, standardize data if necessary
        untransformed_data = untransformed_data.apply(lambda x: (x - np.mean(x)) / np.std(x), axis=1)
        
        # Replace non-finite values with zeros
        untransformed_data.replace([np.inf, -np.inf], np.nan, inplace=True)
        untransformed_data.fillna(0, inplace=True)

        # Calculate correlation matrix
        corr_matrix = untransformed_data.T.corr(method='pearson')
        

        corr_matrix = np.nan_to_num(corr_matrix)
        
        # Convert correlation to distance matrix
        distance_matrix = np.maximum(0, 1 - np.abs(corr_matrix))
        
        # Ensure the diagonal is zero
        np.fill_diagonal(distance_matrix, 0)
        
        # Ensure symmetry
        distance_matrix = (distance_matrix + distance_matrix.T) / 2
        
        # Convert to condensed distance matrix
        condensed_distance = ssd.squareform(distance_matrix)
        
        # Perform linkage
        self.linkage_matrix = linkage(condensed_distance, method='ward')
        print(self.linkage_matrix)

    def plot_seaborn_clustermap(self):

        transformed_data = np.sign(self.data_df) * np.abs(self.data_df) ** (1/2)

        # Load conditions and create column colors
        if self.conditions_file:
            if ".csv" in self.conditions_file:
                conditions_df = pd.read_csv(self.conditions_file)
                condition_map = dict(zip(conditions_df['Sample.Name'], conditions_df['Condition']))
            elif ".xlsx" in self.conditions_file:
                conditions_df = pd.read_excel(self.conditions_file)
                condition_map = dict(zip(conditions_df['Sample.Name'], conditions_df['Condition']))
            
            if not self.colour_map:
                # Create a default color palette for conditions
                condition_palette = sns.color_palette("Set2", len(conditions_df['Condition'].unique()))
                condition_colour_map = dict(zip(conditions_df['Condition'].unique(), condition_palette))
            else:
                condition_colour_map = self.colour_map

            # Generate column colors based on conditions
            col_colors = [condition_colour_map[condition_map[col]] for col in self.data_df.columns]

        else:
            col_colors = None

        # Create clustermap with rows ordered by the linkage matrix
        g = sns.clustermap(
            data=transformed_data,
            row_linkage=self.linkage_matrix,
            col_cluster=True,
            cmap='coolwarm',
            figsize=(19, 12),
            cbar_pos=(0.02, 0.1, 0.01, 1),
            cbar_kws=dict(label='log2(Intensity) z-score', orientation='horizontal'),
            col_colors=col_colors,  # Add column colors to reflect the conditions,
            yticklabels=False
        )

        # Add condition legend
        if self.conditions_file:
            condition_legend_patches = [mpatches.Patch(color=condition_colour_map[cond], label=cond) for cond in condition_colour_map]
            #plt.legend(handles=condition_legend_patches, title='Conditions', loc='upper left', bbox_to_anchor=(1.2, 1), fontsize=8, frameon=False)

            g.ax_heatmap.legend(handles=condition_legend_patches, loc='upper left', bbox_to_anchor=(1, 1), fontsize=12, frameon=False)

        BOTTOM = 0.15
        RIGHT = 0.85

        def shifter(value_list):
            count_list = []
            for word in value_list:
                count = 0
                for char in word:
                    count += 1
                count_list.append(count)
            max_chars = max(count_list)
            shift = float(max_chars-20)/100
            if shift > 0:
                return shift
            else:
                return 0

        BOTTOM += shifter(conditions_df['Sample.Name'].values)
        RIGHT -= shifter(conditions_df['Condition'].values)

        print(f"Bottom: {BOTTOM}")
        print(f"Right: {RIGHT}")

        plt.subplots_adjust(bottom=BOTTOM, right=RIGHT)
        g.figure.subplots_adjust(top=0.9)  # Adjust the top to give space for the title
        g.ax_cbar.set_position([RIGHT+0.01, 0.825, 0.1, 0.04])

        # Adjust font size for the row and column labels
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=12)  # Column labels font size

        # Rotate column labels to avoid overlap
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right')

        if self.col_suffix == 'with_imputed':
            g.figure.suptitle("Hierarchical Clustering of Protein Intensities After Normalization and Imputation", fontsize=16, ha='left', x=0.05)
        if self.col_suffix == 'normed':
            g.figure.suptitle("Hierarchical Clustering of Protein Intensities After Normalization", fontsize=16, ha='left', x=0.05)
        else:
            g.figure.suptitle("Hierarchical Clustering of Protein Intensities Before Normalization", fontsize=16, ha='left', x=0.05)

        # Add a subtitle underneath the main title
        g.figure.text(
            0.05,  # Centered horizontally
            0.93,  # Adjusted vertical position
            "Dendrogram constructed with Ward's Method, Distance = 1 - Pearson Correlation.",
            ha='left',  # Align text to center
            fontsize=14
        )
        #plt.show()
        plt.savefig(self.op_file)


    def run(self):
        self.load_and_prepare_data()
        self.perform_cluster()
        self.plot_seaborn_clustermap()

if __name__ == "__main__":
    DATA_DF = r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G42_MiguelPrado\DAPARv1.34.6_DA_Sheet1_output.tsv"
    CONDITIONS_FILE = r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G42_MiguelPrado\metadata_for_DAPAR.xlsx"
    OP_FILE = "op_file.png"
    COL_SUFFIX = "with_imputed"

    my_cluster= HierarchicalCluster(
        data_df = DATA_DF,
        conditions_file=CONDITIONS_FILE,
        op_file = OP_FILE,
        col_suffix = COL_SUFFIX
    )
    my_cluster.run()


