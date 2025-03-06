from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import scipy.spatial.distance as ssd
import seaborn as sns
import pandas as pd
import matplotlib.patches as mpatches
from sklearn.preprocessing import PowerTransformer
from yellowbrick.cluster import SilhouetteVisualizer

from cluster_gene_lines import ClusterLinePlotter


class KMeansClustering:
    def __init__(self,
                 data,
                 conditions_file,
                 op_file=r"TESTKMEANS.xlsx",
                 op_dir = None,
                 table_op_dir = None,
                 intermediate_op_dir = None,
                 colour_map=None,
                 prot_ids=None,
                 col_suffix='with_imputed',
                 n_clusters=5,
                 linkage_method='ward',
                 id_col='Genes'):
        self.data = data
        self.op_file = op_file
        self.op_dir = op_dir
        if op_dir:
            self.op_dir=op_dir
            if not table_op_dir:
                self.table_op_dir = self.op_dir
            else:
                self.table_op_dir = table_op_dir
            if not intermediate_op_dir:
                self.intermediate_op_dir = self.op_dir
            else:
                self.intermediate_op_dir = intermediate_op_dir
        else:
            self.op_dir = Path(op_file).parent
            self.intermediate_op_dir = self.op_dir
            self.table_op_dir = self.op_dir

        self.conditions_file = conditions_file
        self.colour_map = colour_map
        self.prot_ids = prot_ids
        self.col_suffix = col_suffix
        self.n_clusters = n_clusters
        self.linkage_method = linkage_method
        self.id_col = id_col

        self.data_df = None
        self.conditions_df = None
        self.cluster_labels = None
        self.linkage_matrix = None

    def load_and_prepare_data(self):
        # Load data and conditions
        incoming_data = None
        if isinstance(self.data, str):
            if '.xlsx' in self.data:
                incoming_data = pd.read_excel(self.data).set_index(self.id_col)
            elif '.csv' in self.data:
                incoming_data = pd.read_csv(self.data).set_index(self.id_col)
            elif ('.tsv' in self.data or '.txt' in self.data):
                incoming_data = pd.read_csv(self.data, sep='\t').set_index(self.id_col) 
        elif isinstance(self.data, pd.DataFrame):
            incoming_data = self.data.set_index(self.id_col)
        else:
            return
        print(incoming_data.columns)
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

    def perform_clustering(self):
        # Hierarchical clustering to establish the structure
        
        # Calculate correlation matrix
        corr_matrix = self.data_df.T.corr(method='pearson')
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

        # Determine number of clusters using silhouette score
        max_score = -1
        best_n_clusters = 2
        for n_clusters in range(4, 10):  # Adjust range as needed
            kmeans = KMeans(n_clusters=n_clusters, n_init=1, random_state=42)
            cluster_labels = kmeans.fit_predict(self.data_df)
            #cluster_labels = fcluster(self.linkage_matrix, t=n_clusters, criterion='maxclust')
            score = silhouette_score(self.data_df, cluster_labels)
            print(f"Silhouette score for {n_clusters} clusters: {score}")
            if score > max_score:
                max_score = score
                best_n_clusters = n_clusters

        self.n_clusters = best_n_clusters
        print(f"Optimal number of clusters based on silhouette score: {self.n_clusters}")

    def perform_hkmeans(self):
        # Initial clusters from hierarchical clustering
        initial_cluster_labels = fcluster(self.linkage_matrix, t=self.n_clusters, criterion='maxclust')
        print(initial_cluster_labels)
        print(f"The linkage matrix: {self.linkage_matrix}")
        data_cols = self.data_df.columns
        self.data_df['InitialCluster'] = initial_cluster_labels
        self.data_df.sort_values(by="InitialCluster", inplace=True)
        print(self.data_df)

        # Calculate initial centroids
        initial_centroids = []
        for i in range(1, self.n_clusters + 1):
            cluster_data = self.data_df[self.data_df['InitialCluster'] == i][data_cols]
            print(cluster_data)
            initial_centroids.append(cluster_data.mean(axis=0))
        initial_centroids = np.array(initial_centroids)
        print("Got here!")

        # K-means clustering using the initial centroids
        kmeans = KMeans(n_clusters=self.n_clusters, init=initial_centroids, n_init=1, random_state=42)
        self.cluster_labels = kmeans.fit_predict(self.data_df[data_cols])
        self.data_df['Cluster'] = self.cluster_labels + 1

    def plot_silhouette(self):
        # Fit KMeans model to the data
        model = KMeans(n_clusters=self.n_clusters, random_state=42)
        model.fit(self.data_df)

        # Create a new figure explicitly for the silhouette plot
        fig, ax = plt.subplots(figsize=(8, 6))  # Adjust size as needed
        visualizer = SilhouetteVisualizer(model, colors='yellowbrick', ax=ax)
        
        # Fit the visualizer
        visualizer.fit(self.data_df)
        
        # Finalize the figure (this is required if we don't use show())
        visualizer.finalize()

        # Save to file
        silhouette_stem = Path(self.op_file).stem
        silhouette_path = Path(self.op_dir) / f'{silhouette_stem}_silhouette.png'
        fig.savefig(silhouette_path, bbox_inches='tight')  # Use fig.savefig to save the figure
        plt.close(fig)  # Close the figure explicitly


    def plot_seaborn_clustermap(self):
        # Get the dendrogram order from the linkage matrix
        corr_matrix = self.data_df.T.corr(method='pearson')
        distance_matrix = 1 - corr_matrix
        distance_matrix = (distance_matrix + distance_matrix.T) / 2
        condensed_distance = ssd.squareform(distance_matrix)
        new_linkage_matrix = linkage(condensed_distance, method=self.linkage_method)
        
        # Define colors for each row cluster
        unique_clusters = np.unique(self.cluster_labels)
        cluster_colors = sns.color_palette("tab10", len(unique_clusters))
        row_colors = [cluster_colors[label] for label in self.cluster_labels]

        plot_data = self.data_df.drop(columns=['Cluster', 'InitialCluster'])

        # Define a color palette and map conditions to colors
        condition_map = dict(zip(self.conditions_df['Sample.Name'], self.conditions_df['Condition']))
        if not self.colour_map:
            # Map the column conditions to their respective colors based on the conditions file
            
            condition_palette = sns.color_palette("Set2", len(self.conditions_df['Condition'].unique()))
            condition_colour_map = dict(zip(self.conditions_df['Condition'].unique(), condition_palette))
        else:
            condition_colour_map = self.colour_map

        # Create column colors in the same order as the columns in self.data_df
        self.col_colors = [condition_colour_map[condition_map[col]] for col in plot_data.columns]

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

        BOTTOM += shifter(self.conditions_df['Sample.Name'].values)
        RIGHT -= shifter(self.conditions_df['Condition'].values)

        print(f"Bottom: {BOTTOM}")
        print(f"Right: {RIGHT}")

        # Create clustermap with rows ordered by the linkage matrix
        g = sns.clustermap(
            data=plot_data,
            row_linkage=new_linkage_matrix,
            col_cluster=False,
            cmap='coolwarm',
            figsize=(19, 12),
            row_colors=row_colors,
            col_colors=self.col_colors,  # Column colors based on conditions
            cbar_pos=([RIGHT+0.01, 0.825, 0.1, 0.04]),
            cbar_kws=dict(label='log2(Intensity) z-score', orientation='horizontal'),
            vmin=-2.5,
            vmax=2.5,
            yticklabels=False
        )

        g.figure.subplots_adjust(top=1)  # Adjust the top to give space for the title
        

        # Adjust the heatmap labels
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=10, rotation=45, ha='right')
        g.ax_heatmap.set_ylabel("")

        # Creating the cluster legend
        cluster_counts = {cluster: sum(self.cluster_labels == cluster) for cluster in unique_clusters}
        cluster_patches = [mpatches.Patch(color=cluster_colors[i], label=f'Cluster {i+1} (n={cluster_counts[i]})') for i in range(len(unique_clusters))]


        # Creating the condition legend
        condition_legend_patches = [mpatches.Patch(color=condition_colour_map[cond], label=cond) for cond in condition_colour_map]

        blank_patch = mpatches.Patch(color='white', label='')

        # Combine both legends (clusters + blank patch + conditions)
        legend_patches = cluster_patches + [blank_patch] + condition_legend_patches

        plt.subplots_adjust(bottom=BOTTOM, right=RIGHT)
        g.ax_cbar.set_position([RIGHT+0.01, 0.925, 0.1, 0.04])

        # Add combined legend for clusters and conditions
        g.ax_heatmap.legend(handles=legend_patches, loc='upper left', bbox_to_anchor=(1, 1), fontsize=12, frameon=False)

        # Show the plot
        #plt.show()

        #col_order = g.dendrogram_col.reordered_ind
        #print("Column order after clustering:", col_order)

        # Save the reordered columns for later use
        self.ordered_columns = plot_data.columns

        g.figure.suptitle(x=0.05, y=0.98, ha='left', t=f"Hierarchical KMeans Clustering of Differentially Abundant Protein Intensities (n={int(len(plot_data))})", fontsize=16)
        g.figure.text(
            0.05, 0.93,
            "Initial clustering with Ward's Method, distance = 1 - Pearson Correlation; KMeans cluster number (4-10) determined by best silhouette score",
            ha='left',  # Align text to center
            fontsize=12
        )

        hcluster_stem = Path(self.op_file).stem
        hcluster_path = Path(self.op_dir) / f'{hcluster_stem}.png'
        plt.savefig(hcluster_path)
        plt.close()
        #plt.show()

    def write_excel_output(self):
        # Write the ordered columns to Excel
        table_file = Path(self.table_op_dir) / f'{Path(self.op_file).stem}.xlsx'
        reordered_data_df = self.data_df[self.ordered_columns]
        reordered_data_df['Cluster'] = self.data_df['Cluster']
        reordered_data_df.to_excel(table_file, index=True)

    def make_whisker_plot(self):
        hcluster_stem = Path(self.op_file).stem
        hcluster_path = Path(self.op_dir) / f'{hcluster_stem}.png'
        table_file = Path(self.table_op_dir) / f'{Path(self.op_file).stem}.xlsx'
        cluster_line_plot = ClusterLinePlotter(table_file, hcluster_path, self.conditions_file, colour_map=self.colour_map)
        cluster_line_plot.run()

    def run(self):
        self.load_and_prepare_data()
        if self.data_df.empty:
            print("Data is empty after loading and preparing. Exiting.")
            return
        self.perform_clustering()  # Perform clustering
        self.perform_hkmeans()
        self.plot_seaborn_clustermap()
        self.plot_silhouette()  # Plot silhouette using Yellowbrick
        self.write_excel_output()
        self.make_whisker_plot()

if __name__ == "__main__":
    np.random.seed(42)

    # Generate some random data
    gene_ids = [f"Gene_{i}" for i in range(1, 101)]
    intensity_samples = [f"Intensity_sample_{i}_with_imputed" for i in range(1, 11)]
    data = np.random.normal(size=(100, 10))

    # Create a DataFrame to simulate gene expression data
    demo_data = pd.DataFrame(data, columns=intensity_samples)
    demo_data[0, 0] = 1000  # Extreme large value
    demo_data[1, 1] = 0  # Extreme small value
    demo_data[2, 2] = -500  # Another large value
    demo_data[3, 3] = -10  # Another small value
    demo_data['Genes'] = gene_ids
    print(demo_data)

    # Optionally save to file or pass directly to the KMeansClustering class
    demo_data.to_csv('demo_data.csv',index=True)

    DATA = r'demo_data.csv'

    samples = [f"sample_{i}" for i in range(1, 11)]

    # Corresponding conditions (you can adjust these based on your specific needs)
    conditions = ['Control', 'Treatment', 'Control', 'Treatment', 'Control',
                'Treatment', 'Control', 'Treatment', 'Control', 'Treatment']

    # Additional metadata fields (if needed, can add more later)
    time_points = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1]
    replicates = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5]

    # Create a DataFrame for the conditions
    conditions_df = pd.DataFrame({
        'Sample.Name': samples,  # Matches the column names in the data (before stripping prefix/suffix)
        'Condition': conditions,
        'TimePoint': time_points,  # You can remove this or replace with relevant metadata
        'Replicate': replicates
    })

    # Save to CSV or Excel for use
    conditions_df.to_excel('demo_conditions.xlsx', index=False)

    # Output path (for your reference)

    CONDITIONS_FILE = r'demo_conditions.xlsx'


    clustering = KMeansClustering(
        data=DATA,
        conditions_file=CONDITIONS_FILE,
        #prot_ids = ["Spen", "S100g", "Rbp1", "Timp3", "Pilra", "Sdcbp2", "Serpine1", "Irs2", "Chil4", "Pcsk9", "Acta1", "Txndc11", "Fam120c", "Tln2", "Dock6", "Aspa", "Rcbtb2", "Slc38a1", "Cox19", "Entpd6", "Adgrg2", "Uimc1", "Ripk1", "Smc1b", "Tns2", "Fscn1", "Sh2b1", "Tmlhe", "Bcl9l", "Phka2", "Ibtk", "Maob", "Itgb6", "Mgst1", "Bard1", "Nt5dc2", "Fkbp11", "Brca1", "Igfbp3", "Haspin", "Hebp1", "Acbd6", "Pol", "Ccnb2", "Psca", "Bap1", "Fv4", "Cdh16", "Pfkfb2", "Oga", "Prom1", "Cdh17", "Lgals4", "St3gal4", "Pik3ap1", "Pwwp2a", "Col18a1", "Sprr1a", "Ivl", "Lyn", "Peg3", "Pax2", "Stat1", "Myo7b", "Osbpl1a", "Eif2s3y", "Mettl7b", "Aldob", "Mpp5", "Igf2bp3", "Ddah1", "Stxbp4", "Bicc1", "Misp3", "Acad10", "Lrp2", "Palm3", "Abcc3", "Cdhr2", "Faah", "Aldh1a7", "Agr2", "Igf2bp1", "Ddc", "Mbp", "Sparc", "Calml4", "Mpv17", "Aldh1a1", "Hmgcs2", "Pkp1", "Dab2", "Cobl", "Prkg2", "Fkbp10", "Cyp4b1", "Soga3", "Rrm2b", "Cpm", "Bcat1", "Skp2"],
        col_suffix='with_imputed',  # Column suffix to select relevant data
        n_clusters=4,  # Initial number of clusters (can be updated based on silhouette score)
        linkage_method='ward',  # Use 'ward' for faster, optimized linkage
    )
    clustering.run()