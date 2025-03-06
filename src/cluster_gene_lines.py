from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colorbar import Colorbar
import matplotlib.cm as cm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches



class ClusterLinePlotter():
    def __init__(self, df, op_file, condition_file, colour_map=None):
        self.df = str(df)
        self.op_file = str(op_file)
        self.condition_file = str(condition_file)
        self.colour_map = colour_map
        or_rd_r_colors = [
            "#fff7ec", "#fee8c8", "#fdd49e", "#fdbb84",
            "#fc8d59", "#ef6548", "#d7301f", "#b30000", "#7f0000"
        ]
        self.or_rd_r_cmap = mcolors.ListedColormap(or_rd_r_colors)
    
    def load_data(self):
        if isinstance(self.df, pd.DataFrame):
            return self.df
        elif isinstance(self.df, str):
            file_path = Path(self.df)
            if file_path.exists():
                if file_path.suffix in ('.tsv', '.txt'):
                    self.df = pd.read_csv(file_path, sep='\t')
                elif file_path.suffix == '.csv':
                    self.df = pd.read_csv(file_path)
                elif file_path.suffix == '.xlsx':
                    self.df = pd.read_excel(file_path)
                else:
                    print(f"Can't read input file type: {file_path.suffix}")
                    return None
            else:
                print(f"{file_path} doesn't exist!")
                return None
        else:
            print("No valid data could be loaded")
            return None
        return self.df

    def load_conditions(self):
        """Loads the condition file and creates a mapping of samples to conditions."""
        conditions_df = pd.read_excel(self.condition_file)
        condition_map = dict(zip(conditions_df['Sample.Name'], conditions_df['Condition']))
        unique_conditions = sorted(conditions_df['Condition'].unique())
        if self.colour_map:
        # Use the provided color map, ensuring all conditions have a color
            condition_colors = {condition: self.colour_map.get(condition, 'gray') for condition in unique_conditions}
        else:
            cmap = plt.get_cmap('tab10', len(unique_conditions))
            condition_colors = {condition: cmap(i) for i, condition in enumerate(unique_conditions)}
        return condition_map, condition_colors, unique_conditions

    def plot_boxplots(self, ax, cluster_df, samples, sample_means, norm_distances):
        """Plots the boxplots and the scatter plot of gene intensities."""
        ax.boxplot([cluster_df[col] for col in samples],
                   positions=range(len(samples)), widths=0.5,
                   patch_artist=True,
                   boxprops=dict(facecolor='none', color='black', linewidth=0.5),
                   flierprops=dict(markeredgecolor='none', linewidth=0.5),
                   medianprops=dict(color='black', linewidth=0.5),
                   whiskerprops=dict(color='black', linewidth=0.5),
                   capprops=dict(color='black', linewidth=0.5))

        for j, (_, row_data) in enumerate(cluster_df.iterrows()):
            gene_intensity = [row_data[col] for col in samples]
            colors = self.or_rd_r_cmap(norm_distances[j, :])
            ax.scatter(samples, gene_intensity, color=colors, marker='o', s=10, alpha=0.8)

        ax.plot(samples, sample_means, color='black', marker='o', linestyle='None', markersize=2)

    def plot_condition_means(self, ax, samples, sample_means, conditions, condition_colors, unique_conditions):
        """Plots thick bars representing the mean intensity of each condition with more precise alignment."""
        for condition in unique_conditions:
            # Calculate the mean and standard deviation for the samples in the current condition
            condition_means = [sample_means[i] for i, cond in enumerate(conditions) if cond == condition]
            condition_mean = np.mean(condition_means)
            condition_std = np.std(condition_means)
            
            # Get the indices of the samples belonging to the current condition
            sample_indices = [i for i, cond in enumerate(conditions) if cond == condition]
            
            if sample_indices:
                # Get the x-values (which are the sample indices) for the condition
                start_idx = min(sample_indices)
                end_idx = max(sample_indices)
                
                # Plot the condition mean as a horizontal line aligned with the boxplots
                ax.hlines(y=condition_mean, xmin=start_idx, xmax=end_idx, color=condition_colors[condition], linewidth=1)
                
                # Optionally, plot vertical lines at the start and end to better show the condition boundaries
                ax.vlines(x=start_idx, ymin=condition_mean - 0.05, ymax=condition_mean + 0.05, color=condition_colors[condition], linewidth=0.5)
                ax.vlines(x=end_idx, ymin=condition_mean - 0.05, ymax=condition_mean + 0.05, color=condition_colors[condition], linewidth=0.5)
                
                # Add a rectangle representing the mean ± standard deviation
                rect = patches.Rectangle(
                    (start_idx, condition_mean - condition_std),  # Bottom-left corner
                    end_idx - start_idx,  # Width (difference in x-values)
                    2 * condition_std,  # Height (± standard deviation)
                    linewidth=0, edgecolor=condition_colors[condition], facecolor=condition_colors[condition], alpha=0.3
                )
                ax.add_patch(rect)

    def plot_gene_trajectories(self):
        """Plots mass spec intensity data for all clusters in subplots."""
        clusters = sorted(self.df['Cluster'].unique())
        num_clusters = len(clusters)
        rows = int(np.ceil(np.sqrt(num_clusters)))
        cols = int(np.ceil(num_clusters / rows))

        fig, axes = plt.subplots(rows, cols, figsize=(15, 10), squeeze=True)
        fig.suptitle('Intensities (Z-score) Across Samples by Cluster', fontsize=12, ha='left', x=0.05)

        condition_map, condition_colors, unique_conditions = self.load_conditions()

        for i, cluster in enumerate(clusters):
            row = i // cols
            col = i % cols
            ax = axes[row, col]
            cluster_df = self.df[self.df['Cluster'] == cluster]
            intensity_cols = [col for col in self.df.columns if col not in ['Genes', 'InitialCluster', 'Cluster']]
            intensities = cluster_df[intensity_cols].values
            sample_means = np.mean(intensities, axis=0)
            distance_from_mean = np.abs(intensities - sample_means)
            norm_distances = (distance_from_mean - np.min(distance_from_mean)) / (np.max(distance_from_mean) - np.min(distance_from_mean))
            samples = intensity_cols
            conditions = [condition_map[sample] for sample in samples]

            self.plot_condition_means(ax, samples, sample_means, conditions, condition_colors, unique_conditions)
            self.plot_boxplots(ax, cluster_df, samples, sample_means, norm_distances)
            

            ax.set_title(f'Cluster {cluster}', fontsize=8)
            ax.set_ylabel('Intensity (Z-score)', fontsize=8)
            ax.set_xlabel('Sample', fontsize=8)
            ax.set_xticks(range(len(samples)))
            ax.set_xticklabels(samples, rotation=90, fontsize=6)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        sm = plt.cm.ScalarMappable(cmap=self.or_rd_r_cmap, norm=plt.Normalize(vmin=0, vmax=1))
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.01, pad=0.03)
        cbar.set_label('Normalized Distance from Sample Mean', fontsize=8)
        fig.subplots_adjust(top=0.88, bottom=0.3, left=0.07, right=0.87, hspace=0.3, wspace=0.2)

        # Add condition legend
        handles = [Line2D([0], [0], color=condition_colors[cond], lw=2) for cond in unique_conditions]
        fig.legend(handles=handles, labels=unique_conditions, loc='upper right', fontsize=8)
        plt.show()

    def plot_individual_gene_trajectories(self):
        """Plots mass spec intensity data for each cluster in individual plots."""
        clusters = sorted(self.df['Cluster'].unique())
        condition_map, condition_colors, unique_conditions = self.load_conditions()

        for cluster in clusters:
            fig, ax = plt.subplots(figsize=(19, 12))
            fig.suptitle(f'Intensities (Z-score) Across Samples for Cluster {cluster}', fontsize=16, ha='left', x=0.05)
            cluster_df = self.df[self.df['Cluster'] == cluster]
            intensity_cols = [col for col in self.df.columns if col not in ['Genes', 'InitialCluster', 'Cluster']]
            intensities = cluster_df[intensity_cols].values
            sample_means = np.mean(intensities, axis=0)
            distance_from_mean = np.abs(intensities - sample_means)
            norm_distances = (distance_from_mean - np.min(distance_from_mean)) / (np.max(distance_from_mean) - np.min(distance_from_mean))
            samples = intensity_cols
            conditions = [condition_map[sample] for sample in samples]

            self.plot_condition_means(ax, samples, sample_means, conditions, condition_colors, unique_conditions)
            self.plot_boxplots(ax, cluster_df, samples, sample_means, norm_distances)
            

            #ax.set_title(f'Cluster {cluster}', fontsize=8)
            ax.set_ylabel('Intensity (Z-score)', fontsize=16)
            ax.set_xlabel('Sample', fontsize=16)
            ax.set_xticks(range(len(samples)))
            ax.set_xticklabels(samples, rotation=90, fontsize=12)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

            sm = plt.cm.ScalarMappable(cmap=self.or_rd_r_cmap, norm=plt.Normalize(vmin=0, vmax=1))
            sm.set_array([])
            # Define the colorbar location manually using `add_axes()`
            cbar_ax = fig.add_axes([0.85, 0.6, 0.01, 0.2])  # [left, bottom, width, height] in figure coordinates

            # Create the colorbar at the new location
            cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')

            # Set the label for the colorbar
            cbar.set_label('Normalized Distance\nfrom Sample Mean', fontsize=12)

            # Add condition legend
            handles = [Line2D([0], [0], color=condition_colors[cond], lw=2) for cond in unique_conditions]
            labels = [f"Mean of {cond} (+/- STD)" for cond in unique_conditions]
            ax.legend(handles=handles, labels=labels, loc='lower left', fontsize=10, bbox_to_anchor=(1,0), frameon=False)
            plt.subplots_adjust(
                top=0.88,
                bottom=0.3,
                left=0.07,
                right=0.79,
                hspace=0.2,
                wspace=0.2)
            plt.savefig(str(self.op_file).replace('.png',f'_cluster_{cluster}.png'))

    def run(self):
        self.load_data()
        #self.plot_gene_trajectories()
        self.plot_individual_gene_trajectories()


if __name__ == "__main__":
    DF = r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\Analysis_Sheet1_cluster_assignment.xlsx"
    #DF = r"clustering.xlsx"
    CONDITION_FILE = r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\LOESS\metadata_for_DAPAR.xlsx"

    plotter = ClusterLinePlotter(df=DF, op_file="item.png", condition_file=CONDITION_FILE)
    plotter.run()

