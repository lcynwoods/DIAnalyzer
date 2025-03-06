
import textwrap 
import pandas as pd
import numpy as np
import re
from statsmodels.distributions.empirical_distribution import ECDF
from statsmodels.stats.multitest import multipletests
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as io
import statistics
from pathlib import Path

class GetDA:
        
    def __init__(self,
            bkgd_file,
            metadata_file,
            op_folder,
            comparison_columns,
            pval_cutoff,
            ppushprop_th,
            filterprop_th,
            condition_colour_map,
            table_op_folder = None,
            DA_directions=["up","down"],
            set_ids=None,
            bait_gene=None,
            mark_contams=True,
            Log2FC_cutoff_operation="None", Log2FC_cutoff=None,
            slope = 0.01,
            contaminants_fasta = r"C:\Users\lwoods\Documents\LW_Projects_folder\general\20210914_525_CONTAMINANTS_UNIPROT_FORMAT.fasta",
            extra_POIs=None,
            POI_colour_map=None,
            sheet=None,
            remove_sheet = None,
            remove_list = None,):
        self.bkgd_file = bkgd_file
        self.metadata_file = metadata_file
        self.op_folder = op_folder
        if table_op_folder:
            self.table_op_folder = table_op_folder
        else:
            self.table_op_folder = self.op_folder
        self.mark_contams = mark_contams
        self.DA_directions = DA_directions
        self.bait_gene = bait_gene
        self.POI_colour_map = POI_colour_map
        self.comparison_columns = comparison_columns
        self.set_ids = set_ids
        self.pval_cutoff = pval_cutoff
        self.ppushprop_th = ppushprop_th
        self.filterprop_th = filterprop_th
        self.Log2FC_cutoff_operation = Log2FC_cutoff_operation
        self.Log2FC_cutoff = Log2FC_cutoff
        self.slope = slope
        self.contaminants_fasta = contaminants_fasta
        if sheet:
            self.sheet = sheet
        else:
            self.sheet = "Analysis" 
        self.remove_sheet = remove_sheet
        self.remove_list = remove_list
        self.condition_colour_map =  condition_colour_map

        if extra_POIs:
            self.extra_POIs = [POI.upper() for POI in extra_POIs]
        else:
            self.extra_POIs = []
        self.ordered_dataframe = None
        self.contaminant_ids = None
        self.logFC_cutoff_dict = None

    def load_file_inputs(self):
        self.ordered_dataframe = pd.read_csv(self.bkgd_file, sep="\t")
        if self.contaminants_fasta:
            with open(self.contaminants_fasta, 'r') as f:
                lines = f.readlines()

            # Extract protein IDs (assuming uniprot ids)   
            self.contaminant_ids = [line.strip().lstrip('>').split("|")[1] for line in lines if line.startswith('>')]
           
    def make_op_folder(self):
        (Path(self.op_folder) / f'{self.sheet}').mkdir(parents=True, exist_ok=True)
        self.op_folder = Path(self.op_folder) / f'{self.sheet}'
        (Path(self.table_op_folder) / f'{self.sheet}').mkdir(parents=True, exist_ok=True)
        self.table_op_folder = Path(self.table_op_folder) / f'{self.sheet}'

    def mark_contaminants_dataframe(self):
        def mark_contaminants(row, contaminants_list):
            # Split the Protein_Group by ";"
            proteins = row['Protein_Group'].split(";")

            # Initialize is_contaminant as "No"
            row['is_contaminant'] = "No"
            # Update the Protein_Group and is_contaminant if any protein is a contaminant
            for protein in proteins:
                if contaminants_list and protein in contaminants_list:
                    row['Protein_Group'] = "CON__" + row['Protein_Group']
                    row['is_contaminant'] = "Yes"
                    break  # No need to check further if a contaminant is found
                elif "CON__" in protein:
                    row['is_contaminant'] = "Yes"
                    break  # No need to check further if a contaminant is found

            return row

        if self.mark_contams == True:
            self.ordered_dataframe = self.ordered_dataframe.apply(mark_contaminants, args=(self.contaminant_ids,), axis=1)

    def mark_removed_prots(self):

        def mark_removed_by_row(row, prot_list):
            proteins = row['Protein_Group'].split(";")
            # Initialize "Remove?" as "No"
            row['Removed?'] = "No"
            if len(prot_list) > 0:
                for protein in proteins:
                    if protein in prot_list:
                        row['Removed?'] = "Yes"
            return row

        remove_prots_list = []
        if self.remove_sheet:
            if Path(self.remove_sheet).suffix == ".tsv":
                with open(self.remove_sheet, 'r') as rm_handler:
                    remove_prots_list = [line.strip('\n') for line in rm_handler if line.strip('\n') != "Protein_Groups"]
            elif Path(self.remove_sheet).suffix == ".txt":
                with open(self.remove_sheet, 'r') as rm_handler:
                    remove_prots_list = [line.strip('\n') for line in rm_handler if line.strip('\n') != "Protein_Groups"]
            else:
                print("Please supply list of proteins to be removed as single column txt or tsv list.")
                print("Carrying on without removing proteins")
        if self.remove_list:
            remove_prots_list.extend(self.remove_list)
        print(remove_prots_list)
        self.ordered_dataframe = self.ordered_dataframe.apply(mark_removed_by_row, args=(remove_prots_list,), axis=1)

    def find_logFC_thresh(self):
        self.logFC_cutoff_dict = {}

        # Case 1: If Log2FC_cutoff is a list of numbers
        if isinstance(self.Log2FC_cutoff, list):
            for col, logFC in zip(self.comparison_columns, self.Log2FC_cutoff):
                self.logFC_cutoff_dict[col] = logFC

        # Case 2: If Log2FC_cutoff is a string and "auto", calculate the ECDF-derived logFC
        elif isinstance(self.Log2FC_cutoff, str) and self.Log2FC_cutoff.lower() == "auto":
            for col in self.comparison_columns:
                log2FC_col = f"{col}_logFC"
                positive_logFC_values = self.ordered_dataframe[log2FC_col].abs()
                ecdf_positive = ECDF(positive_logFC_values)

                x_ecdf_positive = np.linspace(min(positive_logFC_values), max(positive_logFC_values))
                y_ecdf_positive = ecdf_positive(x_ecdf_positive)

                delta_y_ecdf_positive = np.diff(y_ecdf_positive)

                # Find the steepest point in the ECDF curve
                steepest_point_idx = np.argmax(delta_y_ecdf_positive)

                # Find the change point after the steepest point
                change_point_idx_ecdf_positive = np.where(delta_y_ecdf_positive[steepest_point_idx:] <= self.slope)[0][0] + steepest_point_idx
                change_point_value_ecdf_positive = x_ecdf_positive[change_point_idx_ecdf_positive]

                # Plot the ECDF
                fig = go.Figure()
                fig.add_trace(go.Scatter(x=x_ecdf_positive, y=y_ecdf_positive, mode='lines', name='Positive ECDF', line=dict(color='black')))
                
                fig.add_annotation(
                    x=change_point_value_ecdf_positive,
                    y=y_ecdf_positive[change_point_idx_ecdf_positive],
                    text=f'Change Point: {change_point_value_ecdf_positive:.2f}',
                    showarrow=True,
                    arrowhead=4,
                    ax=-40,
                    ay=-40,
                    xanchor="right",
                    standoff=10,
                    font=dict(color="#FF5A00")
                )

                fig.add_shape(type="line", x0=change_point_value_ecdf_positive,
                            x1=change_point_value_ecdf_positive, y0=0, y1=1,
                            line=dict(color="#FF5A00", width=2, dash="dash"),
                            name=f'Change Point at {self.slope}: {change_point_value_ecdf_positive:.2f}',
                            )

                fig.update_layout(
                    xaxis_title=f'log2 fold change', 
                    yaxis_title=f'ECDF of positive values for {col}',
                    plot_bgcolor='#EEEAE8', 
                    title=f"ECDF for {col.replace('_', ' ')}<br><br><sup>Log2FC Threshold: Post-Maximum ECDF Rise ≤ {self.slope}.<sup>",
                    font=dict(color='black'),
                    title_font=dict(color='black'),
                    legend_title_font=dict(color='black')
                )

                config = {
                    'toImageButtonOptions': {
                        'format': 'png',  # one of png, svg, jpeg, webp
                        'filename': 'ECDF_plot',
                        'scale': 2  # Multiply title/legend/axis/canvas sizes by this factor
                    }
                }

                # Save the plot
                fig.write_html(f"{self.op_folder}/{col}_ECDF_plot.html", config=config)
                self.logFC_cutoff_dict[col] = change_point_value_ecdf_positive

        # Case 3: If Log2FC_cutoff is a number, apply it directly
        elif isinstance(self.Log2FC_cutoff, (float, int)):
            for col in self.comparison_columns:
                self.logFC_cutoff_dict[col] = self.Log2FC_cutoff

        # If no valid cutoff is provided, print an error message
        else:
            print("Must indicate logFC cutoff preference.")

    def find_logFC_thresh_old(self):
        if not isinstance(self.Log2FC_cutoff,float):
            self.logFC_cutoff_dict = {}
            if isinstance(self.Log2FC_cutoff, list):
                for col,logFC in zip(self.comparison_columns,self.Log2FC_cutoff):
                    self.logFC_cutoff_dict[col] = logFC
        elif isinstance(self.Log2FC_cutoff,float):
            for col,logFC in zip(self.comparison_columns,self.Log2FC_cutoff):
                self.logFC_cutoff_dict[col] = logFC
        elif isinstance((self.Log2FC_cutoff,str) and (str(self.Log2FC_cutoff).lower()=="auto")):
            self.logFC_cutoff_dict = {}
            for col in self.comparison_columns:
                log2FC_col = f"{col}_logFC"
                # Calculate the ECDF
                # Two options: symmetrical and priority
                # Priority
                # positive_logFC_values = self.ordered_dataframe[(self.ordered_dataframe[log2FC_col] > 0)][log2FC_col]
                # Symmetrical
                positive_logFC_values = self.ordered_dataframe[log2FC_col].abs()
                ecdf_positive = ECDF(positive_logFC_values)

                # Generate ECDF values for the positive inverted data
                x_ecdf_positive = np.linspace(min(positive_logFC_values), max(positive_logFC_values))
                y_ecdf_positive = ecdf_positive(x_ecdf_positive)

                # Calculate the change in ECDF values for the positive inverted data
                delta_y_ecdf_positive = np.diff(y_ecdf_positive)

                # Find the steepest point in the curve
                steepest_point_idx = np.argmax(delta_y_ecdf_positive)

                # Now look for change points only after this steepest point
                change_point_idx_ecdf_positive = np.where(delta_y_ecdf_positive[steepest_point_idx:] <= self.slope)[0][0] + steepest_point_idx
                change_point_value_ecdf_positive = x_ecdf_positive[change_point_idx_ecdf_positive]

                # Plot the ECDF for the positive inverted data with change points highlighted

                fig = go.Figure()

                fig.add_trace(go.Scatter(x=x_ecdf_positive, y=y_ecdf_positive, mode='lines', name='Positive ECDF', line=dict(color='black')))
                
                # Annotation for change point
                fig.add_annotation(
                    x=change_point_value_ecdf_positive,
                    y=y_ecdf_positive[change_point_idx_ecdf_positive],
                    text=f'Change Point: {change_point_value_ecdf_positive:.2f}',
                    showarrow=True,
                    arrowhead=4,
                    ax=-40,
                    ay=-40,
                    xanchor="right",   # Ensure the x-coordinate is the right edge of the text box
                    standoff=10,       # Offset the arrow tip by 10 pixels
                    font=dict(color="#FF5A00")
                )

                fig.add_shape(type="line", x0=change_point_value_ecdf_positive,
                            x1=change_point_value_ecdf_positive, y0=0, y1=1,
                            line=dict(color="#FF5A00", width=2, dash="dash"),
                            name=f'Change Point at {self.slope}: {change_point_value_ecdf_positive:.2f}',
                            )     
                
                fig.update_layout(
                            xaxis_title=f'log2 fold change', yaxis_title=f'ECDF of positive values for {col}',
                            template='plotly_white',
                            title=f"ECDF for {col.replace('_',' ')}<br><br><sup>Log2FC Threshold: Post-Maximum ECDF Rise ≤ {self.slope}.<sup>",
                            font=dict(color='black'),
                            title_font=dict(color='black'),
                            legend_title_font=dict(color='black')
                        )
                
                
                config = {
                  'toImageButtonOptions': {
                    'format': 'png', # one of png, svg, jpeg, webp
                    'filename': 'ECDF_plot',
                    'scale': 2 # Multiply title/legend/axis/canvas sizes by this factor
                  }
                }
                #fig.show(config=config)
                #fig.write_image(f"{self.op_folder}/{col}_ECDF_plot.png", width=800, height=400, scale=6)
                fig.write_html(f"{self.op_folder}/{col}_ECDF_plot.html",config=config)
                self.logFC_cutoff_dict[col] = change_point_value_ecdf_positive
        else:
            print("Must indicate logFC cutoff preference.")


    def calculate_alternative_FDR(self, p_val_list, p_val_thresh):
        # Adjust p-values using the BH method
        _, padj, _, _ = multipletests(p_val_list, alpha=p_val_thresh, method='fdr_bh')
        
        # Identify original p-values that are below the significance threshold
        significant_pvals_mask = np.array(p_val_list) <= p_val_thresh
        
        # From those, find the largest p-value that is considered significant
        significant_pvals = np.array(p_val_list)[significant_pvals_mask]
        if significant_pvals.size > 0:
            largest_significant_pval = np.max(significant_pvals)
            # Find the corresponding q-value for this largest significant p-value
            index_of_largest_significant_pval = np.where(np.array(p_val_list) == largest_significant_pval)[0][0]
            BH_fdr = padj[index_of_largest_significant_pval]
        else:
            BH_fdr = np.nan  # Handle case with no significant items

        return BH_fdr
    
            
    def filter_and_save_interactors(self):
        self.interactors_dict = {}  # Dictionary to store filtered interactors dataframes
        self.right_dict = {}
        self.left_dict = {}
        self.BH_FDR_dict = {}
        self.consolidated_df = self.ordered_dataframe.copy()
        columns_to_remove = [col for col in self.consolidated_df if "_logFC" in col or "_pval" in col]
        self.consolidated_df = self.consolidated_df.drop(columns=columns_to_remove)
        print("HERE")
        print(self.consolidated_df.columns)
        print(self.comparison_columns)
        for col in self.comparison_columns:
            right_cond = col.split('_vs_')[0]
            left_cond = col.split('_vs_')[1]
            
            Log2FC_cutoff = self.logFC_cutoff_dict[col]
            log2FC_col = f"{col}_logFC"
            pval_col = f"{col}_pval"

            compare_df_cols = ['Protein_Group', 'Genes']
            if self.mark_contams == True:
                compare_df_cols.append('is_contaminant')

            compare_dataframe = self.ordered_dataframe[compare_df_cols + [log2FC_col, pval_col, "Removed?"]]

            pvals = compare_dataframe[pval_col].values
            #reject, pvals_corrected, _, _ = multipletests(pvals, alpha=self.pval_cutoff, method='fdr_bh')

            # Merge FDR-corrected p-values back to the main dataframe
            # Currently, performing p-value correction in R
            compare_dataframe = compare_dataframe.assign(fdr_corrected_pvals=pvals)
            FDR_df = compare_dataframe[abs(compare_dataframe[log2FC_col]) >= Log2FC_cutoff]
            print(FDR_df)
            BH_FDR = self.calculate_alternative_FDR(FDR_df[pval_col].values.tolist(),self.pval_cutoff)
            self.BH_FDR_dict[col]=BH_FDR

            # Apply filtering criteria based on Log2FC and p-values
            right_df = compare_dataframe[
                (compare_dataframe[log2FC_col] >= Log2FC_cutoff) & 
                (compare_dataframe[pval_col] <= self.pval_cutoff) &
                (compare_dataframe['Removed?'] == "No")
                ]

            left_df = compare_dataframe[
                (compare_dataframe[log2FC_col] <= -Log2FC_cutoff) & 
                (compare_dataframe[pval_col] <= self.pval_cutoff) &
                (compare_dataframe['Removed?'] == "No")]

            if self.mark_contams == True:
                right_df = right_df[right_df['is_contaminant'] == "No"]
                left_df = left_df[left_df['is_contaminant'] == "No"]
        

            right_df[f"{col}_differential_abundance (FDR={(BH_FDR*100):.2f}%)"] = right_cond                
            left_df[f"{col}_differential_abundance (FDR={(BH_FDR*100):.2f}%)"] = left_cond

            all_sig_df = pd.concat([right_df,left_df,], axis = 0)

            print("HERE")
            print(all_sig_df)

            self.consolidated_df = pd.merge(self.consolidated_df, self.ordered_dataframe[
                ['Protein_Group', log2FC_col, pval_col]],
                on='Protein_Group', how='left')
            
            print("HERE")
            print(self.consolidated_df.columns)
            
            self.consolidated_df = pd.merge(self.consolidated_df, all_sig_df[
                ['Protein_Group', f"{col}_differential_abundance (FDR={(BH_FDR*100):.2f}%)"]
                ],
                on='Protein_Group', how='left')
            
            print("HERE")
            print(self.consolidated_df.columns)
            
            self.consolidated_df[f"{col}_differential_abundance (FDR={(BH_FDR*100):.2f}%)"].fillna("No", inplace=True)

            self.right_dict[col] = right_df
            self.left_dict[col] = left_df
            self.interactors_dict[col] = all_sig_df
            
            interactors_tsv_path = f'{self.op_folder}/{col}_logFC_compare_results.tsv'
            interactors_excel_path = f'{self.op_folder}/{col}_logFC_compare_results.xlsx'
            
            all_sig_df[['Protein_Group', 'Genes'] + [log2FC_col, pval_col]].to_csv(interactors_tsv_path, sep="\t", index=False)
            all_sig_df[['Protein_Group', 'Genes'] + [log2FC_col, pval_col]].to_excel(interactors_excel_path, index=False)

            print("Differentially abundant proteins saved to:", interactors_tsv_path, "and", interactors_excel_path)
        
        print(self.interactors_dict)
        #orig_df = pd.read_csv(self.bkgd_file, sep="\t")
        
        full_consolidated_excel_path = f'{self.table_op_folder}/consolidated_results.xlsx'
        self.consolidated_df.to_excel(full_consolidated_excel_path, index=False)
        
    def create_volcano_plot(self):

        def convert_to_percentage(condition):
            # Use regular expression to find the numerical part
            match = re.match(r"([<>=!]=?)(-?\d*\.?\d+)", condition)
            if match:
                operator = match.group(1)
                number = float(match.group(2))

                # Convert number to percentage
                percentage = number * 100

                # Format back into a string with the operator
                result = f"{operator}{percentage:.0f}%"
                return result
            else:
                # Return the original condition if it doesn't match the expected pattern
                return condition

        for comparison in self.comparison_columns:

            right_cond = comparison.split('_vs_')[0]
            left_cond = comparison.split('_vs_')[1]
            
            Log2FC_cutoff = self.logFC_cutoff_dict[comparison]


            self.ordered_dataframe['Genes'] = self.ordered_dataframe['Genes'].fillna('gene')
            num_DA = len(self.interactors_dict[comparison])
            print(num_DA)
            num_right_DA = len(self.right_dict[comparison])
            num_left_DA = len(self.left_dict[comparison])

            # Merge with the dataframe from interactors_dict
            merged_df = self.ordered_dataframe.merge(self.interactors_dict[comparison][['Protein_Group', 'fdr_corrected_pvals']], 
                                                     on='Protein_Group', how='left')
            
            # Fill NaN values in the 'fdr_corrected_pvals' column
            merged_df['fdr_corrected_pvals'] = merged_df['fdr_corrected_pvals'].fillna(np.nan)
            
            # Check if the corrected p-values meet the FDR threshold
            merged_df['meets_FDR'] = merged_df['fdr_corrected_pvals'] <= self.pval_cutoff
            valid_df = merged_df.dropna(axis=0, subset=[f'{comparison}_logFC'])
            valid_prots = [prot.capitalize() for prot in valid_df['Genes'].tolist()]
            print(valid_prots[0:20])
            num_all = len(valid_prots)

            #print ordered_dataframe
            def generate_hover_text(row):
                if self.mark_contams == True:
                    data = [
                        f"<b>Gene:</b> {row['Genes'].capitalize()}",
                        f"<b>Protein ID</b>: {row['Protein_Group']}",
                        f"<b>Contaminant?</b>: {row['is_contaminant']}",
                        f"<b>{comparison} log2FC</b>: {row[f'{comparison}_logFC']:.2f}",
                        f"<b>p-value</b>: {row[f'{comparison}_pval']:.2e}",
                        f"<b>Num. valid samples in {right_cond}:</b> {row[f'{right_cond}_quantified_samples']}",
                        f"<b>Average peptide count in {right_cond}:</b> {row[f'{right_cond}_Avg_Peptide_Count']:.2f}",
                        f"<b>Num. valid samples in {left_cond}:</b> {row[f'{left_cond}_quantified_samples']}",
                        f"<b>Average peptide count in {left_cond}:</b> {row[f'{left_cond}_Avg_Peptide_Count']:.2f}",
                        "<br><b>Click to link to Uniprot</b>."
                    ]
                else:
                    data = [
                        f"<b>Gene:</b> {row['Genes']}",
                        f"<b>Protein ID:</b> {row['Protein_Group']}",
                        f"<b>{comparison} log2FC:</b> {row[f'{comparison}_logFC']:.2e}",
                        f"<b>p-value:</b> {row[f'{comparison}_pval']:.2e}",
                        f"<b># quantified samples in {right_cond}:</b> {row[f'{right_cond}_quantified_samples']}",
                        f"<b># quantified samples in {left_cond}:</b> {row[f'{left_cond}_quantified_samples']}",
                        "<br>Click to link to Uniprot</b>."
                    ]
                return '<br>'.join(data)


            # Identify POIs and contaminants
            POIs = merged_df[merged_df['Genes'].str.upper().isin(self.extra_POIs)]
            contaminants = merged_df[merged_df['Protein_Group'].str.contains("CON__")]

            # Move POIs that are contaminants to the contaminants df
            POI_contaminants = POIs[POIs['Protein_Group'].str.contains("CON__")]
            POIs = POIs[~POIs.index.isin(POI_contaminants.index)]
            contaminants = pd.concat([contaminants, POI_contaminants]).drop_duplicates()

            # Identify non-POIs
            non_POI_df = merged_df[~merged_df['Genes'].str.upper().isin(self.extra_POIs)]

            # Subset non-POIs
            non_DA = non_POI_df[~non_POI_df['meets_FDR']]
            print("SHAPE")
            print(non_DA.shape)
            right_DA = non_POI_df[non_POI_df['Protein_Group'].isin(self.right_dict[comparison]['Protein_Group'].tolist())]
            left_DA = non_POI_df[non_POI_df['Protein_Group'].isin(self.left_dict[comparison]['Protein_Group'].tolist())]


            fig = make_subplots(rows=1, cols=1)

                        # Non-DA and Contaminants traces
            fig.add_trace(
                go.Scatter(
                    x=contaminants[f'{comparison}_logFC'],
                    y=-np.log10(contaminants[f'{comparison}_pval']),
                    mode='markers',
                    opacity=0.4,
                    marker=dict(color='grey', size=6, symbol='x'),
                    name='Contaminants',
                    hoverinfo='skip'
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=non_DA[f'{comparison}_logFC'],
                    y=-np.log10(non_DA[f'{comparison}_pval']),
                    mode='markers',
                    opacity=0.4,
                    marker=dict(color='grey', size=4),
                    name='Non-DA Proteins',
                    hoverinfo='skip'
                )
            )

            # Right-side DA Proteins
            fig.add_trace(
                go.Scatter(
                    x=right_DA[f'{comparison}_logFC'],
                    y=-np.log10(right_DA[f'{comparison}_pval']),
                    mode='markers',
                    opacity=0.6,
                    marker=dict(color=self.condition_colour_map[right_cond], size=8,
                                line=dict(color='white', width=1)
                                ),
                    text=right_DA.apply(generate_hover_text, axis=1),
                    hoverinfo='text',
                    customdata=right_DA['Protein_Group'].apply(lambda x: x.split(';')[0]),
                    name=f'Increased in {right_cond}'
                )
            )

            # Left-side DA Proteins
            fig.add_trace(
                go.Scatter(
                    x=left_DA[f'{comparison}_logFC'],
                    y=-np.log10(left_DA[f'{comparison}_pval']),
                    mode='markers',
                    opacity=0.6,
                    marker=dict(color=self.condition_colour_map[left_cond], size=8,
                                line=dict(color='white', width=1)
                                ),
                    text=left_DA.apply(generate_hover_text, axis=1),
                    hoverinfo='text',
                    customdata=left_DA['Protein_Group'].apply(lambda x: x.split(';')[0]),
                    name=f'Decreased in {right_cond}'
                )
            )

            # Separate traces for each POI
            if self.extra_POIs:
                for POI in self.extra_POIs:
                    poi_df = POIs[POIs['Genes'].str.upper() == POI.upper()]
                    if not poi_df.empty:
                        fig.add_trace(
                            go.Scatter(
                                x=poi_df[f'{comparison}_logFC'],
                                y=-np.log10(poi_df[f'{comparison}_pval']),
                                mode='markers',
                                marker=dict(color=self.POI_colour_map[POI.upper()], size=12,
                                            line=dict(color='white', width=1)
                                            ),
                                text=poi_df.apply(generate_hover_text, axis=1),
                                hoverinfo='text',
                                customdata=poi_df['Protein_Group'].apply(lambda x: x.split(';')[0]),
                                name=f'POI: {POI}'
                            )
                        )
        
            filter_string=convert_to_percentage(self.filterprop_th)
            pp_string=convert_to_percentage(self.ppushprop_th)

            # Update layout and axis labels
            fig.update_layout(
                title=f"Volcano plot of DA proteins: {right_cond} vs {left_cond}"
                    f"<br><sub>N<sub>total</sub>={num_all}; N<sub>DA</sub>={num_DA}; FDR={(100 * float(self.BH_FDR_dict[comparison])):.2f}%"
                    f"<br>logFC threshold = {Log2FC_cutoff:.2f}; -log10(p-value) corresponding to {self.pval_cutoff}."
                    f"<br>Only showing proteins verified in {filter_string} of samples in >=1 condition."
                    f"<br>P-values valid where {pp_string} of the relevant condition is verified, else pushed to 1 (NS).",
                xaxis_title="log2 fold change",
                yaxis_title="-log10 p-value",
                #plot_bgcolor='#EEEAE8',
                template="plotly_white",
                margin=dict(t=200)
            )

            # Wrap text for annotations
            right_text =f'DA proteins increased: <b>{num_right_DA}</b>'
            print(right_text)
            fig.add_annotation(text=right_text,
                               x=0.9, y=0.9,
                               xref="paper",
                               yref="paper",
                               showarrow=False
                               )

            left_text =f'DA proteins decreased: <b>{num_left_DA}</b>'
            fig.add_annotation(text=left_text,
                               x=0.1, y=0.9,
                               xref="paper",
                               yref="paper",
                               showarrow=False
                               )

            fig.update_layout(
            font=dict(color='black'),
            title_font=dict(color='black'),
            legend_title_font=dict(color='black')
            )
            config = {
              'toImageButtonOptions': {
                'format': 'png', # one of png, svg, jpeg, webp
                'filename': 'Volcano_plot',
                'scale': 3 # Multiply title/legend/axis/canvas sizes by this factor
              }
            }

            # Show the plot
            #fig.show(config=config)
            #fig.write_image(f"{self.op_folder}/{comparison}_volcano_plot.png", width=1200, height=500, scale=3)
            fig.write_html(f"{self.op_folder}/{comparison}_volcano_plot.html", config=config, include_plotlyjs='cdn', full_html=True, div_id='volcano-plot')

            # Append JavaScript to keep hover visible for 3 seconds on click
            with open(f"{self.op_folder}/{comparison}_volcano_plot.html", "a") as html_file:
                html_file.write("""
<script>
  document.addEventListener('DOMContentLoaded', function() {
    // Function to open a link based on the point's data
    function openLinkOnClick(eventData) {
    var plot = document.getElementById('volcano-plot');
      const point = eventData.points[0];
      const proteinId = point.customdata;

      if (proteinId) {
        const url = 'https://www.uniprot.org/uniprot/' + proteinId;
        console.log('Opening URL:', url);
        window.open(url, '_blank');
      }
    }

    // Adding click event listener to the plot
    const plotElement = document.getElementById('volcano-plot');
    if (plotElement) {
      plotElement.on('plotly_click', openLinkOnClick);
    } else {
      console.error('Plot element not found. Check if the ID is correct.');
    }
  });
</script>
"""
                )
    
    def capture_interactors(self):
        results_dict={}
        if not self.set_ids:
            self.set_ids=[str(i) for i in range(len(self.comparison_columns)*2)]
        for col,DA_direction,set_id in zip(self.comparison_columns,self.DA_directions, self.set_ids):
            if DA_direction.lower() ==  "up":
                results_dict[str(set_id)]=self.right_dict[col]['Genes'].values.tolist() 
            if DA_direction.lower() ==  "down":
                results_dict[str(set_id)]=self.left_dict[col]['Genes'].values.tolist()
            elif DA_direction.lower() == "both":
                for dic,sub_set_id in zip([self.right_dict,self.left_dict],str(set_id).split(";")):
                    results_dict[sub_set_id]=dic[col]['Genes'].values.tolist()
        return results_dict

    def run_analysis(self):
        self.make_op_folder()
        self.load_file_inputs()
        self.mark_contaminants_dataframe()
        self.mark_removed_prots()
        self.find_logFC_thresh()
        self.filter_and_save_interactors()
        self.create_volcano_plot()
        return self.capture_interactors()


if __name__ == "__main__":

    condition_colour_map = {
        "DMSO": "red",
        "413": "blue"
    }

    da = GetDA(
        bkgd_file=r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\LOESS\DAPARv1.30.6_DA_Sheet1_output.tsv",
        metadata_file=r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\LOESS\metadata_for_DAPAR.xlsx",
        op_folder=r"C:\Users\lwoods\Documents\LW_Projects_folder\TESTS\DIA_summarizer",
        comparison_columns=["413_vs_DMSO"],
        pval_cutoff=0.01,
        ppushprop_th=">=0.7",
        filterprop_th=">=0.7",
        condition_colour_map=condition_colour_map,
        Log2FC_cutoff=1.0,
        DA_directions=["up","down"],
        set_ids=None,
        bait_gene=None,
        mark_contams=True,
        Log2FC_cutoff_operation="None",
        slope = 0.01,
        contaminants_fasta = r"C:\Users\lwoods\Documents\LW_Projects_folder\general\20210914_525_CONTAMINANTS_UNIPROT_FORMAT.fasta",
        extra_POIs=['HMOX1', 'FAKE'],
        POI_colour_map={'HMOX1':'green', 'FAKE':'yellow'},
        remove_sheet = None, remove_list = None
    )
    da.run_analysis()    

