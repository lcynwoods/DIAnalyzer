import numpy as np
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

class PrepareQuantificationData:
    def __init__(self,
            filtered_df,
            conditions_df,
            sorted_conditions,
            colour_map,
            op_folder,
            intermediate_op_folder = None,
            POIs=None,
            POI_colour_map=None,
            use_DIANN_maxLFQ=False,
            biomolecule_level = "Protein",
            is_IP=False,
            is_POC=False,
            report_data=None,
            log_file = None
            ):
        self.spec_full_df = filtered_df
        self.conditions_df = conditions_df
        self.sorted_conditions = sorted_conditions
        self.colour_map = colour_map
        self.POIs = POIs
        self.POI_colour_map = POI_colour_map
        self.use_DIANN_maxLFQ = use_DIANN_maxLFQ
        self.biomolecule_level = biomolecule_level
        self.is_IP = is_IP
        self.is_POC = is_POC
        self.op_folder = op_folder
        if intermediate_op_folder:
            self.intermediate_op_folder = intermediate_op_folder
        else:
            self.intermediate_op_folder = self.op_folder
        self.report_data=report_data,
        self.log_file = log_file

        # Generated variables
        self.intensity_per_prot_df = None
        self.R_input_table = None
        self.average_pept_cols = None

        # Define column names as variables
        # from original report
        self.combined_name_col = 'CombinedName'
        #protein_accessions_col = 'PG.ProteinAccessions'
        self.protein_accessions_col = 'Reattributed.Protein'
        self.genes_col = 'Reattributed.Genes'
        self.species_col = 'Reattributed.Species'
        # Intensity cols
        self.pep_quantity_col = 'Precursor.Quantity'
        self.precursor_id_col = 'EG.PrecursorId'
        self.maxlfq_col = 'PG.MaxLFQ'
        self.pg_quantity_col = 'PG.Quantity'
        self.protein_intensity_col = 'ProteinIntensity'
        self.protein_names_col = 'Protein.Names'
        # From conditions df
        self.sample_name_col = 'Sample.Name'
        self.condition_col = 'Condition'
        self.bio_rep_col = 'Bio.Rep'
        self.tech_rep_col = 'Tech.Rep'
        self.bait_col = 'Bait'
        self.marker_col = 'Marker'
        self.POIs_col = 'POIs'

    
    def convert_polars_to_pandas(self):
        self.spec_full_df = self.spec_full_df.to_pandas()
        self.conditions_df = self.conditions_df.to_pandas()

    def setup_intensity_table(self):
        print("Starting quantification")
        print(self.sorted_conditions)

        average_pept_cols = [f'{condition}_Avg_Peptide_Count' for condition in self.sorted_conditions]

        # Define the list of columns to group by
        group_by_cols = [
            self.combined_name_col,
            self.protein_accessions_col,
            self.genes_col,
            self.species_col,
            self.condition_col,
            self.bio_rep_col,
            self.tech_rep_col,
            self.bait_col,
            self.marker_col,
            self.POIs_col,
            self.protein_names_col,
        ] + average_pept_cols

        # Verify all columns are present
        missing_cols = [col for col in group_by_cols if col not in self.spec_full_df.columns]

        # Working copy of DataFrame
        spec_full_df = self.spec_full_df.copy()

        spec_full_df[group_by_cols] = spec_full_df[group_by_cols].fillna('Unknown')

        # Ensure species column is treated as string
        spec_full_df[self.species_col] = spec_full_df[self.species_col].astype(str)

        unknown_species_df = spec_full_df[
            spec_full_df[self.species_col].isna() |  # Check for NaN
            (spec_full_df[self.species_col].astype(str).str.strip() == "")  # Check for empty string
        ]

        print("Problematic Rows (Species Missing):")
        print(unknown_species_df)



        # Check if DIANN maxLFQ is not used
        print(self.use_DIANN_maxLFQ)
        if not self.use_DIANN_maxLFQ:
            print("SUMMING")
            # Sum peptide quantities grouped by group_by_cols
            intensity_per_prot_df = spec_full_df.groupby(group_by_cols, as_index=False).agg({self.pep_quantity_col: 'sum'})

        elif self.maxlfq_col in spec_full_df.columns:
            print("USING MAXLFQ")
            # Select specific columns and use maxLFQ for intensity
            intensity_per_prot_df = spec_full_df[group_by_cols + [self.maxlfq_col]].copy()
        elif self.pg_quantity_col in spec_full_df.columns:
            # Select specific columns and use PG quantity for intensity
            intensity_per_prot_df = spec_full_df[group_by_cols + [self.pg_quantity_col]].copy()
        else:
            raise ValueError("Neither maxLFQ nor PG.Quantity column is present for intensity calculation.")

        # Rename the intensity column to a common name
        intensity_per_prot_df.columns = group_by_cols + [self.protein_intensity_col]

        # Replace 0 with np.nan for log transformation readiness
        intensity_per_prot_df[self.protein_intensity_col] = intensity_per_prot_df[self.protein_intensity_col].replace(0, np.nan)

        # Debug print to check if maximum intensity is being calculated correctly
        print("Maximum intensity:", intensity_per_prot_df[self.protein_intensity_col].max())

        self.intensity_per_prot_df = intensity_per_prot_df
        self.average_pept_cols = average_pept_cols

    
    def analyze_biomolecules(self):
        if self.biomolecule_level == "Protein":
            self.analyze_proteins()
        elif self.biomolecule_level == "Peptide":
            print("Sorry, not doing peptides yet!")
            return
        
    def analyze_proteins(self):
        print(f'Original DataFrame length: {len(self.intensity_per_prot_df)}')
        intensity_grouped = self.intensity_per_prot_df.groupby([self.protein_accessions_col, self.genes_col, self.protein_names_col]+self.average_pept_cols+[self.species_col,self.combined_name_col])[self.protein_intensity_col].mean().reset_index()
        print(f'Grouped with average pept cols length: {len(intensity_grouped)}')

        pivoted = intensity_grouped.pivot_table(index=[self.protein_accessions_col, self.genes_col, self.protein_names_col]+self.average_pept_cols+[self.species_col], columns=self.combined_name_col, values=self.protein_intensity_col).reset_index()
        pivoted.columns.name=None
        print(f'Pivoted with pept avg length: {len(pivoted)}')

        test_intensity_grouped = self.intensity_per_prot_df.groupby([self.protein_accessions_col,self.combined_name_col])[self.protein_intensity_col].mean().reset_index()
        print(f'Grouped without average pept cols length: {len(test_intensity_grouped)}')
        test_pivoted = test_intensity_grouped.pivot_table(index=[self.protein_accessions_col], columns=self.combined_name_col, values=self.protein_intensity_col).reset_index()
        test_pivoted.columns.name=None
        print(f'Pivoted without pept avg length: {len(test_pivoted)}')

        merge_cols = [self.protein_accessions_col,self.genes_col,self.protein_names_col, self.species_col]
        merge_cols.extend(self.average_pept_cols)

        pivoted = test_pivoted.merge(pivoted[merge_cols], on=self.protein_accessions_col, how="left")
        print(f'After merge length: {len(pivoted)}')

        print(f'with pept avg: {len(pivoted)}')
        print(f'without pept avg: {len(test_pivoted)}') 

        pivoted.insert(loc=0, column='artifical_id', value=[i for i in range(1,len(pivoted)+1)])
        # Renaming columns for DAPAR input file
        pivoted=pivoted.rename(columns={self.protein_accessions_col:'Protein_Group',
                                    self.genes_col:'Genes',
                                    self.species_col : 'Species',
                                    self.protein_names_col: 'Protein_Names'
                                    })

        DAPAR_cols=['Protein_Group','Genes','Protein_Names'] + self.average_pept_cols + ['Species']
        runs = self.conditions_df[self.combined_name_col].values.tolist()


        #########
        print("HERE")
        print(runs)
        print(pivoted.columns)


        for run_col in runs:
            pivoted['Intensity_' + run_col]=pivoted[run_col]
            DAPAR_cols.append('Intensity_' + run_col)
        
        print(f'Length following renaming: {len(pivoted)}')


        R_input_table=pivoted
        R_input_filename = str(Path(self.intermediate_op_folder) / f'{self.biomolecule_level}_quantification_for_R.xlsx')

        R_input_table.to_excel(R_input_filename,columns=DAPAR_cols,index=False)
        R_input_table.to_csv(R_input_filename.replace('xlsx','tsv'),columns=DAPAR_cols,sep='\t',index=False)
        
    def analyze_peptides(self):
        peptide_grouped = self.spec_full_df.groupby([self.precursor_id_col, self.protein_accessions_col, self.genes_col, self.combined_name_col])[self.pep_quantity_col].mean().reset_index()
        pivoted = peptide_grouped.pivot_table(index=[self.precursor_id_col, self.protein_accessions_col, self.genes_col], columns=self.combined_name_col, values=self.pep_quantity_col).reset_index()
        pivoted.insert(loc=0, column='artifical_id', value=[i for i in range(1,len(pivoted)+1)])
        pivoted=pivoted.rename(columns={self.precursor_id_col:'Precursor.Id',
                                        self.protein_accessions_col:'Protein_Group',
                                    self.genes_col:'Genes',
                                    self.protein_names_col:'Protein_Names'
                                    })
        #DAPAR_cols=['Protein_Group','Genes','Precursor.Id']
        #R_input_table=pivoted
    
    def save_plot(self, fig, title, y_title, file_name):
        # HTML Configuration
        config = {
            f'toImageButtonOptions': {
                'format': 'png',  # one of png, svg, jpeg, webp
                'filename': f'{file_name}.png'},
                'scale': 2  # Multiply title/legend/axis/canvas sizes by this factor
        }

        fig.update_layout(
            title=title,
            template='plotly_white',
            title_font_color='black',
            xaxis=dict(
                title="Sample",
                titlefont_size=12,
                color='black'
            ),
            yaxis=dict(
                title=y_title,
                titlefont_size=12,
                color='black'
            ),
            legend=dict(
                font=dict(
                    color='black'
                )
            )
        )
        plot_folder = Path(self.op_folder)
        plot_folder.mkdir(exist_ok=True)
        output_path = plot_folder / f'{file_name}.html'
        fig.write_html(output_path, config=config)

        #if self.report_data:
        #    self.report_data['plots'].append({
        #        'name': f'{title}',
        #        'fig': fig.to_html()
        #    })
  
    def plot_runwise_intensities(self):

        fig4 = go.Figure()

        for condition in self.sorted_conditions:
            fig4.add_trace(go.Violin(
                x=self.intensity_per_prot_df[self.intensity_per_prot_df[self.condition_col] == condition][self.combined_name_col],
                y=np.log2(self.intensity_per_prot_df[self.intensity_per_prot_df[self.condition_col] == condition][self.protein_intensity_col] + 1),
                name=condition,
                line=dict(
                    color=self.colour_map[condition],  # Set the line (border) color
                    width=2  # Border width
                ),
                fillcolor=self.colour_map[condition],  # Set the fill color
                opacity=0.5,
                hoverinfo='skip',
                points=False,
                pointpos=0  # Center the violins on the categories
            ))


        all_samples = self.intensity_per_prot_df[self.combined_name_col].unique()
        # Plot each protein of interest
        intensity_type = "maxLFQ Intensity" if self.use_DIANN_maxLFQ else "Intensity"
        if self.POIs:
            for protein in self.POIs:
                fig5 = None
                POI_data = self.intensity_per_prot_df[self.intensity_per_prot_df[self.genes_col].str.upper() == protein.upper()]
                print(POI_data)
                missing_samples = set(all_samples) - set(POI_data[self.combined_name_col])
                missing_data = pd.DataFrame({
                    self.combined_name_col: list(missing_samples),
                    self.protein_intensity_col: [0] * len(missing_samples),
                    self.genes_col: [protein.lower()] * len(missing_samples),
                    self.protein_accessions_col: [f'{protein}_prot'] * len(missing_samples)
                }).merge(self.conditions_df, on=self.combined_name_col)
                POI_data = pd.concat([POI_data, missing_data]).drop_duplicates()
                POI_data = POI_data[POI_data['Condition'].isin(self.sorted_conditions)]

                hover_texts = [
                    f"<b>Gene:</b> {gene}</b><br><b>{intensity_type}:</b> {intensity}"
                    for gene, intensity in zip(POI_data[self.genes_col], POI_data[self.protein_intensity_col])
                ]

                # Add scatter to intensity plot for each protein of interest
                fig4.add_trace(
                    go.Scatter(
                        x=POI_data[self.combined_name_col],
                        y=np.log2(POI_data[self.protein_intensity_col] + 1),
                        mode='markers',
                        marker=dict(
                            color=self.POI_colour_map[protein],  # Use consistent color for each protein
                            size=10,
                            opacity=[1 if val > 0 else 0 for val in POI_data[self.protein_intensity_col]],
                            line=dict(
                                color='white',
                                width=[1 if val > 0 else 0 for val in POI_data[self.protein_intensity_col]]
                            ),
                        ),
                        name=protein.capitalize(),  # Use the current protein name in the legend
                        hoverinfo='text',
                        hovertext=hover_texts
                    )
                )
                print(POI_data)
        
                fig5 = px.bar(
                    POI_data[POI_data['Condition'].isin(self.sorted_conditions)],
                    x='CombinedName',
                    y='ProteinIntensity',
                    color='Condition',  # Ensure this is set correctly
                    color_discrete_map=self.colour_map,
                    category_orders={'Condition': self.sorted_conditions}
                    )
                
                self.save_plot(fig5, title=f"{protein} {intensity_type} by sample<br><sub>Before normalization and/or imputation.", y_title=f'{protein} {intensity_type}', file_name=f'{protein.capitalize()}_absolute_intensity')
        

        self.save_plot(fig4, title=f'{intensity_type} by sample<br><sub>Before normalization and/or imputation.', y_title=f'log2({intensity_type})', file_name=f'Intensity_by_run')

        
    def run(self):
        self.convert_polars_to_pandas()
        self.setup_intensity_table()
        self.analyze_biomolecules()
        self.plot_runwise_intensities()
        return self.spec_full_df