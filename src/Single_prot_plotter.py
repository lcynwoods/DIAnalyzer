from pathlib import Path
import pandas as pd
import plotly.express as px


class POI_plotter():
    def __init__(self, protein, full_data, sorted_conditions, conditions_df, colour_map, file_name, title, y_title, reference_condition=None, x_title="Sample", col_suffix="", require_untransform=False):
        self.protein = protein
        self.full_data = full_data
        self.file_name = file_name
        self.title = title
        self.y_title = y_title
        self.x_title = x_title
        self.conditions_df = conditions_df
        self.sorted_conditions = sorted_conditions
        self.reference_condition = reference_condition
        self.colour_map = colour_map
        self.col_suffix = col_suffix
        self.require_untransform = require_untransform

        self.genes_col = "Genes"
        self.POI_data = None

    def load_data(self):
        def load_anything(data_item):
            # If self.df is already a DataFrame, do nothing
            if isinstance(data_item, pd.DataFrame):
                return data_item

            # If self.df is a string (a file path)
            elif isinstance(data_item, str):
                file_path = Path(data_item)
                
                # Check if the file exists
                if file_path.exists():
                    # Load the file based on its extension
                    if file_path.suffix in ('.tsv', '.txt'):
                        data_item = pd.read_csv(file_path, sep='\t')
                    elif file_path.suffix == '.csv':
                        data_item = pd.read_csv(file_path)
                    elif file_path.suffix == '.xlsx':
                        data_item = pd.read_excel(file_path)
                    else:
                        print(f"Can't read input file type: {file_path.suffix}")
                        return None
                else:
                    print(f"{file_path} doesn't exist!")
                    return None
            
            else:
                print("No valid data could be loaded")
                return None
                
            return data_item

        self.full_data = load_anything(self.full_data)
        self.conditions_df = load_anything(self.conditions_df)
        # If self.df is neither a DataFrame nor a file path string
        

    def save_plot(self, fig, title, y_title, file_name, x_title="Sample"):
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
                title=x_title,
                titlefont_size=16,
                color='black'
            ),
            yaxis=dict(
                title=y_title,
                titlefont_size=16,
                color='black'
            ),
            legend=dict(
                font=dict(
                    color='black'
                )
            )
        )
        #fig.show()
        fig.write_html(self.file_name, config=config)

    def plot_relative_intensities(self):
        intensity_columns = [col for col in self.full_data.columns if (col.startswith('Intensity_') and col.endswith(self.col_suffix))]

        # Step 2: Extract the data with 'Genes' and intensity columns
        chopped_data = self.full_data[['Genes'] + intensity_columns]

        for col in intensity_columns:
            if self.require_untransform == True:
                chopped_data[col] = chopped_data[col].apply(lambda x: 2 ** x)

        # Step 3: Remove 'Intensity_' prefix from sample names in intensity columns
        chopped_data.columns = ['Genes'] + [col.replace('Intensity_', '').replace(self.col_suffix, "") for col in intensity_columns]


        # Step 3: Pivot the chopped data from wide to long format
        pivoted_data = pd.melt(
            chopped_data, 
            id_vars=['Genes'],  # Keep the 'Genes' column
            var_name='Sample.Name',  # Use the cleaned-up sample names
            value_name='Intensity'  # Intensity values for each sample
        )

        

        # Step 4: Filter the data to only include the specified protein (e.g., Protein A)
        protein_data = pivoted_data[pivoted_data['Genes'].str.lower() == self.protein.lower()]
        print("This is what you want")
        print(protein_data)
        
        # Step 5: Merge with conditions_df to add the condition information to each sample
        merged_data = pd.merge(protein_data, self.conditions_df, on='Sample.Name')

        # Calculate relative intensities if reference condition is provided
        if self.reference_condition:
            reference_value = merged_data[merged_data['Condition'] == self.reference_condition]['Intensity'].mean()
            merged_data['Relative_intensity'] = merged_data['Intensity'] / reference_value * 100
            y_col = 'Relative_intensity'
        else:
            y_col = 'Intensity'

        # Create the bar plot using plotly
        fig = px.bar(
            merged_data,
            x='Sample.Name',  # Each sample gets a bar
            y=y_col,  # The y-values are either the intensity or relative intensity
            color='Condition',  # Color bars by condition
            color_discrete_map=self.colour_map,
            category_orders={'Condition': self.sorted_conditions},
            #barmode='group'  # Group the bars by sample
        )

        # Add the dashed line for the reference condition's average intensity if applicable
        if self.reference_condition:
            fig.add_shape(
                type="line",
                x0=-0.5, x1=len(merged_data['Sample.Name'].unique()) - 0.5, y0=100, y1=100,  # Line across all samples
                line=dict(color="black", width=2, dash="dash"),
                xref='x', yref='y'
            )


        # Show the plot
        #fig.show()

        # Save the plot
        self.save_plot(fig, title=self.title, y_title=self.y_title, file_name=self.file_name)


    def run(self):
        self.load_data()
        self.plot_relative_intensities()

    
if __name__ == "__main__":
    import pandas as pd

    # Example data (full_data)
    # Example full_data
    full_data = pd.DataFrame({
        'Genes': ['ProteinA', 'ProteinB', 'ProteinC'],
        'Intensity_Treatment A': [100, 90, 85],
        'Intensity_Treatment B': [110, 95, 80],
        'Intensity_Treatment C': [105, 98, 88],
        'Intensity_Treatment D': [105, 98, 88],
        'Intensity_Treatment E': [66, 77, 88],
        'Intensity_Treatment F': [105, 98, 88]
    })

    conditions_df = pd.DataFrame({
    'Sample.Name': ['Treatment A', 'Treatment B', 'Treatment C', 'Treatment D', 'Treatment E', 'Treatment F'],
    'Condition': ['Control', 'Treatment1', 'Treatment2','Control', 'Treatment1', 'Treatment2']
    })

    # Variables to fill into your class
    PROTEIN = 'ProteinA'  # Protein of interest
    #FULL_DATA = r"C:\Users\lwoods\Documents\LW_Projects_folder\09_EXT\G39_Gemma_Fabrias_Bernat_Crosas\R01_MireiaCasasempere\P01_CERT_Oncology\myE01_CERTACs\GPA_analysis\LOESS\DAPARv1.30.6_DA_Sheet1_output.tsv"  # The DataFrame with data
    FULL_DATA = full_data
    SORTED_CONDITIONS = ['Control', 'Treatment1', 'Treatment2']  # Custom condition order
    CONDITIONS_DF = conditions_df
    COLOUR_MAP = {
        'Control': 'blue',
        'Treatment1': 'green',
        'Treatment2': 'red'
    }  # Color map for plot
    FILE_NAME = 'protein_plot.html'  # File to save the plot
    TITLE = 'Protein A Intensity Across Conditions'  # Title of the plot
    Y_TITLE = 'Relative Intensity (%)'  # Y-axis title
    X_TITLE = 'Sample'  # X-axis title
    REFERENCE_CONDITION = 'Control'  # Reference condition
    COL_SUFFIX = ''  # Optional column suffix

    test=POI_plotter(protein=PROTEIN,
                full_data=FULL_DATA,
                sorted_conditions=SORTED_CONDITIONS,
                conditions_df=CONDITIONS_DF,
                colour_map=COLOUR_MAP,
                file_name=FILE_NAME,
                title=TITLE,
                y_title=Y_TITLE,
                reference_condition=REFERENCE_CONDITION,
                x_title=X_TITLE,
                col_suffix=COL_SUFFIX)
    test.run()



