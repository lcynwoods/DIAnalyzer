import polars as pl
from pathlib import Path
import logging
import fastexcel
from Bio import SeqIO

class PrepareData:
    def __init__(self, spec_report, conditions_file, column_mapping, fasta_path=None, reattribute_peptides=False, peptides_to_remove=None, remove_outliers=False, min_pepts=0, log_file='./log', report_data=None):

        self.spec_report = spec_report
        self.conditions_file = conditions_file
        self.column_mapping = column_mapping
        self.fasta_path = fasta_path
        self.reattribute_peptides = reattribute_peptides
        self.peptides_to_remove = peptides_to_remove
        self.remove_outliers = remove_outliers
        self.min_pepts = min_pepts
        self.log_file = log_file
        self.logging_setup()
        self.report_data = report_data

        # Initialize returned vars
        self.spec_full_df = None
        self.conditions_df = None
        self.fasta_mapping = None

        # Column values
        self.file_name_col = None
        self.protein_group_col = None
        self.protein_ids_col = None
        self.genes_col = None
        self.precursor_id_col = None
        self.precursor_quantity_col = None
        self.modified_sequence_col = None
        self.protein_names_col = None
        self.run_col = 'Run'
        self.species_col = 'Species'

    def logging_setup(self):
        logging.basicConfig(filename=self.log_file, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s', force=True)

    def make_fasta_mapping(self):
        if self.fasta_path:
            mapping = {}
            with open(self.fasta_path, "r") as fasta_file:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    header = record.description
                    # Parse the header for required information
                    parts = header.split("|")
                    protein_id = parts[1]
                    protein_name = header.split(" ", 2)[2].split(" OS=")[0]
                    organism = header.split("OS=")[1].split(" OX=")[0]
                    gene_name = header.split("GN=")[1].split(" ")[0] if "GN=" in header else "Unknown"
                    mapping[protein_id] = (protein_name, gene_name, organism)
            self.fasta_mapping = mapping

    def load_report(self):
        spec_report_path = Path(self.spec_report)
        if spec_report_path.suffix.lower() == '.parquet':
            spec_full_df = pl.read_parquet(spec_report_path)
        else:
            spec_full_df = pl.read_csv(spec_report_path, separator='\t', encoding='latin-1',
                dtypes={'Decoy.Evidence': float,
                        'PTM.Specific': float,
                        'PTM.Localising':float,
                        'PTM.Site.Confidence':float,
                        'Lib.PTM.Site.Confidence':float,
                        'IM':float,
                        'iIM':float,
                        'Predicted.IM':float,
                        'Predicted.iIM':float,
                        })

        print("Headers before mapping:")
        print(spec_full_df.columns)

        # Mapping columns
        self.file_name_col = self.column_mapping['file_name_col']
        self.protein_group_col = self.column_mapping['protein_group_col']
        self.protein_ids_col = self.column_mapping['protein_ids_col']
        self.genes_col = self.column_mapping['genes_col']
        self.precursor_id_col = self.column_mapping['precursor_id_col']
        self.precursor_quantity_col = self.column_mapping['precursor_quantity_col']
        self.modified_sequence_col = self.column_mapping['modified_sequence_col']
        self.protein_names_col = self.column_mapping['protein_names_col']

        # Ensure that 'Run' is created correctly if it doesn't exist
        if self.run_col not in spec_full_df.columns:
            spec_full_df = spec_full_df.with_columns(
                (pl.col(self.file_name_col)
                .str.split('\\')
                .arr.get(-1)  # Use get(-1) to access the last element
                .str.replace('.raw', '')
                .str.replace('.mzML', '')
                ).alias(self.run_col)  # Apply alias at the end to name the resulting column
            )

        print("Headers after mapping:")
        print(spec_full_df.columns)

        spec_full_df = spec_full_df.with_columns(
            pl.col(self.genes_col).fill_null("gene")
        )
        replaced_rows = spec_full_df.filter(pl.col(self.genes_col) == "gene")

        # Step 5: Print the rows that had null values replaced
        print("Rows where null values have been replaced with 'gene':")
        print(replaced_rows)
        self.null_count=0
        # Apply any necessary transformations
        def fix_protein_accessions(s):
            if not s:
                self.null_count +=1
                return f"protein_{self.null_count}"
            return ';'.join(sorted(set(protein.strip() for protein in s.split(','))))

        def fix_genes(s):
            if not s:
                return "gene"
            return ';'.join(sorted(set(gene.strip().split('_')[0] for gene in s.split(','))))

        def split_for_species(s):
            if not s:
                return "Unknown"
            return s.split('_')[-1].split(";")[0] if '_' in s else "Unknown"

        # Convert columns to lists, apply the transformation, and then convert back to Series
        spec_full_df = spec_full_df.with_columns([
            pl.Series(self.protein_group_col, [fix_protein_accessions(x) for x in spec_full_df[self.protein_group_col].to_list()]).alias(self.protein_group_col),
            pl.Series(self.genes_col, [fix_genes(x) for x in spec_full_df[self.genes_col].to_list()]).alias(self.genes_col),
            pl.Series(self.protein_names_col, [split_for_species(x) for x in spec_full_df[self.protein_names_col].to_list()]).alias(self.species_col)
        ])

        print(f'2. The length is: {len(spec_full_df)} with headers: {spec_full_df.columns} after mapping and fixing columns')
        print(spec_full_df.head())
        self.spec_full_df = spec_full_df
        replaced_rows = spec_full_df.filter(pl.col(self.genes_col) == "gene")

        # Step 5: Print the rows that had null values replaced
        print("Rows where null values have been replaced with 'gene':")
        print(replaced_rows)


    def make_conditions_df(self):
        # Load the conditions DataFrame from an Excel file
        conditions_df = pl.read_excel(self.conditions_file, sheet_name="Sheet1")

        # Assigning technical replicate numbers based on Condition and Bio.Rep
        conditions_df = conditions_df.with_columns([pl.col('Bio.Rep').cum_count().over('Condition') + 1,(pl.col('Condition') + '_BR' + pl.col('Bio.Rep').cast(pl.Utf8) + '_TR' + pl.col('Tech.Rep').cast(pl.Utf8)).alias('CombinedName')])

        # Define required columns and their default values
        required_columns_defaults = {
            'Bait': 'None',
            'Marker': 'None',
            'POIs': 'None',
        }

        # Add any missing columns with their default values
        for col, default_val in required_columns_defaults.items():
            if col not in conditions_df.columns:
                print("missing COLS")
                print(col)
                conditions_df = conditions_df.with_columns(pl.lit(default_val).alias(col))
            else:
                print(f'col here {col}')
                conditions_df = conditions_df.with_columns(pl.col(col).fill_null(default_val))

        # Save the modified DataFrame to the instance variable
        self.conditions_df = conditions_df
        print(conditions_df)

    def merge_report_and_conditions(self):
        print(self.spec_full_df.columns)
        print(self.spec_full_df[self.run_col].unique().to_list())
        print(self.conditions_df.columns)
        print(self.conditions_df['Sample.Name'])
        self.spec_full_df = self.spec_full_df.join(self.conditions_df, left_on=self.run_col, right_on='Sample.Name', how='right')
        self.spec_full_df = self.spec_full_df.with_columns(pl.col('Sample.Name').alias(self.run_col))
        print(self.spec_full_df[self.run_col].unique().to_list())

        # Diagnostic: Check for any missing values post-merge
        missing_values = self.spec_full_df.filter(pl.col(self.run_col).is_null())
        if len(missing_values) > 0:
            print("Warning: Missing values detected in 'Run' column after merging:")
            print(missing_values)
        else:
            print("No missing values in 'Run' column after merging.")

        print(self.spec_full_df['Condition'].unique().to_list())

        # Print sample rows to verify merge success
        print("Sample rows from merged DataFrame:")
        print(self.spec_full_df.head(10))

    def filter_DIANN_data(self):
        print("Filtering")
        num_pept_before = len(self.spec_full_df)

        def reattribute_signal(df):
            # Convert to LazyFrame
            ldf = df.lazy()  

            # Step 1: Initialize new columns as copies of the old ones
            ldf = ldf.with_columns([
                pl.col(self.protein_group_col).alias('Reattributed.Protein'),
                pl.col(self.protein_ids_col).alias('Reordered.Protein.Ids'),
                pl.col(self.protein_names_col).alias('Reattributed.Names'),
                pl.col(self.genes_col).alias('Reattributed.Genes'),
                pl.col(self.species_col).alias('Reattributed.Species')
            ])

            if (self.fasta_mapping and self.reattribute_peptides):

                # Step 2: Add logic to reattribute based on presence of priority proteins
                ldf = ldf.with_columns([
                    # Create a boolean mask for rows that need reattribution
                    (pl.col('POIs').ne("None") | pl.col('Bait').ne("None") | pl.col('Marker').ne("None")).alias('needs_reattribution'),
                    
                    # Define the new protein order using Polars expressions
                    pl.when(pl.col('needs_reattribution'))
                    .then(
                        # Create an expression that combines proteins of interest into a reordered list
                        pl.concat_str([
                            pl.col('Bait').filter(pl.col('Bait').ne("None")),
                            pl.col('Marker').filter(pl.col('Marker').ne("None")),
                            pl.col('POIs').str.split(";").explode().unique().sort(),
                            pl.col(self.protein_ids_col).str.split(";").explode().unique().sort()
                        ], separator=";")
                    )
                    .otherwise(pl.col(self.protein_ids_col))
                    .alias('Reordered.Protein.Ids')
                ])

                # Step 3: Map back protein names, genes, and species using FASTA mapping
                mapping_df = pl.DataFrame([
                    {"protein_id": k, "name": v[0], "gene": v[1], "species": v[2]}
                    for k, v in self.fasta_mapping.items()
                ])

                # Join with the mapping DataFrame to update names, genes, and species
                ldf = ldf.join(mapping_df, left_on='Reattributed.Protein', right_on='protein_id', how='left')

            # Step 4: Collect the final result, triggering the execution
            df_result = ldf.collect()


            return df_result

        def perform_outlier_removal(df):
            unique_proteins = df.unique(subset=['PG.ProteinAccessions', 'PG.Genes'])
            pep_grouped = df.groupby(['EG.PrecursorId', 'PG.ProteinAccessions', 'Run']).agg(pl.mean('PEP.Quantity').alias('mean_quantity'))
            pep_pivot = pep_grouped.pivot(index=['EG.PrecursorId', 'PG.ProteinAccessions'], columns='Run', values='mean_quantity')
            run_cols = [col for col in pep_pivot.columns if col not in ['EG.PrecursorId', 'PG.ProteinAccessions']]

            outlier_corrected_list = []
            for run in run_cols:
                run_df = pep_pivot.select(['EG.PrecursorId', 'PG.ProteinAccessions', run])
                res = run_df.groupby('PG.ProteinAccessions').agg([
                    pl.col(run).quantile(0.05).alias('q1'),
                    pl.col(run).median().alias('median'),
                    pl.col(run).quantile(0.95).alias('q3')
                ])
                res = res.with_columns([
                    (pl.col('q3') - pl.col('q1')).alias('iqr'),
                    (pl.col('q1') - (pl.col('iqr') * 1.5)).alias('fence_low'),
                    (pl.col('q3') + (pl.col('iqr') * 1.5)).alias('fence_high')
                ])
                w_iqr = run_df.join(res, on='PG.ProteinAccessions', how='left')
                w_iqr = w_iqr.with_columns([
                    pl.when(pl.col(run) < pl.col('fence_low')).then(pl.col('fence_low')).otherwise(
                        pl.when(pl.col(run) > pl.col('fence_high')).then(pl.col('fence_high')).otherwise(pl.col(run))
                    ).alias('Corrected_intensity')
                ])
                outlier_corrected_list.append(w_iqr)

            outlier_corrected_df = pl.concat(outlier_corrected_list)
            df = df.join(outlier_corrected_df.select(['EG.PrecursorId', 'Run', 'Corrected_intensity']), on=['EG.PrecursorId', 'Run'], how='left')
            df = df.with_columns([
                pl.col('Corrected_intensity').alias('PEP.Quantity')
            ])
        
        def remove_low_pept_proteins(df, min_pepts):
            min_pepts = int(min_pepts)
            
            # Count peptides for each CombinedName and Reattributed.Protein group
            peptide_counts = (
                df.group_by(['CombinedName', 'Reattributed.Protein'])
                .agg(pl.count().alias('Peptide_Count'))
            )
            
            # Filter out proteins with fewer than the minimum peptide count
            valid_proteins = peptide_counts.filter(
                pl.col('Peptide_Count') >= min_pepts
            ).select(['CombinedName', 'Reattributed.Protein'])
            
            # Keep only rows in the original dataframe that match valid proteins
            df = df.join(
                valid_proteins, 
                on=['CombinedName', 'Reattributed.Protein'], 
                how='inner'
            )
            
            return df

        # DIANN filters
        # Ensure columns are of type float
        spec_full_df = self.spec_full_df.with_columns(
            [
                pl.col('Lib.PG.Q.Value').cast(pl.Float64).round(3),
                pl.col('Protein.Q.Value').cast(pl.Float64).round(3),
                pl.col('PG.Q.Value').cast(pl.Float64).round(3),
            ]
        )
        print("A")
        print(self.spec_full_df['Condition'].unique().to_list())

        # Print initial state
        print("Before filtering:")
        print(spec_full_df.head(10))  # Show first 10 rows for initial inspection
        print("B")
        print(spec_full_df['Condition'].unique().to_list())

        # Apply filtering
        spec_full_df = spec_full_df.filter(
            (pl.col('Lib.PG.Q.Value') < 0.01) &
            (pl.col('Protein.Q.Value') < 0.01) &
            (pl.col('PG.Q.Value') < 0.01)
        )

        if self.peptides_to_remove:
            for peptide in self.peptides_to_remove:
                print("Removing peptide:")
                print(peptide)
                spec_full_df = spec_full_df.filter(pl.col('Stripped.Sequence') != peptide)
        
        if self.remove_outliers:
            print("Removing outliers")
            perform_outlier_removal(spec_full_df)
        
        
        print('5')
        print(spec_full_df['CombinedName'].unique().to_list())
        
        print("Reattributing signal")
        spec_full_df = reattribute_signal(spec_full_df)
        
        if self.min_pepts > 1:
            spec_full_df = remove_low_pept_proteins(spec_full_df, self.min_pepts)
        
        num_pept_after = len(spec_full_df)        

        logging.info(f'Number of peptides before filtering: {num_pept_before}')
        logging.info(f'Number of peptides after filtering: {num_pept_after} ({num_pept_before-num_pept_after} peptides filtered out)')
        logging.info('Note that modified peptides are considered individual peptides.')

        self.spec_full_df = spec_full_df
        print(self.spec_full_df)
        print(self.spec_full_df.columns)
        print('6')
        print(self.spec_full_df['CombinedName'].unique().to_list())


        self.report_data['data_filtering']['peptides_before'] = num_pept_before
        self.report_data['data_filtering']['peptides_after'] = num_pept_after
    
    def find_average_peptides_per_condition(self):
        peptide_counts = (
            self.spec_full_df
            .group_by(['CombinedName', 'Reattributed.Protein'])
            .agg(pl.count().alias('Peptide_Count'))
        )
        print(peptide_counts)
        peptide_counts = peptide_counts.join(
            self.spec_full_df.select(['CombinedName', 'Condition']).unique(), on='CombinedName', how='left'
        )
        print(2)
        print(self.spec_full_df.columns)
        average_peptide_counts = (
            peptide_counts
            .group_by(['Condition', 'Reattributed.Protein'])
            .agg(pl.col('Peptide_Count').mean().alias('Peptide_Count'))
        )
        print(3)
        print(self.spec_full_df.columns)
        average_peptide_pivot = (
            average_peptide_counts
            .pivot(index='Reattributed.Protein', columns='Condition', values='Peptide_Count')
            .fill_null(0)
        )
        average_pept_cols = [f'{col}_Avg_Peptide_Count' for col in average_peptide_pivot.columns if col != 'Reattributed.Protein']
        average_peptide_pivot.columns = ['Reattributed.Protein'] + average_pept_cols
        print(average_pept_cols)
        print(4)
        print(self.spec_full_df.columns)

        # Add averages to report df
        self.spec_full_df = self.spec_full_df.join(average_peptide_pivot, on='Reattributed.Protein', how='left')
        replaced_rows = self.spec_full_df.filter(pl.col(self.genes_col) == "gene")

        # Step 5: Print the rows that had null values replaced
        print("Rows where null values have been replaced with 'gene':")
        print(replaced_rows)

    def write_output(self):
        self.spec_full_df.write_parquet(str(Path(r"C:\Users\lwoods\Documents\LW_Projects_folder\TESTS\DIA_summarizer\after_reattribute_and_filter.parquet")))

    def run(self):
        logging.info('Loading MS report')
        self.load_report()
        logging.info('Successfully loaded MS report')
        logging.info('Organizing conditions')
        self.make_conditions_df()
        logging.info('Organized conditions')
        logging.info('Merging conditions to MS report')
        self.merge_report_and_conditions()
        logging.info('Merged data')
        logging.info('Filtering data')
        self.filter_DIANN_data()
        logging.info('Filtered data')
        logging.info('Adding some new columns...')
        self.find_average_peptides_per_condition()
        logging.info('Calculated average peptides')
        #self.write_output()
        return self.spec_full_df, self.conditions_df
