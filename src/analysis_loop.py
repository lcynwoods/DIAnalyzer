from common_classes import *

class AnalysisLoop():

    def __init__(self,
            analyses_file,
            data_file,
            conditions_df,
            op_folder,
            r_exe_path,
            r_script_path,
            colour_map,
            sorted_conditions,
            reference_condition=None,
            POIs=None,
            POI_colour_map=None,
            species_id=9606,
            go_label="biological_process", #For STRING--currently broken
            use_DIANN_maxLFQ=None,
            biomolecule_level = 'protein', #Hope to have "peptide" in the future
            POC = False,
            mark_contaminants=False,
            contaminants_fasta=None,
            gene_sets =['GO_Biological_Process_2023', 'MSigDB_Hallmark_2020'],
        ):
        
        self.analyses_file = analyses_file
        self.data_file = data_file
        self.conditions_df = conditions_df
        self.op_folder = op_folder
        self.r_exe_path = r_exe_path
        self.r_script_path = r_script_path
        self.colour_map = colour_map
        self.sorted_conditions = sorted_conditions
        self.reference_condition = reference_condition
        self.POIs = POIs
        self.POI_colour_map = POI_colour_map
        self.species_id = species_id
        self.go_label = go_label
        self.use_DIANN_maxLFQ = use_DIANN_maxLFQ

        self.MS_data_type = "DIA"
        self.biomolecule_level = biomolecule_level
        self.POC = POC
        self.contaminants_fasta = contaminants_fasta
        self.mark_contaminants = mark_contaminants
        self.batch_col = "No"
        self.metadata_file = None
        self.DAPAR_output_file = None
        self.gene_sets = gene_sets

    def setup_subdirs(self):
        job_name = Path(self.op_folder).stem
        plot_folder = Path(self.op_folder) /f'{job_name}_data' / 'plots'
        plot_folder.mkdir(exist_ok=True, parents=True)
        self.summary_data_folder = plot_folder / 'summary_data'
        self.summary_data_folder.mkdir(exist_ok=True)
        self.DA_data_folder = plot_folder / 'DA_data'
        self.DA_data_folder.mkdir(exist_ok=True)
        self.POI_data_folder = plot_folder / 'POI_data'
        self.POI_data_folder.mkdir(exist_ok=True)
        self.cluster_data_folder = plot_folder / 'clustering_data'
        self.cluster_data_folder.mkdir(exist_ok=True)
        self.ORA_data_folder = plot_folder / 'ORA_data'
        self.ORA_data_folder.mkdir(exist_ok=True)
        self.GSEA_data_folder = plot_folder / 'GSEA_data'
        self.GSEA_data_folder.mkdir(exist_ok=True)
        table_folder = Path(self.op_folder) /f'{job_name}_data' / 'tables'
        table_folder.mkdir(exist_ok=True, parents=True)
        intermediate_folder = Path(self.op_folder) /f'{job_name}_data' / 'intermediate_files'
        intermediate_folder.mkdir(exist_ok=True, parents=True)
        self.DA_table_folder = table_folder / 'DA_data'
        self.DA_table_folder.mkdir(parents=True, exist_ok=True)
        self.cluster_table_folder = table_folder / 'clustering_data'
        self.cluster_table_folder.mkdir(exist_ok=True)
        self.ORA_table_folder = table_folder / 'ORA_data'
        self.ORA_table_folder.mkdir(exist_ok=True)
        self.ORA_intermediate_folder = intermediate_folder / 'ORA_data'
        self.ORA_intermediate_folder.mkdir(exist_ok=True)
        self.GSEA_table_folder = table_folder / 'GSEA_data'
        self.GSEA_table_folder.mkdir(exist_ok=True)
        self.GSEA_intermediate_folder = intermediate_folder / 'GSEA_data'
        self.GSEA_intermediate_folder.mkdir(exist_ok=True)
    
    def prepare_DAPAR_metadata(self):
        try:
            metadata_cols = ['CombinedName', 'Condition', 'Bio.Rep', 'Tech.Rep']
            for col in self.conditions_df.columns:
                if str(col).upper() == "BATCH":
                    self.batch_col = "Yes"
                    metadata_cols.append(self.batch_col)
            metadata_df = pd.DataFrame(self.conditions_df[metadata_cols])
            metadata_df.columns = ['Sample.Name', 'Condition', 'Bio.Rep', 'Tech.Rep']
            self.metadata_file = str(Path(self.op_folder) / 'metadata_for_DAPAR.xlsx')
            metadata_df.to_excel(self.metadata_file, index=False)
            print(f"Metadata file created at: {self.metadata_file}")
        except Exception as e:
            print(f"Failed to prepare metadata: {e}")
            raise
    
    def run(self):
        self.setup_subdirs()
        self.prepare_DAPAR_metadata()
        sheets = pd.read_excel(self.analyses_file, sheet_name=None)
        print(self.metadata_file)
        all_comparison_sets = {}
        for sheet_name, df in sheets.items():
            comparison_sets = {}
            # Should integrate "comparisons" early
            if 'Right_condition' in df.columns:
                df['Comparisons'] = df['Right_condition'].str.cat(df[['Left_condition']], sep='_vs_')
            raw_comparisons = df['Comparisons'].values.tolist()
            comparisons_list = [comp.strip() for comp in raw_comparisons]
            cleaned_comparisons = ';'.join(comparisons_list)
            DAPAR_desc = f"DAPARv1.34.6_DA_{sheet_name}_output" # TODO: determine version
            DAPAR_output_file = str(Path(self.op_folder) / f'{DAPAR_desc}.tsv')
            user_normalization = df['Normalization'][0]
            if self.use_DIANN_maxLFQ:
                user_normalization = 'DIANN_maxLFQ'
            normalization_type = df['Normalization.Type'][0]
            normalization_type = normalization_type.lower()
            valid_values_filter_thresh = df['Valid.Values.Filter.Thresh'][0]
            valid_values_pp_thresh = df['Valid.Values.PP.Thresh'][0]
            logFC_thresh = df['LogFC.Thresh'][0]
            pval_thresh = df['Pval.Thresh'][0]
            sides =df['Side'].values.tolist()
            result_set_ids = df['ResultSet'].values.tolist()
            ECDF_slope_thresh = 0.01


            R_analysis_subprocess = R_normalization_and_analysis(
                self.r_exe_path,
                self.r_script_path,
                DAPAR_output_file,
                self.data_file,
                self.metadata_file,
                comparisons=cleaned_comparisons,
                MS_data_type=self.MS_data_type,
                biomolecule_level=self.biomolecule_level.lower(),
                DAPAR_norm_type = user_normalization,
                DAPAR_normalization_within = normalization_type,
                th_filter = valid_values_filter_thresh,
                th_push = valid_values_pp_thresh,
                batch_col = self.batch_col
            )
            if not self.POC:
                R_analysis_subprocess.run_DAPAR()
                main_data_file = DAPAR_output_file
                require_log = False
            else:
                main_data_file = self.data_file
                require_log = True

            # Show normalization plots
            # For now, only one plot for maxLFQ available
            if self.use_DIANN_maxLFQ and not self.POC:
                require_log = True
                normalization_method = 'DIANN_maxLFQ'
                normed_intensity_plotter = PlotNormCurve(
                    self.summary_data_folder,
                    main_data_file,
                    self.metadata_file,
                    self.colour_map,
                    normalization_method=user_normalization,
                    require_log_transform=require_log)
                normed_intensity_plotter.run()

            else:
                # For before and after normalization
                normalization_method = "None"      
                raw_intensity_plotter = PlotNormCurve(
                    self.summary_data_folder,
                    main_data_file,
                    self.metadata_file,
                    self.colour_map,
                    normalization_method=normalization_method,
                    require_log_transform=True)
                raw_intensity_plotter.run()

            if not self.POC:
                if user_normalization != "None":
                    normed_intensity_plotter = PlotNormCurve(
                        self.summary_data_folder,
                        main_data_file,
                        self.metadata_file,
                        self.colour_map,
                        normalization_method=user_normalization,
                        require_log_transform=require_log)
                    normed_intensity_plotter.run()
            
            #UMAP plots
            if not self.POC:
                after_PCA_file = str(Path(self.summary_data_folder) / 'PCA_normed_intensities.html')
                after_PCA = PCAComparison(
                    data_file=DAPAR_output_file,
                    op_file= after_PCA_file,
                    title = "PCA Projection of Intensities After Normalization (Log2-transformed)",
                    col_suffix="normed",
                    conditions_file=self.metadata_file,
                    colour_map=self.colour_map,
                    )
                after_PCA.run()

                after_UMAP_file = str(Path(self.summary_data_folder) / 'UMAP_normed_intensities.html')
                after_UMAP = UMAPComparison(
                    data_file=DAPAR_output_file,
                    op_file= after_UMAP_file,
                    title = "UMAP Projection of Intensities After Normalization (Log2-transformed)",
                    col_suffix="normed",
                    conditions_file=self.metadata_file,
                    colour_map=self.colour_map,
                    )
                after_UMAP.run()

                #Hierarchical Clustering (after)
                after_cluster_file = str(Path(self.summary_data_folder) / 'Cluster_normed_intensities.png')
                after_cluster = HierarchicalCluster(
                    data_df=DAPAR_output_file,
                    op_file=after_cluster_file,
                    col_suffix="normed",
                    conditions_file=self.metadata_file,
                    colour_map=self.colour_map                
                    )
                after_cluster.run()

                raw_UMAP_file = str(Path(self.summary_data_folder) / 'UMAP_intensities.html')
                raw_UMAP = UMAPComparison(
                    data_file=DAPAR_output_file,
                    op_file= raw_UMAP_file,
                    title = "UMAP Projection of Intensities Before Normalization (Log2-transformed)",
                    col_suffix="log_transformed",
                    conditions_file=self.metadata_file,
                    colour_map=self.colour_map,
                    )
                #raw_UMAP.run()

                #Hierarchical Clustering (initial)
                raw_cluster_file = str(Path(self.summary_data_folder) / 'Cluster_raw_intensities.png')
                raw_cluster = HierarchicalCluster(
                    data_df=DAPAR_output_file,
                    op_file=raw_cluster_file,
                    col_suffix="log_transformed",
                    conditions_file=self.metadata_file,
                    colour_map=self.colour_map                
                    )
                raw_cluster.run()

            if not self.POC:
                col_suffix = "_with_imputed"
                require_untransform = True
            else:
                col_suffix = ""
                require_untransform = False
            
            if self.POIs:
                for protein in self.POIs:                

                    file_name=str(Path(self.POI_data_folder) / f'{protein.capitalize()}_relative_intensity.html')
                    title = f'{protein.capitalize()} relative intensity<br><sub>Following normalization and imputation; Reference: average across {self.reference_condition} samples</sub>'
                    y_title = f'Relative intensity (%)'
                                
                    plot_normed_POIs = POI_plotter(
                        protein=protein,
                        full_data=main_data_file,
                        reference_condition=self.reference_condition,
                        conditions_df=self.metadata_file,
                        sorted_conditions=self.sorted_conditions,
                        colour_map=self.colour_map,
                        file_name=file_name,
                        title=title,
                        y_title=y_title,
                        x_title="Sample",
                        col_suffix=col_suffix,
                        require_untransform=require_untransform
                        )
                    plot_normed_POIs.run()

            if not self.POC:
                str_log = str(logFC_thresh)
                if not str_log.lower() == "auto":
                    logFC_thresh = float(logFC_thresh)

                DA_finder=GetDA(
                    bkgd_file=DAPAR_output_file,
                    metadata_file=self.metadata_file,
                    op_folder=self.DA_data_folder,
                    table_op_folder=self.DA_table_folder,
                    comparison_columns=comparisons_list,
                    pval_cutoff=pval_thresh,
                    ppushprop_th=valid_values_pp_thresh,
                    filterprop_th=valid_values_filter_thresh,
                    condition_colour_map=self.colour_map,
                    DA_directions=sides,
                    set_ids=result_set_ids,
                    mark_contams=True,
                    Log2FC_cutoff=logFC_thresh,
                    slope = ECDF_slope_thresh,
                    extra_POIs=self.POIs,
                    POI_colour_map=self.POI_colour_map,
                    remove_sheet = None,
                    remove_list = None,
                    sheet = sheet_name,
                )

                comparison_sets=DA_finder.run_analysis()
                all_comparison_sets.update(comparison_sets)
                print(all_comparison_sets)

                # Prepare for overrepresentation analysis
                # Should be performed wth all IDd prots as bkgd
                # Note: Using all IDs rather than all quantified
                id_col = 'Genes'
                col_suffix='with_imputed'

                bkgd_df = pd.read_csv(DAPAR_output_file, sep='\t')

                print(bkgd_df)
                final_data_cols = [col for col in bkgd_df.columns if col_suffix in col]
                print(final_data_cols)
                bkgd_df.dropna(axis=0, how='any', subset=final_data_cols, inplace=True)
                print(bkgd_df)
                bkgd_ids = bkgd_df[id_col].values.tolist()
                print(bkgd_ids)
                bkgd_ids = [str(bg_id).split(';')[0] for bg_id in bkgd_ids]
                print(bkgd_ids)
                DA_genes = []

                for comparison,genes in comparison_sets.items():
                    print(comparison)
                    print(genes)
                    DA_genes.extend(genes)
                    string_gene_list = [str(gene_id).split(';', maxsplit=1)[0] for gene_id in genes]
                    # Having trouble communicating with String ATM
                    #op_file = str(Path(self.op_folder)/f'Analysis_{sheet_name}_{comparison}_network.html')
                    #self.plot_network(string_gene_list, bkgd_ids, op_file)
                    subdir = f'Analysis_{sheet_name}_{comparison}_ORA'
                    ORA_op_dir = Path(self.ORA_data_folder)/f'{subdir}'
                    ORA_table_op_dir = Path(self.ORA_table_folder)/f'{subdir}'
                    ORA_intermediate_op_dir = Path(self.ORA_intermediate_folder)/f'{subdir}'
                    for folder in [ORA_op_dir, ORA_table_op_dir, ORA_intermediate_op_dir]:
                        Path(folder).mkdir(parents=True, exist_ok=True)
                    analysis_group = comparison.replace('_',' ')
                    self.perform_ORA(string_gene_list, bkgd_ids, ORA_op_dir,
                        table_op_dir=ORA_table_op_dir, intermediate_op_dir=ORA_intermediate_op_dir,
                        analysis_group=analysis_group, taxon_id=self.species_id, gene_sets=self.gene_sets)
                
                for comparison in comparisons_list:
                    subdir = f'Analysis_{sheet_name}_{comparison}_GSEA'
                    GSEA_op_dir = Path(self.GSEA_data_folder)/f'{subdir}'
                    GSEA_table_op_dir = Path(self.GSEA_table_folder)/f'{subdir}'
                    GSEA_intermediate_op_dir = Path(self.GSEA_intermediate_folder)/f'{subdir}'
                    GSEA_table_op_dir.mkdir(parents=True, exist_ok=True)
                    GSEA_intermediate_op_dir.mkdir(parents=True, exist_ok=True)
                    GSEA_op_dir.mkdir(parents=True, exist_ok=True)
                    self.perform_GSEA_prerank(logFC_data_file=DAPAR_output_file, analysis_group=comparison,
                        op_dir=GSEA_op_dir, table_op_dir=GSEA_table_op_dir, intermediate_op_dir=GSEA_intermediate_op_dir,
                        taxon_id = 9606, gene_sets=self.gene_sets)

                # Collect all sign-regulated genes from DA analyses
                bkgd_file=DAPAR_output_file
                cluster_op = str(Path(self.cluster_table_folder)/f'Analysis_{sheet_name}_cluster_assignment.xlsx')
                cluster=KMeansClustering(
                    data=bkgd_file,
                    op_file=cluster_op,
                    op_dir=self.cluster_data_folder,
                    table_op_dir=self.cluster_table_folder,
                    conditions_file=self.metadata_file,
                    colour_map=self.colour_map,
                    prot_ids=DA_genes,
                    id_col=id_col,
                    col_suffix=col_suffix
                )
                cluster.run()

                cluster_df = pd.read_excel(cluster_op).set_index(id_col)
                cluster_dict = cluster_df.groupby('Cluster').apply(lambda x: list(x.index)).to_dict()

                for cluster_number,genes in cluster_dict.items():
                    string_gene_list = [str(gene_id).split(';', maxsplit=1)[0] for gene_id in genes]
                    #for go_label in ["molecular_function","biological_process","cellular_function"]:
                        #self.run_enrichment(genes, bkgd_ids, op_file, go_label=go_label)
                        #op_file = str(Path(self.op_folder)/f'Analysis_{sheet_name}_Cluster{cluster_number}_KClustering_{go_label}_enrichment.tsv')                        
                    #op_file = str(Path(self.op_folder)/f'Analysis_{sheet_name}_Cluster{cluster_number}_KClustering_STRING_GO_enrichment.tsv')
                    subfolder = f'Analysis_{sheet_name}_KCluster{cluster_number}_ORA'
                    ORA_op_dir = str(Path(self.ORA_data_folder)/f'{subfolder}')
                    ORA_table_op_dir = str(Path(self.ORA_table_folder)/f'{subfolder}')
                    ORA_intermediate_op_dir = str(Path(self.ORA_table_folder)/f'{subfolder}')
                    for folder in [ORA_op_dir, ORA_table_op_dir, ORA_intermediate_op_dir]:
                        Path(folder).mkdir(parents=True, exist_ok=True)
                    #self.plot_network(string_gene_list, bkgd_ids, op_file)
                    analysis_group = f"Cluster {cluster_number}"
                    self.perform_ORA(string_gene_list,
                        bkgd_ids, ORA_op_dir, table_op_dir=ORA_table_op_dir,
                        intermediate_op_dir=ORA_intermediate_op_dir, analysis_group=analysis_group, 
                        taxon_id=self.species_id, gene_sets=self.gene_sets)
        
            elif self.POC:
                self.plot_MA(main_data_file, comparisons_list)
                self.plot_Rank(main_data_file, comparisons_list)
                # Plot MA
                # Plot Rank

    def plot_network(self, string_gene_list, bkgd_ids, op_file):
        string_network = String_network_builder(
            gene_list = string_gene_list,
            background_gene_list=bkgd_ids,
            species = self.species_id,
            score = 400,
            fdr_enrichment_cutoff = 0.01,
            categories_of_interest = ["COMPARTMENTS", "Process", "Component", "Function", "TISSUES", "Keyword", "KEGG", "SMART", "InterPro", "PMID", "RCTM", "WikiPathways", "HPO", "NetworkNeighborAL"],
            show_n_enrichments = 5,
            colour_strategy = 'top_fdr',
            outfile=op_file
        )
        string_network.run()

    def perform_ORA(self,
        gene_list,
        bkd_ids,
        op_dir,
        analysis_group,
        table_op_dir=None,
        intermediate_op_dir=None,
        taxon_id=9606,
        gene_sets=None):

        enrichr = EnrichrORA(
            gene_list=gene_list,
            op_dir=op_dir,
            table_op_dir=table_op_dir,
            intermediate_op_dir=intermediate_op_dir,
            background=bkd_ids,
            analysis_group=analysis_group,
            taxon_id=taxon_id,
            gene_sets =gene_sets,
            )

        enrichr.run()
    
    def perform_GSEA_prerank(self, logFC_data_file, analysis_group,
            op_dir, table_op_dir=None, intermediate_op_dir=None,
            taxon_id = 9606,
            gene_sets = None,
            ):
        prerank = GSEAPy_Prerank(logFC_data_file, 
            analysis_group, 
            op_dir=op_dir,
            table_op_dir=table_op_dir,
            intermediate_op_dir=intermediate_op_dir,
            taxon_id=taxon_id,
            gene_sets=gene_sets,
            )
        prerank.run()

    
    def plot_MA(self, bkgd_file, comparison_columns):
        MA_plotter = PlotMA(
            bkgd_file=bkgd_file,
            metadata_file=self.metadata_file,
            op_folder=self.op_folder,
            comparison_columns=comparison_columns,
            condition_colour_map=self.colour_map,
            poi_dict=self.POI_colour_map,
            mark_contams=self.mark_contaminants,
            contaminants_fasta=self.contaminants_fasta
        )
        MA_plotter.run_analysis()

    def plot_Rank(self, bkgd_file, comparison_columns):
        rank_plotter = PlotRank(
            bkgd_file=bkgd_file,
            metadata_file=self.metadata_file,
            op_folder=self.op_folder,
            comparison_columns=comparison_columns,
            condition_colour_map=self.colour_map,
            poi_dict=self.POI_colour_map,
            mark_contams=self.mark_contaminants,
            contaminants_fasta=self.contaminants_fasta
        )
        rank_plotter.run_analysis()
