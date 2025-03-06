import requests
import math
import networkx as nx
import plotly.graph_objects as go
import json
import pandas as pd
import textwrap
from itertools import combinations
import plotly.colors as pc

class String_network_builder:
    def __init__(self,
        gene_list,
        background_gene_list = None,
        species = 9606,
        score = 400,
        fdr_enrichment_cutoff = 0.05,
        categories_of_interest = ["COMPARTMENTS", "Process", "Component", "Function", "TISSUES", "Keyword", "KEGG", "SMART", "InterPro", "PMID", "RCTM", "WikiPathways", "HPO", "NetworkNeighborAL"]
,
        show_n_enrichments = 5,
        colour_strategy = 'top_fdr',
        filter_redundant = True,
        outfile='network.tsv'):

        self.gene_list = gene_list
        self.background_gene_list = background_gene_list
        self.categories_of_interest = categories_of_interest
        self.species = species
        self.score = score
        self.fdr_enrichment_cutoff = fdr_enrichment_cutoff
        self.categories_of_interest = categories_of_interest
        self.show_n_enrichments = show_n_enrichments
        self.colour_strategy = colour_strategy
        self.outfile = outfile

        self.data = None
        self.G = None
        self.enrichment_df = None
        self.filtered_enrichment_df = None
        # Import the color palettes
        self.pastel2_colours = [
        'rgb(241,226,204)',#7
         'rgb(255,242,174)',#6
         'rgb(179,226,205)', #1
         'rgb(253,205,172)', #2
         'rgb(203,213,232)', #3
         'rgb(244,202,228)', #4
         'rgb(230,245,201)', #5
         '#F6B7BD'
         ]
        self.set2_colours = [
         'rgb(229,196,148)',
         'rgb(255,217,47)',
         'rgb(102,194,165)',
          'rgb(252,141,98)',
          'rgb(141,160,203)',
          'rgb(231,138,195)',
          'rgb(166,216,84)',
          '#F07F8A'
          ]
        print(self.pastel2_colours)
        print(self.set2_colours)
        self.base_url = "https://version-12-0.string-db.org/api"
        #EB5160


    def get_string_data(self):
        output_format = "json"
        method = "network"
        
        # For POST requests, send the gene list as data instead of embedding in URL
        request_url = f"{self.base_url}/{output_format}/{method}"
        proteins = "\r".join(self.gene_list)
        
        request_data = {
            "identifiers": proteins,
            "species": self.species,
            "required_score": self.score,
            "show_query_node_labels": 1
        }
        
        enrich_method = "enrichment"
        enrich_request_url = f"{self.base_url}/{output_format}/{enrich_method}"
        enrich_data = {
            "identifiers": proteins,
            "species": self.species,
            "show_query_node_labels": 1
        }

        # Add background data if available
        if self.background_gene_list:
            bkgd_prots = "\r".join(self.get_background())
            request_data["background_string_identifiers"] = bkgd_prots
            enrich_data["background_string_identifiers"] = bkgd_prots

        # Perform POST requests with data instead of GET
        response = requests.post(request_url, data=request_data)
        response.raise_for_status()
        self.data = response.json()

        enrich_response = requests.post(enrich_request_url, data=enrich_data)
        enrich_response.raise_for_status()
        enrich_data = enrich_response.json()

        self.enrichment_df = pd.json_normalize(enrich_data)
        if self.enrichment_df is not None:
            enrichment_data_file = str(self.outfile)
            self.enrichment_df.to_csv(enrichment_data_file, sep="\t", index=False)

    def get_background(self):
        print("DOING IT")
        output_format = "json"
        request_url = f"https://string-db.org/api/{output_format}/get_string_ids"
        proteins = "\r".join(self.background_gene_list)
        
        request_data = {
            "identifiers": proteins,
            "species": self.species
        }

        response = requests.post(request_url, data=request_data)
        response.raise_for_status()
        data = response.json()
        response_df = pd.json_normalize(data)
        string_id_list = response_df['stringId'].values.tolist()
        
        return string_id_list

    def create_networkx_graph(self):
        G = nx.Graph()
        for protein in self.gene_list:
            G.add_node(protein)

        self.scores = []
        for interaction in self.data:
            protein1 = interaction['preferredName_A']
            protein2 = interaction['preferredName_B']
            score = interaction['score']
            G.add_edge(protein1, protein2, weight=score)
            self.scores.append(score)
        self.G = G
    
    def filter_enrichments_for_plotting(self):
        if not self.enrichment_df.empty:
            self.filtered_enrichment_df = self.enrichment_df.copy()
            self.filtered_enrichment_df.sort_values(by='fdr', ascending=True, inplace=True)
            if self.categories_of_interest:
                self.filtered_enrichment_df = self.filtered_enrichment_df[self.filtered_enrichment_df['category'].isin(self.categories_of_interest)]

                if self.colour_strategy == 'top_in_category':
                    category_list = self.filtered_enrichment_df['category'].unique().tolist()
                    cat_dfs = {category: self.filtered_enrichment_df[self.filtered_enrichment_df['category'] == category].sort_values(by='fdr', ascending=True) for category in category_list}
                    collected_terms = []

                    i=0
                    for category in cat_dfs:
                        while i < self.show_n_enrichments:
                            for row_index in range(self.show_n_enrichments):
                                row = cat_dfs[category].iloc[row_index]
                                collected_terms.append(row)
                                i+=1
                    print(i)
                            

                    # Create a DataFrame from the sorted rows
                    new_df = pd.DataFrame(collected_terms).head(self.show_n_enrichments)

                    self.filtered_enrichment_df = new_df
                else:
                    if self.colour_strategy == 'top_fdr':
                        self.filtered_enrichment_df= self.filtered_enrichment_df[self.filtered_enrichment_df['category'].isin(self.categories_of_interest)]
                        self.filtered_enrichment_df = self.filtered_enrichment_df.sort_values(by='fdr',ascending=True).head(self.show_n_enrichments)

            else:
                if self.colour_strategy == 'top_fdr':
                    self.filtered_enrichment_df = self.filtered_enrichment_df.sort_values(by='fdr',ascending=True).head(self.show_n_enrichments)
                elif self.colour_strategy == 'most_proteins':
                    # Extract all proteins in enrichments
                    def get_unique_proteins(rows):
                        proteins = set()
                        for row in rows:
                            proteins.update(row['preferredNames'])
                        return proteins

                    best_combination = None
                    max_proteins = 0
                    # Generate all possible combinations of rows
                    for rows in combinations(self.filtered_enrichment_df.to_dict('records'), self.show_n_enrichments):
                        unique_proteins = get_unique_proteins(rows)
                        if len(unique_proteins) > max_proteins:
                            max_proteins = len(unique_proteins)
                            best_combination = rows
                    
                    self.filtered_enrichment_df = pd.DataFrame(best_combination)

            # Expand the enrichment table so that there is an enrichment row for each protein
            # Redundancy filtering based on protein sets
            # Group by the set of proteins and keep the one with the lowest FDR

            print("HERE")
            print(len(self.filtered_enrichment_df))

            def get_protein_set(row):
                return frozenset(row['preferredNames'])

            # Add a column to identify protein sets
            self.filtered_enrichment_df['protein_set'] = self.filtered_enrichment_df.apply(get_protein_set, axis=1)

            # Sort by FDR and drop duplicates, keeping the lowest FDR for each protein set
            self.filtered_enrichment_df.sort_values(by='fdr', ascending=True, inplace=True)
            self.filtered_enrichment_df.drop_duplicates(subset=['protein_set'], keep='first', inplace=True)

            # Drop the temporary column used for filtering
            self.filtered_enrichment_df.drop(columns=['protein_set'], inplace=True)

            # Expand the enrichment table so that there is an enrichment row for each protein
            expanded_rows = []
            for _, row in self.filtered_enrichment_df.iterrows():
                proteins = row['preferredNames']
                for protein in proteins:
                    new_row = row.copy()
                    new_row['preferredName'] = protein
                    expanded_rows.append(new_row)

            self.filtered_enrichment_df = pd.DataFrame(expanded_rows)
            self.filtered_enrichment_df.sort_values(by='fdr', ascending=False, inplace=True)
            print(self.filtered_enrichment_df)
        
        else:
            print("No significant enrichments found")
           
    def plot_networkx_graph(self):
        # Identify nodes with edges and singletons
        nodes_with_edges = {n for edge in self.G.edges() for n in edge}
        singleton_nodes = list(set(self.G.nodes()) - nodes_with_edges)

        # Compute positions for nodes with edges
        pos={}
        
        pos_with_edges = nx.spring_layout(self.G.subgraph(nodes_with_edges), k=1)
        if len(pos_with_edges) > 0:
            pos = {**pos_with_edges}
        
        # Determine the boundaries of the connected network
            network_x = [x for x, y in pos_with_edges.values()]
            network_y = [y for x, y in pos_with_edges.values()]
            min_x, max_x = min(network_x), max(network_x)
            min_y, max_y = min(network_y), max(network_y)

            # Calculate width and height of the network
            network_width = max_x - min_x
            network_height = max_y - min_y

            # Position the grid relative to network boundaries
            grid_offset_x = min_x  # Offset to the left of the network
            grid_offset_y = min_y - (network_height / 2)  # Adjust grid vertically based on network height
        
        else:
            grid_offset_x = 0
            grid_offset_y = 0
        
        grid_spacing = 0.5  # Grid spacing factor

        # Define grid layout for singletons
        grid_size = math.ceil(math.sqrt(len(singleton_nodes)))
        singleton_pos = {}
        for idx, node in enumerate(singleton_nodes):
            x_grid = (idx % grid_size) * grid_spacing + grid_offset_x
            y_grid = -(idx // grid_size) * (grid_spacing*0.75) + grid_offset_y
            singleton_pos[node] = (x_grid, y_grid)

        # Combine both positions
        pos.update(singleton_pos)

        unique_scores = sorted(set(nx.get_edge_attributes(self.G, 'weight').values()))
        edge_traces = []

        for score in unique_scores:
            edge_x = []
            edge_y = []
            edge_customdata = []
            edge_mid_x = []
            edge_mid_y = []
            edge_hover_text = []

            for edge in self.G.edges(data=True):
                if edge[2]['weight'] == score:
                    x0, y0 = pos[edge[0]]
                    x1, y1 = pos[edge[1]]
                    # Calculate and store edge endpoints
                    edge_x.extend([x0, x1, None])
                    edge_y.extend([y0, y1, None])
                    # Calculate midpoint
                    mid_x = (x0 + x1) / 2
                    mid_y = (y0 + y1) / 2
                    edge_mid_x.append(mid_x)
                    edge_mid_y.append(mid_y)
                    # Hover text
                    url = f"https://string-db.org/cgi/network?networkId={edge[0]}.{edge[1]}"
                    edge_hover_text.append(url)

            if len(pos_with_edges) > 0:
            # Create a separate trace for edge hover
                edge_hover_trace = go.Scatter(
                    x=edge_mid_x,
                    y=edge_mid_y,
                    mode='markers',
                    marker=dict(size=0.1, color='rgba(0,0,0,0)'),
                    text=edge_hover_text,
                    showlegend=False,
                    hoverinfo='none',  # Disable hover info for edge line itself
                    customdata=edge_hover_text,  # Store URLs in customdata for click events
                )

                # Create the main edge trace
                edge_trace = go.Scatter(
                    x=edge_x, y=edge_y,
                    line=dict(width=score * 2, color='#888'),  # Scale weight for better visualization
                    mode='lines',
                    hoverinfo='none',  # Disable hover info for edge line itself
                    showlegend=False
                )

                edge_traces.append(edge_trace)
                edge_traces.append(edge_hover_trace)  # Add hover trace after processing all edges for this score

                node_x = []
                node_y = []
                node_texts = []
                for node in nodes_with_edges:
                    x, y = pos[node]
                    node_x.append(x)
                    node_y.append(y)
                    node_texts.append(f'{node}')

                node_trace = go.Scatter(
                    x=node_x, y=node_y,
                    mode='markers+text',
                    textposition='top center',
                    textfont=dict(color='black'),
                    hoverinfo='text',
                    text=node_texts,
                    marker=dict(
                        size=20,  # Size for connected nodes
                        color='#dfe0e2',
                        line=dict(width=3, color="#aaaaaa")
                    ),
                    showlegend=False
                )

        # Singleton trace with smaller markers
        singleton_x = [pos[node][0] for node in singleton_nodes]
        singleton_y = [pos[node][1] for node in singleton_nodes]
        singleton_texts = [f'{node}' for node in singleton_nodes]

        singleton_trace = go.Scatter(
            x=singleton_x, y=singleton_y,
            mode='markers+text',
            textposition='top center',
            textfont=dict(color='black'),
            hoverinfo='text',
            text=singleton_texts,
            marker=dict(
                size=20,  # Smaller size for singletons
                color='#dfe0e2',
                line=dict(width=3, color="#aaaaaa")
            ),
            showlegend=False
        )

        term_traces = []
        if self.filtered_enrichment_df is not None:
            unique_terms = self.filtered_enrichment_df['term'].unique().tolist() if len(self.filtered_enrichment_df) > 0 else None
        else:
            unique_terms = None 
        print(type(unique_terms))
        print(unique_terms)

        if isinstance(unique_terms, list):
            term_colors = {term: colour for term, colour in zip(unique_terms, self.pastel2_colours)}
            term_border_colours = {term: colour for term, colour in zip(unique_terms, self.set2_colours)}

            for term in unique_terms:
                desc = textwrap.shorten(self.filtered_enrichment_df[self.filtered_enrichment_df['term'] == term]['description'].values[0], width=40, placeholder="...")
                fdr = self.filtered_enrichment_df[self.filtered_enrichment_df['term'] == term]['fdr'].values[0]
                print(desc)
                term_x = []
                term_y = []
                term_texts = []
                for node in self.G.nodes():
                    if node in pos:
                        x, y = pos[node]
                        node_data = self.filtered_enrichment_df[(self.filtered_enrichment_df['preferredName'] == node) & (self.filtered_enrichment_df['term'] == term)]
                        if not node_data.empty:
                            term_x.append(x)
                            term_y.append(y)
                            term_texts.append(f'<b>ID: </b>{node}<br><b>Enrichment term: </b>{desc} ({term})<br><b>FDR:</b>{fdr:.2e}')

                term_trace = go.Scatter(
                    x=term_x,
                    y=term_y,
                    mode='markers',
                    marker=dict(
                        size=20,
                        color=term_colors[term],
                        line=dict(width=3, color=term_border_colours[term]),
                    ),
                    text=term_texts,
                    opacity=1,
                    hoverinfo='text',
                    name=f'{desc} ({term})<br><sup>FDR={fdr:.2e}',
                    legendgroup=f'{desc} ({term})<sup><br>FDR={fdr:.2e}'
                )
                term_traces.append(term_trace)
        if len(edge_traces):
            data_cols = edge_traces + [node_trace] + [singleton_trace] + term_traces
        else:
            data_cols = [singleton_trace] + term_traces
        fig = go.Figure(data=data_cols,
                    layout=go.Layout(
                        title='Protein Interaction Network with Enrichment Terms',
                        titlefont_size=16,
                        font=dict(color='black'),
                        plot_bgcolor='white',
                        showlegend=True,
                        hovermode='closest',
                        margin=dict(b=20, l=5, r=5, t=40),
                        xaxis=dict(showgrid=False, zeroline=False, visible=False),
                        yaxis=dict(showgrid=False, zeroline=False, visible=False),
                        uirevision='constant'
                    )
                    )
                    # Add this trace to your figure

     # Save the figure as HTML
        fig.write_html(self.outfile.replace('.tsv', '.html'), config={'toImageButtonOptions': {'format': 'png', 'filename': self.outfile.replace('.tsv', '.png'), 'scale': 2}})

    def run(self):
        self.get_string_data()
        self.create_networkx_graph()
        self.filter_enrichments_for_plotting()
        self.plot_networkx_graph()

# Example usage
gene_list = ['O60884', 'O95816', 'P05141', 'P11142', 'P12236', 'P17066', 'P31689', 'P34931', 'Q14244', 'Q3KQU3', 'Q92598', 'Q99615']
#GENE_LIST = ['A0A087WQ89', 'A2TJV2', 'B2RX12', 'O08914', 'O88338', 'P04370', 'P19258', 'P32114', 'P39876', 'P54869', 'Q3URU2', 'Q5NBX1', 'Q64462', 'Q6NZL0', 'Q80V42', 'Q8K370', 'Q91WQ9', 'Q91Y74', 'Q91Y97', 'Q99MZ6', 'Q9CWS0', 'Q9DD20', 'Q9JLB2', 'Q9Z0Z3']
GENE_LIST = ["Bap1", "Dock6", "Smc1b", "Entpd6", "Fkbp11", "Uimc1", "Cox19", "Phka2", "Maob", "Slc38a1", "Tns2", "Haspin", "Ibtk", "Tln2", "Aspa", "Hebp1", "Mgst1", "Adgrg2", "Sh2b1", "Psca", "Tmlhe", "Pol", "Fv4", "Brca1", "Bard1", "Ccnb2"]
#GENE_LIST = ["Q62504", "Q8C3F2", "P57096", "P68134", "P81122", "Q91VS7", "P97816", "Q99LJ7", "Q3U0P5", "Q99JZ0", "Q5U5Q9", "Q60855", "Q920F6", "Q61553", "Q91ZM2", "Q9Z0T9", "Q91ZE0", "Q67FY2", "Q8CJ12", "Q6ZPR6", "Q71LX4", "Q80W65", "Q8BW75", "Q8BWJ3", "P48754", "P47878", "Q8K0C8", "Q8CGB6", "Q9D1M7", "Q9D061", "P10400", "Q9R257", "P11370", "P22777", "Q8K2P7", "Q8K2W3", "Q99PU7", "P30276", "Q8VDR9", "Q9Z0R0", "A0A2I3BR81", "Q8R3P0", "O70445"]

SPECIES = 10090
CATEGORIES_OF_INTEREST = None
BACKGROUND_GENE_LIST = GENE_LIST
if __name__ == "__main__":
    string_network_builder = String_network_builder(
        gene_list= GENE_LIST,
        #background_gene_list=BACKGROUND_GENE_LIST,
        species= SPECIES,
        categories_of_interest=CATEGORIES_OF_INTEREST
        )
    string_network_builder.run()

