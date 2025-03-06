import pandas as pd
import polars as pl
import plotly.express as px
import colorsys
import re

class ConditionSort():
    def __init__(self, conditions_df, map_strategy='as_is'):
        self.conditions_df = conditions_df
        self.map_strategy = map_strategy
        
        self.sorted_conditions = None
        print(self.map_strategy)
        print(type(self.conditions_df))
        if isinstance(self.conditions_df, pl.DataFrame):
            self.conditions_df = self.conditions_df.to_pandas()

    def sort_as_is(self):
        self.sorted_conditions = self.conditions_df['Condition'].unique()
        print(self.sorted_conditions)

    def sort_alphabetical(self):
        self.sorted_conditions = self.conditions_df.sort_values(by='Condition', inplace=True)
    
    def sort_controls_first(self):
        self.conditions_df.sort_values(by='Condition', inplace=True)
        self.conditions_df['sort_key'] = self.conditions_df['Bait_Ctrl'].apply(lambda x: 0 if x == 'Ctrl' else 1)
        self.conditions_df.sort_values(by='sort_key', inplace=True)
        self.sorted_conditions = self.conditions_df['Condition'].unique()
    
    def sort_controls_last(self):
        self.conditions_df.sort_values(by='Condition', inplace=True)
        self.conditions_df['sort_key'] = self.conditions_df['Bait_Ctrl'].apply(lambda x: 1 if x == 'Ctrl' else 0)
        self.conditions_df.sort_values(by='sort_key', inplace=True)
        self.sorted_conditions = self.conditions_df['Condition'].unique()
    
    def sort(self):
        if self.map_strategy == 'as_is':
            self.sort_as_is()
        elif self.map_strategy == 'alphabetical':
            self.sort_alphabetical()
        elif self.map_strategy == 'controls_first':
            self.sort_controls_first()
        elif self.map_strategy == 'controls_last':
            self.sort_controls_last()

        return self.sorted_conditions
    
    def fetch_POI_data(self):
        if (any(value != "None" for value in self.conditions_df['Bait']) or
            any(value != "None" for value in self.conditions_df['Marker']) or
            any(value != "None" for value in self.conditions_df['POIs'])):

            def get_proteins_of_interest(row):
                proteins = []
                if pd.notna(row['Bait']) and row['Bait'] != "None":
                    proteins.extend(row['Bait'].split(';'))
                if pd.notna(row['Marker']) and row['Marker'] != "None":
                    proteins.extend(row['Marker'].split(';'))
                if pd.notna(row['POIs']) and row['POIs'] != "None":
                    proteins.extend(row['POIs'].split(';'))
                return proteins

            proteins_of_interest = self.conditions_df.apply(get_proteins_of_interest, axis=1).explode().dropna().unique()
            proteins_of_interest = [prot.upper() for prot in proteins_of_interest]

            # Start with a large predefined color palette
            base_colors = (
                px.colors.qualitative.Dark2 +
                px.colors.qualitative.Set3 +
                px.colors.qualitative.Bold +
                px.colors.qualitative.Pastel +
                px.colors.qualitative.Vivid
            )

            def convert_to_hex(color):
                """Convert Plotly RGB or RGBA color to hex format."""
                if color.startswith("#"):  # Already hex
                    return color
                match = re.match(r"rgba?\((\d+),\s*(\d+),\s*(\d+)", color)
                if match:
                    r, g, b = map(int, match.groups())
                    return f'#{r:02x}{g:02x}{b:02x}'
                return "#000000"  # Default to black if invalid color

            # Convert all colors to hex format
            base_colors = [convert_to_hex(c) for c in base_colors]

            # Function to adjust brightness in the HSL space
            def adjust_brightness(color, factor):
                """Darken or lighten a hex color using HSL conversion."""
                color = color.lstrip('#')
                r, g, b = tuple(int(color[i:i+2], 16) for i in (0, 2, 4))
                h, l, s = colorsys.rgb_to_hls(r/255, g/255, b/255)
                new_l = min(1, max(0, l * factor))  # Ensure within 0-1 range
                new_r, new_g, new_b = colorsys.hls_to_rgb(h, new_l, s)
                return f'#{int(new_r*255):02x}{int(new_g*255):02x}{int(new_b*255):02x}'

            # Generate the color map
            protein_colour_map = {}
            num_colors = len(base_colors)
            for i, protein in enumerate(proteins_of_interest):
                base_color = base_colors[i % num_colors]  # Cycle through base colors
                brightness_factor = 1 - ((i // num_colors) * 0.15)  # Adjust brightness every full cycle
                protein_colour_map[protein] = adjust_brightness(base_color, brightness_factor)

            return proteins_of_interest, protein_colour_map
        else:
            return None, None