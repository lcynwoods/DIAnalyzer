

def sort_data_and_assign_colours(conditions_df, reference_condition, colour_palette):
    condition_sorter = ConditionSort(conditions_df, map_strategy='as_is')
    sorted_conditions = condition_sorter.sort()
    if not reference_condition:
        reference_condition = sorted_conditions[0]
    colour_mapper = {condition: color for condition, color in zip(sorted_conditions, colour_palette)}
    POIs, POI_colour_map = condition_sorter.fetch_POI_data()
    return sorted_conditions, reference_condition, colour_mapper, POIs, POI_colour_map