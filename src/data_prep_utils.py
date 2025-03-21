from pathlib import Path
import logging
import sys
sys.path.append(str(Path(__file__).resolve().parent.parent))
from src.common_classes import *

def perform_data_prep(spec_report, conditions_file, column_mapping, reattribute_peptides, peptides_to_remove, remove_outliers, min_pepts, log_file, report_data):
    PEPTIDES_TO_REMOVE = None
    MIN_PEPTS = None

    if peptides_to_remove:
        PEPTIDES_TO_REMOVE = peptides_to_remove.strip().split(',')
    if min_pepts:
        try:
            MIN_PEPTS = int(min_pepts)
            if MIN_PEPTS <= 0:
                MIN_PEPTS = 2
        except ValueError:
            MIN_PEPTS = 2

    prepare_data_obj = PrepareData(
        spec_report=spec_report,
        conditions_file=conditions_file,
        column_mapping=column_mapping,
        reattribute_peptides=reattribute_peptides,
        peptides_to_remove=PEPTIDES_TO_REMOVE,
        remove_outliers=remove_outliers,
        min_pepts=MIN_PEPTS,
        log_file=log_file,
        report_data=report_data
    )
    return prepare_data_obj.run()

def sort_data_and_assign_colours(conditions_df, reference_condition, colour_palette):
    condition_sorter = ConditionSort(conditions_df, map_strategy='as_is')
    sorted_conditions = condition_sorter.sort()
    if not reference_condition:
        reference_condition = sorted_conditions[0]
    colour_mapper = {condition: color for condition, color in zip(sorted_conditions, colour_palette)}
    POIs, POI_colour_map = condition_sorter.fetch_POI_data()
    return sorted_conditions, reference_condition, colour_mapper, POIs, POI_colour_map