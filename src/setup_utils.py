import os
import logging
from pathlib import Path

def check_inputs(spec_report, conditions_file, analyses_file, op_folder, r_exe_path, r_DAPAR_path, r_rrvgo_path, species_id, gh_token):
    inputs = {
        'spec_report': spec_report,
        'conditions_file': conditions_file,
        'analyses_file': analyses_file,
        'op_folder': op_folder,
        'r_exe_path': r_exe_path,
        'r_DAPAR_path': r_DAPAR_path,
        'r_rrvgo_path': r_rrvgo_path
    }
    not_found = []
    empty_paths = []
    invalid_extensions = []

    for name, path in inputs.items():
        if not path:
            empty_paths.append(name)
        elif not Path(path).exists():
            not_found.append(path)
        elif name in ['spec_report', 'conditions_file', 'analyses_file'] and not path.endswith(('.xlsx', '.parquet', '.csv', '.tsv')):
            invalid_extensions.append(path)
    
    if empty_paths:
        empty_paths_str = ', '.join(empty_paths)
        print(f'The following inputs have empty paths: {empty_paths_str}')
        logging.error(f'The following inputs have empty paths: {empty_paths_str}')
    
    if not_found:
        not_found_str = ', '.join(not_found)
        print(f'The following inputs do not exist: {not_found_str}')
        logging.error(f'The following inputs do not exist: {not_found_str}')
    
    if invalid_extensions:
        invalid_extensions_str = ', '.join(invalid_extensions)
        print(f'The following inputs have invalid extensions: {invalid_extensions_str}')
        logging.error(f'The following inputs have invalid extensions: {invalid_extensions_str}')
    
    if empty_paths or not_found or invalid_extensions:
        raise ValueError('Invalid input files. Please check the log for details.')

    # Check for valid species ID
    if species_id not in [9606, 10090]:
        print(f'Invalid species ID: {species_id}')
        logging.error(f'Invalid species ID: {species_id}')
        raise ValueError(f'Invalid species ID: {species_id}')

    # Check for required environment variables
    if not gh_token:
        print('GitHub token (GH_TOKEN) is not set.')
        logging.error('GitHub token (GH_TOKEN) is not set.')
        raise ValueError('GitHub token (GH_TOKEN) is not set.')

def create_folders_and_logger(op_folder, log_file=None):
    if not Path(op_folder).exists():
        Path(op_folder).mkdir(parents=True)
    job_name = Path(op_folder).stem
    plot_folder = Path(op_folder) / f'{job_name}_data' / 'plots'
    summary_data_folder = plot_folder / 'summary_data'
    for folder in [plot_folder, summary_data_folder]:
        Path(folder).mkdir(parents=True, exist_ok=True)

    if not log_file:
        log_file = Path(op_folder) / 'log.log'
    logging.basicConfig(filename=log_file,
                        level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        force=True)