# PUMA project
# October 2018

import os
import time
from process_input_data import build_config, load_json
from metabolite_prediction import predict_metabolite_from_config
from pathway_prediction import pathway_prediction_from_modelConfig
import sys

def main():
    
    # specify location of input testcase files, and location for output data
    base_dir = '/code'
    log_file_path = os.path.join(base_dir, 'output/results_CHO_cell')
    sys.stdout = open(log_file_path, 'w')

    os.environ['PUMA_INPUT_DATA'] = os.path.join(base_dir, 'CHO_cell_data/input_data')
    os.environ['PUMA_OUTPUT_DATA'] = os.path.join(base_dir, 'CHO_cell_data/output_data')


    # building json file
    model_config = build_config(serialize=True, landa=0.5, gamma=0.9, mu=0.5)
    
    # or reload the model configuration
    # model_config = load_json()

    # call pathway activity prediction
    time_start_pathway = time.time()
    samples = pathway_prediction_from_modelConfig(model_config)
    time_end_pathway = time.time()
    elapsed_time_pathway = time_end_pathway - time_start_pathway

    # call metabolite prediction
    time_start_metabolite = time.time()
    predict_metabolite_from_config(model_config, samples)
    time_end_metabolite = time.time()
    elapsed_time_metabolite = time_end_metabolite - time_start_metabolite


if __name__ == "__main__":
    main()
