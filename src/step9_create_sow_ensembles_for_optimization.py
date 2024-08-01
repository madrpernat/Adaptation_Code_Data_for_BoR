import os
import pandas as pd
import shutil

from src.utils import ids


def main():

    # Set working directory
    current_dir = os.path.dirname(__file__)
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    os.chdir(parent_dir)

    five_hundred_sow_info = pd.read_csv(
        filepath_or_buffer='output/python_output/500_sow_info.csv'
    )

    src_dir = 'output/r_output/crss_input_files'
    output_dir = 'output/python_output/sow_ensembles_for_optimization/'

    # Experiments
    # 1) Random sample cluster 0
    # 2) Random sample cluster 6
    # 3) Random sample cluster 8
    # 4) Random sample cluster 11
    # Add other experiments if desired
    experiments = pd.DataFrame([
        {ids.EXPERIMENT: 1, ids.CLUSTER: 0},
        {ids.EXPERIMENT: 2, ids.CLUSTER: 6},
        {ids.EXPERIMENT: 3, ids.CLUSTER: 7},
        {ids.EXPERIMENT: 4, ids.CLUSTER: 11}
    ])

    n_sows = 8

    for i, row in experiments.iterrows():

        # Create a directory for the experiment
        experiment_dir = output_dir + f'experiment_{row[ids.EXPERIMENT]}'
        os.makedirs(experiment_dir, exist_ok=True)

        # Randomly sample sows (# = n_sows) from cluster
        sows = five_hundred_sow_info[
            five_hundred_sow_info[ids.CLUSTER] == row[ids.CLUSTER]
        ].sample(n_sows).sort_values(by=ids.SOW)

        # Save csv with SOW info to experiment directory
        sows.to_csv(
            path_or_buf=experiment_dir + f'/experiment_{row[ids.EXPERIMENT]}_sows.csv',
            index=False
        )

        # Within the experiment directory, create subdirectories for FlowInput and SystemConditionInput. These
        # subdirectories can be copied and pasted into the experiment's borg run directory.
        sow_ids = sows[ids.SOW].tolist()

        destination_flow_dir = experiment_dir + '/FlowInput'
        destination_scd_dir = experiment_dir + '/SystemConditionInput'
        os.makedirs(destination_flow_dir)
        os.makedirs(destination_scd_dir)

        for j in range(len(sows)):

            sow = sow_ids[j]

            sow_flow_dir = destination_flow_dir + f'/trace{j + 1}'
            sow_scd_dir = destination_scd_dir + f'/trace{j + 1}'
            os.makedirs(sow_flow_dir)
            os.makedirs(sow_scd_dir)

            src_flow_dir = src_dir + f'/FlowInput/trace{sow}'
            src_scd_dir = src_dir + f'/SystemConditionInput/trace{sow}'

            for filename in os.listdir(src_flow_dir):
                # Construct full file paths for source and destination
                src_file = os.path.join(src_flow_dir, filename)
                destination_file = os.path.join(sow_flow_dir, filename)

                shutil.copy(src_file, destination_file)

            for filename in os.listdir(src_scd_dir):
                src_file = os.path.join(src_scd_dir, filename)
                destination_file = os.path.join(sow_scd_dir, filename)
                shutil.copy(src_file, destination_file)


if __name__ == '__main__':
    main()
