import os
import pandas as pd

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
    experiments = pd.DataFrame(
        {ids.EXPERIMENT: 1, ids.CLUSTER: 0},
        {ids.EXPERIMENT: 2, ids.CLUSTER: 6},
        {ids.EXPERIMENT: 3, ids.CLUSTER: 7},
        {ids.EXPERIMENT: 4, ids.CLUSTER: 11}
    )

    n_sows = 8

    for i, row in experiments.iterrows():

        experiment_dir = output_dir + f'experiment_{row[ids.EXPERIMENT]}'
        os.makedirs(experiment_dir)

        sows = five_hundred_sow_info[
            five_hundred_sow_info[ids.CLUSTER] == row[ids.CLUSTER]
        ].sample(n_sows).sort_values(by=ids.SOW)

        # Save csv with SOW info
        sows.to_csv(
            path_or_buf=experiment_dir,
            index=False
        )

        # Look to the 500 SOW ensemble CRSS files, find the files for the experiment's sampled SOWs, and copy over to
        # the experiment_dir. You can then copy and paste directly into the borg directory.
        sow_ids = sows[ids.SOW]

        flow_dir = experiment_dir + '/FlowInput'
        scd_dir = experiment_dir + '/SystemConditionInput'
        os.makedirs(flow_dir)
        os.makedirs(scd_dir)

        for j in range(len(sows)):
            sow_flow_dir = flow_dir + f'/trace{j + 1}'
            sow_scd_dir = scd_dir + f'/trace{j + 1}'
            os.makedirs(sow_flow_dir)
            os.makedirs(sow_scd_dir)



    sow_exp1 = five_hundred_sow_info[five_hundred_sow_info[ids.CLUSTER] == 0].sample(8).sort_values(by=ids.SOW)
    sow_exp2 = five_hundred_sow_info[five_hundred_sow_info[ids.CLUSTER] == 6].sample(8).sort_values(by=ids.SOW)
    sow_exp3 = five_hundred_sow_info[five_hundred_sow_info[ids.CLUSTER] == 8].sample(8).sort_values(by=ids.SOW)
    sow_exp4 = five_hundred_sow_info[five_hundred_sow_info[ids.CLUSTER] == 11].sample(8).sort_values(by=ids.SOW)

    sow_exp1.to_csv(
        path_or_buf='output/python_output/sow_ensembles/experiment_1.csv',
        index=False
    )
    sow_exp2.to_csv(
        path_or_buf='output/python_output/sow_ensembles/experiment_2.csv',
        index=False
    )
    sow_exp3.to_csv(
        path_or_buf='output/python_output/sow_ensembles/experiment_3.csv',
        index=False
    )
    sow_exp4.to_csv(
        path_or_buf='output/python_output/sow_ensembles/experiment_4.csv',
        index=False
    )


if __name__ == '__main__':
    main()
