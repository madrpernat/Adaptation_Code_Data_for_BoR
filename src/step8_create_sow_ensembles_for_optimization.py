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

    # Experiments
    # 1) Random sample cluster 0
    # 2) Random sample cluster 6
    # 3) Random sample cluster 8
    # 4) Random sample cluster 11

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
