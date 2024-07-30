import distinctipy
import os
import pandas as pd

from src.utils.functions_library import create_cumulative_timeseries_som_view


# Following step 6, data/neuron_clusters.csv was manually made.


def main():

    # Set working directory
    current_dir = os.path.dirname(__file__)
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    os.chdir(parent_dir)

    five_hundred_sow_info = pd.read_csv(
        filepath_or_buffer='output/python_output/500_sow_info.csv'
    )
    five_hundred_sow_cumulative_timeseries = pd.read_csv(
        filepath_or_buffer='output/python_output/500_sow_cumulative_timeseries.csv'
    )
    neuron_cluster_mapping = pd.read_csv(
        filepath_or_buffer='data/neuron_clusters.csv'
    )

    # CUMULATIVE TIMESERIES VISUALIZATION COLORED BY CLUSTER

    cluster_colors = distinctipy.get_colors(12)  # The clustering configuration has 12 clusters
    neuron_colors = [cluster_colors[i] for i in neuron_cluster_mapping['Cluster']]

    fig = create_cumulative_timeseries_som_view(
        five_hundred_sow_info=five_hundred_sow_info,
        five_hundred_sow_cumulative_timeseries=five_hundred_sow_cumulative_timeseries,
        colors=neuron_colors
    )

    fig.savefig(
        fname='output/python_output/figs/hclusters_som_cumulative_timeseries_view.png',
        dpi=400,
        bbox_inches='tight'
    )


if __name__ == '__main__':
    main()
