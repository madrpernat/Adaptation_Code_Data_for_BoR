import matplotlib
import os
import pandas as pd

from src.utils import ids, plot_configs
from src.utils.functions_library import create_cumulative_timeseries_som_view, create_som_characteristic_view

matplotlib.use('Qt5Agg')


def main():

    # Set working directory
    current_dir = os.path.dirname(__file__)
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    os.chdir(parent_dir)

    # Create figures to visualize and interpret the SOM. Figures saved to python_output/figs/.

    five_hundred_sow_info = pd.read_csv(
        filepath_or_buffer='output/python_output/500_sow_info.csv'
    )
    five_hundred_sow_cumulative_timeseries = pd.read_csv(
        filepath_or_buffer='output/python_output/500_sow_cumulative_timeseries.csv'
    )
    som_neuron_coordinates = pd.read_csv(
        filepath_or_buffer='output/r_output/som_neuron_coordinates.csv'
    )

    # CUMULATIVE TIMESERIES VISUALIZATION

    fig = create_cumulative_timeseries_som_view(
        five_hundred_sow_info=five_hundred_sow_info,
        five_hundred_sow_cumulative_timeseries=five_hundred_sow_cumulative_timeseries,
        colors=['blue'] * 52
    )

    fig.savefig(
        fname='output/python_output/figs/som_cumulative_timeseries_view.png',
        dpi=400,
        bbox_inches='tight'
    )

    # SINGULAR CHARACTERISTIC VIEWS OF SOM

    characteristics = [
        ids.MEDIAN_FLOW, ids.MAX_ANNUAL_FLOW, ids.MIN_ANNUAL_FLOW, ids.IQR_FLOW, ids.DRIEST_10_YEAR_FLOW,
        ids.WETTEST_10_YEAR_FLOW, ids.DEMAND, ids.INIT_STORAGE
    ]

    for characteristic in characteristics:

        ## Calculate average characteristic value for each neuron
        avg_value_per_neuron = five_hundred_sow_info.groupby(ids.NEURON)[characteristic].mean().tolist()

        fig = create_som_characteristic_view(
            characteristic_title=plot_configs.PLOT_TITLES.get(characteristic),
            neuron_values=avg_value_per_neuron,
            neuron_coordinates=som_neuron_coordinates,
            color_scheme=plot_configs.COLOR_SCHEMES.get(characteristic),
            inverse_colorbar=plot_configs.INVERSE_COLORBAR.get(characteristic)
        )

        fig.savefig(
            fname='output/python_output/figs/som_characteristic_view_' + characteristic + '.png',
            dpi=400,
            bbox_inches='tight'
        )


if __name__ == '__main__':
    main()
