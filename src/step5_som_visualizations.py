import matplotlib
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

from src.utils import ids, plot_configs
from src.utils.functions_library import (get_subplot_coordinates, plot_cumulative_timeseries,
                                         plot_som_characteristic_view)

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
    sim_years = np.arange(2027, 2057)

    fig = plt.figure(figsize=(19, 9.5))
    gs = GridSpec(
        nrows=4,
        ncols=28,  # double number of columns to allow for shifting of rows
        figure=fig
    )

    ## Plotting loop
    neurons = np.unique(five_hundred_sow_info[ids.NEURON])

    for neuron in neurons:
        x, y = get_subplot_coordinates(neuron=neuron)
        ax = fig.add_subplot(gs[x, y:y + 2])

        ## Determine which SOWs are in neuron and which are not
        indices = five_hundred_sow_info.index[five_hundred_sow_info["Neuron"] == neuron]
        non_indices = five_hundred_sow_info.index[five_hundred_sow_info["Neuron"] != neuron]

        plot_cumulative_timeseries(
            sim_years=sim_years,
            cumulative_timeseries=five_hundred_sow_cumulative_timeseries,
            indices=indices,
            non_indices=non_indices,
            ax=ax,
            color='blue'
        )
        ax.set_title(
            label="Neuron " + str(neuron),
            y=1.0,
            pad=-14,
            fontsize=11,
            fontweight="bold"
        )

        ## For plots on the bottom, add x-axis tick labels
        if x == 3:
            ax.tick_params(
                labelbottom=True,
                bottom=True,
                length=0
            )
            ax.set_xticks(
                ticks=[2027, 2042, 2056]
            )
            ax.set_xticklabels(
                labels=[2027, 2042, 2056],
                fontsize=13,
                rotation=45
            )
            ax.set_xlabel(
                xlabel="Year",
                fontsize=16
            )

        ## For plots on the left, add y-axis tick labels
        if y == 0 or y == 1:
            ax.tick_params(
                left=True,
                labelleft=True,
                length=0
            )
            ax.set_yticks(
                ticks=[0, 200, 400]
            )
            ax.set_yticklabels(
                labels=[0, 200, 400],
                fontsize=13
            )
            ax.set_ylabel(
                ylabel="Cumulative\nTimeseries (MAF)",
                fontsize=14
            )

    fig.subplots_adjust(
        wspace=0.6,
        hspace=0.1,
        left=0.05,
        right=0.95,
        top=0.95,
        bottom=0.07
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

        fig = plot_som_characteristic_view(
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
