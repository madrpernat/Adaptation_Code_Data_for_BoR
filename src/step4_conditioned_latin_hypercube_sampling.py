import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

from src.utils.functions_library import variable_distribution_comparison_plot
from src.utils.clhs import clhs
from src.utils import ids

matplotlib.use('Qt5Agg')


def main():

    # Set working directory
    current_dir = os.path.dirname(__file__)
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    os.chdir(parent_dir)

    # Perform Uniform Conditioned Latin Hypercube Sampling to decrease the full-factorial SOW ensemble down to 500 SOWs

    ## Read in data
    full_factorial_sow_info = pd.read_csv(
        filepath_or_buffer='output/python_output/full_factorial_sow_info.csv'
    )
    cumulative_timeseries_array = np.load(
        file='output/python_output/full_factorial_cumulative_timeseries.npy'
    )
    cumulative_timeseries_array_scaled = pd.read_csv(
        filepath_or_buffer='output/python_output/full_factorial_cumulative_timeseries_scaled.csv',
        header=None
    ).to_numpy().T

    ## Create df of only the SOW characteristics you want to stratify. Including 'Neuron' as a characteristic ensures
    ## that there is equally sampling among the SOM's neurons.
    variables = [
        ids.MEDIAN_FLOW, ids.MAX_ANNUAL_FLOW, ids.MIN_ANNUAL_FLOW, ids.IQR_FLOW, ids.DRIEST_10_YEAR_FLOW,
        ids.WETTEST_10_YEAR_FLOW, ids.DEMAND, ids.INIT_STORAGE, ids.NEURON
    ]
    clhs_df = full_factorial_sow_info[variables]

    ## Generate uniform cLHS samples. The objective function consists of three components:
    ## 1) Spread of continuous variables,
    ## 2) Spread of the discrete variable ('Neuron' is the only discrete variable),
    ## 3) Preservation of correlation among variables.
    ## The 'weights' parameter adjusts the importance of these components.
    ## Due to the different magnitudes of variables and the higher number of continuous variables (see clhs.py for
    ## details), we experimented with various weightings. We compared histograms of variable distributions between the
    ## full factorial ensemble and the sampled ensemble. The weights below provided a uniform sampling across neurons
    ## and a relatively even distribution of continuous variables.
    clhs_samples = clhs(
        predictors=clhs_df,
        num_samples=500,
        random_state=3,
        max_iterations=len(full_factorial_sow_info),
        weights=[4, 5000, 1]
    )

    # Extract the info and cumulative timeseries of the 500 sampled SOWs
    sample_idx = clhs_samples['sample_indices'].tolist()

    sampled_sows_info = full_factorial_sow_info.loc[sample_idx].sort_index(axis=0)
    sampled_sows_cumulative_timeseries = cumulative_timeseries_array[:, sorted(sample_idx)]
    sampled_sows_cumulative_timeseries_scaled = cumulative_timeseries_array_scaled[:, sorted(sample_idx)]

    # We now have three dataframes describing the 500 SOW ensemble. One contains general info for each SOW, and the
    # others are the scaled and unscaled cumulative timeseries for each SOW. Save these data as csv.
    sampled_sows_info.to_csv(
        path_or_buf='output/python_output/500_sow_info.csv',
        index=False
    )
    pd.DataFrame(sampled_sows_cumulative_timeseries).to_csv(
        path_or_buf='output/python_output/500_sow_cumulative_timeseries.csv',
        index=False
    )
    pd.DataFrame(sampled_sows_cumulative_timeseries_scaled).to_csv(
        path_or_buf='output/python_output/500_sow_cumulative_timeseries_scaled.csv',
        index=False
    )

    # Create figure that compares variable distribution between full factorial ensemble and sampled SOWs
    fig, axes = plt.subplots(
        nrows=3,
        ncols=3,
        figsize=(19, 9.5)
    )

    for i, variable in enumerate(variables):

        if variable == ids.NEURON:
            plot_type = 'hist'
        else:
            plot_type = 'kde'

        variable_distribution_comparison_plot(
            df1=full_factorial_sow_info,
            df2=sampled_sows_info,
            color1='blue',
            color2='red',
            variable=variable,
            xlabel=variable,
            ax=axes[i//3, i % 3],
            plot_type=plot_type
        )

    # Figure adjustments
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.subplots_adjust(
        bottom=0.2,
        hspace=0.42,
        wspace=0.27
    )

    # Add custom legend
    legend_elements = [
        Patch(facecolor='blue', edgecolor='black', alpha=0.4, label='Full Factorial SOW Ensemble'),
        Patch(facecolor='red', edgecolor='black', alpha=0.4, label='500 SOW Ensemble')
    ]
    axes[2, 1].legend(
        handles=legend_elements,
        loc='lower center',
        ncol=2,
        bbox_to_anchor=(0.5, -0.9),
        frameon=False
    )

    fig.savefig(
        fname='output/python_output/figs/full_factorial_vs_500_ensemble_variable_distributions.png',
        dpi=400,
        bbox_inches='tight'
    )


if __name__ == '__main__':
    main()
