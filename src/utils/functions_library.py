# Functions library
import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns

from src.utils import ids


def calculate_trace_metrics(
        traces_timeseries: pd.DataFrame,
        n_traces: int,
        drought_threshold: float
) -> pd.DataFrame:

    """Calculates various hydrological metrics for each trace in traces_timeseries"""

    metrics = {
        "Median": [], "Max": [], "Min": [], "SD": [],
        "IQR": [], "MaxDroughtLength": [], "Driest10yr": [], "Wettest10yr": []
    }

    for trace_id in range(1, n_traces + 1):
        trace_data = traces_timeseries[traces_timeseries[ids.TRACE] == trace_id][ids.LF_ANNUAL].to_numpy()

        metrics["Median"].append(np.median(trace_data))
        metrics["Max"].append(np.max(trace_data))
        metrics["Min"].append(np.min(trace_data))
        metrics["SD"].append(np.std(trace_data))
        metrics["IQR"].append(np.subtract(*np.percentile(trace_data, [75, 25])))
        metrics["MaxDroughtLength"].append(calculate_max_drought_length(trace_data, drought_threshold))

        driest_10yr, wettest_10yr = calculate_driest_and_wettest_decade_avg(trace_data)
        metrics["Driest10yr"].append(driest_10yr)
        metrics["Wettest10yr"].append(wettest_10yr)

    return pd.DataFrame(metrics)


def calculate_max_drought_length(
        trace_data: np.ndarray,
        drought_threshold: float
) -> int:

    """Calculates the maximum drought length for a given trace"""

    drought_lengths = []
    current_drought_length = 0

    for value in trace_data:
        if value < drought_threshold:
            current_drought_length += 1
        else:
            if current_drought_length > 0:
                drought_lengths.append(current_drought_length)
                current_drought_length = 0

    # Add the last drought length if the series ends in a drought
    if current_drought_length > 0:
        drought_lengths.append(current_drought_length)

    return max(drought_lengths, default=0)


def calculate_driest_and_wettest_decade_avg(
        trace_data: np.ndarray
) -> [float, float]:

    """Calculates the driest and wettest 10-year flow averages"""

    if len(trace_data) < 10:
        raise ValueError("Trace data must be at least 10 years long to calculate 10-year periods.")

    # Initialize values
    driest_10yr = float('inf')
    wettest_10yr = float('-inf')

    # Iterate over the trace data to calculate sums of 10-year periods
    for i in range(len(trace_data) - 9):

        ten_year_sum = np.sum(trace_data[i:i + 10])

        # Update driest 10-yr sum if the current sum is smaller
        driest_10yr = min(driest_10yr, ten_year_sum)

        # Update wettest 10-yr sum if the current sum is smaller
        wettest_10yr = max(wettest_10yr, ten_year_sum)

    return driest_10yr / 10, wettest_10yr / 10


def variable_distribution_comparison_plot(
        df1: pd.DataFrame,
        df2: pd.DataFrame,
        color1: str,
        color2: str,
        variable: str,
        xlabel: str,
        ax: matplotlib.axes.Axes,
        plot_type: str = 'kde'
):
    """Overlaid plots (kde or histogram) for a specified variable from two different datasets"""

    plot_funcs = {
        'kde': sns.kdeplot,
        'hist': sns.histplot
    }
    common_params = {
        'x': variable,
        'ax': ax,
        'alpha': 0.4,
    }

    # Set style of plot
    sns.set_theme(
        style='whitegrid',
        rc={'axes.facecolor': '.9', 'grid.color': '.8'}
    )
    sns.set_context(
        context='talk',
        font_scale=1
    )

    # Plot the data
    if plot_type == 'kde':
        plot_funcs['kde'](
            data=df1,
            fill=True,
            color=color1,
            **common_params
        )
        plot_funcs['kde'](
            data=df2,
            fill=True,
            color=color2,
            **common_params
        )
    elif plot_type == 'hist':

        n_bins = len(np.unique(df1[variable]))

        plot_funcs['hist'](
            data=df1,
            color=color1,
            kde=False,
            bins=n_bins,
            stat='density',
            **common_params
        )
        plot_funcs['hist'](
            data=df2,
            color=color2,
            kde=False,
            bins=n_bins,
            stat='density',
            **common_params
        )
    else:
        raise ValueError("Invalid plot type. Use 'kde' or 'hist'.")

    # Add labels and titles
    ax.set_xlabel(xlabel, fontsize=17)
    ax.set_ylabel('Density', fontsize=17)

    # Enlarge axis tick labels
    ax.tick_params(axis='both', which='major', labelsize=14.5)
