# Functions library
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
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
    """Creates overlaid plots (kde or histogram) for a specified variable from two different datasets"""

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


def get_subplot_coordinates(
        neuron: int,
        max_neurons: int = 52,
        rows: int = 4,
        cols: int = 13
) -> [int, int]:

    """Calculate the row and column positions for a subplot based on the neuron index in a Self-Organizing Map"""

    scale = rows / 4

    x = int(np.floor((max_neurons - neuron) / cols))
    y = int(((neuron - 1) % cols)) * 2

    if x == 1 or x == 3:  # determines if the row should be shifted, if so, add 1 to implement the shift
        y += 1

    return int(x * scale), int(y)


def plot_cumulative_timeseries(
        sim_years: np.array,
        cumulative_timeseries: pd.DataFrame,
        indices: [int],
        non_indices: [int],
        ax: matplotlib.axes.Axes,
        color: str
):
    """
    Adds many cumulative timeseries onto a single plot. Timeseries in 'non_indices' (i.e., those not in the neuron being
    plotted) are colored gray and timeseries in 'indices' (i.e., those in the neuron being plotted) are colored in the
    specified color.
    """

    # Plot non_indices timeseries
    for idx in non_indices:
        ax.plot(
            sim_years,
            cumulative_timeseries[str(idx)],
            color='gray',
            linewidth=0.5
        )

    # Plot indices timeseries
    for idx in indices:
        ax.plot(
            sim_years,
            cumulative_timeseries[str(idx)],
            color=color,
            linewidth=0.5
        )

    # Set axis limits and tick params
    ax.set_xlim(
        left=sim_years[0],
        right=sim_years[-1]
    )
    ax.set_ylim(
        bottom=0,
        top=470
    )
    ax.tick_params(
        left=False,
        right=False,
        labelleft=False,
        labelbottom=False,
        bottom=False
    )


def plot_som_characteristic_view(
        characteristic_title: str,
        neuron_values: [float],
        neuron_coordinates: pd.DataFrame,
        color_scheme: str,
        inverse_colorbar: bool
):

    value_range = (
        np.floor(np.min(neuron_values)),
        np.ceil(np.max(neuron_values))
    )

    # Create figure object
    fig, ax = plt.subplots(figsize=(19, 9.5))
    ax.set_aspect('equal')

    # Set figure parameters
    grid_width, grid_height = 13, 4
    hex_radius = 0.52

    # Set color scheme
    colors = get_color_scheme(color_scheme)
    cmap = mcolors.LinearSegmentedColormap.from_list(name='cmap', colors=colors)

    for idx, row in neuron_coordinates.iterrows():

        value = neuron_values[idx]
        neuron_face_color = map_value_to_color(
            value=value,
            value_range=value_range,
            cmap=cmap
        )

        neuron_text_color = get_text_color_for_neuron(neuron_face_color)

        x, y = row['x'], row['y']

        hexagon = patches.RegularPolygon(
            xy=(x, y),
            numVertices=6,
            radius=hex_radius,
            orientation=np.radians(0),
            facecolor=neuron_face_color,
            edgecolor='black'
        )
        ax.add_patch(hexagon)

        ax.text(
            x=x,
            y=y,
            s=round(value, 1),
            ha='center',
            va='center',
            fontsize=27,
            color=neuron_text_color
        )

    # Colorbar
    norm = plt.Normalize(
        vmax=value_range[1],
        vmin=value_range[0]
    )
    sm = plt.cm.ScalarMappable(
        cmap=cmap,
        norm=norm
    )
    sm.set_array([])

    ## Add colorbar to fig
    cbar_ax = fig.add_axes([0.11, 0.2, 0.775, 0.03])
    cbar = fig.colorbar(
        mappable=sm,
        cax=cbar_ax,
        orientation='horizontal'
    )
    cbar.set_label(
        label='Color/Text: Neuron Average ' + characteristic_title,
        fontsize=26,
        labelpad=10
    )
    cbar.ax.tick_params(labelsize=24)

    ## Inverse colorbar to have the largest value on the left
    if inverse_colorbar:
        cbar.ax.invert_xaxis()

    # Set limits and turn off axes
    #ax.set_title('SOW Self-Organizing Map', fontsize=35, pad=0.1)
    ax.set_xlim([0, 2 * hex_radius * grid_width + 2])
    ax.set_ylim([0, 2 * hex_radius * grid_height + 1])
    ax.axis('off')

    fig.subplots_adjust(
        left=0.048,
        right=0.997,
        top=0.965,
        bottom=0.151
    )
    fig.suptitle('SOW Self-Organizing Map', fontsize=35, y=0.835)
    return fig


def get_color_scheme(scheme_name):
    schemes = {"wet_dry": [(0, '#8B0000'), (0.25, '#FFD700'), (0.5, '#008000'), (0.75, "#add8e6"), (1, "#0000FF")],
               "warm_cool": [(0, '#FF4500'), (0.2, '#FFA500'), (0.4, '#FFFF00'), (0.6, '#FFC0CB'), (0.8, "#40E0D0"),
                             (1, "#800080")]}

    return schemes.get(scheme_name, None)


def map_value_to_color(value, value_range, cmap):
    # Normalize the value to the range [0, 1]
    normalized_value = (value - value_range[0]) / (value_range[1] - value_range[0])
    # Map the normalized value to a color in the gradient
    return cmap(normalized_value)


def get_text_color_for_neuron(rgba):
    """
    Determines whether the text color should be white or black based on the luminance of the background color.

    Parameters:
    rgba (tuple): 4-value tuple representing the RGBA values.

    Returns:
    str: 'black' if the background color is light, 'white' if the background color is dark.
    """
    # Convert RGBA values from 0-1 range to 0-255 range
    r, g, b, a = rgba
    r = int(r * 255)
    g = int(g * 255)
    b = int(b * 255)

    # Calculate luminance
    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255

    return 'black' if luminance > 0.28 else 'white'


