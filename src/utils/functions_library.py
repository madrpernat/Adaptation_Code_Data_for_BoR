from typing import List, Tuple, Union

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram

from src.utils import ids


def calculate_trace_metrics(
        traces_timeseries: pd.DataFrame,
        drought_threshold: float
) -> pd.DataFrame:
    """
    Calculates various hydrological metrics for each of the provided traces.

    Args:
        traces_timeseries (pd.DataFrame): Tidy DataFrame containing annual flow data. Each row represents a single
                                          year's data for a specific trace. At a minimum, the DataFrame should include
                                          the following columns:
                                          - 'Trace': Unique identifier for each trace
                                          - 'LF_Annual': Annual flow value for each year of the trace
        drought_threshold (float): The threshold value below which a drought is considered to occur.

    Returns:
        pd.DataFrame: DataFrame containing the calculated metrics for each trace, including:
                      - Median: Median annual flow for each trace
                      - Max: Maximum annual flow for each trace
                      - Min: Minimum annual flow for each trace
                      - SD: Standard deviation of annual flows for each trace
                      - IQR: Inter-quartile range of annual flows for each trace
                      - MaxDroughtLength: Maximum drought period length for each trace
                      - Driest10yr: Average annual flow for the driest consecutive 10-year period for each trace
                      - Wettest10yr: Average annual flow for the wettest consecutive 10-year period for each trace
    """

    n_traces = max(traces_timeseries[ids.TRACE])

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
    """
    Calculates the maximum drought length for a single trace.

    Args:
        trace_data (np.ndarray): A 1-dimensional numpy array containing the annual flow values for a single trace.
        drought_threshold (float): The threshold value below which a drought is considered to occur.

    Returns:
        int: The length of the longest drought period. Returns 0 if there are no drought periods.
    """

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
) -> Tuple[float, float]:
    """
    Calculates the driest and wettest 10-year flow averages for a single trace.

    Args:
        trace_data (np.ndarray): A 1-dimensional numpy array containing the annual flow values for a single trace.

    Returns:
        Tuple[float, float]: A tuple containing two float values:
                             - The average flow for the trace's driest consecutive 10-year period
                             - The average flow for the trace's wettest consecutive 10-year period
    Raises:
        ValueError: If the length of trace_data is less than 10 years.
    """

    if len(trace_data) < 10:
        raise ValueError("Trace data must be at least 10 years long to calculate 10-year periods.")

    # Initialize values
    driest_10yr = float('inf')
    wettest_10yr = float('-inf')

    # Iterate over the trace data to calculate sums of 10-year periods
    for i in range(len(trace_data) - 9):

        ten_year_sum = np.sum(trace_data[i:i + 10])

        # Update the driest 10-yr sum if the current sum is smaller
        driest_10yr = min(driest_10yr, ten_year_sum)

        # Update the wettest 10-yr sum if the current sum is smaller
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
) -> None:
    """
    Creates overlaid plots (kernel density estimate or histogram) for a specified variable from two different datasets.

    Args:
        df1 (pd.DataFrame): The first DataFrame containing the data to be plotted.
        df2 (pd.DataFrame): The second DataFrame containing the data to be plotted.
        color1 (str): The color for the plot of the first DataFrame.
        color2 (str): The color for the plot of the second DataFrame.
        variable (str): The column name of the variable to be plotted.
        xlabel (str): The label for the x-axis.
        ax (matplotlib.axes.Axes): The Matplotlib Axes object on which the plot will be drawn.
        plot_type (str, optional): The type of plot to generate. Can be 'kde' for kernel density estimate or 'hist' for
                                   histogram. Default is 'kde'.

    Raises:
        ValueError: If 'plot_type' is not 'kde' or 'hist'.

    Returns:
        None: The function does not return any value. It directly modifies the provided Axes object.
    """

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
        sns.kdeplot(
            data=df1,
            fill=True,
            color=color1,
            **common_params
        )
        sns.kdeplot(
            data=df2,
            fill=True,
            color=color2,
            **common_params
        )
    elif plot_type == 'hist':

        n_bins = len(np.unique(df1[variable]))

        sns.histplot(
            data=df1,
            color=color1,
            kde=False,
            bins=n_bins,
            stat='density',
            **common_params
        )
        sns.histplot(
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


def create_cumulative_timeseries_som_view(
        five_hundred_sow_info: pd.DataFrame,
        five_hundred_sow_cumulative_timeseries: pd.DataFrame,
        colors: List[str],
        n_neurons: int = 52,
        rows: int = 4,
        cols: int = 13
) -> plt.figure:
    """
    Creates a visual representation of a hexagonal Self-Organizing Map (SOM) where each neuron is displayed as a
    cumulative timeseries plot. The cumulative timeseries for States of the World (SOWs) within the neuron are
    highlighted in a specified color, while those for SOWs outside the neuron are shown in gray.

    Args:
        five_hundred_sow_info (pd.DataFrame): DataFrame containing information about the SOWs, where each row
                                              corresponds to a SOW. The required column is:
                                              - 'Neuron': The neuron number (starting from 1) in the SOM to which the
                                                          SOW is assigned.
        five_hundred_sow_cumulative_timeseries (pd.DataFrame): DataFrame containing the cumulative timeseries data for
                                                               each SOW. Each column corresponds to a SOW, and each row
                                                               represents a year.
        colors (List[str]): List of colors used to highlight the timeseries for each neuron.
        n_neurons (int, optional): The total number of neurons in the SOM. Default is 52.
        rows (int, optional): The number of rows in the subplot grid. Default is 4.
        cols (int, optional): The number of columns in the subplot grid. Default is 13.

    Returns:
        plt.figure: A Matplotlib figure object containing the cumulative timeseries SOM view.
    """

    sim_years = np.arange(2027, 2057)

    fig = plt.figure(figsize=(19, 9.5))
    gs = GridSpec(
        nrows=rows,
        ncols=2 * cols + 1,  # double number of columns and add 1 to allow for shifting of rows
        figure=fig
    )

    ## Plotting loop
    neurons = np.unique(five_hundred_sow_info[ids.NEURON])

    for neuron in neurons:
        x, y = get_subplot_coordinates(
            neuron=neuron,
            n_neurons=n_neurons,
            rows=rows,
            cols=cols
        )
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
            color=colors[neuron - 1]
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

    return fig


def get_subplot_coordinates(
        neuron: int,
        n_neurons: int,
        rows: int,
        cols: int
) -> [int, int]:
    """
    Calculate the GridSpec row and column positions for a subplot based on the neuron index in a Self-Organizing Map
    (SOM).

    Assumptions:
        - The SOM configuration is hexagonal, meaning that alternating rows are shifted. Rows shifted to the right are
          those starting with the bottom row and every other row going up.
        - The GridSpec configuration for plotting the subplot has twice the number of SOM columns plus one (e.g., if the
          SOM has 13 columns, the GridSpec configuration has 27 columns: 13 * 2 + 1).
        - Neurons are numbered starting from 1 in the bottom-left corner, increasing along the bottom row. Numbering
          continues from left to right in the next row up, and so on, with the last neuron being in the top-right corner

    Args:
        neuron (int): The index of the neuron (starting at 1) for which to calculate the subplot coordinates.
        n_neurons (int): The total number of neurons in the SOM.
        rows (int): The number of rows in the SOM.
        cols (int): The number of columns in the SOM.

    Returns:
        Tuple[int, int]: A tuple containing the row and column GridSpec coordinates for the neuron subplot.

    """

    x = int(np.floor((n_neurons - neuron) / cols))
    y = int(((neuron - 1) % cols)) * 2

    if x % 2 != 0:  # check if x is odd, if so, add 1 to implement the row shift
        y += 1

    return int(x), int(y)


def plot_cumulative_timeseries(
        sim_years: np.ndarray,
        cumulative_timeseries: pd.DataFrame,
        indices: List[int],
        non_indices: List[int],
        ax: matplotlib.axes.Axes,
        color: str
) -> None:
    """
    Plots cumulative timeseries onto a single matplotlib Axes object. Timeseries specified by 'non_indices' (i.e., those
    not in the neuron being plotted) are colored gray, and timeseries specified by 'indices' (i.e., those in the neuron
    being plotted) are colored in the specified color.

    Args:
        sim_years (np.ndarray): A 1-dimensional numpy array containing the simulation years.
        cumulative_timeseries (pd.DataFrame): A DataFrame containing the cumulative timeseries data. Each column
                                              represents a timeseries.
        indices (List[int]): A list of indices representing the timeseries to be highlighted in the specified color.
        non_indices (List[int]): A list of indices representing the timeseries to be colored gray.
        ax (matplotlib.axes.Axes): The matplotlib Axes object on which the timeseries will be plotted.
        color (str): The color used to highlight the timeseries specified by 'indices'.

    Returns:
        None: This function does not return any value. It directly modifies the provided Axes object.
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


def create_condensed_som_figure(
        title: str,
        colorbar_label: str,
        neuron_values: List[float],
        neuron_coordinates: pd.DataFrame,
        color_scheme: str,
        n_digits: int,
        annotation_size: int,
        inverse_colorbar: bool
) -> plt.figure:
    """
    Creates a visual representation of a Self-Organizing Map (SOM) where each neuron is depicted as a hexagon with a
    fill color corresponding to its specified value and the specified color scheme.

    Args:
        title (str): The title of the figure.
        colorbar_label (str): The label for the colorbar.
        neuron_values (List[float]): A list of values corresponding to each neuron. The length must be the same as the
                                     number of entries in neuron_coordinates.
        neuron_coordinates (pd.DataFrame): A DataFrame containing the x and y coordinates of each neuron. Required
                                           columns:
                                           - 'x': The x-coordinate of each neuron.
                                           - 'y': The y-coordinate of each neuron.
        color_scheme (str): The name of the color scheme to be used for the hexagons.
        n_digits (int): The number of digits to round the neuron values to for display.
        annotation_size (int): The font size for the neuron value annotations.
        inverse_colorbar (bool): Whether to invert the colorbar to have the largest value on the left.

    Returns:
        plt.figure: A Matplotlib figure object containing the condensed (hexagonal) SOM.
    """

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
            s=round(value, n_digits) if n_digits > 0 else round(value),
            ha='center',
            va='center',
            fontsize=annotation_size,
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
        label='Color/Text: ' + colorbar_label,
        fontsize=26,
        labelpad=10
    )
    cbar.ax.tick_params(labelsize=24)

    ## Inverse colorbar to have the largest value on the left
    if inverse_colorbar:
        cbar.ax.invert_xaxis()

    # Set limits and turn off axes
    ax.set_xlim([0, 2 * hex_radius * grid_width + 2])
    ax.set_ylim([0, 2 * hex_radius * grid_height + 1])
    ax.axis('off')

    fig.subplots_adjust(
        left=0.048,
        right=0.997,
        top=0.965,
        bottom=0.151
    )

    title_y = 0.88 if '\n' in title else 0.835
    fig.suptitle(title, fontsize=35, y=title_y)

    return fig


def get_color_scheme(
        scheme_name: str
) -> List[Tuple[float, str]]:
    """
    Retrieves a predefined color scheme based on the given scheme name.

    Args:
        scheme_name (str): The name of the color scheme to retrieve. Available schemes are:
            - 'wet_dry': A gradient from dark red to blue, passing through yellow and green.
            - 'warm_cool': A gradient from red to purple, passing through orange, yellow, pink, and turquoise.
            - 'good_bad': A gradient from green to red, passing through yellow.
            - 'good_bad_reverse': A gradient from red to green, passing through yellow.

    Returns:
        List[Tuple[float, str]]: A list of tuples representing the color scheme.

    Raises:
        ValueError: If the provided color scheme is not one of the existing predefined schemes.
    """

    schemes = {
        'wet_dry': [
            (0, '#8B0000'), (0.25, '#FFD700'), (0.5, '#008000'), (0.75, '#add8e6'), (1, '#0000FF')
        ],
        'warm_cool': [
            (0, '#FF4500'), (0.2, '#FFA500'), (0.4, '#FFFF00'), (0.6, '#FFC0CB'), (0.8, '#40E0D0'), (1, '#800080')
        ],
        'good_bad': [
            (0, '#008000'), (0.5, '#FFFF00'), (1, '#FF0000')
        ],
        'good_bad_reverse': [
            (0, '#FF0000'), (0.5, '#FFFF00'), (1, '#008000')
        ]
    }

    if scheme_name not in schemes:
        raise ValueError(
            "Please provide an existing color scheme ('wet_dry', 'warm_cool', 'good_bad', 'good_bad_reverse')"
        )

    return schemes.get(scheme_name)


def map_value_to_color(
        value: float,
        value_range: Tuple[float, float],
        cmap: mcolors.LinearSegmentedColormap
) -> Tuple[float, float, float, float]:
    """
    Maps a given value to a color in a specified colormap based on the value's position within a defined range.

    Args:
        value (float): The value to be mapped to a color.
        value_range (Tuple[float, float]): A tuple containing the minimum and maximum values of the range.
        cmap (mcolors.LinearSegmentedColormap): A Matplotlib colormap used to map the normalized value to a color.

    Returns:
        Tuple[float, float, float, float]: A tuple representing the RGBA color mapped from the input value.
    """

    # Normalize the value to the range [0, 1]
    normalized_value = (value - value_range[0]) / (value_range[1] - value_range[0])

    # Map the normalized value to a color in the gradient
    return cmap(normalized_value)


def get_text_color_for_neuron(
        rgba: Tuple[float, float, float, float]
) -> str:
    """
    Determines whether the text color should be white or black based on the luminance of the background color.

    Args:
        rgba (Tuple[float, float, float, float]): 4-value tuple representing the RGBA values.

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


def create_and_save_dendrogram(
        linkage_name: str,
        z: np.ndarray,
        neuron_labels: np.ndarray,
        filename: str
) -> None:
    """
    Creates and saves a dendrogram figure based on a linkage matrix.

    Args:
        linkage_name (str): The type of linkage used, to be used for the plot title.
        z (np.ndarray): The hierarchical clustering linkage matrix (derived via scipy.cluster.hierarchy.linkage)
        neuron_labels (np.ndarray): An array of labels for the neurons.
        filename (str): The path and filename where the dendrogram plot will be saved.

    Returns:
        None: This function does not return any value. It creates and saves the dendrogram plot to the specified file.
    """
    # Create a new figure
    plt.figure(figsize=(19, 9.5))

    # Plot the dendrogram
    dendrogram(
        Z=z,
        labels=neuron_labels
    )

    # Add labels
    plt.title(linkage_name)
    plt.xlabel('Neuron')

    # Save the figure
    plt.savefig(
        fname=filename,
        dpi=400,
        bbox_inches='tight'
    )

    # Close the figure
    plt.close()


def get_policy_reevaluation_data(
        reevaluation_data: pd.DataFrame,
        five_hundred_sow_info: pd.DataFrame,
        experiment: str,
        policy: int
) -> pd.DataFrame:
    """
    Filters the reevaluation data for a given experiment and policy, merges it with SOW neuron information, and pivots
    the data to a wide format where each objective is a column.

    Args:
        reevaluation_data (pd.DataFrame): DataFrame containing the reevaluation data with columns:
                                          - 'Experiment': The name of the experiment.
                                          - 'Policy': The policy ID.
                                          - 'SOW': The state of the world identifier.
                                          - 'Objective': The objective being measured.
                                          - 'Value': The value of the objective.
        five_hundred_sow_info (pd.DataFrame): DataFrame containing SOW information including neuron assignments.
                                              Necessary columns:
                                              - 'SOW': The numerical state of the world identifier.
                                              - 'Neuron': The neuron number in the SOM to which the SOW is assigned.
        experiment (str): The name of the experiment to filter the data by.
        policy (int): The policy ID to filter the data by.

    Returns:
        pd.DataFrame: A DataFrame in wide format for the specified experiment/policy where rows represent SOWs and
                      each objective is a separate column.
    """

    policy_data = reevaluation_data[
        (reevaluation_data[ids.EXPERIMENT] == experiment) &
        (reevaluation_data[ids.POLICY] == policy)
        ]
    policy_data = pd.merge(
        left=policy_data,
        right=five_hundred_sow_info[[ids.SOW, ids.NEURON]],
        on=ids.SOW,
        how='left'
    )
    policy_data_wide = policy_data.pivot(
        index=[ids.EXPERIMENT, ids.POLICY, ids.SOW, ids.NEURON],
        columns=ids.OBJECTIVE,
        values='Value'
    ).reset_index()

    return policy_data_wide
