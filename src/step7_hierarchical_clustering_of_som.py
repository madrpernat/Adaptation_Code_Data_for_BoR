import distinctipy
import matplotlib
import numpy as np
import os
import pandas as pd
from scipy.cluster.hierarchy import cut_tree, linkage

from src.utils import animation_configs
from src.utils.animation_functions import create_animation
from src.utils.functions_library import create_and_save_dendrogram

matplotlib.use('Qt5Agg')


def main():

    # Set working directory
    current_dir = os.path.dirname(__file__)
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    os.chdir(parent_dir)

    som_codes_scaled = pd.read_csv(
        filepath_or_buffer='output/r_output/som_codes_scaled.csv'
    )
    neuron_coordinates = pd.read_csv(
        filepath_or_buffer='output/r_output/som_neuron_coordinates.csv'
    )

    # Perform hierarchical clustering on the neuron weights using different linkage methods
    z_avg_link = linkage(
        y=som_codes_scaled,
        method='average',
        metric='euclidean'
    )
    z_complete_link = linkage(
        y=som_codes_scaled,
        method='complete',
        metric='euclidean'
    )
    z_centroid_link = linkage(
        y=som_codes_scaled,
        method='centroid',
        metric='euclidean'
    )
    z_ward_link = linkage(
        y=som_codes_scaled,
        method='ward',
        metric='euclidean'
    )

    # Create dendrogram for each linkage
    neuron_labels = np.arange(1, 53)

    create_and_save_dendrogram(
        linkage_name="Average Linkage",
        z=z_avg_link,
        neuron_labels=neuron_labels,
        filename='output/python_output/figs/average_linkage_dendrogram.png'
    )
    create_and_save_dendrogram(
        linkage_name="Complete Linkage",
        z=z_complete_link,
        neuron_labels=neuron_labels,
        filename='output/python_output/figs/complete_linkage_dendrogram.png'
    )
    create_and_save_dendrogram(
        linkage_name="Centroid Linkage",
        z=z_centroid_link,
        neuron_labels=neuron_labels,
        filename='output/python_output/figs/centroid_linkage_dendrogram.png'
    )
    create_and_save_dendrogram(
        linkage_name="Ward Linkage",
        z=z_ward_link,
        neuron_labels=neuron_labels,
        filename='output/python_output/figs/ward_linkage_dendrogram.png'
    )

    # Create animation for each linkage method
    colors = distinctipy.get_colors(52)

    create_animation(
        config=animation_configs.avg_link,
        z=z_avg_link,
        neuron_coordinates=neuron_coordinates,
        colors=colors
    )
    create_animation(
        config=animation_configs.complete_link,
        z=z_complete_link,
        neuron_coordinates=neuron_coordinates,
        colors=colors
    )
    create_animation(
        config=animation_configs.centroid_link,
        z=z_centroid_link,
        neuron_coordinates=neuron_coordinates,
        colors=colors
    )
    create_animation(
        config=animation_configs.ward_link,
        z=z_ward_link,
        neuron_coordinates=neuron_coordinates,
        colors=colors
    )

    # The Complete Linkage method was chosen and two different dendrogram cuts were made. Uncomment the below code to
    # make the cuts and write out the neuron-->cluster mapping information.

    # cut_1 = pd.DataFrame(cut_tree(
    #         Z=z_complete_link,
    #         n_clusters=19
    # ))
    # neuron_clusters_1 = pd.DataFrame(
    #     {'Neuron': np.arange(1, 53), 'Cluster': cut_1.loc[:, 0]}
    # )
    # neuron_clusters_1.to_csv(
    #     path_or_buf="output/python_output/neuron_clusters_19.csv",
    #     index=False
    # )
    #
    # cut_2 = pd.DataFrame(cut_tree(
    #     Z=z_complete_link,
    #     n_clusters=8
    # ))
    # neuron_clusters_2 = pd.DataFrame(
    #     {"Neuron": np.arange(1, 53), "Cluster": cut_2.loc[:, 0]}
    # )
    # neuron_clusters_2.to_csv(
    #     path_or_buf="output/python_output/neuron_clusters_8.csv",
    #     index=False
    # )


if __name__ == '__main__':
    main()


