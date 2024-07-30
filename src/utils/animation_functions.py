import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from scipy.cluster.hierarchy import cut_tree, dendrogram
import matplotlib.patches as patches


def create_hexagon_grid(ax, neuron_coordinates, grid_width, grid_height, hex_radius, colors):

    hexagons = []
    texts = []

    for idx, row in neuron_coordinates.iterrows():

        x, y = row['x'], row['y']

        hexagon = patches.RegularPolygon(
            xy=(x, y),
            numVertices=6,
            radius=hex_radius,
            orientation=np.radians(0),
            facecolor=colors[idx],
            edgecolor='black'
        )

        ax.add_patch(hexagon)
        hexagons.append(hexagon)

        text = ax.text(
            x=x,
            y=y,
            s=str(idx + 1),
            ha='center',
            va="center",
            fontsize=14,
            color="black",
            weight="bold"
        )
        texts.append(text)

    ax.set_aspect('equal')
    ax.set_xlim([0, 2 * hex_radius * grid_width + 2])
    ax.set_ylim([0, 2 * hex_radius * grid_height + 1])
    ax.axis('off')

    return hexagons, texts


# Function to update colors based on dendrogram cutoff
def update_hexagon_colors(config, cutoff):
    cut = cut_tree(config['z'], height=cutoff)

    for i in range(len(config['hexagons'])):
        config['hexagons'][i].set_facecolor(config['colors'][cut[i][0]])


# Function to update animation
def update(frame, config):
    # Move the cutoff line up
    config['cutoff_line'].set_ydata([frame, frame])

    # Update hexagon colors based on new cutoff
    update_hexagon_colors(config, frame)

    return [config['cutoff_line']] + config['hexagons'] + config['texts']


def create_animation(config, z, neuron_coordinates, colors):

    config['z'] = z

    fig, (ax_dendro, ax_hex) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(19, 9.5)
    )

    dendrogram(
        Z=config['z'],
        ax=ax_dendro
    )

    ax_dendro.set_title(config['title'], fontsize=23)
    ax_dendro.set_xlabel('Neuron', fontsize=20)
    ax_dendro.set_ylabel(config['y_label'], fontsize=20)
    ax_dendro.tick_params(axis='both', which='major', labelsize=17)

    # Add initial cutoff line
    cutoff_line, = ax_dendro.plot([-15, 600], [0.5, 0.5], 'black', linewidth=2)
    config['cutoff_line'] = cutoff_line

    hexagons, texts = create_hexagon_grid(
        ax=ax_hex,
        neuron_coordinates=neuron_coordinates,
        grid_width=13,
        grid_height=4,
        hex_radius=0.53,
        colors=colors
    )

    config['hexagons'] = hexagons
    config['texts'] = texts
    config['colors'] = colors

    ani = FuncAnimation(
        fig=fig,
        func=update,
        fargs=(config,),
        frames=np.linspace(config['initial_cutoff'], config['final_cutoff'], 190),
        blit=True
    )

    ani.save(
        filename=f"output/python_output/figs/{config['name']}_animation.mp4",
        writer='ffmpeg',
        fps=1
    )
