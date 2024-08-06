import hiplot as hip
import os
import pandas as pd

from src.configs import pc_configs


def main():

    # Set working directory
    current_dir = os.path.dirname(__file__)
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    os.chdir(parent_dir)

    archive_file_name = 'condensed_archive.csv'
    exp_id_file_name = 'experiment_id.txt'
    configs = pc_configs.configs

    for config in configs:

        # Retrieve the list of objectives to plot and initialize an empty df to store policy objective values
        objectives = config['objectives']
        objective_data = pd.DataFrame(
            columns=['Experiment', 'Policy'] + objectives
        )

        for exp_directory in config['exp_directories']:

            # Read the experiment's ID file. This ID will be used to label each policy's objective data to keep track
            # of which policy came from which experiment in the final HiPlot visual.
            with open(os.path.join(exp_directory, exp_id_file_name), 'r') as file:
                exp_id = file.read().strip()

            # Read the objective data for the experiment. Filter to only include the specified objectives.
            exp_objectives = pd.read_csv(
                filepath_or_buffer=os.path.join(exp_directory, archive_file_name)
            )[objectives]

            # Insert the experiment ID and policy ID into the df
            exp_objectives.insert(
                loc=0,
                column='Experiment',
                value=[exp_id] * len(exp_objectives)
            )
            exp_objectives.insert(
                loc=1,
                column='Policy',
                value=range(1, len(exp_objectives) + 1)
            )

            # Append to "all experiments" objective_data df
            objective_data = pd.concat(
                objs=[objective_data, exp_objectives],
                ignore_index=True
            )

        # Create HiPlot visual

        ## HiPlot reverses the order of columns by default, so reverse them here so that they show up in the right order
        pc = hip.Experiment.from_dataframe(objective_data[objective_data.columns[::-1]])
        pc.display_data(hip.Displays.TABLE).update({"hide": ["uid", "from_uid"]})

        ## Hide 'Policy' axis from PC plot (but keep in table)
        pc.display_data(hip.Displays.PARALLEL_PLOT).update({'hide': ["Policy"]})

        ## Order table columns to match order in PC plot
        ordered_columns = ['Experiment', 'Policy'] + objectives
        pc.display_data(hip.Displays.TABLE).update({
            'order': ordered_columns
        })

        ## Save as html file
        pc.to_html(os.path.join(config['output'], f'pc_plot__{config["name"]}.html'))


if __name__ == '__main__':
    main()
