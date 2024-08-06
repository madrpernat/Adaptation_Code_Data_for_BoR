import matplotlib
import os
import pandas as pd
import re

from src.configs import policy_performance_configs
from src.settings.objective_settings import objective_settings
from src.utils import ids
from src.utils.functions_library import create_condensed_som_figure, get_policy_reevaluation_data

matplotlib.use('Qt5Agg')


def main():

    # Set working directory
    current_dir = os.path.dirname(__file__)
    parent_dir = os.path.abspath(os.path.join(current_dir, os.pardir))
    os.chdir(parent_dir)

    output_dir = 'output/figs/policy_performance'

    # Read in data
    five_hundred_sow_info = pd.read_csv(
        filepath_or_buffer='output/500_sow_info.csv'
    )
    som_neuron_coordinates = pd.read_csv(
        filepath_or_buffer='output/som_neuron_coordinates.csv'
    )
    reevaluation_data = pd.read_csv(
        filepath_or_buffer='output/new_experiments_reevaluation_data.csv'
    )

    # Plotting configurations
    objective_view_configs = policy_performance_configs.objective_view_configs
    acceptability_view_configs = policy_performance_configs.acceptability_view_configs

    for config in objective_view_configs:

        experiment = config['experiment']
        policy = config['policy']

        policy_data = get_policy_reevaluation_data(
            reevaluation_data=reevaluation_data,
            five_hundred_sow_info=five_hundred_sow_info,
            experiment=experiment,
            policy=policy
        )

        for objective in config['objectives']:

            obj_settings = objective_settings[objective]
            label = obj_settings['label']
            label_w_unit = obj_settings['label_w_unit']
            scale = obj_settings['scale']
            n_digits = obj_settings['n_digits']
            annotation_size = obj_settings['annot_size']

            objective_data = policy_data[[ids.NEURON, objective]]

            avg_objective_per_neuron = objective_data.groupby(ids.NEURON).mean().reset_index()[objective]
            avg_objective_per_neuron = avg_objective_per_neuron / scale

            fig = create_condensed_som_figure(
                title=f'Experiment {experiment}, Policy {policy}\nObjective View: {label}',
                colorbar_label=f'Neuron Average {label_w_unit}',
                neuron_values=avg_objective_per_neuron,
                neuron_coordinates=som_neuron_coordinates,
                color_scheme='good_bad' if scale > 0 else 'good_bad_reverse',
                n_digits=n_digits,
                annotation_size=annotation_size,
                inverse_colorbar=False if scale > 0 else True
            )

            fig.savefig(
                fname=f'{output_dir}/exp{experiment}_policy{policy}__{objective}.png',
                dpi=400,
                bbox_inches='tight'
            )

    for config in acceptability_view_configs:

        experiment = config['experiment']
        policy = config['policy']

        policy_data = get_policy_reevaluation_data(
            reevaluation_data=reevaluation_data,
            five_hundred_sow_info=five_hundred_sow_info,
            experiment=experiment,
            policy=policy
        )

        thresholds_met = pd.DataFrame()
        objectives = config['objectives']
        thresholds = config['thresholds']

        # We will build the acceptability definition as we iterate through each of the objectives/thresholds below
        acceptability_def = ""

        # For each objective, check whether each SOW met or didn't meet threshold
        for objective, threshold in zip(objectives, thresholds):

            # Get objective settings to create string acceptability defn
            obj_settings = objective_settings[objective]
            label = obj_settings['label']
            scale = obj_settings['scale']
            unit = obj_settings['unit']
            n_digits = obj_settings['n_digits']

            objective_values = policy_data[objective]
            threshold_met = objective_values < threshold
            thresholds_met[objective] = threshold_met

            op = '<' if scale > 0 else '>'
            threshold_label = threshold / scale
            threshold_label = round(threshold_label, n_digits) if n_digits > 0 else round(threshold_label)
            acceptability_def = acceptability_def + f'{label} {op} {threshold_label} {unit}, '

        # All_True only if ALL objective thresholds were met, otherwise False
        thresholds_met['AllTrue'] = thresholds_met[objectives].all(axis=1)

        # Calculate % of sows acceptable within each neuron
        thresholds_met[ids.NEURON] = policy_data[ids.NEURON]
        pct_acceptable_per_neuron = thresholds_met.groupby(ids.NEURON)['AllTrue'].mean().reset_index()['AllTrue'] * 100

        fig = create_condensed_som_figure(
            title=f'Experiment {experiment}, Policy {policy}\nAcceptability: {acceptability_def[:-2]}',
            colorbar_label='Percent SOW Acceptable',
            neuron_values=pct_acceptable_per_neuron,
            neuron_coordinates=som_neuron_coordinates,
            color_scheme='good_bad_reverse',
            n_digits=1,
            annotation_size=27,
            inverse_colorbar=True
        )

        acceptability_def_wo_sym = re.sub(
            pattern="[><= .,]",
            repl="",
            string=acceptability_def[:-2]
        )

        fig.savefig(
            fname=f'{output_dir}/exp{experiment}_policy{policy}__{acceptability_def_wo_sym}.png',
            dpi=400,
            bbox_inches='tight'
        )


if __name__ == '__main__':
    main()
