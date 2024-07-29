import numpy as np
import pandas as pd
from sklearn.preprocessing import scale

from src.utils import ids
from src.utils.functions_library import calculate_trace_metrics


def main():

    # Constants
    traces_timeseries_file = "../data/traces_timeseries.csv"
    scd_file = "../data/scd_csv.csv"
    output_dir = "python_output"
    start_year = 2027
    end_year = 2057
    drought_threshold = 14.65

    # Prepare hydrologic trace data
    ## Set simulation years
    sim_years = np.arange(start_year, end_year)
    n_years = len(sim_years)

    ## Read traces timeseries data and filter to simulation years
    traces_timeseries = pd.read_csv(traces_timeseries_file)
    traces_timeseries = traces_timeseries[
        (traces_timeseries[ids.YEAR] >= start_year) & (traces_timeseries[ids.YEAR] <= end_year)
        ]
    n_traces = max(traces_timeseries[ids.TRACE])

    ## Initialize another dataframe to store general trace info (e.g., median flow)
    traces_info = traces_timeseries[
        [ids.TRACE, ids.SCENARIO, ids.TRACE_NUMBER]
    ].drop_duplicates().reset_index(drop=True)

    ## Calculate and concatenate trace metrics
    traces_metrics = calculate_trace_metrics(
        traces_timeseries=traces_timeseries,
        n_traces=n_traces,
        drought_threshold=drought_threshold
    )
    traces_info = pd.concat(
        objs=[traces_info, traces_metrics],
        axis=1
    )

    # Read System Condition and Demand file
    scd_samples = pd.read_csv(scd_file)
    n_scd = len(scd_samples)

    # Generate full-factorial ensemble
    ## Repeat each trace n_scd times
    full_factorial_ensemble = traces_info.loc[traces_info.index.repeat(n_scd)].reset_index(drop=True)

    ## Repeat scd samples n_traces times
    repeated_scd_samples = pd.concat(
        objs=[scd_samples] * n_traces,
        ignore_index=True
    )

    ## Concatenate dataframes, creating all unique combinations of trace + scd sample
    full_factorial_ensemble = pd.concat(
        objs=[full_factorial_ensemble, repeated_scd_samples],
        axis=1
    )

    # Calculate and store Cumulative Timeseries for each SOW in full factorial ensemble
    ## Convert dataframes to NumPy for efficient calculations
    traces_timeseries_array = traces_timeseries[[ids.TRACE, ids.LF_ANNUAL]].to_numpy()
    full_factorial_ensemble_array = full_factorial_ensemble[[ids.TRACE, ids.DEMAND, ids.INIT_STORAGE]].to_numpy()

    ## NUMERIC indices to call columns (need this for working with NumPy, as opposed to the string ids used above)
    TRACE_IDX = 0
    DEMAND_IDX = 1
    INIT_STORAGE_IDX = 2
    LF_ANNUAL_IDX = 1

    ## Initialize array for storing cumulative timeseries
    n_sow = full_factorial_ensemble_array.shape[0]
    cumulative_timeseries_array = np.zeros((n_years, n_sow))

    # Vectorized calculation of cumulative timeseries for each SOW
    for i in range(n_sow):
        trace = full_factorial_ensemble_array[i, TRACE_IDX]
        demand = full_factorial_ensemble_array[i, DEMAND_IDX]
        init_storage = full_factorial_ensemble_array[i, INIT_STORAGE_IDX]

        ## Filter LF_Annual values corresponding to current trace
        lf_annual = traces_timeseries_array[traces_timeseries_array[:, TRACE_IDX] == trace, LF_ANNUAL_IDX]

        ## Calculate cumulative timeseries
        flow_minus_demand_timeseries = lf_annual - demand
        cumulative_timeseries = np.cumsum(flow_minus_demand_timeseries) + init_storage

        # Store cumulative timeseries
        cumulative_timeseries_array[:, i] = cumulative_timeseries

    # Save:
    ## 1) info of each SOW in full factorial ensemble, and;
    full_factorial_ensemble.to_csv(
        path_or_buf=f'{output_dir}/full_factorial_sow_info.csv',
        index=False
    )
    ## 2) cumulative timeseries of each SOW in full factorial ensemble
    np.save(
        file=f'{output_dir}/full_factorial_cumulative_timeseries.npy',
        arr=cumulative_timeseries_array
    )

    # Scale cumulative timeseries for SOM fitting and save
    cumulative_timeseries_array_scaled = scale(cumulative_timeseries_array.T)
    np.savetxt(
        fname=f'{output_dir}/full_factorial_cumulative_timeseries_scaled.csv',
        X=cumulative_timeseries_array_scaled,
        delimiter=','
    )

    ## Save scaling factors separately
    means = np.mean(cumulative_timeseries_array, axis=1)
    stds = np.std(cumulative_timeseries_array, axis=1)
    scale_factors = pd.DataFrame({'Mean': means, 'Std': stds})
    scale_factors.to_csv(
        path_or_buf=f'{output_dir}/full_factorial_cumulative_timeseries_scaling_factors.csv',
        index=False
    )


if __name__ == '__main__':
    main()
