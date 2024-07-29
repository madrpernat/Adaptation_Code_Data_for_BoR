import matplotlib
import pandas as pd


def main():

    five_hundred_sow_info = pd.read_csv(
        filepath_or_buffer='output/python_output/500_sow_info.csv'
    )
    five_hundred_sow_cumulative_timeseries = pd.read_csv(
        filepath_or_buffer='output/python_output/500_sow_cumulative_timeseries.csv'
    )
    full_factorial_sow_info = pd.read_csv(
        filepath_or_buffer='output/python_output/full_factorial_sow_info.csv'
    )

    #



if __name__ == '__main__':
    main()
