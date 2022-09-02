#!/usr/bin/env python3

from argparse import ArgumentParser
from configparser import ConfigParser
from datetime import datetime, timedelta
from os.path import join as pjoin
import numpy as np


def get_result(result_folder, column, start_time, time_step):
    tss_file = pjoin(result_folder, "run.tss")
    to_skip = None

    # determine the lines to skip
    with open(tss_file) as f:
        f.readline()
        to_skip = int(f.readline()) + 2

    # read data
    data = np.loadtxt(tss_file, skiprows=to_skip)

    # compute timestamps for data rows
    time_list = np.empty(data.shape[0])
    for i in range(data.shape[0]):
        time_list[i] = (start_time + i * time_step).timestamp()

    return np.column_stack((time_list, data[:, column]))


def to_float_nan(v):
    try:
        return float(v)
    except ValueError:
        return np.nan


def get_gauge(file_name):
    time_list = []
    value_list = []

    with open(file_name) as f:
        for line in f:
            line = line.split(",")

            time_list.append(
                datetime.strptime(line[0], "%d.%m.%Y %H:%M:%S").timestamp()
            )
            value_list.append(to_float_nan(line[1]))

    return np.column_stack((time_list, value_list))


def parse_config(config_file):
    config = ConfigParser()
    config.read(config_file)

    start_time = datetime.strptime(config["run"]["starttime"], "%Y-%m-%d %H:%M:%S")
    time_step = timedelta(seconds=int(config["run"]["timestepsecs"]))

    return start_time, time_step


def select_evaluation_window(data, start_time, end_time):
    indexes = np.where(
        (data[:, 0] >= start_time.timestamp()) & (data[:, 0] <= end_time.timestamp())
    )[0]
    return data[indexes]


def calc_nse(model_data, reference_data):
    return (
        1
        - ((model_data - reference_data) ** 2).sum()
        / ((model_data - reference_data.mean()) ** 2).sum()
    )


def main():
    parser = ArgumentParser(description="Evaluate a wflow run")
    parser.add_argument("config_file", help="simulation configuration file")
    parser.add_argument("result_folder", help="folder containing the results")
    parser.add_argument("result_column", type=int, help="column to compare in run.tss")
    parser.add_argument("gauge_file", help="gauge data")
    parser.add_argument("start_time", help="start date of the evaluation")
    parser.add_argument("end_time", help="end date of the evaluation")
    args = parser.parse_args()

    args.start_time = datetime.strptime(args.start_time, "%Y-%m-%dT%H:%M:%S")
    args.end_time = datetime.strptime(args.end_time, "%Y-%m-%dT%H:%M:%S")

    (model_start, time_step) = parse_config(args.config_file)

    results = get_result(args.result_folder, args.result_column, model_start, time_step)
    results = select_evaluation_window(
        results,
        args.start_time,
        args.end_time,
    )

    gauge = get_gauge(args.gauge_file)
    gauge = select_evaluation_window(gauge, args.start_time, args.end_time)

    nse = calc_nse(results[:, 1], gauge[:, 1])
    print(f"nse: {nse}")


if __name__ == "__main__":
    main()
