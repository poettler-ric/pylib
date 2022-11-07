#!/usr/bin/env python3

"""
Display data from a wflow run.tss file.
"""

__author__ = "Richard Pöttler"
__copyright__ = "Copyright (c) 2022 Richard Pöttler"
__license__ = "MIT"
__email__ = "richard.poettler@gmail.com"


from argparse import ArgumentParser
from configparser import ConfigParser
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_tss(filename):
    """Parse a given tss file."""

    skip_lines = 0
    with open(filename, encoding="utf8") as stream:
        stream.readline()
        skip_lines = int(stream.readline()) + 2
    return np.loadtxt(filename, skiprows=skip_lines)


def read_measured(filename, skip_lines):
    """Read real measurement data."""

    return pd.read_csv(
        filename,
        encoding="latin1",
        sep=";",
        decimal=",",
        skiprows=skip_lines,
        names=["date", "q_m3s"],
        index_col="date",
        parse_dates=["date"],
        dayfirst=True,
    )


def main():
    """Main routine."""

    parser = ArgumentParser(description="Prepare wflow files")
    parser.add_argument("config", help="configuration file destination")
    parser.add_argument("tss", help="tss file to display")
    args = parser.parse_args()

    config = ConfigParser()
    config.read(args.config)

    from_date = datetime.strptime(
        config["DEFAULT"]["from_date"], "%Y-%m-%dT%H:%M:%S"
    ).strftime("%Y-%m-%d %H:%M:%S")
    to_date = datetime.strptime(
        config["DEFAULT"]["to_date"], "%Y-%m-%dT%H:%M:%S"
    ).strftime("%Y-%m-%d %H:%M:%S")

    model_data = parse_tss(args.tss)
    axes = plt.subplot()

    if config["DEFAULT"]["y_max"]:
        plt.ylim(top=config.getint("DEFAULT", "y_max"))

    for key, key_config in config.items():
        if key.startswith("level."):
            name = key[len("level.") :]
            measured = read_measured(
                key_config["filename"], key_config.getint("skip_lines")
            )
            axes.plot(
                measured.loc[from_date:to_date].index,
                model_data[:, key_config.getint("tss_id")],
                label=f"model {name}",
            )
            axes.plot(measured.loc[from_date:to_date].q_m3s, label=f"measured {name}")

    axes.legend()
    plt.show()


if __name__ == "__main__":
    main()
