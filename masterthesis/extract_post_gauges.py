#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 09:22:17 2022

@author: iwbworkstation
"""

import datetime
import os
import os.path as path
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import numpy.typing as typing
import pandas as pd
import pytz
from osgeo import gdal

TEXTWITH_IN = 4.7747
SIMULATION_START_DATE = datetime.datetime(2013, 1, 1, 0, 0)
SIMULATION_END_DATE = datetime.datetime(2017, 12, 30, 0, 0)
CALIBRATION_START_DATE = datetime.datetime(2014, 1, 1, 0, 0)
OFFICIAL_END_DATE = datetime.datetime(2018, 1, 1, 0, 0)

CATCHMENT_FILE = (
    "/data/home/richi/master_thesis/model_MAR/staticmaps/wflow_subcatch.map"
)

LINE_WIDTH = 1.0
LAYOUT = "tight"


def nse(measured, simulated) -> float:
    return (
        1
        - ((simulated - measured) ** 2).sum()
        / ((measured - measured.mean()) ** 2).sum()
    )


def getCatchmentMask(catchment_file: str) -> typing.NDArray:
    catchment_mask = getRaster(catchment_file)
    catchment_mask[catchment_mask == 1] = 0
    catchment_mask[catchment_mask == np.nan] = 1
    return catchment_mask


def nce(measured, simulated) -> float:
    return 1.0 - abs(measured - simulated).sum() / measured.sum()


def kge(measured, simulated) -> float:
    m1, m2 = np.mean(measured), np.mean(simulated)
    r = np.sum((measured - m1) * (simulated - m2)) / (
        np.sqrt(np.sum((measured - m1) ** 2)) * np.sqrt(np.sum((simulated - m2) ** 2))
    )
    beta = m2 / m1
    gamma = (np.std(simulated) / m2) / (np.std(measured) / m1)
    return 1 - np.sqrt((r - 1) ** 2 + (beta - 1) ** 2 + (gamma - 1) ** 2)


def getRaster(filename, no_data=np.nan, mask=None):
    ds = gdal.Open(filename)
    raster = ds.GetRasterBand(1).ReadAsArray()
    no_data_value = ds.GetRasterBand(1).GetNoDataValue()
    raster[raster == no_data_value] = no_data
    if mask is not None:
        raster = ma.masked_array(raster, mask=mask)
    return raster


def stat_outmap_max(filenames, mask=None):
    result = np.zeros((len(filenames)))
    for i, filename in enumerate(filenames):
        raster = getRaster(filename)
        if mask is not None:
            raster = ma.masked_array(raster, mask=mask)
        result[i] = np.nanmax(raster)
    return result


def stat_outmap_sum(filenames, mask=None):
    result = np.zeros((len(filenames)))
    for i, filename in enumerate(filenames):
        raster = getRaster(filename)
        if mask is not None:
            raster = ma.masked_array(raster, mask=mask)
        result[i] = np.nansum(raster)
    return result


def stat_outmap_average(filenames, mask=None):
    result = np.zeros((len(filenames)))
    for i, filename in enumerate(filenames):
        raster = getRaster(filename)
        if mask is not None:
            raster = ma.masked_array(raster, mask=mask)
        result[i] = np.nanmean(raster)
    return result


def sumMaps(filenames, mask=None):
    result = getRaster(filenames[0])
    for filename in filenames[1:]:
        result += getRaster(filename)
    return ma.masked_array(result, mask)


def generateYearlyAverageEvaporation(directory: str) -> typing.NDArray:
    _, _, filenames = next(os.walk(directory))
    filenames = sorted(filenames)

    warumup_duration = CALIBRATION_START_DATE - SIMULATION_START_DATE
    warm_up_hours = int(divmod(warumup_duration.total_seconds(), 3600)[0])
    calibration_duration = SIMULATION_END_DATE - CALIBRATION_START_DATE
    calibration_hours = int(divmod(calibration_duration.total_seconds(), 3600)[0])

    prefixedFilenames = [
        path.join(directory, f)
        for f in filenames[warm_up_hours : warm_up_hours + calibration_hours + 1]
    ]
    catchmentMask = getCatchmentMask(CATCHMENT_FILE)
    result = sumMaps(prefixedFilenames, mask=catchmentMask)
    hoursInYear = 24 * 365.25
    result = result / (len(prefixedFilenames) / hoursInYear)
    return result


def statPotentialEvaporation(directory: str) -> typing.NDArray:
    _, _, filenames = next(os.walk(directory))
    filenames = sorted(filenames)

    warumup_duration = CALIBRATION_START_DATE - SIMULATION_START_DATE
    warm_up_hours = int(divmod(warumup_duration.total_seconds(), 3600)[0])
    calibration_duration = SIMULATION_END_DATE - CALIBRATION_START_DATE
    calibration_hours = int(divmod(calibration_duration.total_seconds(), 3600)[0])

    prefixedFilenames = [
        path.join(directory, f)
        for f in filenames[warm_up_hours : warm_up_hours + calibration_hours + 1]
    ]
    catchmentMask = getCatchmentMask(CATCHMENT_FILE)
    return stat_outmap_average(prefixedFilenames, mask=catchmentMask)


def stat_outmaps(directory):
    # Outmaps are generated from start datetime of the simulation until the end datetime
    # INCLUSIVE.
    # This means that start "2010-01-01 00:00" until "2010-01-02 00:00" will generate 25
    # maps.
    _, _, filenames = next(os.walk(directory))
    filenames = sorted(filenames)
    prefixed_files = {}
    for filename in filenames:
        prefix = filename[:3]
        if prefix not in prefixed_files:
            prefixed_files[prefix] = []
        prefixed_files[prefix].append(path.join(directory, filename))

    warumup_duration = CALIBRATION_START_DATE - SIMULATION_START_DATE
    warm_up_hours = int(divmod(warumup_duration.total_seconds(), 3600)[0])

    return {
        "precipitation_sum": stat_outmap_sum(prefixed_files["pre"][warm_up_hours:]),
        "precipitation_max": stat_outmap_max(prefixed_files["pre"][warm_up_hours:]),
        "precipitation_average": stat_outmap_average(
            prefixed_files["pre"][warm_up_hours:]
        ),
        "evaporation_average": stat_outmap_average(
            prefixed_files["pot"][warm_up_hours:]
        ),
        "temperature_average": stat_outmap_average(
            prefixed_files["tem"][warm_up_hours:]
        ),
    }


def eval_peak_distance(meas_data, sim_data, hq1, window):
    # find indices in meas data larger than hq1
    indices = np.where(meas_data > hq1)

    # create chuncks

    chunks = []

    chunk = []

    for i in range(0, len(indices[0]) - 1):
        # print(indices[0][i])
        if len(chunk) == 0:
            chunk.append(indices[0][i])
        if abs(indices[0][i] - indices[0][i + 1]) < window:
            chunk.append(indices[0][i + 1])
        else:
            chunks.append(chunk)
            chunk = []

    # reinit array blow up chuncks to measurement array
    chuncks_blown = []

    for chunk in chunks:
        min_val = min(chunk) - window
        max_val = max(chunk) + window

        range_new = list(range(min_val, max_val))

        chuncks_blown.append(range_new)

    # compute distances

    distances = []

    for chunck in chuncks_blown:
        val_meas_i = np.asarray(meas_data[chunck])
        val_sim_i = np.asarray(sim_data[chunck])

        max_scale_values = max(max(val_meas_i), max(val_sim_i))
        min_scale_values = min(min(val_meas_i), min(val_sim_i))

        # normalize
        time_norm = np.linspace(0, 1, len(chunck))
        val_meas_norm = (val_meas_i - min_scale_values) / (
            max_scale_values - min_scale_values
        )
        val_sim_norm = (val_sim_i - min_scale_values) / (
            max_scale_values - min_scale_values
        )

        # compute distance
        x_time_sim = time_norm[np.where(val_sim_norm == max(val_sim_norm))[0][0]]
        y_val_sim = val_sim_norm[np.where(val_sim_norm == max(val_sim_norm))[0][0]]
        x_time_meas = time_norm[np.where(val_meas_norm == max(val_meas_norm))[0][0]]
        y_val_meas = val_meas_norm[np.where(val_meas_norm == max(val_meas_norm))[0][0]]

        dist = np.sqrt(
            (x_time_meas - x_time_sim) ** 2.0 + (y_val_meas - y_val_sim) ** 2.0
        )
        distances.append(dist)

    eval_metric = np.mean(np.asarray(distances))

    return eval_metric


def read_data(timezone, filepath, lower, upper):
    # .timestamp() to convert to unix epoch
    # convert back to timestamp -> datetime.datetime.fromtimestamp(t0)
    time_list = []
    value_list = []
    with open(filepath, encoding="iso-8859-1") as f:
        for line in f:
            try:
                line = line.strip()
                line = line.split(";")
                date = line[0].strip(" ")
                value = line[1].strip(" ").replace(",", ".")

                date_obj = datetime.datetime.strptime(date, "%d.%m.%Y %H:%M:%S")

                time_localized = timezone.localize(date_obj)

                time_unix = time_localized.timestamp()

                if time_unix >= lower and time_unix <= upper:
                    time_list.append(time_unix)

                    if value.replace(".", "0").isdigit():
                        value_list.append(float(value))
                    else:
                        value_list.append(np.nan)

            except:
                continue

    value_list = np.asarray(value_list)
    nans, x = np.isnan(value_list), lambda z: z.nonzero()[0]
    value_list[nans] = np.interp(x(nans), x(~nans), value_list[~nans])

    gauge_data = np.column_stack((time_list, value_list))

    return gauge_data


def read_sim(filepath, lower, eval_lower, timestep, higher):
    # determine header
    f = open(filepath, "r")
    data = f.read()
    data = data.split("\n")
    header = int(data[1]) + 2
    f.close()

    discharge_all = np.genfromtxt(filepath, skip_header=header)

    time_list = []
    discharge_list = []

    for i in range(0, len(discharge_all)):
        timei = lower + i * timestep

        if timei >= eval_lower and timei <= higher:
            time_list.append(timei)
            discharge_list.append(discharge_all[i, 1:])

    time_list_np = np.asarray(time_list)
    discharge_np = np.asarray(discharge_list)

    return time_list_np, discharge_np


def get_calibration_column_names(df):
    result = []
    for c in df.columns:
        if c.startswith("calibrated_"):
            result.append(c)
    return result


def compute_data_frame():
    # gauges direcotry
    list_gauges = {
        "1": {
            "name": "kapfenberg/Diemlach-Muerz QOW3082",
            "path": r"/data/home/richi/master_thesis/model_MAR/gauges/Qow3082.csv",
            "sCatch": 1,
            "Index": 4,
        },
    }

    timezone = pytz.timezone("Etc/GMT-1")

    lower_obj = datetime.datetime.strptime("2013-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")
    # lower localized
    lower_localized = timezone.localize(lower_obj)
    # lwoer timestamp
    lower = lower_localized.timestamp()

    eval_lower_obj = datetime.datetime.strptime(
        "2014-01-01 00:00:00", "%Y-%m-%d %H:%M:%S"
    )
    # lower localized
    eval_lower_localized = timezone.localize(eval_lower_obj)
    # lwoer timestamp
    eval_lower = eval_lower_localized.timestamp()

    higher_obj = datetime.datetime.strptime("2017-12-28 00:00:00", "%Y-%m-%d %H:%M:%S")
    # lower localized
    higher_localized = timezone.localize(higher_obj)
    # lwoer timestamp
    higher = higher_localized.timestamp()

    timestep = 3600

    res_path = r"/data/home/richi/master_thesis/model_MAR/results_INCA/run.tss"
    res_path2 = r"/data/home/richi/master_thesis/model_MAR/results_ERA5_NSE/run.tss"
    res_path3 = r"/data/home/richi/master_thesis/model_MAR/results_ERA5_KGE/run.tss"

    # read INCA simulation
    time_sim, discharge_sim = read_sim(res_path, lower, eval_lower, timestep, higher)
    time_sim2, discharge_sim2 = read_sim(res_path2, lower, eval_lower, timestep, higher)
    time_sim3, discharge_sim3 = read_sim(res_path3, lower, eval_lower, timestep, higher)

    gauge_path = list_gauges["1"]["path"]
    gauge_index = int(list_gauges["1"]["Index"]) - 1
    gauge_data = read_data(timezone, gauge_path, eval_lower, higher)

    if (
        not np.array_equal(time_sim, time_sim2)
        or not np.array_equal(time_sim, time_sim3)
        or not np.array_equal(time_sim, gauge_data[:, 0])
    ):
        raise ValueError("time columns must be equal")

    date = gauge_data[:, 0]
    to_datetime = np.vectorize(datetime.datetime.fromtimestamp)
    datetimes = to_datetime(date)
    measured = gauge_data[:, 1]

    inca_stats = stat_outmaps(
        "/data/home/richi/master_thesis/model_MAR/results_INCA/outmaps"
    )
    era5_stats = stat_outmaps(
        "/data/home/richi/master_thesis/model_MAR/results_ERA5_NSE/outmaps"
    )
    era5_evaporation = statPotentialEvaporation("/data/home/richi/tmp/warped")

    return pd.DataFrame(
        {
            "date": datetimes,
            "measured": measured,
            "calibrated_inca_nse": discharge_sim[:, gauge_index],
            "calibrated_era5_nse": discharge_sim2[:, gauge_index],
            "calibrated_era5_kge": discharge_sim3[:, gauge_index],
            "inca_precipitation_sum": inca_stats["precipitation_sum"][: len(datetimes)],
            "inca_precipitation_max": inca_stats["precipitation_max"][: len(datetimes)],
            "inca_precipitation_average": inca_stats["precipitation_average"][
                : len(datetimes)
            ],
            "inca_evaporation_average": inca_stats["evaporation_average"][
                : len(datetimes)
            ],
            "inca_temperature_average": inca_stats["temperature_average"][
                : len(datetimes)
            ],
            "era5_precipitation_sum": era5_stats["precipitation_sum"][: len(datetimes)],
            "era5_precipitation_max": era5_stats["precipitation_max"][: len(datetimes)],
            "era5_precipitation_average": era5_stats["precipitation_average"][
                : len(datetimes)
            ],
            "era5_evaporation_average": era5_stats["evaporation_average"][
                : len(datetimes)
            ],
            "era5_evaporation_average_download": era5_evaporation[: len(datetimes)],
            "era5_temperature_average": era5_stats["temperature_average"][
                : len(datetimes)
            ],
        }
    )


def plot_df(df):
    column_names = get_calibration_column_names(df)
    fig, axs = plt.subplots(len(column_names), layout=LAYOUT, dpi=300)
    for index, name in enumerate(column_names):
        axs[index].plot(
            df["date"], df["measured"], "k--", linewidth=LINE_WIDTH, label="measured"
        )
        axs[index].plot(
            df["date"], df[name], "b", linewidth=LINE_WIDTH, label="simulated"
        )
        axs[index].grid()
        axs[index].legend()
    plt.show()
    plt.close(fig)


def plotDetails(df: pd.DataFrame, outfile: str = "") -> None:
    DETAIL_START_DATE = datetime.datetime(2016, 7, 1, 0, 0)
    DETAIL_END_DATE = datetime.datetime(2016, 8, 1, 0, 0)
    detail = df[(df["date"] >= DETAIL_START_DATE) & (df["date"] <= DETAIL_END_DATE)]

    fig, ax = plt.subplots(
        nrows=2, layout=LAYOUT, gridspec_kw={"height_ratios": [1, 3]}, dpi=300
    )
    fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN)

    ax[1].set_xlabel("Date")
    ax[1].set_xlim(left=DETAIL_START_DATE, right=DETAIL_END_DATE)
    plt.xticks(rotation=45, ha="right")

    ax[1].set_ylabel("Discharge [\\si{\\cubic\\meter\\per\\second}]")
    ax[1].plot(
        detail["date"],
        detail["measured"],
        linewidth=LINE_WIDTH,
        color="black",
        label="Measured",
    )
    ax[1].plot(
        detail["date"],
        detail["calibrated_inca_nse"],
        linewidth=LINE_WIDTH,
        color="blue",
        label="Discharge INCA",
    )
    ax[1].plot(
        detail["date"],
        detail["calibrated_era5_nse"],
        linewidth=LINE_WIDTH,
        color="red",
        label="Discharge ERA5",
    )
    ax[1].legend()

    ax[0].set_ylabel("Precipitation [\\si{\\milli\\meter}]")
    ax[0].set_xticks([])
    ax[0].set_xlim(left=DETAIL_START_DATE, right=DETAIL_END_DATE)
    ax[0].plot(
        detail["date"],
        detail["inca_precipitation_average"],
        linewidth=LINE_WIDTH,
        color="blue",
        label="Avg. Precipitation INCA",
    )
    ax[0].plot(
        detail["date"],
        detail["era5_precipitation_average"],
        linewidth=LINE_WIDTH,
        color="red",
        label="Avg. Precipitation ERA5",
    )
    ax[0].legend()

    if outfile:
        if outfile.endswith(".pgf"):
            plt.savefig(outfile, backend="pgf")
        else:
            plt.savefig(outfile)

    else:
        plt.show()
    plt.close(fig)


def plotWeather(df):
    fig, axs = plt.subplots(3, layout=LAYOUT, dpi=300)
    axs[0].set_title("Precipitation Sum")
    axs[0].plot(
        df["date"],
        df["inca_precipitation_sum"],
        "k--",
        linewidth=LINE_WIDTH,
        label="inca",
    )
    axs[0].plot(
        df["date"],
        df["era5_precipitation_sum"],
        "b",
        linewidth=LINE_WIDTH,
        label="era5",
    )
    axs[0].grid()
    axs[0].legend()
    axs[1].set_title("Temperature Average")
    axs[1].plot(
        df["date"],
        df["inca_temperature_average"],
        "k--",
        linewidth=LINE_WIDTH,
        label="inca",
    )
    axs[1].plot(
        df["date"],
        df["era5_temperature_average"],
        "b",
        linewidth=LINE_WIDTH,
        label="era5",
    )
    axs[1].grid()
    axs[1].legend()
    axs[2].set_title("Evaporation Average")
    axs[2].plot(
        df["date"],
        df["inca_evaporation_average"],
        "k--",
        linewidth=LINE_WIDTH,
        label="inca",
    )
    axs[2].plot(
        df["date"],
        df["era5_evaporation_average"],
        "b",
        linewidth=LINE_WIDTH,
        label="era5",
    )
    axs[2].grid()
    axs[2].legend()
    plt.show()
    plt.close(fig)


def parseDate(s):
    return datetime.datetime.strptime(s, "%Y-%m-%d %H:%M:%S")


def export_pgf(df, folder):
    plt.rcParams.update(
        {
            "pgf.preamble": r"\usepackage{siunitx}",
            # "pgf.texsystem": "pdflatex",
            "font.family": "serif",
            # "text.usetex": True,
            "pgf.rcfonts": False,
        }
    )
    for name in get_calibration_column_names(df):
        fig, axs = plt.subplots(layout=LAYOUT, dpi=300)
        fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN)
        plt.xlabel("Date")
        axs.set_xlim(left=df["date"].iloc[0], right=OFFICIAL_END_DATE)
        plt.xticks(rotation=45, ha="right")
        plt.ylabel("Discharge [\\si{\\cubic\\meter\\per\\second}]")
        axs.plot(df["date"], df["measured"], linewidth=LINE_WIDTH, label="Measured")
        axs.plot(df["date"], df[name], linewidth=LINE_WIDTH, label="Simulated")
        axs.grid()
        axs.legend()
        fig.savefig(path.join(folder, f"{name}.pgf"), backend="pgf")
        plt.close(fig)

    fig, axs = plt.subplots(layout=LAYOUT, dpi=300)
    fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN)
    plt.xlabel("Date")
    axs.set_xlim(left=df["date"].iloc[0], right=OFFICIAL_END_DATE)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Precipitation Sum [\\si{\\milli\\meter\\per\\hour}]")
    axs.plot(
        df["date"], df["inca_precipitation_sum"], linewidth=LINE_WIDTH, label="INCA"
    )
    axs.plot(
        df["date"], df["era5_precipitation_sum"], linewidth=LINE_WIDTH, label="ERA5"
    )
    axs.grid()
    axs.legend()
    fig.savefig(path.join(folder, "precipitation_sum.pgf"), backend="pgf")
    plt.close(fig)

    fig, axs = plt.subplots(layout=LAYOUT, dpi=300)
    fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN)
    plt.xlabel("Date")
    axs.set_xlim(left=df["date"].iloc[0], right=OFFICIAL_END_DATE)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("Average Precipitation [\\si{\\milli\\meter\\per\\hour}]")
    axs.plot(
        df["date"], df["inca_precipitation_average"], linewidth=LINE_WIDTH, label="INCA"
    )
    axs.plot(
        df["date"], df["era5_precipitation_average"], linewidth=LINE_WIDTH, label="ERA5"
    )
    axs.grid()
    axs.legend()
    fig.savefig(path.join(folder, "precipitation_average.pgf"), backend="pgf")
    plt.close(fig)

    fig, axs = plt.subplots(2, layout=LAYOUT, dpi=300)
    fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN * 1.5)

    axs[0].set_title("Average hourly temperature INCA")
    axs[0].set_xlabel("Date")
    axs[0].tick_params(labelrotation=45)
    axs[0].set_xlim(left=df["date"].iloc[0], right=OFFICIAL_END_DATE)
    axs[0].set_ylabel("Temperature [\\si{\\degreeCelsius}]")
    axs[0].plot(
        df["date"], df["inca_temperature_average"], linewidth=LINE_WIDTH, label="INCA"
    )
    axs[0].grid()

    axs[1].set_title("Average hourly temperature ERA5")
    axs[1].set_xlabel("Date")
    axs[1].tick_params(labelrotation=45)
    axs[1].set_xlim(left=df["date"].iloc[0], right=OFFICIAL_END_DATE)
    axs[1].set_ylabel("Temperature [\\si{\\degreeCelsius}]")
    axs[1].plot(
        df["date"], df["era5_temperature_average"], linewidth=LINE_WIDTH, label="ERA5"
    )
    axs[1].grid()
    fig.savefig(path.join(folder, "temperature_average.pgf"), backend="pgf")
    plt.close(fig)

    plotDetails(df, path.join(folder, "detail_201607.pgf"))

    fig, ax = plt.subplots(layout=LAYOUT, dpi=300)
    fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN)
    ax.axis("equal")
    ax.set_xlabel("Precipitation Sum INCA [mm]")
    ax.set_ylabel("Precipitation Sum ERA5 [mm]")
    minimum, maximum = scatterMinMax(
        df["inca_precipitation_sum"], df["era5_precipitation_sum"]
    )
    ax.plot([minimum, maximum], [minimum, maximum], color="black")
    ax.scatter(
        df["inca_precipitation_sum"],
        df["era5_precipitation_sum"],
        marker=".",
        color="blue",
    )
    fig.savefig(path.join(folder, "scatter_precipitation_sum.png"))
    plt.close(fig)

    fig, ax = plt.subplots(layout=LAYOUT, dpi=300)
    fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN)
    ax.axis("equal")
    ax.set_xlabel("Precipitation Maximum INCA [mm]")
    ax.set_ylabel("Precipitation Maximum ERA5 [mm]")
    minimum, maximum = scatterMinMax(
        df["inca_precipitation_max"], df["era5_precipitation_max"]
    )
    ax.plot([minimum, maximum], [minimum, maximum], color="black")
    ax.scatter(
        df["inca_precipitation_max"],
        df["era5_precipitation_max"],
        marker=".",
        color="blue",
    )
    fig.savefig(path.join(folder, "scatter_precipitation_max.png"))
    plt.close(fig)

    fig, ax = plt.subplots(layout=LAYOUT, dpi=300)
    fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN)
    ax.axis("equal")
    ax.set_xlabel("Precipitation Average INCA [mm]")
    ax.set_ylabel("Precipitation Average ERA5 [mm]")
    minimum, maximum = scatterMinMax(
        df["inca_precipitation_average"], df["era5_precipitation_average"]
    )
    ax.plot([minimum, maximum], [minimum, maximum], color="black")
    ax.scatter(
        df["inca_precipitation_average"],
        df["era5_precipitation_average"],
        marker=".",
        color="blue",
    )
    fig.savefig(path.join(folder, "scatter_precipitation_average.png"))
    plt.close(fig)

    fig, ax = plt.subplots(layout=LAYOUT, dpi=300)
    fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN)
    ax.axis("equal")
    ax.set_xlabel("Temperature average INCA [$^\\circ$C]")
    ax.set_ylabel("Temperature average ERA5 [$^\\circ$C]")
    minimum, maximum = scatterMinMax(
        df["inca_temperature_average"], df["era5_temperature_average"]
    )
    ax.plot([minimum, maximum], [minimum, maximum], color="black")
    ax.scatter(
        df["inca_temperature_average"],
        df["era5_temperature_average"],
        marker=".",
        color="blue",
    )
    fig.savefig(path.join(folder, "scatter_temperature_average.png"))
    plt.close(fig)

    fig, ax = plt.subplots(layout=LAYOUT, dpi=300)
    fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN)
    ax.axis("equal")
    ax.set_xlabel("Potential evaporation average INCA [mm]")
    ax.set_ylabel("Potential evaporation average ERA5 [mm]")
    minimum, maximum = scatterMinMax(
        df["inca_evaporation_average"], df["era5_evaporation_average"]
    )
    ax.plot([minimum, maximum], [minimum, maximum], color="black")
    ax.scatter(
        df["inca_evaporation_average"],
        df["era5_evaporation_average"],
        marker=".",
        color="blue",
    )
    fig.savefig(path.join(folder, "scatter_evaporation_average.png"))
    plt.close(fig)

    fig, ax = plt.subplots(layout=LAYOUT, dpi=300)
    fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN)
    ax.axis("equal")
    ax.set_xlabel("Potential evaporation average ERA5 calculated [mm]")
    ax.set_ylabel("Potential evaporation average ERA5 downloaded [mm]")
    minimum, maximum = scatterMinMax(
        df["era5_evaporation_average"], df["era5_evaporation_average_download"]
    )
    ax.plot([minimum, maximum], [minimum, maximum], color="black")
    ax.scatter(
        df["era5_evaporation_average"],
        df["era5_evaporation_average_download"],
        marker=".",
        color="blue",
    )
    fig.savefig(path.join(folder, "scatter_evaporation_average_era5_era5_download.png"))
    plt.close(fig)

    # plot precipitation grid difference
    df["prod"] = df["inca_precipitation_sum"] * df["era5_precipitation_sum"]
    id = df["prod"].idxmax()
    duration = df.at[id, "date"] - SIMULATION_START_DATE
    hours = divmod(duration.total_seconds(), 60**2)[0]
    filename = f"pre{hours/1000:09.3f}"

    inca_prec_raster = getRaster(
        path.join(
            "/data/home/richi/master_thesis/model_MAR/results_INCA/outmaps", filename
        )
    )
    era5_prec_raster = getRaster(
        path.join(
            "/data/home/richi/master_thesis/model_MAR/results_ERA5_NSE/outmaps",
            filename,
        )
    )
    maximum = max(np.nanmax(inca_prec_raster), np.nanmax(era5_prec_raster))

    fig, ax = plt.subplots(2, layout=LAYOUT, dpi=300)
    fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN * 1.5)

    ax[0].set_title("ERA5 Precipitation")
    ax[0].tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)
    dem_ax = ax[0].imshow(era5_prec_raster, interpolation="none", vmin=0, vmax=maximum)
    fig.colorbar(
        dem_ax,
        ax=ax[0],
        label="Precipitation [\\si{\\milli\\meter}]",
        shrink=0.7,
    )

    ax[1].set_title("INCA Precipitation")
    ax[1].tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)
    dem_ax = ax[1].imshow(inca_prec_raster, interpolation="none", vmin=0, vmax=maximum)
    fig.colorbar(
        dem_ax,
        ax=ax[1],
        label="Precipitation [\\si{\\milli\\meter}]",
        shrink=0.7,
    )

    fig.savefig(path.join(folder, "precipitation_grid_era5_inca.pgf"), backend="pgf")
    plt.close(fig)


def scatterMinMax(x: np.ndarray, y: np.ndarray):
    return min(x.min(), y.min()), max(x.max(), y.max())


def print_stats(df):
    for name in get_calibration_column_names(df):
        print(f"= {name}")
        print(f"NSE = {nse(df['measured'], df[name]):.2f}")
        print(f"NCE = {nce(df['measured'], df[name]):.2f}")
        print(f"KGE = {kge(df['measured'], df[name]):.2f}")

    print(
        "highest precipitation sum: "
        + f"INCA: {df['inca_precipitation_sum'].max():.2f} "
        + f"ERA5: {df['era5_precipitation_sum'].max():.2f}"
    )
    print(
        "highest precipitation max: "
        + f"INCA: {df['inca_precipitation_max'].max():.2f} "
        + f"ERA5: {df['era5_precipitation_max'].max():.2f}"
    )
    print(
        "highest precipitation average: "
        + f"INCA: {df['inca_precipitation_average'].max():.2f} "
        + f"ERA5: {df['era5_precipitation_average'].max():.2f}"
    )
    print(
        "highest evaporation average: "
        + f"INCA: {df['inca_evaporation_average'].max():.2f} "
        + f"ERA5: {df['era5_evaporation_average'].max():.2f} "
        + f"ERA5 download: {df['era5_evaporation_average_download'].max():.2f}"
    )

    precipitation_sum = df[["inca_precipitation_sum", "era5_precipitation_sum"]]
    precipitation_sum_corr = precipitation_sum.corr().loc[
        "inca_precipitation_sum", "era5_precipitation_sum"
    ]
    print(f"correlation precipitation sum: {precipitation_sum_corr:.2f}")

    precipitation_max = df[["inca_precipitation_max", "era5_precipitation_max"]]
    precipitation_max_corr = precipitation_max.corr().loc[
        "inca_precipitation_max", "era5_precipitation_max"
    ]
    print(f"correlation precipitation max: {precipitation_max_corr:.2f}")

    precipitation_average = df[
        ["inca_precipitation_average", "era5_precipitation_average"]
    ]
    precipitation_average_corr = precipitation_average.corr().loc[
        "inca_precipitation_average", "era5_precipitation_average"
    ]
    print(f"correlation precipitation average: {precipitation_average_corr:.2f}")

    temperature_average = df[["inca_temperature_average", "era5_temperature_average"]]
    temperature_average_corr = temperature_average.corr().loc[
        "inca_temperature_average", "era5_temperature_average"
    ]
    print(f"correlation temperature average: {temperature_average_corr:.2f}")

    evaporation_average = df[["inca_evaporation_average", "era5_evaporation_average"]]
    evaporation_average_corr = evaporation_average.corr().loc[
        "inca_evaporation_average", "era5_evaporation_average"
    ]
    print(f"correlation evaporation average: {evaporation_average_corr:.2f}")

    evaporation_average = df[
        ["inca_evaporation_average", "era5_evaporation_average_download"]
    ]
    evaporation_average_corr = evaporation_average.corr().loc[
        "inca_evaporation_average", "era5_evaporation_average_download"
    ]
    print(
        f"correlation evaporation average (inca, era5 download): {evaporation_average_corr:.2f}"
    )

    evaporation_average = df[
        ["era5_evaporation_average", "era5_evaporation_average_download"]
    ]
    evaporation_average_corr = evaporation_average.corr().loc[
        "era5_evaporation_average", "era5_evaporation_average_download"
    ]
    print(
        f"correlation evaporation average (era5, era5 download): {evaporation_average_corr:.2f}"
    )

    df["prod"] = df["inca_precipitation_sum"] * df["era5_precipitation_sum"]
    id = df["prod"].idxmax()
    print(f"date where most sources have the most precipitation: {df.at[id, 'date']}")


def writeData(df: pd.DataFrame, filename: str) -> None:
    if filename.endswith(".csv"):
        df.to_csv(filename, index=False)
    df.to_pickle(filename)


def readData(filename: str) -> pd.DataFrame:
    if filename.endswith(".csv"):
        return pd.read_csv(filename, parse_dates=["date"], date_parser=parseDate)
    return pd.read_pickle(filename)


def main():
    gdal.UseExceptions()

    parser = ArgumentParser()
    parser.add_argument("-o", "--output", metavar="file", help="write data to file")
    parser.add_argument("-p", "--plot", action="store_true", help="show plot")
    parser.add_argument("-s", "--stats", action="store_true", help="print stats")
    parser.add_argument(
        "-l", "--latex", metavar="folder", help="write pgf files to folder"
    )
    parser.add_argument("-i", "--input", metavar="file", help="use data from file")
    parser.add_argument(
        "-w", "--weather", action="store_true", help="plot weather data"
    )
    parser.add_argument("-d", "--details", action="store_true", help="plot details")
    parser.add_argument("-t", "--test", action="store_true", help="debugging")
    args = parser.parse_args()

    df = compute_data_frame() if not args.input else readData(args.input)

    if args.plot:
        plot_df(df)
    if args.latex:
        export_pgf(df, args.latex)
    if args.output:
        writeData(df, args.output)
    if args.stats:
        print_stats(df)
    if args.weather:
        plotWeather(df)
    if args.details:
        plotDetails(df)
    if args.test:
        from icecream import ic

        # test = getRaster("/data/home/richi/tmp/warped/Muerz_PEV_20170506_2100.tiff")
        #
        # ic(np.nanmin(test), np.nanmax(test))
        # return

        fig, ax = plt.subplots(layout=LAYOUT, dpi=300)
        fig.set_size_inches(w=TEXTWITH_IN, h=TEXTWITH_IN)
        ax.tick_params(bottom=False, labelbottom=False, left=False, labelleft=False)

        eva_raster = generateYearlyAverageEvaporation("/data/home/richi/tmp/warped")
        ic(np.nanmin(eva_raster), np.nanmax(eva_raster))

        dem_ax = ax.imshow(eva_raster, interpolation="none")
        label = "Average yearly evaporation [\\si{\\milli\\meter}]"
        fig.colorbar(
            dem_ax,
            ax=ax,
            label=label,
            orientation="horizontal",
        )

        plt.show()
        # fig.savefig(path.join(folder, "precipitation_grid_inca.pgf"), backend="pgf")
        plt.close(fig)


if __name__ == "__main__":
    main()
