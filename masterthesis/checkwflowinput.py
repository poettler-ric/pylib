#!/usr/bin/env python3

__author__ = "Richard Pöttler"
__copyright__ = "Copyright (c) 2022 Richard Pöttler"
__license__ = "MIT"
__email__ = "richard.poettler@gmail.com"


from datetime import datetime, timedelta
from logging import info, debug
from os import listdir
import logging
import numpy as np
import os.path
import pcraster as pcr
import re

TEMPERATURE_MIN = -10  # oC
TEMPERATURE_MAX = 35  # oC

PRECIPITATION_SUM = 1200  # mm/y
PRECIPITATION_THRESHOLD = 100  # mm/y
EVAPOTRANSPIRATION_SUM = 600  # mm/y
EVAPOTRANSPIRATION_THRESHOLD = 50  # mm/y

INMAP_FOLDER = "/data/home/richi/master_thesis/wflow_model/inmaps_calibration"


START_TIME = "2013-09-01 00:00:00"
END_TIME = "2015-02-28 23:00:00"
ONE_YEAR = timedelta(days=365.25)


def get_map_min_max(m):
    map_np = pcr.pcr2numpy(pcr.readmap(m), np.NaN)
    return np.nanmin(map_np), np.nanmax(map_np)


def sum_inmaps(folder, filepattern):
    map_sum = None
    abs_min = np.inf
    abs_max = np.NINF

    for filename in listdir(folder):
        if filepattern.match(filename):
            map_np = pcr.pcr2numpy(pcr.readmap(os.path.join(folder, filename)), np.NaN)
            abs_min = np.nanmin([abs_min, np.nanmin(map_np)])
            abs_max = np.nanmax([abs_max, np.nanmax(map_np)])
            if map_sum is None:
                map_sum = map_np
            else:
                map_sum += map_np
    return map_sum, abs_min, abs_max


def is_temperature_valid():
    is_valid = True
    abs_min = np.inf
    abs_max = np.NINF

    filepattern = re.compile(r"^TEMP\d+.\d+$")
    for filename in listdir(INMAP_FOLDER):
        if filepattern.match(filename):
            temp_min, temp_max = get_map_min_max(os.path.join(INMAP_FOLDER, filename))
            abs_min = np.nanmin([abs_min, temp_min])
            abs_max = np.nanmax([abs_max, temp_max])

            if temp_min < TEMPERATURE_MIN or temp_max > TEMPERATURE_MAX:
                is_valid = False
                debug(f"temperature: {filename} min: {temp_min} max: {temp_max}")

    return is_valid, abs_min, abs_max


def is_precipitation_valid():
    prec_min = np.nan
    prec_max = np.nan
    is_valid = False

    start_time = datetime.strptime(START_TIME, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.strptime(END_TIME, "%Y-%m-%d %H:%M:%S")
    duration = end_time - start_time

    filepattern = re.compile(r"^P\d+.\d+$")
    precipitation_sum, local_min, local_max = sum_inmaps(INMAP_FOLDER, filepattern)
    if precipitation_sum is not None:
        prec_min = ONE_YEAR / duration * np.nanmin(precipitation_sum)
        prec_max = ONE_YEAR / duration * np.nanmax(precipitation_sum)

        is_valid = (
            prec_min >= PRECIPITATION_SUM - PRECIPITATION_THRESHOLD
            and prec_max <= PRECIPITATION_SUM + PRECIPITATION_THRESHOLD
        )

    debug(
        f"precipitation min: {prec_min} max: {prec_max} local_min: {local_min} local_max: {local_max}"
    )
    return is_valid, prec_min, prec_max


def is_evapotranspiration_valid():
    evap_min = np.nan
    evap_max = np.nan
    is_valid = False

    start_time = datetime.strptime(START_TIME, "%Y-%m-%d %H:%M:%S")
    end_time = datetime.strptime(END_TIME, "%Y-%m-%d %H:%M:%S")
    duration = end_time - start_time

    filepattern = re.compile(r"^PET\d+.\d+$")
    evapotranspiration_sum, local_min, local_max = sum_inmaps(INMAP_FOLDER, filepattern)
    if evapotranspiration_sum is not None:
        evap_min = ONE_YEAR / duration * np.nanmin(evapotranspiration_sum)
        evap_max = ONE_YEAR / duration * np.nanmax(evapotranspiration_sum)

        is_valid = (
            evap_min >= EVAPOTRANSPIRATION_SUM - EVAPOTRANSPIRATION_THRESHOLD
            and evap_max <= EVAPOTRANSPIRATION_SUM + EVAPOTRANSPIRATION_THRESHOLD
        )

    debug(
        f"evapotranspiration min: {evap_min} max: {evap_max} local_min: {local_min} local_max: {local_max}"
    )
    return is_valid, evap_min, evap_max


def main():
    logging.basicConfig(level=logging.INFO)
    # logging.basicConfig(level=logging.DEBUG)

    is_valid, abs_min, abs_max = is_temperature_valid()
    if not is_valid:
        info(
            f"temperature not valid. min: {abs_min} ({TEMPERATURE_MIN}) max: {abs_max} ({TEMPERATURE_MAX})"
        )
    is_valid, abs_min, abs_max = is_precipitation_valid()
    if not is_valid:
        info(
            f"precipitaion not valid. min: {abs_min} ({PRECIPITATION_SUM - PRECIPITATION_THRESHOLD}) "
            + f"max: {abs_max} ({PRECIPITATION_SUM + PRECIPITATION_THRESHOLD})"
        )
    is_valid, abs_min, abs_max = is_evapotranspiration_valid()
    if not is_valid:
        info(
            f"evapotranspiration not valid. min: {abs_min} ({EVAPOTRANSPIRATION_SUM - EVAPOTRANSPIRATION_THRESHOLD}) "
            + f"max: {abs_max} ({EVAPOTRANSPIRATION_SUM + EVAPOTRANSPIRATION_THRESHOLD})"
        )


if __name__ == "__main__":
    main()
