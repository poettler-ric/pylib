#!/usr/bin/env python3

from argparse import ArgumentParser
from logging import info
from os.path import join as pjoin
import logging

import cdsapi


ALL_MONTHS = [
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "11",
    "12",
]
ALL_DAYS = [
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "23",
    "24",
    "25",
    "26",
    "27",
    "28",
    "29",
    "30",
    "31",
]
ALL_HOURS = [
    "00:00",
    "01:00",
    "02:00",
    "03:00",
    "04:00",
    "05:00",
    "06:00",
    "07:00",
    "08:00",
    "09:00",
    "10:00",
    "11:00",
    "12:00",
    "13:00",
    "14:00",
    "15:00",
    "16:00",
    "17:00",
    "18:00",
    "19:00",
    "20:00",
    "21:00",
    "22:00",
    "23:00",
]
AUSTRIA = [
    49,
    9,
    46,
    19,
]
VALID_VARIABLES = [
    "2m_temperature",
    "potential_evaporation",
    "total_precipitation",
    "surface_solar_radiation_downwards",
]


def main():
    parser = ArgumentParser(description="Downloads ERA5 grib files")
    parser.add_argument(
        "variable", choices=VALID_VARIABLES, help="variable to download"
    )
    parser.add_argument("year", type=int, help="year to download")
    parser.add_argument("output_folder", help="folder where to put the download")
    args = parser.parse_args()

    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    log_formatter = logging.Formatter("%(levelname)s %(asctime)s: %(message)s")
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)

    filename = pjoin(args.output_folder, f"{args.variable}_{args.year}.grib")

    c = cdsapi.Client()
    info(f"Starting download ({args.variable} for {args.year})")
    c.retrieve(
        "reanalysis-era5-land",
        {
            "format": "grib",
            "variable": args.variable,
            "year": args.year,
            "month": ALL_MONTHS,
            "day": ALL_DAYS,
            "time": ALL_HOURS,
            "area": AUSTRIA,
        },
        filename,
    )
    info("Download finished")


if __name__ == "__main__":
    main()
