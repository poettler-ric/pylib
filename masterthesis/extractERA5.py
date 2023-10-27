#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 09:48:06 2023

@author: sebastian

precipitation in m -> *1000 for mm
temperature in kelvin -> - 273.15 for celsius
radiation in J/m2 -> /timestep for W/m2

temperature = t2m
precipitation = tp
radiation = ssrd

"""

from datetime import datetime, timedelta

import netCDF4
import rasterio

# netcdf files
# first netcdf file
# netcdf_file_path = '/home/sebastian/working_dir/wflow/model_MAR/ForcingData/ERA5/indownload/part1_2010_2013.nc'
# second netcdf file
netcdf_file_path = "/data/home/richi/master_thesis/model_MAR/ForcingData/ERA5/indownload/part2_2014_sep2017.nc"

# outpaths
RR_out = r"/data/home/richi/master_thesis/model_MAR/ForcingData/ERA5/RR_calib"
GL_out = r"/data/home/richi/master_thesis/model_MAR/ForcingData/ERA5/GL_calib"
TT_out = r"/data/home/richi/master_thesis/model_MAR/ForcingData/ERA5/TT_calib"

# process first netcdf
# open dataset
dataset = netCDF4.Dataset(netcdf_file_path, "r")
timestep_seconds = 3600
reference_date = datetime(1900, 1, 1)

# extract time
time = dataset.variables["time"][:]
result_datetime = [
    reference_date + timedelta(seconds=float(i * 3600))
    for i in dataset.variables["time"][:]
]

# loop over datetime
for i in range(0, len(result_datetime)):
    # construct names
    rr_path = RR_out + "/" + result_datetime[i].strftime("Muerz_RR_%Y%m%d_%H%M.tif")
    tt_path = TT_out + "/" + result_datetime[i].strftime("Muerz_TT_%Y%m%d_%H%M.tif")
    gl_path = GL_out + "/" + result_datetime[i].strftime("Muerz_GL_%Y%m%d_%H%M.tif")

    # read lat lon
    lat = dataset.variables["latitude"][:]  # is y
    lon = dataset.variables["longitude"][:]  # is x

    transform = None  # Replace with the appropriate transformation matrix if available
    crs_data = "EPSG:4326"  # Replace with the correct EPSG code

    transform = rasterio.transform.from_origin(
        lon.min() - abs(lon[0] - lon[1]) / 2,
        lat.max() + abs(lat[0] - lat[1]) / 2,  # Upper-left coordinates
        abs(lon[0] - lon[1]),  # Pixel width
        abs(lat[0] - lat[1]),  # Pixel height
    )

    # extract variables
    rr = dataset.variables["tp"][i] * 1000.0
    gl = dataset.variables["ssrd"][i] / timestep_seconds
    tt = dataset.variables["t2m"][i] - 273.15

    rr[rr < 0] = 0
    gl[gl < 0] = 0

    # set valaues to zero if necessary

    # out raseter
    with rasterio.open(
        rr_path,
        "w",
        driver="GTiff",
        height=rr.shape[0],
        width=rr.shape[1],
        count=1,  # Number of bands
        dtype=rr.dtype,
        crs=crs_data,
        transform=transform,
    ) as dst:
        dst.write(rr, 1)

    with rasterio.open(
        tt_path,
        "w",
        driver="GTiff",
        height=tt.shape[0],
        width=tt.shape[1],
        count=1,  # Number of bands
        dtype=rr.dtype,
        crs=crs_data,
        transform=transform,
    ) as dst:
        dst.write(tt, 1)

    with rasterio.open(
        gl_path,
        "w",
        driver="GTiff",
        height=gl.shape[0],
        width=gl.shape[1],
        count=1,  # Number of bands
        dtype=rr.dtype,
        crs=crs_data,
        transform=transform,
    ) as dst:
        dst.write(gl, 1)
