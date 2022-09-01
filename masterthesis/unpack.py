#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:02:48 2022

@author: iwbworkstation

# INCA projection 
# proj4: +proj=lcc +lat_0=47.5 +lon_0=13.33333 +lat_1=49 +lat_2=46 +x_0=400000 +y_0=400000 +datum=WGS84 +units=m +no_defs +type=crs

"""

import netCDF4 as nc
from netCDF4 import date2num, num2date, Dataset
from osgeo import gdal
import numpy as np
import os
import re

# breakout for testing
breakout = None

# standard file format and desired years
fileformat = "incal-hourly_%Y0101T0000_%Y0131T2300.nc"
years = ["2014", "2015"]

filepath = (
    r"/home/richi/master_thesis/tmp_weather/GL_RR_T2M_20140101-20150228"
)

# name in INCA File, Destination path,output format, scaling factor
extract_variables = [
    [
        "RR",
        "/home/richi/master_thesis/tmp_weather/inca_data/RR",
        "Muerz_RR_%Y%m%d_%H%M.tif",
        "0.001",
    ],
    [
        "T2M",
        "/home/richi/master_thesis/tmp_weather/inca_data/T2M",
        "Muerz_TT_%Y%m%d_%H%M.tif",
        "0.01",
    ],
    [
        "GL",
        "/home/richi/master_thesis/tmp_weather/inca_data/GL",
        "Muerz_GL_%Y%m%d_%H%M.tif",
        "0.01",
    ],
]

i = 0
for subdir, dirs, files in os.walk(filepath):
    # loops of files in folder AND SORTES them alphabetically
    # metadata is used from precipitation
    # find occurences of year %Y
    years_ext = [m.start() for m in re.finditer("%Y", fileformat)]
    years_ext = [int(years_ext[i] + i * 2.0) for i in range(0, len(years_ext))]
    # file loop
    for file in sorted(files):
        # check if file contains our date
        # i know, this implementation could cause some problems at certain constellations
        year_list = [file[start : start + 4] for start in years_ext]
        # check if any year is present in list
        if any(year in year_list for year in years):
            # process netcdf file
            # get time of netcdf file
            net_path = subdir + "/" + file
            ds = nc.Dataset(net_path)
            time = ds["time"][:]
            time_array = num2date(ds["time"][:], ds["time"].units)

            ds = None
            # time loop
            for t in range(0, len(time_array)):
                # variable loop
                current_time = time_array[t].strftime(extract_variables[0][2])

                for var in extract_variables:
                    # read variable
                    ds = gdal.Open("NETCDF:{0}:{1}".format(net_path, var[0]))
                    # get respective raster band
                    band = ds.GetRasterBand(t + 1)
                    # retrive data
                    data = band.ReadAsArray(
                        0, 0, ds.RasterXSize, ds.RasterYSize
                    ).astype(np.float)
                    # scale parameters
                    data = data * float(var[3])

                    # create empty GTIFF
                    # build filepath
                    tif_path = var[1] + "/" + current_time
                    # create driver
                    driver = gdal.GetDriverByName("GTiff")
                    # create ds out
                    ds_out = driver.Create(
                        tif_path, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Float32
                    )
                    # get and set geotransform and projection
                    geotrans = ds.GetGeoTransform()
                    proj = ds.GetProjection()
                    ds_out.SetGeoTransform(geotrans)
                    ds_out.SetProjection(proj)

                    # write data to raster
                    ds_out.GetRasterBand(1).WriteArray(data)

                    ds = None
                    ds_out = None

            i += 1

        if i == breakout:
            break
