#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:02:48 2022

@author: iwbworkstation
"""

from osgeo import gdal
from os.path import join as pjoin
import netCDF4 as nc
import os


class Variable:
    def __init__(self, name, nc_path, tiff_path, tiff_filename_pattern, factor):
        self.name = name
        self.nc_path = nc_path
        self.tiff_path = tiff_path
        self.tiff_filename_pattern = tiff_filename_pattern
        self.factor = factor


# name in INCA File, Destination path,output format, scaling factor
extract_variables = [
    Variable(
        "RR",
        "/data/home/richi/data/weather/austria/muerz/inca_data/RR",
        "/data/home/richi/tmp/inca_data/RR",
        "Muerz_RR_%Y%m%d_%H%M.tif",
        0.001,
    ),
    Variable(
        "T2M",
        "/data/home/richi/data/weather/austria/muerz/inca_data/T2M",
        "/data/home/richi/tmp/inca_data/T2M",
        "Muerz_TT_%Y%m%d_%H%M.tif",
        0.01,
    ),
    Variable(
        "GL",
        "/data/home/richi/data/weather/austria/muerz/inca_data/GL",
        "/data/home/richi/tmp/inca_data/GL",
        "Muerz_GL_%Y%m%d_%H%M.tif",
        0.01,
    ),
]

for var in extract_variables:
    for root, dirs, files in os.walk(var.nc_path):
        for f in sorted(files):
            net_path = pjoin(root, f)
            ds = nc.Dataset(net_path)
            time = ds["time"][:]
            time_array = nc.num2date(ds["time"][:], ds["time"].units)

            ds = None
            # time loop
            for i in range(0, len(time_array)):
                ds = gdal.Open(f"NETCDF:{net_path}:{var.name}")
                # get respective raster band
                band = ds.GetRasterBand(i + 1)
                # retrive data
                data = band.ReadAsArray(0, 0, ds.RasterXSize, ds.RasterYSize).astype(
                    float
                )
                # scale parameters
                data = data * var.factor

                # create empty GTIFF
                # build filepath
                tif_path = pjoin(
                    var.tiff_path, time_array[i].strftime(var.tiff_filename_pattern)
                )
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
