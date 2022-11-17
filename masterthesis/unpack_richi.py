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
import multiprocessing as mp


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
        "/home/richi/mnt/zamg_inca/RR",
        "/home/richi/mnt/tiff_tmp/RR",
        "Austria_RR_%Y%m%d_%H%M.tif",
        0.001,
    ),
    Variable(
        "T2M",
        "/home/richi/mnt/zamg_inca/T2M",
        "/home/richi/mnt/tiff_tmp/T2M",
        "Austria_TT_%Y%m%d_%H%M.tif",
        0.01,
    ),
    Variable(
        "GL",
        "/home/richi/mnt/zamg_inca/GL",
        "/home/richi/mnt/tiff_tmp/GL",
        "Austria_GL_%Y%m%d_%H%M.tif",
        0.01,
    ),
]


def unpack_file(net_path, i, var, time_array):
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

def get_time_array(net_path):
    ds = nc.Dataset(net_path)
    return nc.num2date(ds["time"][:], ds["time"].units)


def main():
    with mp.Pool(mp.cpu_count()) as pool:
        for var in extract_variables:
            for root, _, files in os.walk(var.nc_path):
                for f in sorted(files):
                    net_path = pjoin(root, f)
                    time_array = get_time_array(net_path)
                    results = []

                    # time loop
                    for i in range(0, len(time_array)):
                        results.append(pool.apply_async(unpack_file, (net_path, i, var, time_array)))

                    [result.wait() for result in results]


if __name__ == "__main__":
    main()
