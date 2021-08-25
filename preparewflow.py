#!/bin/env python3

from argparse import ArgumentParser
from configparser import ConfigParser, ExtendedInterpolation
from json import loads
from logging import info, debug
from pyproj import Transformer
from scipy.interpolate import NearestNDInterpolator
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

import logging
import math
import numpy as np
import os
import pcraster as pcr
import shapefile
import xarray as xr

# import gdal

DEFAULT_MAX_STREAMORDER = 4

"""Mapping from log level strings to logging levels"""
LOG_LEVEL_MAP = {
    "critical": logging.CRITICAL,
    "error": logging.ERROR,
    "warning": logging.WARNING,
    "info": logging.INFO,
    "debug": logging.DEBUG,
}


# generates a map of all the cell centers of the rasters
def init_cellcenter(rows, cols, cell_size, xmin, ymin):

    cell_centers = np.zeros((rows, cols, 2))

    for i in range(0, rows):
        for j in range(0, cols):
            cell_centers[i][j][0] = xmin + cell_size / 2.0 + cell_size * j
            cell_centers[i][j][1] = ymin - cell_size / 2.0 - cell_size * i

    return cell_centers


def gen_inpoly(shapefile, coord_map, rows, cols):

    polygon = Polygon(shapefile["coordinates"][0])

    raster = np.zeros((rows, cols))

    i = 0
    for i in range(0, rows):
        j = 0
        for j in range(0, cols):
            point = Point(coord_map[i][j][0], coord_map[i][j][1])
            if polygon.contains(point):
                raster[i][j] = 1

    return raster


def generate_river_points(shapefile, cell_size):

    # read shape features
    points_2D = []

    for i in range(0, len(shapefile.shapes())):
        feature = shapefile.shapes()[i].__geo_interface__["coordinates"]

        # resamples to 80% of the raster size
        d = np.diff(feature, axis=0)
        segdists = np.hypot(d[:, 0], d[:, 1])
        divisions = np.ceil(segdists / (cell_size * 0.8))
        points_2D.append((feature[0][0], feature[0][1]))
        for j in range(0, len(feature) - 1):
            x1 = feature[j][0]
            x2 = feature[j + 1][0]
            y1 = feature[j][1]
            y2 = feature[j + 1][1]
            n = int(divisions[j])
            for i in range(1, n):
                a = float(i) / n
                x = (1 - a) * x1 + a * x2
                y = (1 - a) * y1 + a * y2
                points_2D.append((x, y))
            points_2D.append((x2, y2))

    points_2D = np.asarray(points_2D)
    return points_2D


# TODO: delete cell_size
def burn_in_river(cell_centers, rows, cols, riv_points, cell_size):

    river_array = np.empty((rows, cols))
    river_array[:] = np.NaN

    for i in range(0, len(riv_points)):
        riv1 = np.sqrt(
            (cell_centers[:, :, 0] - riv_points[i][0]) ** 2.0
            + (cell_centers[:, :, 1] - riv_points[i][1]) ** 2.0
        )
        minriv = np.amin(abs(riv1))
        index = np.where(riv1 == minriv)
        river_array[index[0][0]][index[1][0]] = 1

    return river_array


# TODO: delete cell_size
def burn_coords(cell_centers, rows, cols, riv_points, cell_size):

    river_array = np.empty((rows, cols))
    river_array[:] = -9999

    for i in range(0, len(riv_points)):
        riv1 = np.sqrt(
            (cell_centers[:, :, 0] - riv_points[i][0]) ** 2.0
            + (cell_centers[:, :, 1] - riv_points[i][1]) ** 2.0
        )
        minriv = np.amin(abs(riv1))
        index = np.where(riv1 == minriv)
        river_array[index[0][0]][index[1][0]] = i + 1

    return river_array


def gen_river_connectivity(river_array, rows, cols):

    river_array_corrected = np.copy(river_array)

    for i in range(0, rows):
        for j in range(0, cols):
            if river_array[i][j] == 1:
                if river_array[i - 1][j - 1] == 1 or river_array[i - 1][j + 1] == 1:
                    if math.isnan(river_array[i - 1][j]) == True:
                        river_array_corrected[i - 1][j] = 1

    return river_array_corrected


def read_soil_to_dict(soils_folder):
    for subdir, dirs, files in os.walk(soils_folder):
        for file in files:
            filepath = os.path.join(subdir, file)
            map = pcr.readmap(filepath)
            map_np = pcr.pcr2numpy(map, 0.0)
            strings = filepath.split("/")
            mapstring = strings[-1]
            namestring = mapstring[:2]
            depthstring = strings[-2]
            dictionary[depthstring][namestring] = map_np

    # Populate variables

    # Populate variables uniform
    i = 0
    for key in dictionary:
        # Create outermost layer
        Cli = dictionary[key]["CL"] / 1000.0
        Cli[Cli == 0] = np.median(Cli[Cli > 0])
        SAi = dictionary[key]["SA"] / 1000.0
        SAi[SAi == 0] = np.median(SAi[SAi > 0])
        SIi = dictionary[key]["SI"] / 1000.0
        SIi[SIi == 0] = np.median(SIi[SIi > 0])
        BDi = dictionary[key]["BD"] * 0.01
        BDi[BDi == 0] = np.median(BDi[BDi > 0])
        OCi = dictionary[key]["OC"] / 10000.0
        OCi[OCi == 0] = np.median(OCi[OCi > 0])
        thetaRi = (
            0.09878
            + 0.002127 * Cli
            - (8.366 * 10 ** -4) * SIi
            - 0.0767 / (OCi + 1)
            + SIi * Cli * (3.853 * 10 ** -5)
            + 0.00233 * Cli / (OCi + 1)
            + 9.498 * 10 ** -4 * SIi / (OCi + 1)
        )
        thetaSi = (
            0.6819
            + 0.06480 / (OCi + 1)
            - 0.119 * BDi ** 2.0
            - 0.02668
            + (8.031 * 10 ** -4) * SIi
            + 0.02312 * BDi ** 2.0 / (OCi + 1.0)
            + Cli * 0.001489
            + 0.01908 * BDi ** 2.0
            - 0.001109 * Cli
            - (2.315 * 10 ** -5) * SIi * Cli
            - 0.001197 * SIi * BDi ** 2.0
            - (1.068 * 10 ** -4) * Cli * BDi ** 2.0
        )
        ksat_veri = 240.19 * np.exp(
            19.52348 * thetaSi
            - 8.96847
            - 0.028212 * Cli
            + 1.8107 * 10 ** -4 * SAi ** 2.0
            - 9.4125 * 10 ** -3 * Cli ** 2.0
            - 8.395215 * thetaSi ** 2.0
            + 0.077718 * SAi * thetaSi
            - 0.00298 * SAi ** 2.0 * thetaSi ** 2.0
            - 0.019492 * Cli ** 2 * thetaSi ** 2.0
            + 1.73 * 10 ** -5 * SAi ** 2.0 * Cli
            + 0.02733 * Cli ** 2 * thetaSi
            + 0.001434 * SAi ** 2.0 * thetaSi
            - 3.5 * 10 ** -6 * Cli ** 2.0 * SAi
        )
        lambda_i = np.exp(
            -0.784
            + 0.018 * SAi
            - 1.062 * thetaSi
            - SAi ** 2.0 * 5 * 10 ** -5
            - 0.003 * Cli ** 2.0
            + 1.111 * thetaSi ** 2.0
            - 0.031 * SAi * thetaSi
            + 3.10 ** -4 * SAi ** 2.0 * thetaSi ** 2.0
            - 0.006 * Cli ** 2.0 * thetaSi ** 2.0
            - 2 * 10 ** -6 * SAi ** 2.0 * Cli
            + 0.008 * Cli ** 2.0 * thetaSi
            - 0.007 * thetaSi ** 2.0 * Cli
        )
        ci = 3 + 2 / lambda_i
        if i == 0:
            thetaR = np.copy(thetaRi)
            thetaS = np.copy(thetaSi)
            ksat_ver = np.copy(ksat_veri)
            c = np.copy(ci)
        if i > 0:
            thetaR = np.dstack((thetaR, thetaRi))
            thetaS = np.dstack((thetaS, thetaSi))
            ksat_ver = np.dstack((ksat_ver, ksat_veri))
            c = np.dstack((c, ci))

        i = i + 1

    return thetaS, thetaR, c, ksat_ver


def create_outlet_map(config, rows, cols, cell_size, cell_centers):
    """Creates the outlet map"""
    info("Generate outlet map")

    outlets = []
    for _, coords in config["Outlets"].items():
        outlets.append(loads(coords))

    ## burn in river and stuff like that
    outlet_array = burn_coords(cell_centers, rows, cols, outlets, cell_size)

    outlet_pcr = pcr.numpy2pcr(pcr.Nominal, outlet_array, 10)
    pcr.report(outlet_pcr, config["Outfiles"]["outlet_map"])


def create_catchment_mask(config, rows, cols, cell_centers):
    """Crete the catchment mask"""
    info("Create catchment mask")

    # reads the catchment shapefile
    shape = shapefile.Reader(config["Shapefiles"]["shape_catchment"])
    feature = shape.shapeRecords()[0]
    # contains shape geometry
    first = feature.shape.__geo_interface__
    # creates a numpy array of the mask
    raster = gen_inpoly(first, cell_centers, rows, cols)
    # write raster out
    mask_raster = pcr.numpy2pcr(pcr.Nominal, raster, 10)
    pcr.report(mask_raster, config["Outfiles"]["catchment_mask"])

    return mask_raster


def create_river_burn(config, rows, cols, cell_size, cell_centers):
    """Create river burnin map"""
    info("Burn in river")

    riv_shape = shapefile.Reader(config["Shapefiles"]["rivershape"])
    riv_points = generate_river_points(riv_shape, cell_size)
    riv_array = burn_in_river(cell_centers, rows, cols, riv_points, cell_size)
    riv_corrected = gen_river_connectivity(riv_array, rows, cols)
    ## turn off correction
    riv_corrected = riv_array
    riv_pcr = pcr.numpy2pcr(pcr.Nominal, riv_corrected, 10)
    pcr.report(riv_pcr, config["Outfiles"]["river_burn"])

    return riv_corrected, riv_pcr


def create_ldd_map(config, dem, riv_corrected):
    """Create local drainage direction map"""
    info("Create local drainage direction")

    # removing nans
    riv_where_nan = np.isnan(riv_corrected)
    riv_corrected[riv_where_nan] = 0.0
    riv_pcr_no_nan = pcr.numpy2pcr(pcr.Scalar, riv_corrected, 10)

    # determine regional slope where the river should run
    # ldddem = pcr.ifthen(pcr.boolean(mask_raster), dem)
    ldddem = pcr.ifthenelse(riv_pcr_no_nan >= 1.0, dem - 1000.0, dem)
    ldd = pcr.lddcreate(ldddem, 10.0e35, 10.0e35, 10.0e35, 10.0e35)
    lddrep = pcr.lddrepair(ldd)
    # lddmasked = pcr.ifthen(pcr.boolean(mask_raster), lddrep)
    pcr.report(lddrep, config["Outfiles"]["ldd_map"])

    ##riv_pcr = pcr.ifthen(pcr.scalar(mask_raster) >= 1, riv_pcr)
    ##disttocatch = pcr.spread(pcr.nominal(mask_raster), 0.0, 1.0)
    ##demmax = pcr.ifthenelse(pcr.scalar(mask_raster) >= 1.0,demmax,demmax + (pcr.celllength() * 100.0) / disttocatch,)

    return ldd


def create_streamorder(config, rows, cols, mask_raster, ldd):
    """Create streamorder map"""
    info("Create stream order map")

    # manually adjust maximum streamorder
    stro = pcr.streamorder(ldd)
    stro_scalar = pcr.scalar(stro)
    stro_np = pcr.pcr2numpy(stro_scalar, 0.0)

    ist_max = np.amax(stro_np)
    factor = ist_max / config.getint(
        "Configuration", "max_stream_order", fallback=DEFAULT_MAX_STREAMORDER
    )

    for i in range(0, rows):
        for j in range(0, cols):
            stro_np[i][j] = np.floor(stro_np[i][j] / factor)
            if stro_np[i][j] == 0.0:
                stro_np[i][j] = 1.0

    stro_corr = pcr.numpy2pcr(pcr.Scalar, stro_np, 10)
    stro_masked = pcr.ifthen(pcr.boolean(mask_raster), stro_corr)
    pcr.report(stro_masked, config["Outfiles"]["streamorder_map"])

    return stro_np


def create_river_width(config, rows, cols, riv_pcr, stro_np):
    """
    Compute width on basis of strahler order
    Downing et al (2012): Global abundace and size distribution of streams and rivers.
    """
    info("Create river width")

    width_np = np.copy(stro_np)

    for i in range(0, rows):
        for j in range(0, cols):
            width_np[i][j] = 0.542 * math.exp(0.842 * width_np[i][j])

    width_pcr = pcr.numpy2pcr(pcr.Scalar, width_np, 10)
    riv_masked = pcr.ifthen(pcr.boolean(riv_pcr), width_pcr)
    pcr.report(riv_masked, config["Outfiles"]["river_width_map"])


def create_soil_maps(config, rows, cols):
    """Create soil maps"""
    info("Create unifrom soil map")

    soil_np = np.ones((rows, cols))
    soil_pcr = pcr.numpy2pcr(pcr.Nominal, soil_np, 10)
    pcr.report(soil_pcr, config["Outfiles"]["river_width_map"])

    # print('Create soil thickness map')
    # soil_thick_np = np.ones((rows,cols)) * soil_thickness
    # soil_thick_pcr = pcr.numpy2pcr(pcr.Scalar,soil_thick_np,10)
    # pcr.report(soil_thick_pcr, working_folder + '/' + soil_thickness_map)
    # pcr.report(soil_thick_pcr, working_folder + '/' + min_soil_thickness_map)
    #
    # thetaS, thetaR, c, ksat_ver = read_soil_to_dict(soils_folder)
    #
    # print('Create thetaS')
    # thetaS_pcr = pcr.numpy2pcr(pcr.Scalar,np.copy(thetaS[:,:,0]),10)
    # out_thetaS = working_folder + '/' + thetaS_file
    # pcr.report(thetaS_pcr, out_thetaS)
    # print('Create thetaR')
    # thetaR_pcr = pcr.numpy2pcr(pcr.Scalar,np.copy(thetaR[:,:,0]),10)
    # out_thetaR = working_folder + '/' + thetaR_file
    # pcr.report(thetaR_pcr, out_thetaR)
    #
    # print('ksatver')
    # ksatver_pcr = pcr.numpy2pcr(pcr.Scalar,np.copy(ksat_ver[:,:,0]),10)
    # out_ksat_ver = working_folder + '/' + ksat_ver_file
    # pcr.report(ksatver_pcr, out_ksat_ver)
    #
    # print('Create M')
    # M = np.zeros((rows,cols))
    # for i in range(0,rows):
    #    for j in range(0,cols):
    #        ks_depth = ksat_ver[i,j,:]
    #        y = ks_depth/ksat_ver[i,j,0]
    #        fit = np.polyfit(soil_depth, np.log(y), 1, w=np.sqrt(y))
    #        f = -fit[0]
    #        M[i][j] = (thetaS[i][j][0]-thetaR[i][j][0])/f
    #
    # M_pcr = pcr.numpy2pcr(pcr.Scalar,M,10)
    # out_ksat_ver = working_folder + '/' + M_file
    # pcr.report(M_pcr, out_ksat_ver)
    #
    # print('Create c')
    #
    # for i in range(0,len(take_c)):
    #    c_pcr = pcr.numpy2pcr(pcr.Scalar,np.copy(c[:,:,take_c[i]]),10)
    #    out_c = working_folder + '/c_' + str(i) + '.map'
    #    pcr.report(c_pcr, out_c)


def generate_landuse_lookup(path):
    """Read landuse lookup and create a dictionary for it."""
    lookup_np = np.loadtxt(path, delimiter=",")
    lookup_dict = {}
    for row in lookup_np:
        lookup_dict[int(row[0])] = {
            "N": row[1],
            "Sl": row[2],
            "Swood": row[3],
            "Kext": row[4],
            "RD": row[5],
        }
    return lookup_dict


def create_land_use(config, rows, cols):
    """Creates land use maps"""
    info("Create landuse maps")

    landuse = pcr.readmap(config["Paths"]["landuse_file"])
    lookup = generate_landuse_lookup(config["Paths"]["landuse_lookup"])

    lan_np = pcr.pcr2numpy(landuse, 0.0)

    N = np.zeros((rows, cols))
    Sl = np.zeros((rows, cols))
    Swood = np.zeros((rows, cols))
    Kext = np.zeros((rows, cols))
    RD = np.zeros((rows, cols))

    for i in range(0, rows):
        for j in range(0, cols):
            row = lookup[int(lan_np[i][j])]
            N[i][j] = row["N"]
            Sl[i][j] = row["Sl"]
            Swood[i][j] = row["Swood"]
            Kext[i][j] = row["Kext"]
            RD[i][j] = row["RD"]

    N_pcr = pcr.numpy2pcr(pcr.Scalar, N, 10)
    pcr.report(N_pcr, config["Outfiles"]["N_file"])
    Sl_pcr = pcr.numpy2pcr(pcr.Scalar, Sl, 10)
    pcr.report(Sl_pcr, config["Outfiles"]["Sl_file"])
    Swood_pcr = pcr.numpy2pcr(pcr.Scalar, Swood, 10)
    pcr.report(Swood_pcr, config["Outfiles"]["Swood_file"])
    Kext_pcr = pcr.numpy2pcr(pcr.Scalar, Kext, 10)
    pcr.report(Kext_pcr, config["Outfiles"]["Kext_file"])
    RD_pcr = pcr.numpy2pcr(pcr.Scalar, RD, 10)
    pcr.report(RD_pcr, config["Outfiles"]["rooting_file"])
    pcr.report(landuse, config["Outfiles"]["landuse_map"])


def get_dem_info(dem):
    """Determines raster infos of the dem"""

    # Get values of the clone
    rows = dem.clone().nrRows()
    cols = dem.clone().nrCols()

    cell_size = dem.clone().cellSize()

    # coordinates are in upper left corner
    xmin = dem.clone().west()
    ymin = dem.clone().north()

    return rows, cols, cell_size, xmin, ymin


def split_timedelta64_ns(td):
    """Splits timedelta into days, hours, minutes and seconds portions."""
    seconds = int(td / (10 ** 9))  # [ns] -> [s]
    minutes = int(seconds / 60)
    seconds %= 60
    hours = int(minutes / 60)
    minutes %= 60
    days = int(hours / 24)
    hours %= 60
    return days, hours, minutes, seconds


def create_inmap_temperature(config, rows, cols, cell_centers):
    """Creates temperature inmaps.

    Values are in Kelvin."""
    info("Create temperature inmaps")

    grib = xr.open_dataset(config["Weatherfiles"]["temperature"], engine="cfgrib")

    # create cell centers in input projection
    xscale = grib.coords["longitude"].data
    yscale = grib.coords["latitude"].data

    input_centers = np.zeros((len(yscale), len(xscale), 2))
    for i, ypos in enumerate(yscale):
        for j, xpos in enumerate(xscale):
            input_centers[i][j][0] = xpos
            input_centers[i][j][1] = ypos

    # project cell centers
    transformer = Transformer.from_crs(
        config["Projections"]["in_temperature"], config["Projections"]["out"]
    )
    input_centers[:, :, 0], input_centers[:, :, 1] = transformer.transform(
        input_centers[:, :, 0], input_centers[:, :, 1]
    )

    # reshapes cell centers for interpolator
    input_centers_flat = input_centers.reshape(
        len(input_centers[:, 0]) * len(input_centers[0, :]), 2
    )
    centers_flat = cell_centers.reshape(rows * cols, 2)

    # loop over timesteps
    is_first = True
    is_second = True
    first_step = None
    max_steps = config.getint("Weatherfiles", "max_steps", fallback=0)
    counter = 0
    for step in grib["t2m"]:
        if np.isnan(step).all() and is_first:
            # skip first empty records
            d = np.datetime_as_string(step.time + step.step, unit="s")
            debug(f"skipping: {d}")
            continue
        elif is_first:
            # print start
            d = np.datetime_as_string(step.time + step.step, unit="s")
            info(f"Recording starts at: {d}")
            first_step = step
            is_first = False
        elif not is_first and is_second:
            # print step
            days, hours, minutes, seconds = split_timedelta64_ns(
                (step.time + step.step) - (first_step.time + first_step.step)
            )
            info(f"Step is: {days} d {hours} h {minutes} m {seconds} s")
            is_second = False
        elif max_steps and counter >= max_steps:
            info("max_steps reached")
            return

        # build interpolator
        input_values_flat = step.data.reshape(len(input_centers_flat))
        (step_rows, step_cols) = np.shape(step)
        assert step_rows == len(yscale), "length of rows doesn't match"
        assert step_cols == len(xscale), "length of columns doesn't match"
        interp = NearestNDInterpolator(input_centers_flat, input_values_flat)

        # create map
        interpolated = interp(centers_flat).reshape(rows, cols)
        pcr.report(
            pcr.numpy2pcr(pcr.Scalar, interpolated, -9999),
            f"{config['Paths']['inmaps']}/TEMP{counter/1000:011.3f}",
        )

        counter += 1


def create_inmap_precipitation(config, rows, cols, cell_centers):
    """Creates precipitation inmaps.

    Values are in meters.
    """
    info("Create precipitation inmaps")

    grib_file = config["Weatherfiles"]["precipitation"]
    grib_projection = config["Projections"]["in_precipitation"]
    grib_variable = "tp"
    file_template = config["Paths"]["inmaps"] + "/P{:011.3f}"
    create_inmap_era5_grib_steps(
        config,
        rows,
        cols,
        cell_centers,
        grib_file,
        grib_projection,
        grib_variable,
        file_template,
    )


def create_inmap_evaporation(config, rows, cols, cell_centers):
    """Creates evaporation inmaps.

    Values are in meters."""
    info("Create evaporation inmaps")

    grib_file = config["Weatherfiles"]["evaporation"]
    grib_projection = config["Projections"]["in_evaporation"]
    grib_variable = "pev"
    file_template = config["Paths"]["inmaps"] + "/PET{:011.3f}"
    create_inmap_era5_grib_steps(
        config,
        rows,
        cols,
        cell_centers,
        grib_file,
        grib_projection,
        grib_variable,
        file_template,
    )


def create_inmap_era5_grib_steps(
    config,
    rows,
    cols,
    cell_centers,
    grib_file,
    grib_projection,
    grib_variable,
    file_template,
):
    """Creates mapstacks from era5 grib files with multiple steps."""

    grib = xr.open_dataset(grib_file, engine="cfgrib")

    # create cell centers in input projection
    xscale = grib.coords["longitude"].data
    yscale = grib.coords["latitude"].data

    input_centers = np.zeros((len(yscale), len(xscale), 2))
    for i, ypos in enumerate(yscale):
        for j, xpos in enumerate(xscale):
            input_centers[i][j][0] = xpos
            input_centers[i][j][1] = ypos

    # project cell centers
    transformer = Transformer.from_crs(grib_projection, config["Projections"]["out"])
    input_centers[:, :, 0], input_centers[:, :, 1] = transformer.transform(
        input_centers[:, :, 0], input_centers[:, :, 1]
    )

    # reshapes cell centers for interpolator
    input_centers_flat = input_centers.reshape(
        len(input_centers[:, 0]) * len(input_centers[0, :]), 2
    )
    centers_flat = cell_centers.reshape(rows * cols, 2)

    # loop over timesteps
    is_first = True
    is_second = True
    first_step = None
    max_steps = config.getint("Weatherfiles", "max_steps", fallback=0)
    counter = 0
    for steps in grib[grib_variable]:
        for step in steps:
            if np.isnan(step).all() and is_first:
                # skip first empty records
                d = np.datetime_as_string(step.time + step.step, unit="s")
                debug(f"skipping: {d}")
                continue
            elif is_first:
                # print start
                d = np.datetime_as_string(step.time + step.step, unit="s")
                info(f"Recording starts at: {d}")
                first_step = step
                is_first = False
            elif not is_first and is_second:
                # print step
                days, hours, minutes, seconds = split_timedelta64_ns(
                    (step.time + step.step) - (first_step.time + first_step.step)
                )
                info(f"Step is: {days} d {hours} h {minutes} m {seconds} s")
                is_second = False
            elif max_steps and counter >= max_steps:
                info("max_steps reached")
                return

            # build interpolator
            input_values_flat = step.data.reshape(len(input_centers_flat))
            (step_rows, step_cols) = np.shape(step)
            assert step_rows == len(yscale), "length of rows doesn't match"
            assert step_cols == len(xscale), "length of columns doesn't match"
            interp = NearestNDInterpolator(input_centers_flat, input_values_flat)

            # create map
            interpolated = interp(centers_flat).reshape(rows, cols)
            pcr.report(
                pcr.numpy2pcr(pcr.Scalar, interpolated, -9999),
                file_template.format(counter / 1000),
            )

            counter += 1


def main():
    """Main function to prepare the files"""

    parser = ArgumentParser(description="Prepare wflow files")
    parser.add_argument("config_file", help="configuration file destination")
    args = parser.parse_args()

    config = ConfigParser(interpolation=ExtendedInterpolation())
    config.read(args.config_file)

    logging.basicConfig(
        level=LOG_LEVEL_MAP.get(
            config.get("Configuration", "log_level", fallback="INFO").lower(),
            logging.INFO,
        ),
        format="%(levelname)s %(asctime)s: %(message)s",
    )

    # Soil stuff
    # soils_folder = "/home/iwbworkstation/Desktop/working_dir/50m_data/2_Soil"
    # dictionary = {
    #     "0-5": {},
    #     "5-15": {},
    #     "15-30": {},
    #     "30-60": {},
    #     "60-100": {},
    #     "100-200": {},
    # }
    # soil_depth = [25, 100, 225, 450, 800, 1500]
    # soil_thickness = 2000.0
    # takes c from which layers
    # take_c = [1, 2, 3, 4]

    pcr.setglobaloption("unitcell")
    pcr.setclone(config["Paths"]["masterdem"])

    rows, cols, cell_size, xmin, ymin = get_dem_info(pcr)
    debug(f"rows: {rows} cols: {cols}")
    debug(f"cell_size: {cell_size}")
    debug(f"xmin: {xmin} ymin: {ymin}")

    cell_centers = init_cellcenter(rows, cols, cell_size, xmin, ymin)

    # resolve dependencies
    need_outlet_map = config.getboolean("Jobs", "outlet_map", fallback=False)
    need_land_use = config.getboolean("Jobs", "land_use_map", fallback=False)
    need_soil_map = config.getboolean("Jobs", "soil_map", fallback=False)
    need_river_width = config.getboolean("Jobs", "river_width", fallback=False)
    need_stream_order = (
        config.getboolean("Jobs", "stream_order", fallback=False) or need_river_width
    )
    need_ldd_map = (
        config.getboolean("Jobs", "ldd_map", fallback=False) or need_stream_order
    )
    need_river_burn = (
        config.getboolean("Jobs", "river_burn", fallback=False)
        or need_ldd_map
        or need_river_width
    )
    need_catchment_mask = (
        config.getboolean("Jobs", "catchment_mask", fallback=False) or need_stream_order
    )
    need_inmap_precipitation = config.getboolean(
        "Jobs", "inmap_precipitation", fallback=False
    )
    need_inmap_temperature = config.getboolean(
        "Jobs", "inmap_temperature", fallback=False
    )
    need_inmap_evaporation = config.getboolean(
        "Jobs", "inmap_evaporation", fallback=False
    )

    # execute tasks
    if need_catchment_mask:
        mask_raster = create_catchment_mask(config, rows, cols, cell_centers)

    if need_outlet_map:
        create_outlet_map(config, rows, cols, cell_size, cell_centers)

    if need_river_burn:
        riv_corrected, riv_pcr = create_river_burn(
            config, rows, cols, cell_size, cell_centers
        )

    if need_ldd_map:
        dem = pcr.readmap(config["Paths"]["masterdem"])
        ldd = create_ldd_map(config, dem, riv_corrected)

    if need_stream_order:
        stro_np = create_streamorder(config, rows, cols, mask_raster, ldd)

    if need_river_width:
        create_river_width(config, rows, cols, riv_pcr, stro_np)

    if need_soil_map:
        create_soil_maps(config, rows, cols)

    if need_land_use:
        create_land_use(config, rows, cols)

    if need_inmap_precipitation:
        create_inmap_precipitation(config, rows, cols, cell_centers)

    if need_inmap_temperature:
        create_inmap_temperature(config, rows, cols, cell_centers)

    if need_inmap_evaporation:
        create_inmap_evaporation(config, rows, cols, cell_centers)

    debug("Tasks complete")


if __name__ == "__main__":
    main()
