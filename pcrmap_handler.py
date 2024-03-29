# python!
"""


"""


#######################################################################
#                           Imports
#######################################################################

import netCDF4 as nc
import pcraster as pcr
import numpy as np
import os
import math
import pyproj
from scipy.interpolate import NearestNDInterpolator

from argparse import ArgumentParser
from configparser import ConfigParser, ExtendedInterpolation
from logging import info, debug
from preparewflow import init_cellcenter, get_dem_info, LOG_LEVEL_MAP
import logging


#######################################################################
#                        Functions
#######################################################################


def init_cellcenters_esrii(metadata):

    ncols = int(metadata[0][1])
    nrows = int(metadata[1][1])
    x_center = float(metadata[2][1])
    y_center = float(metadata[3][1])
    cellsize = float(metadata[4][1])

    cell_centers = np.zeros((nrows, ncols, 2))
    debug(
        f"ncols: {ncols} nrows: {nrows} x_center: {x_center} "
        + f"y_center: {y_center} cellsize: {cellsize}"
    )

    for i in range(0, nrows):
        for j in range(0, ncols):
            cell_centers[i][j][0] = x_center + cellsize * j
            cell_centers[i][j][1] = (y_center + cellsize * (nrows - 1)) - cellsize * i

    return cell_centers


def convert_lat_to_rad(lat):

    x1 = float(lat[0:2])
    x2 = float(lat[3:5])
    x3 = lat[-1]
    if x3 == "N":
        dec = x1 + x2 / 60.0
    if x3 == "S":
        dec = -x1 - x2 / 60.0
    radians = math.pi / 180.0 * dec

    return radians


def main():
    """Prepare mapstacks"""

    parser = ArgumentParser(description="Prepare wflow mapstacks")
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

    #######################################################################
    #                        pcraster blueprint
    #######################################################################

    pcr.setglobaloption("unitcell")
    pcr.setclone(config["Paths"]["masterdem"])

    rows, cols, cell_size, xmin, ymin = get_dem_info(pcr)
    cell_centers = init_cellcenter(rows, cols, cell_size, xmin, ymin)

    info("Original cell centers initialized")

    #######################################################################
    #                        NetCDF settings
    #######################################################################

    # General settings
    rain = r"/media/iwbworkstation/Volume/Dissertation/3_Hydrology/0_Model_Data/4_Precepitation/ForcingOct04-05-06/RR"
    temp = r"/media/iwbworkstation/Volume/Dissertation/3_Hydrology/0_Model_Data/4_Precepitation/ForcingOct04-05-06/TT"
    resultfolder = (
        r"/home/iwbworkstation/Desktop/working_dir/NSGA2_hbv/Model_200m_HBV04-06/inmaps"
    )

    # minutes is accurate enough
    timestep = 15  # in min
    # averaging temperature
    tempaveraging = 1440  # in minutes
    # printouts for progress 'day' or 'month'

    #######################################################################

    # Initialize Projections
    inProj = pyproj.Proj(init=config["projections"]["in_precipitation"])
    outProj = pyproj.Proj(init=config["projections"]["out"])

    # Count all the files in the given folder
    path, dirs, files = next(os.walk(rain))
    count = len(files)

    # allocate rain array
    rain_input = np.zeros((count, rows, cols), dtype="float32")

    # Rainloop
    i = 0
    for subdir, dirs, files in os.walk(rain):
        # loops of files in folder AND SORTES them alphabetically
        for file in sorted(files):
            # read some metadata from the first file and compute cell centers
            # this is only done once since they don't change
            # create file path
            path = rain + "/" + file
            if i == 0:
                # start time metadata
                startyear = file[-17:-13]
                startmonth = file[-13:-11]
                startday = file[-11:-9]
                starthour = file[-9:-7]
                startminute = file[-6:-4]
                # read INCA file in ESRII format, header only
                metadata = np.genfromtxt(path, max_rows=6, dtype=str)
                # compute cell centers in original projection
                centers = init_cellcenters_esrii(metadata)
                centers[:, :, 0], centers[:, :, 1] = pyproj.transform(
                    inProj, outProj, centers[:, :, 0], centers[:, :, 1]
                )
                # reshapes cell centers matrix for interpolator
                centers_flatt = centers.reshape(
                    len(centers[:, 0]) * len(centers[0, :]), 2
                )
                # flatten cell centers from master dem for interpolator
                orig_centers = cell_centers.reshape(rows * cols, 2)

            # read INCA file in ESRII format, precipitation data
            values = np.genfromtxt(path, skip_header=6)
            # reshape values array for interpolator
            values_flatt = values.reshape(len(centers[:, 0]) * len(centers[0, :]))
            # create the nearest neighbor interpolation function
            myInterpolator1 = NearestNDInterpolator(centers_flatt, values_flatt)
            # interpolated precipitation values
            rain_interp = myInterpolator1(orig_centers)
            # reshape them back to be writebale to the netcdf file
            rain_inter2D = rain_interp.reshape(rows, cols)
            # convert to pcr
            rain_inter2D_pcr = pcr.numpy2pcr(pcr.Scalar, rain_inter2D, -9999)
            pcrpath = resultfolder + "/P" + "{:011.3f}".format((i) / 1000.0)
            pcr.report(rain_inter2D_pcr, pcrpath)
            # write the each timestep to global precipitation matrix
            rain_input[i, :, :] = rain_inter2D[:, :]
            # count i up
            i += 1

    info("Precipitation array generated")

    # allocate temperature matrix
    temp_input = np.zeros((count, rows, cols), dtype="float32")
    julian_month = []

    ## Temperature loop also calculate PET
    i = 0
    for subdir, dirs, files in os.walk(temp):
        # loops of files in folder AND SORTES them alphabetically
        # metadata is used from precipitation
        for file in sorted(files):
            path = temp + "/" + file
            j_month = file[-13:-11].lstrip("0")
            julian_month.append(j_month)
            # create and flatten temperature values
            values = np.genfromtxt(path, skip_header=6)
            values_flatt = values.reshape(len(centers[:, 0]) * len(centers[0, :]))
            # create new temperature interpolator
            myInterpolator2 = NearestNDInterpolator(centers_flatt, values_flatt)
            # interpolate and reshape temperature values
            temp_interp = myInterpolator2(orig_centers)
            temp_inter2D = temp_interp.reshape(rows, cols)
            temp_inter2D_pcr = pcr.numpy2pcr(pcr.Scalar, temp_inter2D, -9999)
            pcrpath = resultfolder + "/TEMP" + "{:08.3f}".format((i) / 1000.0)
            pcr.report(temp_inter2D_pcr, pcrpath)
            # write to global temperature matrix
            temp_input[i, :, :] = temp_inter2D[:, :]

            i += 1

    info("Temperature array generated")

    # how many files or arrays is in the averaging period?
    steps_ = tempaveraging / timestep

    # Extraterrestrial radiation
    S0 = [10.7, 16.3, 23.8, 32.6, 38.9, 41.85, 40.5, 35.15, 27.05, 18.65, 12.0, 9.3]

    # inint PET array copy from temperature
    PET_out = np.copy(temp_input)

    currentstep = 0
    # loop while is smaller, the rest will be taken care of automatically
    while currentstep < temp_input.shape[0]:
        # compute lower index
        lower = int(currentstep)
        # compute higher index
        higher = int(currentstep + steps_)
        # get current month
        curr_month = julian_month[lower]
        # get current S from lookup table
        S0_curr = S0[int(curr_month) - 1]
        # compute tmax, tmin and tmean over averaging period
        tmax = np.amax(temp_input[lower:higher, :, :], axis=0)
        tmin = np.amin(temp_input[lower:higher, :, :], axis=0)
        tmean = np.mean(temp_input[lower:higher, :, :], axis=0)

        # compute evatranspirtation for each field
        ET = 0.0023 * S0_curr * np.sqrt(tmax - tmin) * (tmean + 17.8)
        # now divide equally on averaging period
        ET = ET / float(steps_)
        # now fill the array
        PET_out[lower:higher, :, :] = ET

        currentstep += steps_

    for i in range(0, len(PET_out)):
        PET_temp = PET_out[i, :, :]
        pet_pcr = pcr.numpy2pcr(pcr.Scalar, PET_temp, -9999)
        pcrpath = resultfolder + "/PET" + "{:09.3f}".format((i) / 1000.0)
        pcr.report(pet_pcr, pcrpath)


if __name__ == "__main__":
    main()
