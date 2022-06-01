#!/bin/env python3
"""


"""


#######################################################################
#                           Imports
#######################################################################

import netCDF4 as nc
import pcraster as pcr
import numpy as np
import datetime
import os
import time
import math
import pyproj
from scipy.interpolate import NearestNDInterpolator
from scipy import interpolate
from tqdm import tqdm


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

    for i in range(0, nrows):
        for j in range(0, ncols):
            cell_centers[i][j][0] = x_center + cellsize * j
            cell_centers[i][j][1] = (y_center + cellsize * (nrows - 1)) - cellsize * i

    return cell_centers


def init_cellcenter(rows, cols, cell_size, xmin, ymin):

    cell_centers = np.zeros((rows, cols, 2))

    for i in range(0, rows):
        for j in range(0, cols):
            cell_centers[i][j][0] = xmin + cell_size / 2.0 + cell_size * j
            cell_centers[i][j][1] = ymin - cell_size / 2.0 - cell_size * i

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


#######################################################################
#                        time things
#######################################################################

start_time = time.time()

#######################################################################
#                        pcraster blueprint
#######################################################################

working_folder = r"/home/iwbworkstation/Desktop/working_dir/model_rerun_paper1/Model_200m_HBV/staticmaps"
masterdem = "wflow_dem.map"

pcr.setglobaloption("unitcell")

pcr.setclone(working_folder + "/" + masterdem)

rows = pcr.clone().nrRows()
cols = pcr.clone().nrCols()

cell_size = pcr.clone().cellSize()

# coordinates are in upper left corner
xmin = pcr.clone().west()
ymin = pcr.clone().north()
xmax = xmin + cell_size * cols
ymax = ymin - cell_size * rows

cell_centers = init_cellcenter(rows, cols, cell_size, xmin, ymin)

print("Original cell centers initialized in " + str(time.time() - start_time) + " s")

#######################################################################
#                        NetCDF settings
#######################################################################

rain = r"/media/iwbworkstation/My Passport/Dissertation/3_Hydrology/0_Model_Data/4_Precepitation/Forcing040506/RR/RR"
temp = r"/media/iwbworkstation/My Passport/Dissertation/3_Hydrology/0_Model_Data/4_Precepitation/Forcing040506/TT/TT"
rad = r"/media/iwbworkstation/Volume/Dissertation/3_Hydrology/0_Model_Data/4_Precepitation/Strahlung_Muerz"
resultfolder = r"/home/iwbworkstation/Desktop/working_dir/model_rerun_paper1/NSGA2/Model_200m_HBV/inmaps_040506"

INCA_epsg = "epsg:31258"
model_epsg = "epsg:32633"

#######################################################################

# Initialize Projections
inProj = pyproj.Proj(init=INCA_epsg)
outProj = pyproj.Proj(init=model_epsg)

# Count all the files in the given folder
path, dirs, files = next(os.walk(rain))
count = len(files)

# allocate rain array
rain_input = np.zeros((count, rows, cols), dtype="float32")

current_prog = 0

pbar1 = tqdm(total=count)

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
            # read INCA file in ESRII format, header only
            metadata = np.genfromtxt(path, max_rows=6, dtype=str)
            # compute cell centers in original projection
            centers = init_cellcenters_esrii(metadata)
            centers[:, :, 0], centers[:, :, 1] = pyproj.transform(
                inProj, outProj, centers[:, :, 0], centers[:, :, 1]
            )
            # reshapes cell centers matrix for interpolator
            centers_flatt = centers.reshape(len(centers[:, 0]) * len(centers[0, :]), 2)
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
        pbar1.update(1)

        i += 1


print("Precipitation array generated in " + str(time.time() - start_time) + " s")
#####################################
#     Temperature
#####################################

# allocate temperature matrix
temp_input = np.zeros((count, rows, cols), dtype="float32")

pbar2 = tqdm(total=count)

## Temperature loop also calculate PET
i = 0
for subdir, dirs, files in os.walk(temp):
    # loops of files in folder AND SORTES them alphabetically
    # metadata is used from precipitation
    for file in sorted(files):
        path = temp + "/" + file
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

        # handle PET
        # datestring
        datepart = file.split("_")[2]
        # FIXME: magic string
        path_rad = rad + "/" + "Muerz_GL_" + datepart + "00_00.asc"

        rad_vals = np.genfromtxt(path_rad, skip_header=6)

        # remove potential nan values
        # check if nan in values list
        pot_nans = np.isnan(np.sum(rad_vals))
        # if nans are encountered
        if pot_nans:
            # create mask of nan values
            mask = np.isnan(rad_vals)
            # interpolate values
            rad_vals[mask] = np.interp(
                np.flatnonzero(mask), np.flatnonzero(~mask), rad_vals[~mask]
            )

        rad_flatt = rad_vals.reshape(len(centers[:, 0]) * len(centers[0, :]))
        myInterpolator3 = NearestNDInterpolator(centers_flatt, rad_flatt)
        # interpolate and reshape
        rad_interp = myInterpolator3(orig_centers)
        rad_inter2D = rad_interp.reshape(rows, cols)

        # compute PET
        phi = 0.646 + 0.0006 * temp_inter2D
        lamb = 1000.0 * (2501 - 2.38 * temp_inter2D)
        delta = (6.107 * 7.5 * 273.3 / (273.3 + temp_inter2D) ** 2.0) * np.exp(
            7.5 * temp_inter2D / (273.3 + temp_inter2D)
        )
        rho = 1000.0

        # limit pet to 0
        pet = np.maximum(
            1000 * 0.75 * (delta / (delta + phi)) * (900 * rad_inter2D / lamb / rho),
            0.0,
        )
        pet_pcr = pcr.numpy2pcr(pcr.Scalar, pet, -9999)
        pcrpath_pet = resultfolder + "/PET" + "{:09.3f}".format((i) / 1000.0)
        pcr.report(pet_pcr, pcrpath_pet)

        pbar2.update(1)

        i += 1

# compute parameters required for PET calculation

# inint PET array copy from temperature


#


# while currentstep < temp_input.shape[0]:
#
#
#
#
#    # compute lower index
#    lower = int(currentstep)
#    # compute higher index
#    higher = int(currentstep + steps_)
#    # get current month
#    curr_month = julian_month[lower]
#    # get current S from lookup table
#    S0_curr = S0[int(curr_month)-1]
#    # compute tmax, tmin and tmean over averaging period
#    tmax = np.amax(temp_input[lower:higher,:,:],axis = 0)
#    tmin = np.amin(temp_input[lower:higher,:,:],axis = 0)
#    tmean = np.mean(temp_input[lower:higher,:,:],axis = 0)
#
#    # compute evatranspirtation for each field
#    ET = 0.0023*S0_curr*np.sqrt(tmax-tmin)*(tmean+17.8)
#    # now divide equally on averaging period
#    ET = ET / float(steps_)
#    # now fill the array
#    PET_out[lower:higher,:,:] = ET
#
#    currentstep += steps_
#
# for i in range(0,len(PET_out)):
#    PET_temp = PET_out[i,:,:]
#    pet_pcr = pcr.numpy2pcr(pcr.Scalar,PET_temp,-9999)
#    pcrpath = resultfolder + '/PET' + '{:09.3f}'.format((i)/1000.)
#    pcr.report(pet_pcr,pcrpath)
