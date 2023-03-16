#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
++++++++++++++++++++++++++++ wflowModelC +++++++++++++++++++++++++++++++

    @brief
        Converts gdal readable formats to netcdf
    
    @note
        Currently only works from high time steps to lower ones
        
    @libraries

    @history 30/05/2022 -- Sebastian Gegenleithner
        First version
             02/06/2022 -- Sebastian Gegenleithner
        Improved buffering and overlap
             10/06/2022 -- Sebastian Gegenleithner
        Improved nan handling
             20/06/2022 -- Sebastian Gegenleithner
        buffer set to whole days, fixed a small bug
        
"""

############################################
#               Imports
############################################
import datetime
from osgeo import gdal
import numpy as np
import netCDF4 as nc
import time
import os
import cftime
import osgeo
import osgeo.ogr
import pandas as pd
import sys
from tqdm import trange
from tqdm import tqdm
import shutil

from scipy import interpolate

# from tqdm.auto import tqdm

############################################
#               Functions
############################################


def set_metadata():
    # initiate metadata entries
    metadata = {}
    metadata["title"] = "wflow input mapstack"
    metadata["institution"] = "TU Graz - (c) S. Gegenleithner"
    metadata["source"] = "NETCDFhandler.py"
    metadata["history"] = time.ctime()
    metadata["references"] = "wflow/pcr2netcdf.py - (c) J. Schellekens"
    metadata["Conventions"] = "CF-1.4"

    return metadata


############################################
#               Classes
############################################

# Time series to netcdf
# converts a time series to netcdf format
class Series2nc:
    # init
    def __init__(
        self,
        out_file,
        start_date,
        end_date,
        pathRR,
        pathTT,
        outProj,
        bounds,
        resolution,
        **kwargs
    ):

        # --------------- retrieve optional arguments ------------------
        # path to radiation, default is None
        self.pathGL = kwargs["pathGL"] if "pathGL" in kwargs else None
        # path to potential evapotranspiration, default is None
        self.pathPET = kwargs["pathPET"] if "pathPET" in kwargs else None
        # time step of the output data, if None we automatically take the
        # smallest time step of the input
        self.time_step = kwargs["timestep"] if "timestep" in kwargs else None
        # defined the allowed memory allocation, default is 10 GB
        self.allowed_ram = kwargs["allowed_ram"] if "allowed_ram" in kwargs else 3.0
        # --------------------------------------------------------------

        # --------------- Set required variables ------------------
        # start time as datetime object
        self.start_date = start_date
        # end tiem as datetime object
        self.end_date = end_date
        # output netcdf file path and name
        self.out_file = out_file
        # path to precipitation folder
        self.pathRR = pathRR
        # path to temperature folder
        self.pathTT = pathTT
        # extents of the model, bounding box
        self.boundingbox = bounds
        # resolution of the model
        self.resolution = resolution
        # projection of the output netcdf file
        self.outproj = outProj
        # --------------------------------------------------------------

        # --------------- Dervived variables ------------------
        # init metadata of the netcdf file
        self.metadata = set_metadata()
        # determine numer of rows and cols
        self.rows = int(
            (self.boundingbox[3] - self.boundingbox[1]) / self.resolution
        )
        self.cols = int(
            (self.boundingbox[2] - self.boundingbox[0]) / self.resolution
        )
        # set options for PET computation
        # Default PET option is blainy-griddle (if only) temperature and
        # precipitation is given
        # 0 = blainey criddle (Temperature and precipitation values)
        # 1 = PET is given
        # 2 = Makkink (If TT, RR and GL (downwards radiation) is given)
        self.pet_option = 0
        if self.pathPET:
            print(
                "[%s] -> PET was given, directly used"
                % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
            )
            self.pet_option = 1
            # automaticall determine timestep from files
            if not self.time_step:
                self.time_step = min(
                    float(self.pathTT["timestep_secs"]),
                    float(self.pathRR["timestep_secs"]),
                    float(self.pathPET["timestep_secs"]),
                )

            # determine maximum timestamp
            self.max_time_step = max(
                float(self.pathTT["timestep_secs"]),
                float(self.pathRR["timestep_secs"]),
                float(self.pathPET["timestep_secs"]),
            )

            # allocated large variables
            alloc_vars = 3
        # we use makkink if no pet was given but radiation was given
        elif self.pathGL and not self.pathPET:
            print(
                "[%s] -> GL detected, we use MAKKINK to compute PET"
                % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
            )
            self.pet_option = 2

            # automaticall determine timestep from files
            if not self.time_step:
                self.time_step = min(
                    float(self.pathTT["timestep_secs"]),
                    float(self.pathRR["timestep_secs"]),
                    float(self.pathGL["timestep_secs"]),
                )
            # allocated large variables

            # determine maximum timestamp
            self.max_time_step = max(
                float(self.pathTT["timestep_secs"]),
                float(self.pathRR["timestep_secs"]),
                float(self.pathGL["timestep_secs"]),
            )

            alloc_vars = 4
        else:
            print(
                "[%s] -> TT and RR were found, we use Blainey-Criddle tp compute PET"
                % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
            )
            self.pet_option = 0
            # automaticall determine timestep from files
            if not self.time_step:
                self.time_step = min(
                    float(self.pathTT["timestep_secs"]),
                    float(self.pathRR["timestep_secs"]),
                )

            # determine maximum timestamp
            self.max_time_step = max(
                float(self.pathTT["timestep_secs"]), float(self.pathRR["timestep_secs"])
            )

            alloc_vars = 3

        # estimate buffer size from large variables and size
        # compute the number of bytes per time step
        byte_per_timestep = np.random.rand(self.rows, self.cols).nbytes
        # compute number of bytes per variable
        byte_per_variable = byte_per_timestep * len(self._get_timerange())
        # compute size of all variables
        total_bytes = byte_per_variable * alloc_vars
        # translate to GB
        total_gb = total_bytes / 1000.0 / 1000.0 / 1000.0
        # set buffersize
        self.buffsize = np.ceil(total_gb / self.allowed_ram)
        # self.buffsize= 3

    def _interpolate_spatial(self, data):

        array = np.ma.masked_invalid(data)

        x = np.arange(0, array.shape[1])
        y = np.arange(0, array.shape[0])

        xx, yy = np.meshgrid(x, y)
        # get only the valid values
        x1 = xx[~array.mask]
        y1 = yy[~array.mask]
        newarr = array[~array.mask]

        array_inter = interpolate.griddata(
            (x1, y1), newarr.ravel(), (xx, yy), method="nearest"
        )

        return array_inter

    # compute time range between start and end
    def _get_timerange(self):
        # derive range
        r = int(
            (
                self.end_date
                + datetime.timedelta(seconds=self.time_step)
                - self.start_date
            ).total_seconds()
            / self.time_step
        )

        return [
            self.start_date + datetime.timedelta(seconds=(self.time_step * i))
            for i in range(r)
        ]

    # read a single map, and warp using gdal
    def _read_map(self, filepath, inproj):
        # read map to memory

        ds = gdal.Warp(
            "",
            filepath,
            outputBounds=self.boundingbox,
            dstSRS=self.outproj,
            srcSRS=inproj,
            srcNodata=None,
            dstNodata=None,
            xRes=self.resolution,
            yRes=self.resolution,
            format="MEM",
        )
        # get geotransform object
        geotrans = ds.GetGeoTransform()
        # get origin x and y
        originX = geotrans[0]
        originY = geotrans[3]
        resX = geotrans[1]
        resY = geotrans[5]
        cols = ds.RasterXSize
        rows = ds.RasterYSize

        x = np.linspace(
            originX + resX / 2, originX + resX / 2 + resX * (cols - 1), cols
        )
        y = np.linspace(
            originY + resY / 2, originY + resY / 2 + resY * (rows - 1), rows
        )
        # no data value
        RasterBand = ds.GetRasterBand(1)
        FillVal = RasterBand.GetNoDataValue()
        # read data
        data = RasterBand.ReadAsArray(0, 0, cols, rows)
        # check if the data contains nan values
        if np.isnan(np.sum(data)):
            # close data gaps by interpolation
            print(
                "[%s] -> some nan values found in %s -> interpolating"
                % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"), filepath)
            )
            data = self._interpolate_spatial(data)
            print(data)

        # clear allocation
        ds = None

        return x, y, data.copy(), FillVal

    # prepare an empty netcdf file with dimensions and metadata
    def _prepare_nc(self):
        # prepares netcdf file, format, compression, etc. given
        Format = "NETCDF4"
        complevel = 9
        zlib = True
        calendar = "gregorian"
        # writes the netcdf file
        nc_out = nc.Dataset(
            self.out_file, "w", format=Format, zlib=zlib, complevel=complevel
        )
        # create time dimesnion
        nc_out.createDimension("time", 0)
        # get timerange
        # read first file in directory to init x and y coordinates
        # loop over directory
        for file in os.listdir(self.pathRR["root"]):
            filepath = self.pathRR["root"] + "/" + file
            # catches non maching files
            try:
                # this is used to catch errors in time format
                time_catcher = datetime.datetime.strptime(file, self.pathRR["naming"])
                try:
                    # tries to read the map, throws and error if not possible
                    x, y, data, FillVal = self._read_map(
                        filepath, self.pathRR["projection"]
                    )
                except:
                    print("Something went wrong whilst reading the file")
                # break at first valid file
                break
            except:
                print(
                    "time format %s does not match with %s"
                    % (self.pathRR["naming"], file)
                )
        # compute timerange
        timerange = self._get_timerange()
        # set units of time axis
        epoch = timerange[0]
        units = "seconds since %04d-%02d-%02d %02d:%02d:%02d.0 00:00" % (
            epoch.year,
            epoch.month,
            epoch.day,
            epoch.hour,
            epoch.minute,
            epoch.second,
        )
        # find start and endday
        startDayNr = cftime.date2num(
            timerange[0].replace(tzinfo=None), units=units, calendar=calendar
        )
        endDayNr = cftime.date2num(
            timerange[-1].replace(tzinfo=None), units=units, calendar=calendar
        )
        timeAR = np.linspace(startDayNr, endDayNr, num=len(timerange))
        # add time to netcdf
        DateHour = nc_out.createVariable(
            "time", "f8", ("time",), fill_value=FillVal, zlib=zlib, complevel=complevel
        )
        # Settings for time axis
        DateHour.units = units
        DateHour.calendar = calendar
        DateHour.standard_name = "time"
        DateHour.long_name = "time"
        DateHour.axis = "T"
        DateHour[:] = timeAR
        # handle projection
        srs = osgeo.osr.SpatialReference()
        res = srs.ImportFromEPSG(int(self.outproj[5:]))
        projStr = srs.ExportToProj4()
        # create x and y variables
        nc_out.createDimension("y", len(y))
        nc_out.createDimension("x", len(x))
        y_var = nc_out.createVariable(
            "y", "f4", ("y",), fill_value=FillVal, zlib=zlib, complevel=complevel
        )
        # Settings for spatial variables
        y_var.standard_name = "projection_y_coordinate"
        y_var.long_name = "y-coordinate in Cartesian system"
        y_var.units = "m"
        y_var.axis = "Y"
        x_var = nc_out.createVariable(
            "x", "f4", ("x",), fill_value=FillVal, zlib=zlib, complevel=complevel
        )
        x_var.standard_name = "projection_x_coordinate"
        x_var.long_name = "x-coordinate in Cartesian system"
        x_var.units = "m"
        x_var.axis = "X"
        y_var[:] = y
        x_var[:] = x
        # projection settings
        crs = nc_out.createVariable("crs", "c")
        crs.long_name = self.outproj
        crs.grid_mapping_name = "universal_transverse_mercator"
        crs.utm_zone_number = srs.GetUTMZone()
        crs.semi_major_axis = srs.GetSemiMajor()
        crs.inverse_flattening = srs.GetInvFlattening()
        crs._CoordinateTransformType = "Projection"
        crs._CoordinateAxisTypes = "y x"
        crs.proj4_params = projStr
        # update netcdf
        for attr in self.metadata:
            nc_out.setncattr(attr, self.metadata[attr])
        # sync and close netcdf file
        nc_out.sync()
        nc_out.close()

    # adds a variable to netcdf structure
    def _add_variable(self, name, format, st_name, unit):
        # Format etc. given as above
        Format = "NETCDF4"
        complevel = 9
        zlib = True
        calendar = "gregorian"
        # open netcdf file in append mode
        nc_out = nc.Dataset(
            self.out_file, "a", format=Format, zlib=zlib, complevel=complevel
        )
        # add variable
        nc_var = nc_out.createVariable(
            name,
            format,
            ("time", "y", "x"),
            fill_value=-9999.0,
            zlib=zlib,
            complevel=9,
            least_significant_digit=None,
        )
        # Settings for the added dvariable
        nc_var.units = unit
        nc_var.standard_name = st_name
        # sync and close netcdf
        nc_out.sync()
        nc_out.close()

    # Add data to existing structure with given buffer range
    def _add_data(self, var, data, buff_range):
        # Format etc. given as above
        Format = "NETCDF4"
        complevel = 9
        zlib = True
        calendar = "gregorian"
        # open netcdf file in append mode
        nc_out = nc.Dataset(
            self.out_file, "a", format=Format, zlib=zlib, complevel=complevel
        )
        # give to which variable we are writing
        nc_var = nc_out.variables[var]
        # add data at buffer location
        nc_var[buff_range[0] : buff_range[1] + 1, :, :] = data
        # sync and close
        nc_out.sync()
        nc_out.close()

    # determines range of the buffer
    def _det_buffer_range(self, timerange):
        # allocate buffer list
        buff_list = []

        # enforce first entry to be the beginning of the first buffer
        buff_list.append(0)

        # get indices of buffer seperators
        for i in range(
            int(np.ceil(len(timerange) / self.buffsize)),
            len(timerange),
            int(np.ceil(len(timerange) / self.buffsize)),
        ):
            # enforce buffer to nearest candidate
            j = 0
            buffer_boolean = True
            while buffer_boolean == True:
                # enforce buffer to full days, 86400. seconds in a day
                if (
                    (timerange[i + j] - datetime.datetime(1970, 1, 1)).total_seconds()
                    / 86400.0
                ).is_integer():
                    buffer_boolean = False
                else:
                    j += 1

            buff_list.append(i + j)

        # add to list
        buff_list.append(len(timerange) - 1)
        # splits buffer in index pairs
        buffers = [
            [buff_list[i], buff_list[i + 1]] for i in range(0, len(buff_list) - 1)
        ]

        return buffers

    # Create a list of datetimes and corresponding filenames
    def _index_P(self, timerange):
        # create indexed files precipitation
        p = []
        # loop over precipitation
        for file in sorted(os.listdir(self.pathRR["root"])):
            filepath = self.pathRR["root"] + "/" + file
            # check if files match naming convention, else the file is skipped
            try:
                date_obj = datetime.datetime.strptime(file, self.pathRR["naming"])
                # get timesteps
                if (timerange[0] - date_obj).total_seconds() <= 0.0 and (
                    timerange[-1] - date_obj
                ).total_seconds() >= 0.0:
                    p.append([date_obj, filepath])
            except:
                print(
                    "time format %s does not match with %s"
                    % (self.pathRR["naming"], file)
                )
                continue
        return p

    # Create a list of datetimes and corresponding filenames
    def _index_TEMP(self, timerange):
        # create indexed files temperature
        temp = []
        # loop over temperature
        for file in sorted(os.listdir(self.pathTT["root"])):
            filepath = self.pathTT["root"] + "/" + file
            # check if files match naming convention, else the file is skipped
            try:
                date_obj = datetime.datetime.strptime(file, self.pathTT["naming"])
                # get timesteps
                if (timerange[0] - date_obj).total_seconds() <= 0.0 and (
                    timerange[-1] - date_obj
                ).total_seconds() >= 0.0:
                    temp.append([date_obj, filepath])
            except:
                print(
                    "WARNING: time format %s does not match with %s"
                    % (self.pathTT["naming"], file)
                )
                continue

        return temp

    # Create a list of datetimes and corresponding filenames
    def _index_GL(self, timerange):
        # create indexed files radiation
        gl = []
        # loop over radiation
        for file in sorted(os.listdir(self.pathGL["root"])):
            filepath = self.pathGL["root"] + "/" + file
            # check if files match naming convention, else the file is skipped
            try:
                date_obj = datetime.datetime.strptime(file, self.pathGL["naming"])
                # get timesteps
                if (timerange[0] - date_obj).total_seconds() <= 0.0 and (
                    timerange[-1] - date_obj
                ).total_seconds() >= 0.0:
                    gl.append([date_obj, filepath])
            except:
                print(
                    "WARNING: time format %s does not match with %s"
                    % (self.pathGL["naming"], file)
                )
                continue

        return gl

    # Create a list of datetimes and corresponding filenames
    def _index_PET(self, timerange):
        # create indexed files potential evapotranspiration
        pet = []
        # loop over potential evapotranspiration
        for file in sorted(os.listdir(self.pathPET["root"])):
            filepath = self.pathPET["root"] + "/" + file
            # check if files match naming convention, else the file is skipped
            try:
                date_obj = datetime.datetime.strptime(file, self.pathPET["naming"])
                # get timesteps
                if (timerange[0] - date_obj).total_seconds() <= 0.0 and (
                    timerange[-1] - date_obj
                ).total_seconds() >= 0.0:
                    pet.append([date_obj, filepath])
            except:
                print(
                    "WARNING: time format %s does not match with %s"
                    % (self.pathPET["naming"], file)
                )
                continue
        return pet

    # computes potential evapotranspiration if not given
    def _det_PET(self, temp, radiation):
        # compute PET based on makkink formula
        # latent heat of vaporization water
        lamb = 1000000 * (2.501 - 0.002361 * temp)  # [J / kg]

        # psychometric constant
        # psi = cp * P / epsilon / lambda * 10 **3 [kPa / °C]
        # cp = specific heat of moist air = 1.013 [kJ / kg / °C]
        # P = athmospheric pressure = 101325 Pa
        # epsilon = molecular weight of water vapor = 0.622
        # lambda ) latent heat of vaporization water
        cp = 1013  # [J / kg / °C]
        P = 101325  # [Pa]
        epsilon = 0.622  # [-]

        psi = cp * P / epsilon / lamb  # [Pa / °C]

        # gradient of vapor pressure
        es = 0.6108 * np.exp(17.27 * temp / (237.3 + temp)) * 1000.0  # [Pa]
        delta = 4098.0 * es / (237.3 + temp) ** 2.0  # [Pa / C°]

        # specific weight of water
        rho = 1000  # [kg / m³]

        # makking potential evapotranspiration choisnel
        # 0.75 for whole europe
        pet = np.maximum(
            1000
            * 0.75
            * (delta / (delta + psi))
            * (self.time_step * radiation / lamb / rho),
            0.0,
        )
        return pet

    # More or less the main routine
    def _convert(self):
        # clear progressbar instances to avoid printout issues
        p_bar = None
        temp_bar = None
        gl_bar = None

        # starting converting
        print(
            "[%s] -> Starting to convert INCA to netcdf files"
            % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
        )
        # converts
        # init base structure of the hdf file
        print(
            "[%s] -> Starting to prepare empty netcdf file"
            % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
        )
        self._prepare_nc()
        print(
            "[%s] -> Finished" % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
        )
        # add a variable to the structure
        # variable names
        # add more variables here
        names = ["P", "TEMP", "PET", "GL"]
        formats = ["f4", "f4", "f4", "f4"]
        standard_names = [
            "Precipitation",
            "Temperature",
            "PotEvaporation",
            "ShortWRadiation",
        ]
        units = ["mm", "C", "mm", "Wpm2"]

        # add variables from list
        print(
            "[%s] -> Adding variables to netcdf file"
            % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
        )
        for i in range(0, len(names)):
            self._add_variable(names[i], formats[i], standard_names[i], units[i])
        print(
            "[%s] -> Finished" % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
        )

        # get timerange input
        timerange = self._get_timerange()

        # split timerange in buffersize
        buffers = self._det_buffer_range(timerange)

        print(
            "[%s] -> Buffersize was chosen to be %s "
            % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"), len(buffers))
        )
        # loop over buffers
        for k in range(0, len(buffers)):
            print(
                "[%s] -> Starting buffer %s/%s "
                % (
                    datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"),
                    k + 1,
                    len(buffers),
                )
            )
            # determine overlap, we use 5 except if the buffer elements are too short
            overlap = int(self.max_time_step / self.time_step * 4)
            # set overlap to 0 if we only have one buffer element or if we are
            if k == 0:
                overlap = 0
            # buffer timerange

            timerange_buff = timerange[buffers[k][0] - overlap : buffers[k][1] + 1]

            if self.pet_option == 0:
                # temperature based model
                print(
                    "Sorry, temperature based approach is not implemented yet. However, I would not recommend using this anyway"
                )
                sys.exit(0)

            # pet directly given
            if self.pet_option == 1:
                # not tested yet
                print("Sorry, this option is not implemented yet")
                sys.exit(0)

            # radiation given : makkink
            if self.pet_option == 2:
                # get index files
                print(
                    "[%s] -> Indexing folder structure"
                    % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
                )
                p = self._index_P(timerange_buff)
                temp = self._index_TEMP(timerange_buff)
                gl = self._index_GL(timerange_buff)
                print(
                    "[%s] -> Finished"
                    % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
                )

                print(
                    "[%s] -> Read maps for given buffer (This will take a while)"
                    % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
                )
                # init dataframes
                df_p = pd.DataFrame(
                    index=[
                        pd.to_datetime(datetime.datetime.strftime(i, "%Y-%m-%d %H:%M"))
                        for i in timerange_buff
                    ],
                    columns=np.arange(self.rows * self.cols),
                    dtype="float",
                )
                print("df_p")
                print(df_p)
                df_temp = pd.DataFrame(
                    index=[
                        pd.to_datetime(datetime.datetime.strftime(i, "%Y-%m-%d %H:%M"))
                        for i in timerange_buff
                    ],
                    columns=np.arange(self.rows * self.cols),
                    dtype="float",
                )
                df_gl = pd.DataFrame(
                    index=[
                        pd.to_datetime(datetime.datetime.strftime(i, "%Y-%m-%d %H:%M"))
                        for i in timerange_buff
                    ],
                    columns=np.arange(self.rows * self.cols),
                    dtype="float",
                )

                # fill dataframes
                p_bar = tqdm(total=len(p), position=0, leave=True)
                # t = trange(len(p), desc='Bar desc', leave=True)
                for i in range(0, len(p)):
                    x, y, data, FillVal = self._read_map(
                        p[i][1], self.pathRR["projection"]
                    )
                    # replace unplausible values with 0
                    data[data < 0.0] = 0
                    df_p.loc[
                        pd.to_datetime(
                            datetime.datetime.strftime(p[i][0], "%Y-%m-%d %H:%M")
                        )
                    ] = data.reshape(self.rows * self.cols)
                    p_bar.set_description(
                        "[%s] -> Processing Precipitation buffer %s/%s"
                        % (
                            datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"),
                            k + 1,
                            len(buffers),
                        )
                    )
                    p_bar.update(1)
                p_bar.close()

                print("df_p_filled")
                print(df_p)

                print(
                    "[%s] -> Precipitation: Finished"
                    % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
                )
                temp_bar = tqdm(total=len(temp), position=0, leave=True)
                for i in range(0, len(temp)):
                    x, y, data, FillVal = self._read_map(
                        temp[i][1], self.pathTT["projection"]
                    )
                    df_temp.loc[
                        pd.to_datetime(
                            datetime.datetime.strftime(temp[i][0], "%Y-%m-%d %H:%M")
                        )
                    ] = data.reshape(self.rows * self.cols)
                    temp_bar.set_description(
                        "[%s] -> Processing Temperature buffer %s/%s"
                        % (
                            datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"),
                            k + 1,
                            len(buffers),
                        )
                    )
                    temp_bar.update(1)
                temp_bar.close()
                print(
                    "[%s] -> Temperature: Finished"
                    % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
                )
                gl_bar = tqdm(total=len(gl), position=0, leave=True)
                for i in range(0, len(gl)):
                    x, y, data, FillVal = self._read_map(
                        gl[i][1], self.pathGL["projection"]
                    )
                    # replace unplausible values with 0
                    data[data > 10000.0] = 0
                    df_gl.loc[
                        pd.to_datetime(
                            datetime.datetime.strftime(gl[i][0], "%Y-%m-%d %H:%M")
                        )
                    ] = data.reshape(self.rows * self.cols)
                    gl_bar.set_description(
                        "[%s] -> Processing Radiation buffer %s/%s"
                        % (
                            datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"),
                            k + 1,
                            len(buffers),
                        )
                    )
                    gl_bar.update(1)
                gl_bar.close()

                print(
                    "[%s] -> Radiation: Finished"
                    % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
                )

                print(
                    "[%s] -> Finished"
                    % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
                )

                # resample dataframe if required
                # temperature and radiation is resampled with a mean approach and linear interpolation between NaNs
                # works in both directions upsampling and downsampling

                print(
                    "[%s] -> Resample time structure and interpolate in space"
                    % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
                )
                df_temp = df_temp.interpolate()
                df_gl = df_gl.interpolate()

                df_p = df_p.bfill().interpolate() / (
                    float(self.pathRR["timestep_secs"]) / self.time_step
                )

                print("df_p_res")
                print(df_p)

                print(
                    "[%s] -> Finished"
                    % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
                )

                print(
                    "[%s] -> Compute Potential Evapotranspiration"
                    % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
                )
                df_pet = self._det_PET(df_temp, df_gl)
                print(
                    "[%s] -> Finished"
                    % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
                )

                # slice dataframes to buffsize
                df_temp = df_temp[overlap : len(df_temp)]
                df_p = df_p[overlap : len(df_p)]
                df_pet = df_pet[overlap : len(df_pet)]
                df_gl = df_gl[overlap : len(df_gl)]

            # add data from buffer to netcdf
            print(
                "[%s] -> Adding buffered data to netcdf"
                % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
            )
            self._add_data(
                "P",
                df_p.to_numpy().reshape(len(df_p), self.rows, self.cols),
                buffers[k],
            )
            self._add_data(
                "TEMP",
                df_temp.to_numpy().reshape(len(df_temp), self.rows, self.cols),
                buffers[k],
            )
            self._add_data(
                "GL",
                df_gl.to_numpy().reshape(len(df_gl), self.rows, self.cols),
                buffers[k],
            )
            self._add_data(
                "PET",
                df_pet.to_numpy().reshape(len(df_pet), self.rows, self.cols),
                buffers[k],
            )
            print(
                "[%s] -> Finished"
                % (datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
            )
