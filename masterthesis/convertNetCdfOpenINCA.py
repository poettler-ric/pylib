#!/usr/bin/env python3
"""
++++++++++++++++++++++++++++ wflowModelC +++++++++++++++++++++++++++++++

    @brief
        Converts gdal readable forcings to netcdf with P, TEMP, PET
    
    @note
        
    @libraries

    @history 31/05/2022 -- Sebastian Gegenleithner
        First version
        
"""

############################################
#               Imports
############################################

# Import from model creator libary
from ModelCreator.NETCDFhandler import Series2nc

import datetime

# rain info as dictionary
# naming has to include some kind of time stamp
raindata = {
    "root": r"/media/iwbworkstation/My Passport/Dissertation/3_Hydrology/Muerz/basedata/INCA_DATA/free_Data/netcdf_files/RR",
    "projection": "epsg:31287",
    "naming": "Muerz_RR_%Y%m%d_%H%M.tif",
    "timestep_secs": "3600",
}
tempdata = {
    "root": r"/media/iwbworkstation/My Passport/Dissertation/3_Hydrology/Muerz/basedata/INCA_DATA/free_Data/netcdf_files/TT",
    "projection": "epsg:31287",
    "naming": "Muerz_RR_%Y%m%d_%H%M.tif",
    "timestep_secs": "3600",
}
raddata = {
    "root": r"/media/iwbworkstation/My Passport/Dissertation/3_Hydrology/Muerz/basedata/INCA_DATA/free_Data/netcdf_files/GL",
    "projection": "epsg:31287",
    "naming": "Muerz_RR_%Y%m%d_%H%M.tif",
    "timestep_secs": "3600",
}

# resultfolder = r'/home/iwbworkstation/Desktop/working_dir/model_rerun_paper1/sbm_rerun/Stanzbach/inmaps_040506_GL_corr.nc'
resultfolder = r"/home/iwbworkstation/Desktop/working_dir/model_rerun_paper1/sbm_rerun/Muerz/inmaps_openINCA_new_3years.nc"

outproj = "epsg:32633"
# as datetime object
start_date = datetime.datetime.strptime("01.01.2015 00:00", "%d.%m.%Y %H:%M")
# as datetime object
# end_date = datetime.datetime.strptime('31.08.2005 23:45', '%d.%m.%Y %H:%M')
# end_date = datetime.datetime.strptime('31.12.2006 23:45', '%d.%m.%Y %H:%M')
end_date = datetime.datetime.strptime("01.01.2018 00:00", "%d.%m.%Y %H:%M")

bounds = [497000.0, 5249000.0, 566200.0, 5295000.0]
resolution = 400.0
time_step = 900.0

## bounds
# bounds = [528900.0,5249000.0,547100.0,5261000.0]
# resolution = 200.
#
# time_step = 900.

zamgnc = Series2nc(
    resultfolder,
    start_date,
    end_date,
    raindata,
    tempdata,
    outproj,
    bounds,
    resolution,
    pathGL=raddata,
    timestep=time_step,
)

zamgnc._convert()
