#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NSGA 2 optimizer for wflow hbv

Sebastian Gegenleithner

possible parameters and ranges for a small catchment

parameter explanation

### Snow ####
Cfmax [mm / d /C°], degree day factor for snowmelt, can vary in space and time, values 3 for non forested 0.6 for forested areas
	https://ec-jrc.github.io/lisflood-model/2_04_stdLISFLOOD_snowmelt/
CFR [-], refreezing constant, around 0.05
TT [C°], threshold for snowfall, around 0 c°
TTI [c°], defines how snow can partly all as snow or rain, around 1.0 c°
WHC [-], The model allows for a certain volume of melted water to remain within
	the snowpack, given as a fraction of the corresponding snow water equivalent of the snowpack, about 0.1

### Soil parameters ####
FC [mm], Field capacity, water holding capacity of the soil.
	https://soilgrids.org/ -> Cambisols around 100
BetaSeepage [-], Values of beta vary generally between 1 and 3. Larger values of beta reduce runoff and indicate a higher absorption capacity of the soil
LP [-], soil dependent evaporation factor, around 0.8
K4 [mm/d], recession coefficient of baselow, about 100 days -> 0.01 [mm/d] -> 0.000104166 [mm / 15 min]
KD [mm/d], recession coefficient of deep storage, about 200 days -> 0.005 [mm/d] -> 0.00005208 [mm / 15 min]
!! SET KQUICKFLOW TO 1
K0 [mm/d], recession coefficient of the flood peak, about 0.5 days -> 2 [mm/d] -> 0.020833 [mm / 15 min]
KQuickFlow [mm/d], recession coefficient of the interflow during flooding events, about 2 days -> 0.5 [mm/d] -> 0.005208 [mm / 15 min]
SUZ [mm], threshold over which K0 is used, around 60
PERC [mm / d], perculation from upper to lower zone, about 4 [mm/d] -> 0.041667 [mm / 15 min] 
AlphaPerc [-], percent of what reaches lower and deep storage, around 0.5
Pcorr [-], correction factor for precipitation (including snow), around 1.0
RFCF [-], correction factor for rain fall (excluding snow), around 1.0
SFCF [-], correction factor for snow fall (excluding rain), around 1.0
Cflux [mm / d], maximum capillary rise, aronud 1.0 [mm / d] -> 0.002778 [mm / 15 min]
	file:///tmp/mozilla_iwbworkstation0/water-10-01549.pdf (Improved Process Representation in the Simulation of the Hydrology of a Meso-Scale Semi-Arid Catchment)
ICF [mm], maximum interception storage in areas dependent on landuse, around 0-3 mm
	https://journals.openedition.org/belgeo/17806
CEVPF [-], correction factor for potential evapotranspiration, usually 1.0 or 1.15 for forested areas
EPF [-], exponent of correction factor for potential evapotranspiration, usually 0.0
ECORR [-], correction factor for total evaporation, usually 1.0

# River runoff
!! set NRiverMethod = 1
N [Manning values], surface manning for landuse values
	van osnabrugge, interpolate, simulate assimilate page 165
    
# ranges
# Literature
# https://nhess.copernicus.org/articles/20/3521/2020/
# https://www.researchgate.net/figure/The-HBV-model-parameters-and-their-ranges-used-for-calibration_tbl1_317065756
# https://www.sciencedirect.com/science/article/pii/S0022169403004013?via%3Dihub
# https://www.researchgate.net/figure/HBV-parameter-ranges-upper-band-UB-lower-band-LB-fixed-parameters-have-lower-and_tbl2_319364688

Cfmax [mm / d /C°] ->  0.001 - 10 [mm / d /C°] -> 0.0000104166 - 0.104167 [mm / 15min /C°]
CFR [-] -> 0 - 0.1 [-]
TT [°C] -> -3 - 4 [°C]
TTI [°C] -> 0 - 7 [°C]
WHC [-] -> 0.0 - 0.8 [-]
FC [mm] -> 0.1 - 1000 [mm]
BetaSeepage [-] -> 0.01 - 7.0 [-]
LP [-] -> 0.1 - 1.0 [-]
# recession coefficients from fastest to slowest
K0 [mm/d] -> 0.05 - 2.0 [d] -> 0.5 - 20 [mm/d] -> 0.005208 - 0.208 [mm / 15min]
KQuickFlow [mm/d] -> 2.0 - 30 [d] -> 0.03333 - 0.5 [mm/d] -> 0.0004719 - 0.00521 [mm / 15min]
K4 [mm/d] -> 20 - 120 [d] -> 0.008333 - 0.05 [mm/d] -> 0.0000868021 - 0.0005208333 [mm / 15min]
KD [mm/d] -> 60 - 300 [d] -> 0.003333 - 0.016667 [mm/d] -> 0.0000347188 - 0.000173611 [mm / 15min]
SUZ [mm] -> 5 - 60 [mm]
PERC [mm/d] -> 0.01 - 30 [mm/d] -> 0.00010416667 - 0.3125
AlphaPerc [-] -> 0.1 - 0.9 [-]
Cflux [mm / d] -> 0.0 - 4.0 [mm/d] -> 0.0 - 0.0416667 [mm / 15min]
ICF [mm] -> 0.0 - 5.0 [mm]
CEVPF [-] -> 0.7 - 1.5 [-]
N [s / m^(1/3))] -> 0.01 - 0.8 [s / m^(1/3))]

"""
############################
#       Imports
############################
import subprocess
import sys
import numpy as np
import multiprocessing
from datetime import datetime
import re
import shutil
from pymoo.model.problem import Problem
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.optimize import minimize
from pymoo.factory import get_termination
import pytz
import matplotlib.pyplot as plt

from argparse import ArgumentParser
from logging import error, info, debug
import logging
from os.path import isfile
from os.path import join as pjoin
from configparser import ConfigParser, ExtendedInterpolation

"""Mapping from log level strings to logging levels"""
LOG_LEVEL_MAP = {
    "critical": logging.CRITICAL,
    "error": logging.ERROR,
    "warning": logging.WARNING,
    "info": logging.INFO,
    "debug": logging.DEBUG,
}


############################
#    Helper functions
############################

# cuts data to given time window for evaluation
def cut_data(data, mask_array):
    # timezone
    timezone = pytz.timezone("UTC")
    # lower value
    lower_obj = datetime.strptime(mask_array[0], "%d.%m.%Y %H:%M")
    # lower localized
    lower_localized = timezone.localize(lower_obj)
    # lwoer timestamp
    lower = lower_localized.timestamp()

    # upper value
    upper_obj = datetime.strptime(mask_array[1], "%d.%m.%Y %H:%M")
    # lower localized
    upper_localized = timezone.localize(upper_obj)
    # lwoer timestamp
    upper = upper_localized.timestamp()

    # create mask
    mask = np.where((data[:, 0] >= lower) & (data[:, 0] <= upper))[0]
    # cutout data
    array = data[mask]

    return array


def to_float_nan(v):
    try:
        return float(v)
    except ValueError:
        return np.nan


def read_gauge_data(filepath, timezone):
    time_list = []
    value_list = []

    with open(filepath) as f:
        for line in f:
            line = line.split(",")

            time_list.append(
                timezone.localize(
                    datetime.strptime(line[0], "%d.%m.%Y %H:%M:%S")
                ).timestamp()
            )

            value_list.append(to_float_nan(line[1]))

    return np.column_stack((time_list, value_list))


############################
# wflow model wrapper class
############################
class wflowModel:
    def __init__(
        self, modelpath, inifile, wflow_model, use_col, use_tables, num_cores, timezone
    ):
        # init general model settings
        # path to model folder
        self.modelpath = modelpath
        # name of inifile
        self.inifile = inifile
        # name of modelpath
        self.wflow_model = wflow_model
        # define which coloumn to use for extracting results
        self.use_col = use_col
        # define tables folder
        self.use_tables = use_tables
        # number of cores to utilize
        self.num_cores = num_cores

        # model timezone
        self.timezone = timezone

        # init variables from inifile
        self.init_from_ini()

    # parallel processor
    def run_parallel(self, parameters, ensemble):
        # generate number of tasks, which is equal to the number of samples
        number_of_tasks = ensemble.shape[0]
        # generate a pool of tasks
        # pool = multiprocessing.Pool(number_of_tasks)
        pool = multiprocessing.Pool(self.num_cores)
        # asynch run of model
        info("Spawning processes")
        runs = []
        for i in range(0, number_of_tasks):
            run = pool.apply_async(self.run, (i, parameters, ensemble[i, :]))
            runs.append(run)
        # close and join pooling results
        pool.close()
        pool.join()

        array = None
        for i in range(0, len(runs)):
            run_i = runs[i].get()
            array = run_i if array is None else np.column_stack((array, run_i[:, 1]))

        info("Joining processes")
        return array

    def init_from_ini(self):
        config = ConfigParser()
        config.read(pjoin(self.modelpath, self.inifile))

        self.starttime = self.timezone.localize(
            datetime.strptime(config["run"]["starttime"], "%Y-%m-%d %H:%M:%S")
        ).timestamp()
        self.endtime = self.timezone.localize(
            datetime.strptime(config["run"]["endtime"], "%Y-%m-%d %H:%M:%S")
        ).timestamp()
        self.timestepsecs = float(config["run"]["timestepsecs"])

    def get_results(self, ID):
        # path to results
        respath = self.modelpath + "/" + str(ID) + "/run.tss"

        time_list = []
        value_list = []

        with open(respath) as f:
            i = 0
            for line in f:
                line = line.split()
                try:
                    time_step = self.starttime + i * self.timestepsecs
                    value_step = float(line[self.use_col])
                    value_list.append(value_step)
                    time_list.append(time_step)
                    i = i + 1
                except:
                    continue

        sim_data = np.column_stack((time_list, value_list))

        return sim_data

    # run one model instance
    def run(self, ID, parameters, parameter_values):
        #####################
        #   possible inputs
        #####################
        # -C: casename (directory)
        # -f: force overwrite existing files
        # -T: Set end time of the run: yyyy-mm-dd hh:mm:ss
        # -S: Set start time of the run: yyyy-mm-dd hh:mm:ss
        # -N: No lateral flow, use runoff response function to generate fast runoff
        # -s: Set the model timesteps in seconds
        # -I: re-initialize the initial model conditions with default
        # -i: Set input table directory (default is intbl)
        # -x: run for subcatchment only (e.g. -x 1)
        # -R: set the name runId within the current case
        # -L: set the logfile
        # -c: name of wflow the configuration file (default: Casename/wflow_hbv.ini).
        # -h: print usage information
        # -U: The argument to this option should be a .tss file with measured discharge in
        #     [m^3/s] which the program will use to update the internal state to match
        #     the measured flow. The number of columns in this file should match the
        #     number of gauges in the wflow_gauges.map file.
        # -u: list of gauges/columns to use in update. Format:
        #     -u [1 , 4 ,13]
        #     The above example uses column 1, 4 and 13
        # -P: set parameter change string (e.g: -P "self.FC = self.FC * 1.6") for non-dynamic variables
        # -p: set parameter change string (e.g: -P "self.Precipitation = self.Precipitation * 1.11") for
        #     dynamic variables
        # -l: loglevel (most be one of DEBUG, WARNING, ERROR)
        # -X overwrites the initial values at the end of each timestep

        # generate inputs
        process_list = []

        # mandatory commands
        process_list.append("-C {0}".format(self.modelpath))
        process_list.append("-R {0}".format(ID))
        process_list.append("-c {0}".format(self.inifile))
        process_list.append("-f")
        process_list.append("-i {0}".format(self.use_tables))

        # change variables
        for i in range(0, len(parameters)):
            # genate change string
            change_string = "self.%s = self.%s * %s" % (
                parameters[i],
                parameters[i],
                parameter_values[i],
            )
            # append process
            process_list.append("-P {0}".format(change_string))

        # spawn child process
        debug(f"Starting process {ID}")
        process = subprocess.Popen(
            [sys.executable, self.wflow_model] + process_list,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        # wait for process to finish
        process.communicate()

        # retrieve results
        sim_data = self.get_results(ID)
        return sim_data


############################
#  Optmizer problem class
############################
class OptimizationProblem(Problem):
    def __init__(self, *args, **kwargs):

        self.parameters = kwargs["parameters"]
        self.modelpath = kwargs["modelpath"]
        self.inifile = kwargs["inifile"]
        self.wflow_model = kwargs["wflow_model"]
        self.opt_functions = kwargs["opt_functions"]
        self.const_functions = kwargs["const_functions"]
        self.meas_path = kwargs["gauge_root"]
        self.use_col = kwargs["use_col"]
        self.mask_array = kwargs["mask_array"]
        self.use_tables = kwargs["use_tables"]
        self.num_cores = kwargs["num_cores"]
        self.timezone_gauge = kwargs["timezone_gauge"]
        self.timezone_model = kwargs["timezone_model"]

        # init wflow model
        self.myModel = wflowModel(
            self.modelpath,
            self.inifile,
            self.wflow_model,
            self.use_col,
            self.use_tables,
            self.num_cores,
            self.timezone_model,
        )

        # define time dependent variables
        time_variables = [
            "Cfmax",
            "K0",
            "KQuickFlow",
            "K4",
            "KD",
            "PERC",
            "Cflux",
            "Leakage",
        ]

        # set lookup table for optimizer
        opt_lookup = {
            "Cfmax": "0.001>10",
            "CFR": "0.0>0.1",
            "TT": "-3.0>4.0",
            "TTI": "0.0>7.0",
            "WHC": "0.0>0.8",
            "FC": "0.1>1000.0",
            "BetaSeepage": "0.01>7.0",
            "LP": "0.1>1.0",
            "K0": "0.5>20.0",
            "KQuickFlow": "0.0333333>0.5",
            "K4": "0.00833333>0.05",
            "KD": "0.003333>0.016667",
            "SUZ": "5.0>60.0",
            "PERC": "0.01>30.0",
            "AlphaPerc": "0.1>0.9",
            "Cflux": "0.0>4.0",
            "ICF": "0.0>5.0",
            "CEVPF": "0.7>1.5",
            "N": "0.01>0.8",
            "Leakage": "0.0>2.0",
        }

        # generate lower and higher ranges
        lower_range = []
        upper_range = []

        for param in self.parameters:
            lookup_param = opt_lookup[param].split(">")
            lower_i = float(lookup_param[0])
            higher_i = float(lookup_param[1])
            # parameters including timestep have to be converted
            # 86400. = seconds in day
            if param in time_variables:
                lower_i = lower_i / (86400.0 / self.myModel.timestepsecs)
                higher_i = higher_i / (86400.0 / self.myModel.timestepsecs)
            lower_range.append(lower_i)
            upper_range.append(higher_i)

        self.lower_range = np.asarray(lower_range)
        self.upper_range = np.asarray(upper_range)

        # init measurements
        self.gaugedata = read_gauge_data(self.meas_path, self.timezone_gauge)

        # cut measurement data
        self.cut_measurements = cut_data(self.gaugedata, self.mask_array)

        # init from parent class
        super().__init__(
            n_var=len(self.parameters),
            n_obj=len(self.opt_functions),
            n_constr=len(self.const_functions),
            xl=self.lower_range,
            xu=self.upper_range,
        )

    # evaluate optimizing thingsin local zone
    def _evaluate(self, x, out, *args, **kwargs):

        # perform model runs for this generation
        sim_results = self.myModel.run_parallel(self.parameters, x)

        # cut simulation results
        sim_results_cut = cut_data(sim_results, self.mask_array)

        # evaluate objective functions
        temp = []
        temp_g = []

        for i in range(0, len(x)):
            # loop over objective functions
            obj_functions = []
            const_functions = []
            for obj in self.opt_functions:
                # evaluate nash sutcliffe
                if obj == "NSE":
                    f = (
                        (
                            1
                            - (
                                (
                                    sim_results_cut[:, i + 1]
                                    - self.cut_measurements[:, 1]
                                )
                                ** 2
                            ).sum()
                            / (
                                (
                                    self.cut_measurements[:, 1]
                                    - self.cut_measurements[:, 1].mean()
                                )
                                ** 2
                            ).sum()
                        )
                    ) * -1
                    info("NSE = %s" % f)
                    obj_functions.append(f)
                if obj == "VOL":
                    f = (
                        abs(
                            self.cut_measurements[:, 1] - sim_results_cut[:, i + 1]
                        ).sum()
                        / self.cut_measurements[:, 1].sum()
                    )
                    obj_functions.append(f)

            for const in self.const_functions:
                # evalute max
                if const == "MAXQ":
                    g = (
                        max(sim_results[:, i + 1])
                        - max(self.cut_measurements[:, 1]) * 2.0
                    )
                    const_functions.append(g)

            temp.append(obj_functions)
            temp_g.append(const_functions)

        # append objective functions
        out["F"] = np.column_stack([temp])
        out["G"] = np.column_stack([temp_g])


#################################
#          Main routine
#################################


def main():
    parser = ArgumentParser(description="Prepare wflow files")
    parser.add_argument("config_file", help="configuration file destination")
    args = parser.parse_args()

    if not isfile(args.config_file):
        error(f"Configuration file {args.config_file} doesn't exist")
        exit(1)

    config = ConfigParser(interpolation=ExtendedInterpolation())
    config.read(args.config_file)

    root_logger = logging.getLogger()
    log_formatter = logging.Formatter("%(levelname)s %(asctime)s: %(message)s")
    log_level = LOG_LEVEL_MAP.get(
        config.get("Configuration", "log_level", fallback="INFO").lower(),
        logging.INFO,
    )
    root_logger.setLevel(log_level)

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)

    log_file = config.get("Configuration", "log_file", fallback=None)
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(log_formatter)
        root_logger.addHandler(file_handler)

    #################################
    #          Settings
    #################################

    # wflow model settings
    # path to wflow model to optimize
    modelpath = config["Model"]["modelpath"]
    # name of the steering file
    inifile = config["Model"]["inifile"]
    # path to wflow model, we have to modify some things. For some reason subprocess messes up the strings
    # so for the arguments we have to strip the arguments in the source code like:
    # for o, a in opts:
    # o = o.strip()
    # a = a.strip()
    wflow_model = config["Model"]["wflow_model"]
    # define which tables to use. Parameters which should be optimized should be set to 1 so we get correct
    # values by multiplication
    use_tables = r"intbl_NSGA2"
    # define coloumn for which we evaluate the discharge of the modelling results
    # in run.tss
    use_col = config.getint("Calibration", "run_tss_col")
    # define the number of physical cpu cores to utilize
    num_cores = config.getint("Configuration", "num_cores", fallback=1)

    # measurement settings
    gauge_root = config["Calibration"]["gauge_root"]
    # set timezone of gauge data !!careful timezones are inverted!!
    # so gmt -1 is actually gmt +1
    timezone_gauge = pytz.timezone("Etc/GMT-1")
    timezone_model = pytz.timezone("UTC")

    # optimizer settings
    # define parameters to optimize. The value ranges and time conversions are taken care off
    # automatically
    parameters = [
        "Cfmax",
        "CFR",
        "TT",
        "TTI",
        "WHC",
        "FC",
        "BetaSeepage",
        "LP",
        "K0",
        "KQuickFlow",
        "K4",
        "KD",
        "SUZ",
        "PERC",
        "AlphaPerc",
        "Cflux",
        "ICF",
        "CEVPF",
        "N",
        "Leakage",
    ]
    # define optimizer range for evaluation. This is not the simulation period!!
    # it defines merely between which time stamps the objective functions are
    # evaluated
    # optimizer range, range = [start,end] %d.%m.%Y %H:M%
    mask_array = [config["Calibration"]["from_time"], config["Calibration"]["to_time"]]
    # Define objective functions and constraints
    # options are:
    # NSE: minimize Nash sutcliffe
    # VOL: minimize volume error
    opt_functions = ["NSE", "VOL"]
    # MAXQ: constraint of maximum discharge
    const_functions = ["MAXQ"]
    const_functions = []
    # define population size for NSGA alogrithm
    pop_size = config.getint("Calibration", "pop_size")
    # define number of generations
    n_gen = config.getint("Calibration", "generations")

    #################################
    #    Init and run optimizer
    #################################

    # init our optimization problem
    problem = OptimizationProblem(
        parameters=parameters,
        modelpath=modelpath,
        inifile=inifile,
        wflow_model=wflow_model,
        opt_functions=opt_functions,
        const_functions=const_functions,
        gauge_root=gauge_root,
        use_col=use_col,
        mask_array=mask_array,
        use_tables=use_tables,
        num_cores=num_cores,
        timezone_gauge=timezone_gauge,
        timezone_model=timezone_model,
    )

    # get termination condition, here number of generations
    termination = get_termination("n_gen", n_gen)

    # define optimization alogrithm
    algorithm = NSGA2(
        pop_size,
        n_offsprings=pop_size,
        sampling=get_sampling("real_random"),
        crossover=get_crossover("real_sbx", prob=0.9, eta=15),
        mutation=get_mutation("real_pm", eta=20),
        eliminate_duplicates=True,
    )

    # minimize the problems residuals
    res = minimize(
        problem, algorithm, termination, seed=1, save_history=True, verbose=True
    )

    #################################
    #    Retrieve results for
    #   optimized parameter set
    #################################

    # init wflow object
    reWflow = wflowModel(
        modelpath, inifile, wflow_model, use_col, use_tables, num_cores, timezone_model
    )

    # retrieve the optimized parameter sets
    optpars = res.X
    # remove temp working directories
    for i in range(0, pop_size):
        shutil.rmtree(modelpath + "/{0}".format(i))

    info("Start recomputing optimized results")
    reWflow.run_parallel(parameters, optpars)
    info("Finished optimization")


if __name__ == "__main__":
    main()
