#!/usr/bin/env python3

__author__ = "Richard Pöttler"
__copyright__ = "Copyright (c) 2021 Richard Pöttler"
__license__ = "MIT"
__email__ = "richard.poettler@gmail.com"

# base units
SECOND = 1
METER = 1
KILOGRAM = 1

# magnitudes
MILLI = 10 ** -3
CENTI = 10 ** -2
DECI = 10 ** -1

PERCENT = 10 ** -2
PERMILLE = 10 ** -3

# abbreviations
M = METER
S = SECOND

# combined units
MM = MILLI * METER  # [mm] -> [m]
CM = CENTI * METER  # [cm] -> [m]
DM = DECI * METER  # [dm] -> [m]

CM2 = CM ** 2  # [cm^2] -> [m^2]

MPS = METER / SECOND  # [m/s]
MPS2 = METER / SECOND ** 2  # [m/s^2]
LPS = DM ** 3 / SECOND  # [l/s] -> [m^3/s]

M3PS = METER ** 3 / SECOND  # [m^3/s]

MINUTE = 60 * SECOND  # [minute] -> [s]
HOUR = 60 * MINUTE  # [hour] -> [s]
DAY = 24 * HOUR  # [day] -> [s]

# constants
GRAVITY = 9.80665 * MPS2
DENSITY_WATER = 1 * KILOGRAM / DM ** 3
KINEMATIC_VISCOSITY = 1.3 * 10 ** -6 * M ** 2 / S  # at 10 degree celsius
