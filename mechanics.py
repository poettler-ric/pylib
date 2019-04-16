#!/usr/bin/env python3

from sympy import S
from sympy.physics.units import (
    convert_to,
    Quantity,
    kilo, mega, giga,
    force, length, pressure,
    N, pa, kPa, percent, permille,
    newton, pascal, meter,
    mm, cm, m,
)

One = S.One

torque = force * length

kN = kilonewton = Quantity("kilonewton", abbrev="kN")
kilonewton.set_dimension(force)
kilonewton.set_scale_factor(kilo * newton)

MN = meganewton = Quantity("meganewton", abbrev="MN")
meganewton.set_dimension(force)
meganewton.set_scale_factor(mega * newton)

GN = giganewton = Quantity("giganewton", abbrev="GN")
giganewton.set_dimension(force)
giganewton.set_scale_factor(giga * newton)

MPa = megapascal = Quantity("megapascal", abbrev="MPa")
megapascal.set_dimension(pressure)
megapascal.set_scale_factor(mega * pascal)

GPa = gigapascal = Quantity("gigapascal", abbrev="GPa")
gigapascal.set_dimension(pressure)
gigapascal.set_scale_factor(giga * pascal)

Nm = newtonmeter = Quantity("newtonmeter", abbrev="Nm")
newtonmeter.set_dimension(torque)
newtonmeter.set_scale_factor(newton * meter)

kNm = kilonewtonmeter = Quantity("kilonewtonmeter", abbrev="kNm")
kilonewtonmeter.set_dimension(torque)
kilonewtonmeter.set_scale_factor(kilo * newton * meter)

MNm = meganewtonmeter = Quantity("meganewtonmeter", abbrev="MNm")
meganewtonmeter.set_dimension(torque)
meganewtonmeter.set_scale_factor(mega * newton * meter)

GNm = giganewtonmeter = Quantity("giganewtonmeter", abbrev="GNm")
giganewtonmeter.set_dimension(torque)
giganewtonmeter.set_scale_factor(giga * newton * meter)
