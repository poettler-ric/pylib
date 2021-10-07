#!/usr/bin/env python3

from math import log, sqrt, atan
from scipy.optimize import fsolve
from units import GRAVITY, KINEMATIC_VISCOSITY

__author__ = "Richard Pöttler"
__copyright__ = "Copyright (c) 2021 Richard Pöttler"
__license__ = "MIT"
__email__ = "richard.poettler@gmail.com"


def f_ruehlmann(y):
    """F function for the Ruehlmann calculation."""
    return 1 / 6 * log((y ** 2 + y + 1) / (y - 1) ** 2) + 1 / sqrt(3) * atan(
        (1 + 2 * y) / sqrt(3)
    )


def ruehlmann_rect(h, hn, h0, h_cr, i):
    """Calculates the distance of a waterlevel according to Ruehlmann."""
    y = h / hn
    y0 = h0 / hn
    return (
        hn
        / i
        * (
            h0 / hn
            - h / hn
            + (1 - (h_cr / hn) ** 3) * (f_ruehlmann(y) - f_ruehlmann(y0))
        )
    )


def t_n_rect(discharge, strickler_roughness, inclination, width, start=1):
    """Calculates the normal depth of a rectangular channel."""

    A = lambda w, h: w * h
    U = lambda w, h: w + 2 * h
    R = lambda w, h: A(w, h) / U(w, h)

    return fsolve(
        lambda h: strickler_roughness
        * inclination ** (1 / 2)
        * R(width, h) ** (2 / 3)
        * A(width, h)
        - discharge,
        start,
    )[0]


def i_r_rect(discharge, strickler_roughness, width, t_1, t_2):
    """Calculates the inclination of the energy line based on the strickler value."""
    # calculate area and wetted perimeter on median values
    a_m = width * (t_1 + t_2) / 2
    u_m = width + t_1 + t_2
    r_m = a_m / u_m

    return (discharge / (a_m * strickler_roughness * r_m ** (2 / 3))) ** 2


def l_transition_i_r_rect(discharge, strickler_roughness, width, t_1, t_2, inclination):
    """Calculates the transition lenght based only on the inclination of the energy line"""
    i_r = i_r_rect(discharge, strickler_roughness, width, t_1, t_2)
    v_1 = discharge / (width * t_1)
    v_2 = discharge / (width * t_2)
    return (t_2 + v_2 ** 2 / (2 * GRAVITY) - (t_1 + v_1 ** 2 / (2 * GRAVITY))) / (
        inclination - i_r
    )


def t_crit_rect(discharge, width):
    """Calculates the critical depth of a rectangular channel."""
    return (discharge ** 2 / (width ** 2 * GRAVITY)) ** (1 / 3)


def lambda_turbulent_rough(k, d):
    """Calculates lambda for pipe loss for rough conditions"""
    return (1 / (2 * log(k / d / 3.71, 10))) ** 2


def lambda_turbulent_transition(k, d, re):
    """Calculates lambda for pipe loss for transition (rough->smoth) conditions"""
    return fsolve(
        lambda lam: (1 / (2 * log(2.51 / (re * sqrt(lam)) + k / d / 3.71, 10))) ** 2
        - lam,
        lambda_turbulent_rough(k, d),
    )[0]


def d_hyd(width, height):
    """Calculates hydraulic diameter for rectangular profiles"""
    return 4 * (width * height) / (2 * (width + height))


def reynolds_number(v, d):
    """Calculates the reynolds number"""
    return v * d / KINEMATIC_VISCOSITY


def pipe_loss(l, a, k, d, q):
    """Calculates the pipe loss for a given discharge"""
    v = q / a
    re = reynolds_number(v, d)
    lam = (
        lambda_turbulent_rough(k, d)
        if re * k / d > 1300
        else lambda_turbulent_transition(k, d, re)
    )
    return lam * l / d / (2 * GRAVITY * a ** 2) * q ** 2


def local_loss(nu, a, q):
    """Calculates a local loss for a given discharge"""
    return nu / (2 * GRAVITY * a ** 2) * q ** 2
