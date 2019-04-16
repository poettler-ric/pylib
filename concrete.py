#!/usr/bin/env python3

from math import isclose
from mechanics import (
    convert_to,
    permille,
    One,
    m,
    MN,
    MNm,
    MPa, GPa,
)
from os.path import dirname, join
from pandas import read_csv
from sympy import (
    symbols,
    solve,
)

# Bst 550 A
_default_E_y = 200 * GPa
epsilon_ylim = 4 * permille

alpha_cc = .85

gamma_c = 1.5
gamma_y = 1.15

__CONCRETE_CSV = join(dirname(__file__), 'concrete.csv')
c_properties = read_csv(__CONCRETE_CSV, sep=';')
c_properties.fck = c_properties.fck * MPa
c_properties.fck_cube = c_properties.fck_cube * MPa
c_properties.fcm = c_properties.fcm * MPa
c_properties.fctm = c_properties.fctm * MPa
c_properties.fctk_005 = c_properties.fctk_005 * MPa
c_properties.fctk_095 = c_properties.fctk_095 * MPa
c_properties.Ecm = c_properties.Ecm * GPa
c_properties.epsilon_c1 = c_properties.epsilon_c1 * permille
c_properties.epsilon_cu1 = c_properties.epsilon_cu1 * permille
c_properties.epsilon_c2 = c_properties.epsilon_c2 * permille
c_properties.epsilon_cu2 = c_properties.epsilon_cu2 * permille
c_properties.epsilon_c3 = c_properties.epsilon_c3 * permille
c_properties.epsilon_cu3 = c_properties.epsilon_cu3 * permille


def fcd(fck):
    return alpha_cc * fck / gamma_c


def fyd(fyk):
    return fyk / gamma_y


def plastification_threshold(fyk, E_y=_default_E_y):
    return convert_to(fyd(fyk) / E_y, permille)


def ny(fyk, area, epsilon_y, E_y=_default_E_y):
    if epsilon_y < plastification_threshold(fyk, E_y):
        return area * E_y * epsilon_y
    else:
        return area * fyd(fyk)


def ny_x(fyk, area, epsilon_y, E_y=_default_E_y):
    return area * E_y * epsilon_y


def alpha_r(fck, epsilon_c):
    row = c_properties[c_properties.fck == fck]
    ec2 = row.epsilon_c2.values.item(0)
    ecu2 = row.epsilon_cu2.values.item(0)
    assert 0 <= epsilon_c and epsilon_c <= 3.5 * permille

    plain = epsilon_c / permille
    if 0 <= epsilon_c and epsilon_c < ec2:
        return plain / 2 / plain ** 2 / 12
    elif ec2 <= epsilon_c and epsilon_c <= ecu2:
        return (3 * plain - 2) / (3 * plain)


def k_a(fck, epsilon_c):
    row = c_properties[c_properties.fck == fck]
    ec2 = row.epsilon_c2.values.item(0)
    ecu2 = row.epsilon_cu2.values.item(0)
    assert 0 <= epsilon_c and epsilon_c <= 3.5 * permille

    plain = epsilon_c / permille
    if 0 <= epsilon_c and epsilon_c < ec2:
        return (8 - plain) / (24 - 4 * plain)
    elif ec2 <= epsilon_c and epsilon_c <= ecu2:
        return (3 * plain**2 - 4 * plain + 2) / (6 * plain**2 - 4 * plain)


def nc_block_coefficient(fck, b, epsilon_c, x):
    return alpha_r(fck, epsilon_c) * fcd(fck) * x * b, k_a(fck, epsilon_c) * x


def nc_x_block_coefficient(fck, b, epsilon_c):
    x = symbols('x')
    return alpha_r(fck, epsilon_c) * fcd(fck) * x * b, k_a(fck, epsilon_c) * x


def mn_pure_bending(fck, fyk, b, h, d1, d2, epsilon_c, as1, as2):
    x = symbols('x')
    d = h - d1
    epsilon_y1_x = epsilon_c / x * (d - x)
    epsilon_y2_x = epsilon_c / x * (x - d2)
    nc_x, d_nc_x = nc_x_block_coefficient(fck, b, epsilon_c)

    threshold = plastification_threshold(fyk)

    # assume x = d2
    # -> epsilon_y2 = 0
    nc = convert_to(nc_x.subs(x, d2), MN)
    ny1 = convert_to(ny(fyk, as1, epsilon_y1_x.subs(x, d2)), MN)

    assert not isclose(nc / MN, ny1 / MN)

    # assume x > d2
    # steel 1 plasticized
    # steel 2 elastic
    nc = nc_x
    ny1 = ny(fyk, as1, 4 * permille)
    ny2_x = ny_x(fyk, as2, epsilon_y2_x)

    res = solve(nc + ny2_x - ny1, x)
    to_m = lambda x: convert_to(x, [m, One])
    valid_values = lambda x: 0 * m <= x and x <= h
    res = list(filter(valid_values, map(to_m, res)))
    # we handle only one valid result
    assert len(res) == 1
    res = res.pop()

    # check assumptions
    assert epsilon_y1_x.subs(x, res) > threshold
    assert epsilon_y2_x.subs(x, res) < threshold
    nc = nc_x.subs(x, res)
    ny2 = ny2_x.subs(x, res)
    m_rd = nc * (d - d_nc_x.subs(x, res)) + ny2 * (d - d2)
    return convert_to(m_rd, MNm).evalf(), 0 * MN


def mn_balanced_point(fck, fyk, b, h, d1, d2, epsilon_c, as1, as2):
    x = h / 2
    d = h - d1

    nc, d_nc = nc_block_coefficient(fck, b, epsilon_c, x)

    epsilon_y1 = epsilon_c / x * (d - x)
    epsilon_y2 = epsilon_c / x * (x - d2)

    ny1 = ny(fyk, as1, epsilon_y1)
    ny2 = ny(fyk, as2, epsilon_y2)

    n_rd = ny1 - nc - ny2
    m_rd = n_rd * (h / 2 - d1) + nc * (d - d_nc) + ny2 * (d - d2)
    return convert_to(m_rd, MNm).evalf(), convert_to(n_rd, MN).evalf()
