#!/usr/bin/env python3

from sympy import pi
from concrete import (
    alpha_cc, gamma_c, gamma_y,
    mn_pure_bending,
    mn_balanced_point,
)
from mechanics import (
    convert_to,
    permille,
    mm, cm, m,
    kPa, MPa, GPa,
)

if __name__ == '__main__':
    print("=== start ===")
    import sympy
    sympy.printing.str.StrPrinter._default_settings['abbrev'] = True

    h = .45 * m
    b = .4 * m
    d1 = d2 = .05 * m
    as1 = as2 = 2 * (30 * mm)**2 / 4 * pi

    fck = 30 * MPa
    fyk = 550 * MPa

    res = list()

    epsilon_c = 3.5 * permille
    res.append(mn_pure_bending(fck, fyk, b, h, d1, d2, epsilon_c, as1, as2))
    res.append(mn_balanced_point(fck, fyk, b, h, d1, d2, epsilon_c, as1, as2))
    print(res)
