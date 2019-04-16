#!/usr/bin/env python3

from sympy import pi
from concrete import (
    mn_pure_bending,
    mn_balanced_point,
    mn_decompression,
    mn_pure_compression,
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
    res.append(mn_pure_bending(fck, fyk, b, h, d1, d2, as1, as2))
    res.append(mn_balanced_point(fck, fyk, b, h, d1, d2, as1, as2))
    res.append(mn_decompression(fck, fyk, b, h, d1, d2, as1, as2))
    res.append(mn_pure_compression(fck, fyk, b, h, d1, d2, as1, as2))
    print(res)
