#!/usr/bin/env python3

from sympy import pi
from concrete import (
    mn_centric_pull,
    mn_excentric_pull,
    mn_pure_bending,
    mn_balanced_point,
    mn_decompression,
    mn_pure_compression,
)
from mechanics import (
    convert_to,
    permille,
    mm, cm, m,
    kN,
    kNm,
    kPa, MPa, GPa,
)

if __name__ == '__main__':
    import sympy
    sympy.printing.str.StrPrinter._default_settings['abbrev'] = True

    h = .55 * m
    b = .35 * m
    d1 = d2 = .05 * m
    as1 = as2 = 5 * (16 * mm)**2 / 4 * pi

    fck = 25 * MPa
    fyk = 550 * MPa

    results = list()
    results.append(mn_centric_pull(fck, fyk, b, h, d1, d2, as1, as2))
    results.append(mn_excentric_pull(fck, fyk, b, h, d1, d2, as1, as2))
    results.append(mn_pure_bending(fck, fyk, b, h, d1, d2, as1, as2))
    results.append(mn_balanced_point(fck, fyk, b, h, d1, d2, as1, as2))
    results.append(mn_decompression(fck, fyk, b, h, d1, d2, as1, as2))
    results.append(mn_pure_compression(fck, fyk, b, h, d1, d2, as1, as2))

    names = ("centric pull", "excentric pull", "bending", "balanced point",
             "decompression", "pure preassure")
    _to_kN = lambda x: convert_to(x, kN)
    _to_kNm = lambda x: convert_to(x, kNm)
    for n, r in zip(names, results):
        print("{}: {} {}".format(n, _to_kNm(r[0]), _to_kN(r[1])))
