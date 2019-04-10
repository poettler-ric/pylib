#!/usr/bin/env python3

from math import pi, sin, cos
from scipy.integrate import quad

STUDENT_ID = "1530149"

USE_TABLE_VALUES = True

kPa2MPa = lambda x: x * 10**-3

MPa2kPa = lambda x: x * 10**3
GPa2kPa = lambda x: x * 10**6

A = int(STUDENT_ID[-3])
B = int(STUDENT_ID[-2])
C = int(STUDENT_ID[-1])

alpha = pi / 6  # rad
l1 = 1.9 + (A / 3 * .1)  # m
l2 = 8  # m
l3 = 1  # m
l4 = 3.9 + (C / 3 * .1)  # m

h = .55  # m
b = .35  # m
d1 = d2 = .05  # m
d = h - d1  # m

gk = 26 + B / 4  # kN/m
qk = 40  # kN/m

gd = 1.35 * gk  # kN/m
qd = 1.5 * qk  # kN/m

gamma_c = 1.5
gamma_y = 1.15
alpha_cc = .85

fck = MPa2kPa(25)  # kPa
fcd = alpha_cc * fck / gamma_c  # kPa

fctm = MPa2kPa(.3 * kPa2MPa(fck)**(2 / 3))  # kPa
fctk005 = .7 * fctm  # kPa

fcm = fck + MPa2kPa(8)  # kPa
Ec = GPa2kPa(22 * (kPa2MPa(fcm) / 10)**.3)  # kPa

fyk = MPa2kPa(550)  # kPa
fyd = fyk / gamma_y  # kPa

Ey = MPa2kPa(200_000)  # kPa

Iy = b * h**3 / 12  # m^4

if USE_TABLE_VALUES:
    fctm = MPa2kPa(2.6)
    fctk005 = MPa2kPa(1.8)
    fcm = MPa2kPa(33)
    Ec = GPa2kPa(31)


class Beam:
    pass


class Support:
    pass


class System0:

    def __init__(self, q1, q2, q3):
        self.support = Support()
        self.support.H = 0
        self.support.V = q1 * l1 + q2 * l2 + q3 * l3
        self.support.M = -q1 * l1**2 / 2 + q2 * l2**2 / 2 * cos(alpha) + \
            q3 * l3 * (l2 + l3 / 2) * cos(alpha)

        self.beam1 = Beam()
        self.beam1.length = l1
        self.beam1.N = lambda x: 0
        self.beam1.Q = lambda x: -q1 * x
        self.beam1.M = lambda x: -q1 * x**2 / 2

        self.beam2 = Beam()
        self.beam2.length = l2
        self.beam2.N = lambda x: \
            q2 * sin(alpha) * l2 * (x / l2 - 1) - q3 * sin(alpha) * l3
        self.beam2.Q = lambda x: \
            q2 * cos(alpha) * l2 * (1 - x / l2) + q3 * cos(alpha) * l3
        self.beam2.M = lambda x: \
            q2 * cos(alpha) * l2**2 * (-x**2 / (2 * l2**2) + x / l2 - 1 / 2) + \
            q3 * cos(alpha) * l3**2 * (x / l3 - l2 / l3 - 1 / 2)

        self.beam3 = Beam()
        self.beam3.length = l3
        self.beam3.N = lambda x: q3 * sin(alpha) * l3 * (x / l3 - 1)
        self.beam3.Q = lambda x: q3 * cos(alpha) * l3 * (1 - x / l3)
        self.beam3.M = lambda x: \
            q3 * cos(alpha) * l3**2 * (-x**2 / (2 * l3**2) + x / l3 - 1 / 2)

        self.beam4 = Beam()
        self.beam4.length = l4
        self.beam4.N = lambda x: -self.support.V
        self.beam4.Q = lambda x: -self.support.H
        self.beam4.M = lambda x: self.support.M


class System1:

    def __init__(self):
        self.support = Support()
        self.support.H = 0
        self.support.V = -1
        self.support.M = -1 * l2 * cos(alpha)

        self.beam1 = Beam()
        self.beam1.length = l1
        self.beam1.N = lambda x: 0
        self.beam1.Q = lambda x: 0
        self.beam1.M = lambda x: 0

        self.beam2 = Beam()
        self.beam2.length = l2
        self.beam2.N = lambda x: 1 * sin(alpha)
        self.beam2.Q = lambda x: -1 * cos(alpha)
        self.beam2.M = lambda x: 1 * cos(alpha) * l2 * (1 - x / l2)

        self.beam3 = Beam()
        self.beam3.length = l3
        self.beam3.N = lambda x: 0
        self.beam3.Q = lambda x: 0
        self.beam3.M = lambda x: 0

        self.beam4 = Beam()
        self.beam4.length = l4
        self.beam4.N = lambda x: -self.support.V
        self.beam4.Q = lambda x: -self.support.H
        self.beam4.M = lambda x: self.support.M


class CombinedSystem:

    def __init__(self, sys0, sys1, x1):
        self.support = Support()
        self.support.H = sys0.support.H + x1 * sys1.support.H
        self.support.V = sys0.support.V + x1 * sys1.support.V
        self.support.M = sys0.support.M + x1 * sys1.support.M

        self.beam1 = Beam()
        self.beam1.length = l1
        self.beam1.N = lambda x: sys0.beam1.N(x) + x1 * sys1.beam1.N(x)
        self.beam1.Q = lambda x: sys0.beam1.Q(x) + x1 * sys1.beam1.Q(x)
        self.beam1.M = lambda x: sys0.beam1.M(x) + x1 * sys1.beam1.M(x)

        self.beam2 = Beam()
        self.beam2.length = l2
        self.beam2.N = lambda x: sys0.beam2.N(x) + x1 * sys1.beam2.N(x)
        self.beam2.Q = lambda x: sys0.beam2.Q(x) + x1 * sys1.beam2.Q(x)
        self.beam2.M = lambda x: sys0.beam2.M(x) + x1 * sys1.beam2.M(x)

        self.beam3 = Beam()
        self.beam3.length = l3
        self.beam3.N = lambda x: sys0.beam3.N(x) + x1 * sys1.beam3.N(x)
        self.beam3.Q = lambda x: sys0.beam3.Q(x) + x1 * sys1.beam3.Q(x)
        self.beam3.M = lambda x: sys0.beam3.M(x) + x1 * sys1.beam3.M(x)

        self.beam4 = Beam()
        self.beam4.length = l4
        self.beam4.N = lambda x: sys0.beam4.N(x) + x1 * sys1.beam4.N(x)
        self.beam4.Q = lambda x: sys0.beam4.Q(x) + x1 * sys1.beam4.Q(x)
        self.beam4.M = lambda x: sys0.beam4.M(x) + x1 * sys1.beam4.M(x)


def d_xx(E, I, sys0, sys1):
    m1 = lambda x: sys0.beam1.M(x) * sys1.beam1.M(x)
    m2 = lambda x: sys0.beam2.M(x) * sys1.beam2.M(x)
    m3 = lambda x: sys0.beam3.M(x) * sys1.beam3.M(x)
    m4 = lambda x: sys0.beam4.M(x) * sys1.beam4.M(x)
    return 1 / (E * I) * \
        (quad(m1, 0, sys0.beam1.length)[0] +
         quad(m2, 0, sys0.beam2.length)[0] +
         quad(m3, 0, sys0.beam3.length)[0] +
         quad(m4, 0, sys0.beam4.length)[0])


def energy_theorem(E, I, sys0, sys1):
    d10 = d_xx(E, I, sys1, sys0)
    d11 = d_xx(E, I, sys1, sys1)
    return d10, d11, -d10 / d11


def heading1(heading):
    print()
    print("=" * (12 + len(heading)))
    print(f"===   {heading}   ===")
    print("=" * (12 + len(heading)))


def heading2(heading):
    print()
    print(heading)
    print("=" * len(heading))


def heading3(heading):
    print()
    print(heading)
    print("-" * len(heading))


def print_parameters():
    heading2("Parameters")
    print(f"A: {A} B: {B} C: {C}")
    print(f"l1: {l1}m l2: {l2}m l3: {l3}m l4: {l4}m")
    print(f"d: {d}m")
    print(f"gk: {gk}kN/m gd: {gd}kN/m")
    print(f"qk: {qk}kN/m qd: {qd}kN/m")
    print(f"fck: {kPa2MPa(fck)}MPa fcd: {kPa2MPa(fcd)}MPa")
    print(f"fctk005: {kPa2MPa(fctk005)}MPa fctm: {kPa2MPa(fctm)}MPa")
    print(f"Ec: {kPa2MPa(Ec)}MPa")
    print(f"fyk: {kPa2MPa(fyk)}MPa fyd: {kPa2MPa(fyd)}MPa")
    print(f"Ey: {kPa2MPa(Ey)}MPa")
    print(f"Iy: {Iy}m^4")


def print_support(name, support, unit_force, unit_momentum):
    heading3(f"Support: {name}")
    print(f"Horizontal: {support.H}{unit_force}")
    print(f"Vertical: {support.V}{unit_force}")
    print(f"Horizontal: {support.M}{unit_momentum}")


def print_beam(name, beam, unit_force, unit_momentum):
    heading3(f"Beam: {name}")
    print(f"N(0): {beam.N(0)}{unit_force} N(l): {beam.N(beam.length)}{unit_force}")
    print(f"Q(0): {beam.Q(0)}{unit_force} Q(l): {beam.Q(beam.length)}{unit_force}")
    print(f"M(0): {beam.M(0)}{unit_momentum} "
          + f"M(l/2): {beam.M(beam.length/2)}{unit_momentum} "
          + f"M(l): {beam.M(beam.length)}{unit_momentum}")


def print_system(name, sys, unit_force, unit_momentum):
    heading2(name)
    print_support("support", sys.support, unit_force, unit_momentum)
    print_beam("1", sys.beam1, unit_force, unit_momentum)
    print_beam("2", sys.beam2, unit_force, unit_momentum)
    print_beam("3", sys.beam3, unit_force, unit_momentum)
    print_beam("4", sys.beam4, unit_force, unit_momentum)


def print_energy_theorem(d10, d11, x1):
    heading2("Energy Theorem")
    print(f"d10: {d10} d11: {d11} x1: {x1}")

if __name__ == '__main__':
    unit_force = "kN"
    unit_momentum = "kNm"

    print(f"student#: {STUDENT_ID}")
    print_parameters()

    heading1("Bending")

    sys0 = System0(gd + qd, gd + qd, gd)
    print_system("System 0", sys0, unit_force, unit_momentum)
    sys1 = System1()
    print_system("System 1", sys1, unit_force, unit_momentum)
    d10, d11, x1 = energy_theorem(Ec, Iy, sys0, sys1)
    print_energy_theorem(d10, d11, x1)
    sys = CombinedSystem(sys0, sys1, x1)
    print_system("System Bending", sys, unit_force, unit_momentum)
