#!/usr/bin/env python3

from math import pi, sin, cos
from rmath import simpson as s

A = 1
B = 4
C = 9

l1 = 1 # [m]
l2 = 10 # [m]
l3 = 2.5 + (B/3 * 0.1) # [m]
l4 = 6 + (A/3 *.1) # [m]

b = .35 # [m]
h = .55 # [m]

Iy = b * h**3 / 12 # [m^2]
E_c = 22 * (33/10)**.3 * 10**6 # [kPa]

alpha = 15 / 180 * pi # [rad]

l1s = l1 * cos(alpha) # [m]
l2s = l2 * cos(alpha) # [m]
l1ss = l1 * sin(alpha) # [m]
l2ss = l2 * sin(alpha) # [m]

gk = 15 + C/4 # [kN]
qk = 20 # [kN]

Bh1 = 0
Bv1 = -1
Mb1 = l2s

def loads(gk, qk, field1, field2, field3):
    result = [ 1.35 * gk, 1.35 * gk, 1.35 * gk ]
    if field1:
        result[0] = 1.35 * gk + 1.5 * qk
    if field2:
        result[1] = 1.35 * gk + 1.5 * qk
    if field3:
        result[2] = 1.35 * gk + 1.5 * qk
    return result

def N11(x):
    return 0

def Q11(x):
    return 0

def M11(x):
    return 0

def N21(x):
    return -sin(alpha)

def Q21(x):
    return cos(alpha)

def M21(x):
    return cos(alpha) * x

def N31(x):
    return 0

def Q31(x):
    return 0

def M31(x):
    return 0

def N41(x):
    return -Bv1

def Q41(x):
    return 0

def M41(x):
    return -Mb1


def Bh0(q1, q2, q3):
    return 0

def Bv0(q1, q2, q3):
    return l1*q1 + l2*q2 + l3*q3

def Mb0(q1, q2, q3):
    return q3* l3**2 / 2 - q1*l1*(l2s + l1s / 2) - q2 * l2 *l2s / 2

def N10(x, q1, q2, q3):
    return q1 * sin(alpha) * x

def Q10(x, q1, q2, q3):
    return -q1 * cos(alpha) * x

def M10(x, q1, q2, q3):
    return -q1 * x**2 /2 * cos(alpha)

def N20(x, q1, q2, q3):
    return N10(l1, q1, q2, q3) + q2 * sin(alpha) * x

def Q20(x, q1, q2, q3):
    return Q10(l1, q1, q2, q3) - q2 * cos(alpha) * x

def M20(x, q1, q2, q3):
    return M10(l1, q1, q2, q3) + Q10(l1, q1, q2, q3) * x - q2 * cos(alpha) * x**2 /2

def N30(x, q1, q2, q3):
    return 0

def Q30(x, q1, q2, q3):
    return q3 * x

def M30(x, q1, q2, q3):
    return -q3 * x**2 /2

def N40(x, q1, q2, q3):
    return -Bv0(q1, q2, q3)

def Q40(x, q1, q2, q3):
    return 0

def M40(x, q1, q2, q3):
    return -Mb0(q1, q2, q3)

def print_system1():
    print("=== System 1 ===")
    print("n11(0): {}".format(N11(0)))
    print("n11(l1): {}".format(N11(l1)))
    print("n21(0): {}".format(N21(0)))
    print("n21(l2): {}".format(N21(l2)))
    print("n31(0): {}".format(N31(0)))
    print("n31(l3): {}".format(N31(l3)))
    print("n41(0): {}".format(N41(0)))
    print("n41(l4): {}".format(N41(l4)))
    print("===")
    print("q11(0): {}".format(Q11(0)))
    print("q11(l1): {}".format(Q11(l1)))
    print("q21(0): {}".format(Q21(0)))
    print("q21(l2): {}".format(Q21(l2)))
    print("q31(0): {}".format(Q31(0)))
    print("q31(l3): {}".format(Q31(l3)))
    print("q41(0): {}".format(Q41(0)))
    print("q41(l4): {}".format(Q41(l4)))
    print("===")
    print("m11(0): {}".format(M11(0)))
    print("m11(l1): {}".format(M11(l1)))
    print("m21(0): {}".format(M21(0)))
    print("m21(l2): {}".format(M21(l2)))
    print("m31(0): {}".format(M31(0)))
    print("m31(l3): {}".format(M31(l3)))
    print("m41(0): {}".format(M41(0)))
    print("m41(l4): {}".format(M41(l4)))

def print_system0(q1, q2, q3):
    print("=== System 0 ===")
    print("n10(0): {}".format(N10(0, q1, q2, q3)))
    print("n10(l1): {}".format(N10(l1, q1, q2, q3)))
    print("n20(0): {}".format(N20(0, q1, q2, q3)))
    print("n20(l2): {}".format(N20(l2, q1, q2, q3)))
    print("n30(0): {}".format(N30(0, q1, q2, q3)))
    print("n30(l3): {}".format(N30(l3, q1, q2, q3)))
    print("n40(0): {}".format(N40(0, q1, q2, q3)))
    print("n40(l4): {}".format(N40(l4, q1, q2, q3)))
    print("===")
    print("q10(0): {}".format(Q10(0, q1, q2, q3)))
    print("q10(l1): {}".format(Q10(l1, q1, q2, q3)))
    print("q20(0): {}".format(Q20(0, q1, q2, q3)))
    print("q20(l2): {}".format(Q20(l2, q1, q2, q3)))
    print("q30(0): {}".format(Q30(0, q1, q2, q3)))
    print("q30(l3): {}".format(Q30(l3, q1, q2, q3)))
    print("q40(0): {}".format(Q40(0, q1, q2, q3)))
    print("q40(l4): {}".format(Q40(l4, q1, q2, q3)))
    print("===")
    print("m10(0): {}".format(M10(0, q1, q2, q3)))
    print("m10(l1/2): {}".format(M10(l1/2, q1, q2, q3)))
    print("m10(l1): {}".format(M10(l1, q1, q2, q3)))
    print("m20(0): {}".format(M20(0, q1, q2, q3)))
    print("m20(l2/2): {}".format(M20(l2/2, q1, q2, q3)))
    print("m20(l2): {}".format(M20(l2, q1, q2, q3)))
    print("m30(0): {}".format(M30(0, q1, q2, q3)))
    print("m30(l3/2): {}".format(M30(l3/2, q1, q2, q3)))
    print("m30(l3): {}".format(M30(l3, q1, q2, q3)))
    print("m40(0): {}".format(M40(0, q1, q2, q3)))
    print("m40(l4): {}".format(M40(l4, q1, q2, q3)))

print("E_c: {} kPa".format(E_c))

print("l1: {} m".format(l1))
print("l2: {} m".format(l2))
print("l3: {} m".format(l3))
print("l4: {} m".format(l4))

print("gk: {} kN/m".format(gk))
print("qk: {} kN/m".format(qk))

print_system1()

# field moment 2
print("=== fieldmoment l2")
(q1, q2, q3) = loads(gk, qk, False, True, False)
print("q1: {} kN/m".format(q1))
print("q2: {} kN/m".format(q2))
print("q3: {} kN/m".format(q3))
print_system0(q1, q2, q3)

# energy theorem
print("=== energy theorem")
d10_m2 = lambda x: M20(x, q1, q2, q3) * M21(x)
d10_m4 = lambda x: M40(x, q1, q2, q3) * M41(x)
d10 = 1/(E_c * Iy) * (s(d10_m2, 0, l2) + s(d10_m4, 0, l4))
print("d10: {}".format(d10))
d11_m2 = lambda x: M21(x,) * M21(x)
d11_m4 = lambda x: M41(x,) * M41(x)
d11 = 1/(E_c * Iy) * (s(d11_m2, 0, l2) + s(d11_m4, 0, l4))
print("d11: {}".format(d11))
x1 = -d10 / d11
print("x1: {}".format(x1))

m2 = lambda x: M20(x, q1, q2, q3) + M21(x) * x1
print("max moment beam 2: {} kNm".format(m2(l2/2)))
n2 = lambda x: N20(x, q1, q2, q3) + N21(x) * x1
print("corresponding normal force beam 2: {} kN".format(n2(l2/2)))

m4 = lambda x: M40(x, q1, q2, q3) + M41(x) * x1
print("max moment beam 4: {} kNm".format(m4(l4)))
n4 = lambda x: N40(x, q1, q2, q3) + N41(x) * x1
print("corresponding normal force beam 4: {} kN".format(n4(l4)))

# field moment 3
print("=== fieldmoment l3")
(q1, q2, q3) = loads(gk, qk, False, False, True)
print("q1: {} kN/m".format(q1))
print("q2: {} kN/m".format(q2))
print("q3: {} kN/m".format(q3))
print_system0(q1, q2, q3)

# energy theorem
print("=== energy theorem")
d10_m2 = lambda x: M20(x, q1, q2, q3) * M21(x)
d10_m4 = lambda x: M40(x, q1, q2, q3) * M41(x)
d10 = 1/(E_c * Iy) * (s(d10_m2, 0, l2) + s(d10_m4, 0, l4))
print("d10: {}".format(d10))
d11_m2 = lambda x: M21(x,) * M21(x)
d11_m4 = lambda x: M41(x,) * M41(x)
d11 = 1/(E_c * Iy) * (s(d11_m2, 0, l2) + s(d11_m4, 0, l4))
print("d11: {}".format(d11))
x1 = -d10 / d11
print("x1: {}".format(x1))

m3 = lambda x: M30(x, q1, q2, q3) + M31(x) * x1
print("max moment beam 3: {} kNm".format(m3(l3)))
n3 = lambda x: N30(x, q1, q2, q3) + N31(x) * x1
print("corresponding normal force beam 3: {} kN".format(n3(l3)))
sf3 = lambda x: Q30(x, q1, q2, q3) + Q31(x) * x1
print("corresponding shearing force beam 3: {} kN".format(sf3(l3)))

