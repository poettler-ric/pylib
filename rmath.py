from numpy import arange, nextafter

def simpson(f, a, b, n=3):
    """
    Simple implementation of simpson's rule.

    f is the function to integrate. a is the starting point, b the end point and
    n is the number of points for which to caclulate the function value.
    """
    steps = n - 1
    result = 0
    c = 0 # counter
    dx = (b - a) / steps
    for i in arange(a, nextafter(b, b+dx), dx):
        factor = 1
        if c == 0:
            pass
        elif c == steps:
            pass
        elif c % 2:
            factor = 4
        else:
            factor = 2
        c += 1
        result += dx / 3 * f(i) * factor
    return result
