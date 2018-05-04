from numpy import arange, nextafter

def simpson(f, a, b, n=3):
    """
    Simple implementation of simpson's rule.

    f is the function to integrate. a is the starting point, b the end point and
    n is the number of points for which to caclulate the function value.
    """
    steps = n - 1
    result = 0
    i = 0 # counter
    dx = (b - a) / steps
    for x in arange(a, nextafter(b, b+dx), dx):
        factor = 1
        if i == 0:
            pass
        elif i == steps:
            pass
        elif i % 2:
            factor = 4
        else:
            factor = 2
        i += 1
        result += dx / 3 * f(x) * factor
    return result
