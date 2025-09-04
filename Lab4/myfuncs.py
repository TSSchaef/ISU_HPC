
def factorial(x):
    if x == 0 or x == 1:
        return 1
    s = 1
    for k in range(2, x + 1):
        s = s * k
    return s


def exponentiate(x):
    e = 2.7182818284590451
    x0 = int(round(x))
    z = x - x0

    # e^x  e^x0 (1 + z + (z^2/2!) + z^3/3! ...)

    taylorDepth = 20
    taylorSum = 1 + z
    for k in range(2, taylorDepth):
        taylorSum += (z**k / factorial(k))

    return (e**x0) * (taylorSum)


def logarithm(x):
    NewtonDepth = 100
    s = x

    # Newton Method (repeated tangents)
    for k in range(0, NewtonDepth):
        s = s - 1 + (x * exponentiate(-1 * s))
    return s
    

def sqrt(x):
    if x <= 0:
        return 0

    NewtonDepth = 100
    s = x
    # Newton Method (repeated tangents)
    for k in range(0, NewtonDepth):
        s = (s + (x / s)) / 2
    return s
