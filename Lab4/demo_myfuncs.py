import myfuncs
import numpy

print("My exponential e^2.467: ", myfuncs.exponentiate(2.467))
print("Numpy exponential e^2.467: ", numpy.exp(2.467))
print("My exponential e^12.738: ", myfuncs.exponentiate(12.738))
print("Numpy exponential e^12.738: ", numpy.exp(12.738), "\n")

print("My log 2.467: ", myfuncs.logarithm(2.467))
print("Numpy log 2.467: ", numpy.log(2.467))
print("My log 12.738: ", myfuncs.logarithm(12.738))
print("Numpy log 12.738: ", numpy.log(12.738), "\n")

print("My sqrt 2.467: ", myfuncs.sqrt(2.467))
print("Numpy sqrt 2.467: ", numpy.sqrt(2.467))
print("My sqrt 12.738: ", myfuncs.sqrt(12.738))
print("Numpy sqrt 12.738: ", numpy.sqrt(12.738), "\n")

print("My factorial 6: ", myfuncs.factorial(6))
print("My factorial 9: ", myfuncs.factorial(9))
