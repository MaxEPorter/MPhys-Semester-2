import datetime
from astropy.time import Time
import numpy as np
import pandas as ps

dt = datetime.datetime.now()
tm = Time(dt)
print(tm.mjd)


a = [1, 2, 3, 4]
b = [2, 3, 4, 5]
print(np.multiply(a, b))

ex = ps.read_csv("C:/Users/Max/Documents/Uni/MPhys/semester 2/flux_densities/J2217+5733/pb_30jul2012_J2217.csv")
print(ex)
print(ex.drop([3], axis=0))
print(ex.drop([5], axis=0))
print(ex.drop([5], axis=0).reset_index())

for i in range(10):
    print(i)