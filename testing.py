import datetime
from astropy.time import Time
import numpy as np

dt = datetime.datetime.now()
tm = Time(dt)
print(tm.mjd)


a = [1, 2, 3, 4]
b = [2, 3, 4, 5]
print(np.multiply(a, b))