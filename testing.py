import datetime
from astropy.time import Time

dt = datetime.datetime.now()
tm = Time(dt)
print(tm.mjd)


