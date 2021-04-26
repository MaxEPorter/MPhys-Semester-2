import pandas
import matplotlib.pyplot as plt
import numpy as np
import datetime
from astropy.time import Time


def calcRedX2():
    # calculates reduced chi squared
    chisqrd = 0
    dof = len(x) - 1
    for i in range(0, len(x)):
        E = slope * x[i] + intercept
        O = y[i]
        chisqrd += ((O - E) ** 2) / E

    redchisqrd = chisqrd / dof

    print('Reduced Chi squared:', redchisqrd)



data = pandas.read_csv('test1.csv')
print(data)
x = [Time(datetime.datetime.strptime(i, '%d-%b-%y')).mjd for i in data['date']]
print(x)


