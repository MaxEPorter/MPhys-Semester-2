import pandas
import matplotlib.pyplot as plt
import numpy as np
import datetime
from astropy.time import Time


data = pandas.read_csv('test1.csv')
print(data)
dates = [Time(datetime.datetime.strptime(i, '%d/%m/%Y')).mjd for i in data['date']]
print(dates)
