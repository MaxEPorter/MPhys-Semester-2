import numpy as np
import pandas
import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import datetime


def extract():

    data = []

    for filename in os.listdir('C:/Users/Max/Documents/Uni/MPhys/semester 2/flux_densities/J2217+5733'):
        data.append(pandas.read_csv('C:/Users/Max/Documents/Uni/MPhys/semester 2/flux_densities/J2217+5733/' + filename))

    for filename in os.listdir('C:/Users/Max/Documents/Uni/MPhys/semester 2/flux_densities/J2208+5500'):
        data.append(pandas.read_csv('C:/Users/Max/Documents/Uni/MPhys/semester 2/flux_densities/J2208+5500/' + filename))

    conc = pandas.concat(data)

    return conc


def combine_sources():

    data = extract()

    ra = []
    de = []
    for i, j in zip(data['ra'], data[' dec']):
        ra.append(i)
        de.append(j)

    # 10arsec = 0.002777 deg
    rab = np.arange(330.6, 336.1, 0.15)
    deb = np.arange(54, 58.4, 0.15)

    plt.hist2d(ra, de, bins=[rab, deb], cmap='Greys')
    plt.colorbar()

    print(pandas.cut(data['ra']))




if __name__ == '__main__':
    extract()
    combine_sources()

    plt.show()



