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
    j17 = pandas.concat(data)

    data = []
    for filename in os.listdir('C:/Users/Max/Documents/Uni/MPhys/semester 2/flux_densities/J2208+5500'):
        data.append(pandas.read_csv('C:/Users/Max/Documents/Uni/MPhys/semester 2/flux_densities/J2208+5500/' + filename))
    j08 = pandas.concat(data)

    print('loaded j17 and j08')

    return j17, j08


def combine_sources(data):

    ra = []
    de = []
    for i, j in zip(data['ra'], data[' dec']):
        ra.append(i)
        de.append(j)

    # 10arsec = 0.002777 deg
    rab = np.arange(330.6, 336.1, 0.1)
    deb = np.arange(54, 58.4, 0.1)

    plt.hist2d(ra, de, bins=[rab, deb], cmap='Greys')
    plt.colorbar()

    #print(pandas.cut(data['ra']))
    d1 = data.assign(acut=pandas.cut(data['ra'], bins=rab, labels=False), bcut=pandas.cut(data[' dec'], bins=deb, labels=False))

    d2 = d1.assign(bin=pandas.Categorical(d1.filter(regex='cut').apply(tuple, 1)))

    plt.figure()
    plt.scatter(d2['acut'], d2['bcut'])

    return d2


def isolate_sources(data, tname='C:/Users/Max/Documents/Uni/MPhys/semester 2/flux_densities/J2217+5733/pb_30jul2012_J2217.csv'):

    template = pandas.read_csv(tname)

    print('\nisolating...')
    print('using template {}'.format(tname))
    radius = 0.02

    print('looking {:.1f} arcsec square around template'.format(radius/0.000277))

    sources = []

    for index, row in template.iterrows():
        sources.append(pandas.DataFrame())  # columns=['ra', ' ra_err']))

        # 1deg = 0.000277arsec
        rmin = row['ra'] - radius
        rmax = row['ra'] + radius
        dmin = row[' dec'] - radius
        dmax = row[' dec'] + radius

        for i in range(data.shape[0]):

            ra = data.iloc[i, 0]
            dec = data.iloc[i, 2]

            if rmax > ra > rmin and dmax > dec > dmin:

                sources[index] = sources[index].append(data.iloc[[i]].copy(), ignore_index=True)

    newrow = 0
    for i in sources:
        newrow += i.shape[0]

    print('{} sources'.format(len(sources)))
    print('was {} rows, now {}'.format(data.shape[0], newrow))

    return sources


def remove_bad(soudata):

    print('\nflagging...')
    multi = 10
    print('multiplyer factor {}'.format(multi))
    counter = 0

    for i in range(len(soudata)):
        dev = np.std(soudata[i][' int_flux'])
        ch = multi*dev
        for j in range(soudata[i].shape[0]):
            if soudata[i].iloc[j, 10] > ch:
                soudata[i].drop([j], axis=0)
                counter += 1

    print('removed {} rows'.format(counter))


def test(data):

    stds = []

    done = []
    for index, row in data.iterrows():

        pix = row['bin']

        if pix in done:
            continue

        flux = []
        for index2, row2 in data.iterrows():
            if row2['bin'] == pix:
                # looping over a source now
                flux.append(row2[' int_flux'])

        stds.append(np.std(flux))

        done.append(pix)


    plt.figure()
    plt.hist(stds, bins=20)


if __name__ == '__main__':
    j2217, j2208 = extract()
    #sourced = combine_sources(j2217)
    #test(sourced)


    j17sour = isolate_sources(j2217,
                     tname='C:/Users/Max/Documents/Uni/MPhys/semester 2/flux_densities/J2217+5733/pb_30jul2012_J2217.csv')
    remove_bad(j17sour)


    plt.show()



