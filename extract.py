import numpy as np
import pandas
import os
import matplotlib.pyplot as plt
import datetime
from astropy.time import Time
import json


def extract():

    locs = {}
    with open('source_locations.json', 'r') as f:
        locs = json.load(f)

    data = []

    for filename in os.listdir(locs['j2217']):
        data.append(pandas.read_csv(locs['j2217'] + '/' + filename))

        # get datetime
        day = int(filename[3:5])
        mon = filename[5:8]
        year = int(filename[8:12])
        d = '{}-{}-{}'.format(day, mon, year)
        d = datetime.datetime.strptime(d, '%d-%b-%Y')

        dates = []
        for i in range(data[-1].shape[0]):
            dates.append(Time(d).mjd)

        data[-1]['dates'] = dates

    j17 = pandas.concat(data)

    data = []
    for filename in os.listdir(locs['j2208']):
        data.append(pandas.read_csv(locs['j2208'] + '/' + filename))

        # get datetime
        minus = filename.find('-')
        day = int(filename[minus+1:minus+3])
        mon = filename[minus+3:minus+6]
        year = int(filename[minus+6:minus+8])
        d = '{}-{}-{}'.format(day, mon, year)
        d = datetime.datetime.strptime(d, '%d-%b-%y')

        dates = []
        for i in range(data[-1].shape[0]):
            dates.append(Time(d).mjd)

        data[-1]['dates'] = dates

    j08 = pandas.concat(data)

    print('loaded j17 and j08')

    return j17, j08


def isolate_sources(data, tname):

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

    n = n_entrys(sources)

    print('{} sources'.format(len(sources)))
    print('was {} rows, now {}'.format(data.shape[0], n))

    return sources


def remove_bad(sources, multiplier=10, max_fl=30):

    print('\nflagging...')
    print('multiplyer factor {}'.format(multiplier))

    counter_anom1 = 0
    counter_anom2 = 0

    for i in range(len(sources)):

        dev = np.std(sources[i][' int_flux'])
        mean = np.average(sources[i][' int_flux'])

        ch = multiplier*dev

        j = 0
        length = sources[i].shape[0]
        while j < length:

            try:
                if abs(sources[i].iloc[j, 10] - mean) > ch:
                    sources[i] = sources[i].drop([j], axis=0).reset_index(drop=True)
                    counter_anom1 += 1
                    j = 0

                elif sources[i].iloc[j, 10] > max_fl:
                    sources[i] = sources[i].drop([j], axis=0).reset_index(drop=True)
                    counter_anom2 += 1
                    j = 0

                else:
                    j += 1

            except KeyError:
                continue
            except IndexError:
                if sources[i].shape[0] == 0 or j >= sources[i].shape[0]:
                    break
                print('indexerror')
                print(j)
                print(sources[i])
                exit()

    print('anom1 - {}\nanom2 - {}'.format(counter_anom1, counter_anom2))
    print('{} rows remain'.format(n_entrys(sources)))

    return sources


def n_entrys(sources):
    newrow = 0
    for i in sources:
        newrow += i.shape[0]
    return newrow



if __name__ == '__main__':
    with open('source_locations.json', 'r') as f:
        flocs = json.load(f)

    j2217, j2208 = extract()

    j17sour = isolate_sources(j2217,
                     tname=flocs['j2217_template'])
    flagged = remove_bad(j17sour, multiplier=0, max_fl=0)

    # print(flagged)

    plt.show()



