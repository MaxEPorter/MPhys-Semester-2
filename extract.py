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

        # print(filename)
        # print(data[-1].shape[0])

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


def isolate_sources(data, tname, ignore):

    template = pandas.read_csv(tname)
    bad_region = pandas.read_csv(ignore)

    print('\nisolating...')
    print('using template {}'.format(tname))
    radius = 0.02

    print('looking {:.1f} arcsec square around template'.format(radius/0.000277))

    sources = []
    listlen = -1
    for index, row in template.iterrows():

        # Ignore bad region from file
        skip = False
        for j in range(bad_region.shape[0]):
            if bad_region['ra_max'][j] > row['ra'] > bad_region['ra_min'][j] and bad_region['dec_max'][j] > row[' dec'] > bad_region['dec_min'][j]:
                skip = True
                break

        if skip:
            continue

        # 1deg = 0.000277arsec
        rmin = row['ra'] - radius
        rmax = row['ra'] + radius
        dmin = row[' dec'] - radius
        dmax = row[' dec'] + radius

        added_source = False
        for i in range(data.shape[0]):

            ra = data.iloc[i, 0]
            dec = data.iloc[i, 2]

            if rmax > ra > rmin and dmax > dec > dmin:

                if added_source:
                    pass
                else:
                    sources.append(pandas.DataFrame())  # columns=['ra', ' ra_err']))
                    listlen += 1
                    added_source = True

                sources[listlen] = sources[listlen].append(data.iloc[[i]].copy(), ignore_index=True)


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


def average_sources(sources, write=False):

    s_average = pandas.DataFrame(columns=['ra', 'ra_err', 'dec', 'dec_err', 'smaj', 'smaj_err', 'smin', 'smin_err',
                                          'pa', 'pa_err', 'int_flux', 'int_flux_err', 'pk_flux', 'pk_flux_err', 'n_points'])

    with open('source_locations.json', 'r') as f:
        floc = json.load(f)

    for source in sources:

        if source.shape[0] < 3:
            continue

        raweights = 1./np.power(source[' ra_err'], 2)
        decweights = 1./np.power(source[' dec_err'], 2)
        smajweights = 1./np.power(source[' smaj_err'], 2)
        sminweights = 1./np.power(source[' smin_err'], 2)
        paweights = 1./np.power(source[' pa_err'], 2)
        fluxweights = 1./np.power(source[' int_flux_err'], 2)
        pkweights = 1./np.power(source[' pk_flux_err'], 2)

        raaverage  = np.average(source['ra'], weights=raweights)
        raerr = 1./np.sqrt(np.sum(raweights))

        decaverage = np.average(source[' dec'], weights=decweights)
        decerr = 1./np.sqrt(np.sum(decweights))

        smajaverage = np.average(source[' smaj'], weights=smajweights)
        smajerr = 1./np.sqrt(np.sum(smajweights))

        sminaverage = np.average(source[' smin'], weights=sminweights)
        sminerr = 1./np.sqrt(np.sum(sminweights))

        paaverage = np.average(source[' pa'], weights=paweights)
        paerr = 1./np.sqrt(np.sum(paweights))

        fluxaverage = np.average(source[' int_flux'], weights=fluxweights)
        fluxerr = 1./np.sqrt(np.sum(fluxweights))

        pkaverage = np.average(source[' pk_flux'], weights=pkweights)
        pkerr = 1./np.sqrt(np.sum(pkweights))

        s_average = s_average.append({'ra': raaverage, 'ra_err': raerr, 'dec': decaverage, 'dec_err': decerr,
                                      'smaj': smajaverage, 'smaj_err': smajerr, 'smin': sminaverage,
                                      'smin_err': sminerr, 'pa': paaverage, 'pa_err': paerr, 'int_flux': fluxaverage,
                                      'int_flux_err': fluxerr, 'pk_flux': pkaverage, 'pk_flux_err': pkerr,
                                      'n_points': source.shape[0]}, ignore_index=True)

        if write:

            rastr = '{:.6f}'.format(raaverage).replace('.', '-')
            decstr = '{:.6f}'.format(decaverage).replace('.', '-')
            source.to_csv(floc['sources'] + '/{} {}.csv'.format(rastr, decstr), encoding='utf-8', index=False)

    if write:

        s_average = average_sources(sources)
        s_average.to_csv(floc['sources'] + '/average.csv', encoding='utf-8', index=False)

    print(s_average)
    return s_average


def plot_data(sources):

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax.set_xlabel('dates')
    ax.set_ylabel('int_flux')
    ax2 = fig.add_subplot(212)
    ax2.set_xlabel('dates')
    ax2.set_ylabel(r'$\frac{1}{\sigma^2}$')

    fog = plt.figure()
    ox = fog.add_subplot()
    ox.set_xlabel('ascension')
    ox.set_ylabel('declination')
    for i in sources:
        ax.scatter(i['dates'], i[' int_flux'])
        ax2.scatter(i['dates'], 1/np.power(i[' int_flux_err'], 2))

        ox.scatter(i['ra'], i[' dec'])
    ox.invert_xaxis()


def plot_all(data):

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.scatter(data['ra'], data[' dec'])
    ax.set_xlabel('ascension')
    ax.set_ylabel('declination')
    ax.set_title('All data')
    ax.invert_xaxis()


def n_entrys(sources):
    newrow = 0
    for i in sources:
        newrow += i.shape[0]
    return newrow


if __name__ == '__main__':
    with open('source_locations.json', 'r') as f:
        flocs = json.load(f)

    j2217, j2208 = extract()
    # plot_all(j2217)
    j17sour = isolate_sources(j2217,
                     tname=flocs['j2217_template'], ignore='j2217_regionignore.csv')
    f17_flagged = remove_bad(j17sour, multiplier=10, max_fl=30)
    # plot_data(j17sour)
    # plot_data(f17_flagged)
    # average_sources(f17_flagged, write=False)

    j08_sour = isolate_sources(j2208, tname=flocs['j2208_template'], ignore='j2208_regionignore.csv')
    j08_flag = remove_bad(j08_sour, multiplier=10, max_fl=30)
    plot_all(j2208)
    plot_data(j08_flag)

    plt.show()



