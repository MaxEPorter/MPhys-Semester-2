import numpy as np
import pandas
import os
import matplotlib.pyplot as plt
import datetime
from astropy.time import Time
import json

plt.style.use('seaborn-whitegrid')
plt.rcParams["font.family"] = "serif"


def extract_all():
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
        day = int(filename[minus + 1:minus + 3])
        mon = filename[minus + 3:minus + 6]
        year = int(filename[minus + 6:minus + 8])
        d = '{}-{}-{}'.format(day, mon, year)
        d = datetime.datetime.strptime(d, '%d-%b-%y')

        dates = []
        for i in range(data[-1].shape[0]):
            dates.append(Time(d).mjd)

        data[-1]['dates'] = dates

    j08 = pandas.concat(data)

    print('loaded j17 and j08')

    return j17, j08


def extract():

    locs = {}
    with open('source_locations.json', 'r') as f:
        locs = json.load(f)

    data = []

    for filename in os.listdir(locs['j2217']):
        if filename == 'pb_04may2012_J2217.csv':
            continue

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

        if filename == 'j22-17nov13.csv':
            continue

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
    counter_too_few_row = 0
    counter_duplicates = 0
    counter_too_few_source = 0

    for i in range(len(sources)):

        dev = np.std(sources[i][' int_flux'])
        mean = np.average(sources[i][' int_flux'])

        ch = multiplier*dev

        j = 0
        length = sources[i].shape[0]
        while j < length:

            try:
                # check for outliers
                if abs(sources[i].iloc[j, 10] - mean) > ch:
                    sources[i] = sources[i].drop([j], axis=0).reset_index(drop=True)
                    counter_anom1 += 1
                    j = 0

                # check for very high jy bad data
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

        # check for duplicate datapoints
        done = []
        lenyfv = sources[i].shape[0]
        for j in range(sources[i].shape[0]):

            try:
                date = sources[i].iloc[j, 14]
            except IndexError:
                continue

            if sources[i].iloc[j, 14] not in done:
                done.append(sources[i].iloc[j, 14])
                duplicates = []
                for k in range(sources[i].shape[0]):
                    if sources[i].iloc[k, 14] == done[-1]:
                        duplicates.append(k)
                if len(duplicates) == 1:
                    continue
                else:
                    closest = 100
                    closindx = 0
                    for k in duplicates:
                        dis = abs(mean - sources[i].iloc[k, 10])
                        if dis < closest:
                            closest = dis
                            closindx = k
                    for k in duplicates:
                        if k == closindx:
                            continue
                        else:
                            sources[i] = sources[i].drop([k], axis=0)
                            counter_duplicates += 1

            sources[i] = sources[i].reset_index(drop=True)

        # check number of sources
        if sources[i].shape[0] < 8:
            counter_too_few_source += 1
            for j in range(sources[i].shape[0]):
                sources[i] = sources[i].drop([j], axis=0)
                counter_too_few_row += 1

    print('anom1 - {}\nanom2 - {}\ntoo_few_row - {}\ntoo_few_source - {}\nduplicates - {}'.format(counter_anom1, counter_anom2, counter_too_few_row, counter_too_few_source, counter_duplicates))

    new_sources = []
    for i in sources:
        if i.shape[0] < 5:
            continue
        else:
            new_sources.append(i.copy())

    print('{} rows remain'.format(n_entrys(sources)))

    return new_sources


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


def plot_sources(sources):

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


def plot_by_epoch(data):

    print(data)

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.invert_xaxis()
    #colours = np.linspace(0, 1, 16)
    #counter = 0
    done = []
    for i in range(data.shape[0]):
        ras = []
        decs = []
        date = data.iloc[i, 14]
        if date not in done:
            done.append(date)
            for j in range(data.shape[0]):
                # print(data.iloc[j, 14], date)
                if data.iloc[j, 14] == date:
                    ras.append(data.iloc[j, 0])
                    decs.append(data.iloc[j, 2])


            ax.scatter(ras, decs)#, color=(colours[counter], 0.5, 0))
            #counter += 1
        else:
            continue


def source_by_pos(sources, ra, dec):

    rad = 0.02
    ra_min = ra - rad
    ra_max = ra + rad
    dec_min = dec - rad
    dec_max = dec + rad

    for i in sources:
        if i.shape[0] == 0:
            continue
        if ra_min < i['ra'][0] < ra_max and dec_min < i[' dec'][0] < dec_max:
            return i

    print('Cant find source by position')
    return None


def read_source(ra, dec):

    rad = 0.02
    ra_min = ra - rad
    ra_max = ra + rad
    dec_min = dec - rad
    dec_max = dec + rad

    with open('source_locations.json', 'r') as f:
        locs = json.load(f)

    source = None
    counter = 0
    for filename in os.listdir(locs['sources']):

        r = float('{}.{}'.format(filename[0:4], filename[5:11]))
        d = float('{}.{}'.format(filename[12:14], filename[15:-4]))

        if ra_min < r < ra_max and dec_min < d < dec_max:

            counter += 1
            source = pandas.read_csv(locs['sources'] + '/' + filename)
            print('found {} near {},{}'.format(filename, ra, dec))

    if counter > 1:
        print('WARNING FOUND MULTIPLE SOURCES AT {}, {}'.format(ra, dec))
    elif counter == 0:
        print('FOUND NO SOURCES AT {},{}'.format(ra, dec))

    return source


def read_sources():

    sources = []
    with open('source_locations.json', 'r') as f:
        locs = json.load(f)

    for filename in os.listdir(locs['sources']):
        sources.append(pandas.read_csv(locs['sources'] + '/' + filename))

    return sources


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
    j17sour = isolate_sources(j2217,
                     tname=flocs['j2217_template'], ignore='j2217_regionignore.csv')
    j17_flagged = remove_bad(j17sour, multiplier=10, max_fl=30)
    plot_by_epoch(j2217)

    plot_all(j2217)
    # plot_data(j17sour)
    #plot_sources(j17_flagged)
    # average_sources(j17_flagged, write=False)

    j08_sour = isolate_sources(j2208, tname=flocs['j2208_template'], ignore='j2208_regionignore.csv')
    j08_flag = remove_bad(j08_sour, multiplier=10, max_fl=30)
    # average_sources(j08_flag, write=True)

    all_sources = j17_flagged + j08_flag
    # average_sources(all_sources, write=False)

    # sources = read_sources()
    # plot_sources(sources)


    # plot_all(j2208)
    # plot_sources(j08_flag)

    plt.show()



