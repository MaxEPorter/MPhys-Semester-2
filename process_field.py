import pandas
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

import input_field


def epoch_to_date(df):

    d = []
    count = 0
    for i in df.columns:
        if count < 2:
            count += 1
            continue

        # d.append(datetime.strptime(i, '%d/%m/%y'))
        print(i)
        d.append(i)

    n = []
    for i in d:
        n.append((i - d[0]).total_seconds())

    return d


def plot_fluxes(flx, err, index):

    d = epoch_to_date(flx)

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.errorbar(d, flx.iloc[index, 2:], err.iloc[index, 2:], linestyle=None, fmt='.')

    ax.set_title(flx.iloc[index]['NVSS id'])
    ax.set_xlabel('Date')
    ax.set_ylabel('Janskys')


def hist_dev(flx, err):

    std = []
    for i in range(flx.shape[0]):
        std.append(np.std(flx.iloc[i, 2:].to_numpy()))

    fig = plt.figure()
    ax = fig.add_subplot()
    print(std)
    ax.hist(std, bins=20)
    ax.set_ylabel('freq')
    ax.set_xlabel('std')


if __name__ == '__main__':
    fluxes, errors = input_field.get_sources(input_field.get_fname())
    print('Fluxes')
    print(fluxes)
    print('Errors')
    print(errors)

    # last arg is the index of source in flux table, eg 0 is first source
    plot_fluxes(fluxes, errors, 7)

    # plot standard deviation of all sources in table as histogram
    # hist_dev(fluxes, errors)

    plt.show()
