import extract
import variability
import matplotlib.pyplot as plt
import json
import numpy as np
import pandas

plt.style.use('seaborn-whitegrid')
plt.rcParams["font.family"] = "serif"

'''
1) Non variable graph (w/ chi squared + variability coefficient)
2) bunch of variable graphs (w/ chi squared + variability coefficient)
(w/ source co-ordinates for each graph)
(w/ field too plz)
iphone
chocolate
yacht
'''


def request_1(sources):

    example_coord = [335.83, 56.82]
    coord = [332.216, 54.702]

    selected_source = extract.source_by_pos(all_sources, *coord)
    analysed = variability.analyse([selected_source])

    fog = plt.figure()
    ox = fog.add_subplot()
    for i in sources:
        ox.scatter(i['ra'], i[' dec'], color='mediumseagreen')

    ox.scatter(selected_source['ra'], selected_source[' dec'], color='crimson')
    ox.invert_xaxis()

    fweights = 1. / np.power(selected_source[' int_flux_err'], 2)
    faverage = np.average(selected_source[' int_flux'], weights=fweights)

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.errorbar(selected_source['dates'], selected_source[' int_flux'], yerr=selected_source[' int_flux_err'], fmt='.', color='black')
    ax.plot([min(selected_source['dates']), max(selected_source['dates'])], [faverage, faverage], color='mediumseagreen')
    ax.set_ylabel('Flux density (Jy)')
    ax.set_xlabel('time (MJD)')
    #ax.legend()

    info = 'RA={}\nDec={}' '\n' r'$\xi_\nu={:.3f}\pm {:.4f} $' '\n' r'$ V_\nu = {:.3f}$' '\n' r'$\eta_\nu = {:.3f} $'.format(*coord, analysed['int_flux'][0], analysed['int_flux_err'][0], analysed['variation'][0], analysed['eta'][0])
    #ax.text(56200, 1.2, info)
    ax.annotate(info, xy=(0, 1), xytext=(0, 0.7), xycoords='axes fraction')


def plot_removed():
    with open('source_locations.json', 'r') as f:
        flocs = json.load(f)

    j2217, j2208 = extract.extract()
    j17sour = extract.isolate_sources(j2217,
                                      tname=flocs['j2217_template'], ignore='j2217_regionignore.csv')
    j17_flagged = extract.remove_bad(j17sour, multiplier=4, max_fl=30)

    j08_sour = extract.isolate_sources(j2208, tname=flocs['j2208_template'], ignore='j2208_regionignore.csv')
    j08_flag = extract.remove_bad(j08_sour, multiplier=4, max_fl=30)
    # average_sources(j08_flag, write=True)

    flagged_sources = j17_flagged + j08_flag

    #all_sources = extract.extract_all()
    all_sources = extract.extract()
    all_sources = pandas.concat(all_sources)

    fig = plt.figure()
    ax = fig.add_subplot()


    #ax.scatter(all_sources['ra'], all_sources[' dec'], color='crimson', label='Rejected')
    #for s in flagged_sources:
    #    ax.scatter(s['ra'], s[' dec'], color='mediumseagreen')

    #'''

    f_r = []
    f_d = []
    a_r = []
    a_d = []
    for all_index in range(all_sources.shape[0]):

        found = False
        s_counter = 0
        while s_counter < len(flagged_sources):

            for flag_index in range(flagged_sources[s_counter].shape[0]):
                if flagged_sources[s_counter].iloc[flag_index, 0] == all_sources.iloc[all_index, 0] and flagged_sources[s_counter].iloc[flag_index, 2] == all_sources.iloc[all_index, 2]:
                    a_r.append(all_sources.iloc[all_index, 0])
                    a_d.append(all_sources.iloc[all_index, 2])
                    found = True
                    s_counter = 1000
                    break

            s_counter += 1

        if not found:
            f_r.append(all_sources.iloc[all_index, 0])
            f_d.append(all_sources.iloc[all_index, 2])

    ax.scatter(a_r, a_d, color='mediumseagreen', label='Accepted')
    ax.scatter(f_r, f_d, color='crimson', label='Rejected')
    
    # '''

    ax.set_xlabel('R.A (Deg)')
    ax.set_ylabel('Dec (Deg)')
    ax.invert_xaxis()
    leg = ax.legend(fancybox=True, frameon=True)
    leg.get_frame().set_edgecolor('black')
    leg.get_frame().set_facecolor('white')


if __name__ == '__main__':

    '''
    with open('source_locations.json', 'r') as f:
        flocs = json.load(f)

    j2217, j2208 = extract.extract()

    j17sour = extract.isolate_sources(j2217, tname=flocs['j2217_template'], ignore='j2217_regionignore.csv')
    j17_flagged = extract.remove_bad(j17sour, multiplier=3, max_fl=30)

    j08_sour = extract.isolate_sources(j2208, tname=flocs['j2208_template'], ignore='j2208_regionignore.csv')
    j08_flag = extract.remove_bad(j08_sour, multiplier=3, max_fl=30)

    all_sources = j17_flagged + j08_flag

    # request_1(all_sources)

    '''

    plot_removed()

    plt.show()
