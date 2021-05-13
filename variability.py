import extract
import matplotlib.pyplot as plt
import numpy as np
import scipy.odr as odr
import json
from astropy.coordinates import SkyCoord

plt.style.use('seaborn-whitegrid')
plt.rcParams["font.family"] = "serif"


def f_lin(B, x):
    # linear function, B[0] gradient, b[1] intercept
    return B[0]*x + B[1]


def f_flat(B, x):
    return B[0] + 0*x


def std_hist(sources):

    std = []
    for i in sources:
        temp = np.std(i[' int_flux'])
        if temp>50:
            continue

        std.append(temp)

    plt.hist(std, bins=60)


def analyse(sources):

    print('Analysing...')

    analysed = {
        'ra': [],
        'ra_err': [],
        'dec': [],
        'dec_err': [],
        'smaj': [],
        'smaj_err': [],
        'smin': [],
        'smin_err': [],
        'pa': [],
        'pa_err': [],
        'int_flux': [],
        'int_flux_err': [],
        'pk_flux': [],
        'pk_flux_err': [],
        'n_points': [],

        'eta': [],

        'variation': [],
        'variation_err': [],

        'unusual': []
    }

    for index, source in enumerate(sources):
        # looping over sources

        if source.shape[0] < 4:
            continue

        raweights = 1. / np.power(source[' ra_err'], 2)
        decweights = 1. / np.power(source[' dec_err'], 2)
        smajweights = 1. / np.power(source[' smaj_err'], 2)
        sminweights = 1. / np.power(source[' smin_err'], 2)
        paweights = 1. / np.power(source[' pa_err'], 2)
        fluxweights = 1. / np.power(source[' int_flux_err'], 2)
        pkweights = 1. / np.power(source[' pk_flux_err'], 2)

        analysed['ra'].append(np.average(source['ra'], weights=raweights))
        analysed['ra_err'].append(1. / np.sqrt(np.sum(raweights)))

        analysed['dec'].append(np.average(source[' dec'], weights=decweights))
        analysed['dec_err'].append( 1. / np.sqrt(np.sum(decweights)))

        analysed['smaj'].append(np.average(source[' smaj'], weights=smajweights))
        analysed['smaj_err'].append(1. / np.sqrt(np.sum(smajweights)))

        analysed['smin'].append(np.average(source[' smin'], weights=sminweights))
        analysed['smin_err'].append(1. / np.sqrt(np.sum(sminweights)))

        analysed['pa'].append(np.average(source[' pa'], weights=paweights))
        analysed['pa_err'].append(1. / np.sqrt(np.sum(paweights)))

        analysed['int_flux'].append(np.average(source[' int_flux'], weights=fluxweights))
        analysed['int_flux_err'].append(1. / np.sqrt(np.sum(fluxweights)))

        analysed['pk_flux'].append(np.average(source[' pk_flux'], weights=pkweights))
        analysed['pk_flux_err'].append(1. / np.sqrt(np.sum(pkweights)))

        analysed['n_points'].append(source.shape[0])

        fweights = 1./np.power(source[' int_flux_err'], 2)
        faverage = np.average(source[' int_flux'], weights=fweights)

        """
        # dates is x,  int_flux is y
        linear_mod = odr.Model(f_flat)
        data = odr.Data(source['dates'], source[' int_flux'], we=1./np.power(source[' int_flux_err'], 2))
        od = odr.ODR(data, linear_mod, beta0=[0.5])
        output = od.run()
        # output.pprint()
        """

        # print('odr av -> {}, faverage -> {}'.format(output.beta, faverage))

        # ---------------------  flux density coefficient of variation ------------------
        standard_dev = np.std(source[' int_flux'], ddof=1)
        analysed['variation'].append(standard_dev/faverage)


        # ------------------------ CHISQRD -----------------------
        eta = 0
        dof = source.shape[0] - 1
        for i in range(source.shape[0]):
            # looping over row in source

            eta += np.power(source[' int_flux'][i] - faverage, 2)/np.power(source[' int_flux_err'][i], 2)

        red_eta = eta / dof
        analysed['eta'].append(red_eta)

        if red_eta < 100:
            analysed['unusual'].append(index)
        elif analysed['variation'][-1] < .3:
            analysed['unusual'].append(index)

    print('found {} unusual points'.format(analysed['unusual']))
    for i in analysed['unusual']:
        print('{}, {}'.format(analysed['ra'][i], analysed['dec'][i]))

    return analysed


def plot_eta(analysed):

    plt.figure()
    plt.hist(analysed['eta'], color='mediumpurple', edgecolor='black', bins=10**np.linspace(-1, 3, 40))
    plt.xscale('log')
    plt.xlabel(r'$\eta_\nu$')
    plt.ylabel('Frequency')

    plt.figure()
    plt.hist(analysed['eta'], color='mediumpurple', edgecolor='black', bins=np.linspace(0, 15, 40))
    plt.xlabel(r'$\eta_\nu$')
    plt.ylabel('Frequency')


def plot_variablity(analysed):

    plt.figure()
    plt.hist(analysed['variation'], bins=40, color='mediumpurple', edgecolor='black')
    plt.xlabel('flux density coefficient of variation')
    plt.ylabel('Frequency')


def plot_var_vs_eta(analysed):

    plt.figure()
    plt.scatter(analysed['int_flux'], analysed['variation'])
    plt.xlabel('average flux')
    plt.ylabel('variation')

    plt.figure()
    plt.scatter(analysed['variation'], analysed['eta'])
    plt.xlabel('variation')
    plt.ylabel(r'$\eta_\nu$')
    plt.yscale('log')

    plt.figure()
    plt.scatter(analysed['eta'], analysed['int_flux'])


def light_curve(sources, ra, dec):

    if dec < 56:
        xlim = [330.5, 333.5]
        ylim = [54, 56]
    else:
        xlim = [332.8, 336.3]
        ylim = [56.5, 58.5]

    fig = plt.figure(figsize=(7, 9))
    ox = fig.add_subplot(211)
    for i in sources:
        ox.scatter(i['ra'], i[' dec'], color='mediumseagreen')

    s = extract.source_by_pos(sources, ra, dec)

    ox.scatter(s['ra'], s[' dec'], color='crimson')
    ox.invert_xaxis()
    ox.set_xlim(xlim)
    ox.set_ylim(ylim)

    fweights = 1. / np.power(s[' int_flux_err'], 2)
    faverage = np.average(s[' int_flux'], weights=fweights)

    # fig = plt.figure()
    ax = fig.add_subplot(212)
    ax.errorbar(s['dates'], s[' int_flux']/1000, yerr=s[' int_flux_err']/1000, fmt='.', color='black')
    ax.plot([min(s['dates']), max(s['dates'])], [faverage, faverage], color='mediumseagreen')
    ax.set_ylabel('Flux density (mJy)')
    ax.set_xlabel('time (MJD)')


def plot_unusual(sources, analysed):

    counter = 0
    for i in analysed['unusual']:

        if analysed['dec'][i] < 56:
            xlim = [330.5, 333.5]
            ylim = [54, 56]
            field = 'J2208'
        else:
            xlim = [332.8, 336.3]
            ylim = [56.5, 58.5]
            field = 'J2217'

        fig = plt.figure(figsize=(7, 10))
        ox = fig.add_subplot(211)
        for j in sources:
            ox.scatter(j['ra'], j[' dec'], color='mediumseagreen')

        s = extract.source_by_pos(sources, analysed['ra'][i], analysed['dec'][i])

        ox.scatter(s['ra'], s[' dec'], color='crimson')
        ox.set_xlim(xlim)
        ox.set_ylim(ylim)
        ox.invert_xaxis()

        fweights = 1. / np.power(s[' int_flux_err'], 2)
        faverage = np.average(s[' int_flux'], weights=fweights)

        pos = SkyCoord(analysed['ra'][i], analysed['dec'][i], frame='icrs', unit='deg')
        # pos_str = '{}h{}m{}s +{}d{}m{}s'.format(*pos.ra.hms, pos.dec.day, pos.dec.minute, pos.dec.second)

        # fig = plt.figure()
        ax = fig.add_subplot(212)
        ax.errorbar(s['dates'], s[' int_flux']*1000, yerr=s[' int_flux_err']*1000, fmt='.', color='black')
        ax.plot([min(s['dates']), max(s['dates'])], [faverage*1000, faverage*1000], color='mediumseagreen',
                label='{}' '\n' '{}' '\n' r'$\xi_\nu= ({:.1f}\pm {:.1f} )mJy$' '\n' r'$ V_\nu = {:.3f}$' '\n' r'$\eta_\nu = {:.3f} $'.format(field, pos.to_string('hmsdms'), analysed['int_flux'][i]*1000, analysed['int_flux_err'][i]*1000, analysed['variation'][i], analysed['eta'][i]))
        ax.set_ylabel('Flux density (mJy)')
        ax.set_xlabel('time (MJD)')

        ax.legend()
        leg = ax.legend(handlelength=0, handletextpad=0, fancybox=True, frameon=True)
        leg.get_frame().set_edgecolor('black')
        leg.get_frame().set_facecolor('white')
        for item in leg.legendHandles:
            item.set_visible(False)

        counter += 1
        #plt.savefig('../usual figures/{}.png'.format(counter))

    #for i in analysed['unusual']:
        #light_curve(extract.read_sources(), analysed['ra'][i], analysed['dec'][i])



def test_chisqrd(sources):

    test_index = 2

    linear_mod = odr.Model(f_flat)
    data = odr.Data(sources[test_index]['dates'], sources[test_index][' int_flux'], we=1./np.power(sources[test_index][' int_flux_err'], 2))
    od = odr.ODR(data, linear_mod, beta0=[1])
    output = od.run()
    print('\ntest chisqrd odr output')
    output.pprint()

    plt.figure()
    plt.errorbar(sources[test_index]['dates'], sources[test_index][' int_flux'], yerr=sources[test_index][' int_flux_err'], fmt='.')
    x = np.linspace(min(sources[2]['dates']), max(sources[2]['dates']), 100)
    y = f_flat(output.beta, x)
    plt.plot(x, y)
    plt.xlabel('mjd')
    plt.ylabel('int_flux')
    plt.title('example source')


if __name__ == '__main__':

    with open('source_locations.json', 'r') as f:
        flocs = json.load(f)

    j2217, j2208 = extract.extract()

    j2217_sources = extract.isolate_sources(j2217, tname=flocs['j2217_template'], ignore='j2217_regionignore.csv')
    j2217_sources_flagged = extract.remove_bad(j2217_sources, multiplier=10, max_fl=30)

    j2208_sources = extract.isolate_sources(j2208, tname=flocs['j2208_template'], ignore='j2208_regionignore.csv')
    j2208_sources_flagged = extract.remove_bad(j2208_sources, multiplier=10, max_fl=30)

    all_sources = j2217_sources_flagged + j2208_sources_flagged

    #chi_sqrd(j2217_sources_flagged)
    # analyse(all_sources)

    # test_chisqrd(j2217_sources_flagged)
    light_curve(all_sources, 330.885, 54.505)

    plt.show()




