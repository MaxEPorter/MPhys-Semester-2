import extract
import matplotlib.pyplot as plt
import numpy as np
import scipy.odr as odr
import json

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

    etas = []
    variation = []
    faverages = []
    for source in sources:
        # looping over sources

        if source.shape[0] < 4:
            continue

        # dates is x,  int_flux is y
        """
        linear_mod = odr.Model(f_flat)
        data = odr.Data(source['dates'], source[' int_flux'], we=1./np.power(source[' int_flux_err'], 2))
        od = odr.ODR(data, linear_mod, beta0=[0.5])
        output = od.run()
        # output.pprint()
        """

        fweights = 1./np.power(source[' int_flux_err'], 2)
        faverage = np.average(source[' int_flux'], weights=fweights)
        faverages.append(faverage)

        # print('odr av -> {}, faverage -> {}'.format(output.beta, faverage))

        # ---------------------  flux density coefficient of variation ------------------
        variation.append(np.std(source[' int_flux'],
                                ddof=1)/faverage)

        # ------------------------ CHISQRD -----------------------
        eta = 0
        dof = source.shape[0] - 1
        for i in range(source.shape[0]):
            # looping over row in source

            eta += np.power(source[' int_flux'][i] - faverage, 2)/np.power(source[' int_flux_err'][i], 2)

        red_eta = eta/dof
        etas.append(eta / dof)

    print(etas)
    plt.figure()
    plt.hist(etas, bins=10**np.linspace(0, 5, 40))
    plt.xscale('log')
    plt.xlabel(r'$\eta_\nu$')

    plt.figure()
    plt.hist(variation, bins=40)
    plt.xlabel('flux density coefficient of variation')

    plt.figure()
    plt.scatter(faverages, variation)
    plt.xlabel('average flux')
    plt.ylabel('variation')

    plt.figure()
    plt.scatter(variation, etas)
    plt.xlabel('variation')
    plt.ylabel(r'$\eta_\nu$')
    plt.yscale('log')


def light_curve(sources, ra, dec):

    s = extract.source_by_pos(sources, ra, dec)
    fweights = 1. / np.power(s[' int_flux_err'], 2)
    faverage = np.average(s[' int_flux'], weights=fweights)

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.errorbar(s['dates'], s[' int_flux'], yerr=s[' int_flux_err'], fmt='.', color='black')
    ax.plot([min(s['dates']), max(s['dates'])], [faverage, faverage], color='mediumseagreen')
    ax.set_ylabel('Flux density (Jy)')
    ax.set_xlabel('time (MJD)')


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




