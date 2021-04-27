import extract
import matplotlib.pyplot as plt
import numpy as np
import scipy.odr as odr


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


def chi_sqrd(sources):

    etas = []
    variation = []
    for source in sources:
        # looping over sources

        if source.shape[0] < 5:
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
        fweightssum = np.sum(fweights)

        # ---------------------  flux density coefficient of variation ------------------
        variation.append(np.std(source[' int_flux'],
                                ddof=1)/faverage)

        # ------------------------ CHISQRD -----------------------
        eta = 0
        dof = source.shape[0] - 1
        for i in range(source.shape[0]):
            # looping over row in source

            eta += np.power(source[' int_flux'] - faverage, 2)/np.power(source[' int_flux_err'], 2)

        etas.append(eta / dof)

    plt.figure()
    plt.hist(etas, bins=200)
    plt.xlabel(r'Reduced $\eta_\nu$')

    plt.figure()
    plt.hist(variation, bins=40)
    plt.xlabel('flux density coefficient of variation')


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

    j2217, j2208 = extract.extract()
    j2217_sources = extract.isolate_sources(j2217)
    j2217_sources_flagged = extract.remove_bad(j2217_sources, multiplier=10, max_fl=30)

    chi_sqrd(j2217_sources_flagged)
    test_chisqrd(j2217_sources_flagged)

    plt.show()




