import extract
import variability
import matplotlib.pyplot as plt


sources = extract.read_sources()
analysed = variability.analyse(sources)

extract.plot_sources(sources)

variability.plot_eta(analysed)
#variability.plot_variablity(analysed)
#variability.plot_var_vs_eta(analysed)

variability.plot_unusual(analysed)



plt.show()
