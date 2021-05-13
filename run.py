import extract
import variability
import matplotlib.pyplot as plt
import json


sources = extract.read_sources()
with open('source_locations.json', 'r') as f:
    flocs = json.load(f)

j2217, j2208 = extract.extract()
j17sour = extract.isolate_sources(j2217,
                          tname=flocs['j2217_template'], ignore='j2217_regionignore.csv')
j17_flagged = extract.remove_bad(j17sour, multiplier=4, max_fl=30)

j08_sour = extract.isolate_sources(j2208, tname=flocs['j2208_template'], ignore='j2208_regionignore.csv')
j08_flag = extract.remove_bad(j08_sour, multiplier=4, max_fl=30)
# average_sources(j08_flag, write=True)

all_sources = j17_flagged + j08_flag
#extract.plot_sources(j17_flagged)
#extract.plot_sources(j08_flag)


analysed = variability.analyse(all_sources)

# extract.plot_sources(sources)

variability.plot_eta(analysed)
variability.plot_variablity(analysed)
variability.plot_var_vs_eta(analysed)

#variability.plot_unusual(all_sources, analysed)

#variability.light_curve(sources, 331.9553164560985, 55.58776007997343) # flagged double source

#j08_ana = variability.analyse(j08_flag)
#variability.plot_eta(j08_ana)

#j17_ana = variability.analyse(j17_flagged)
#variability.plot_eta(j17_ana)

plt.show()
