import datetime
from astropy.time import Time
import pandas
####
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
####


def createArrays():
    # sucks up all those juicy data points (slurrp)

    # converts csv file into array
    #xerror.append(float(splitUp[2]))
    #yerror.append(float(splitUp[3]))

    #x.append(float(tm.mjd))

    y.append(float(splitUp[1]))





def plotGraph():
    # plots points and labels
    ax.scatter(x, y, alpha=1, color=pointcolour, edgecolors='white', label=Label)
    #plt.errorbar(x, y, yerr=xerror, xerr=yerror, fmt='.k', color='red',
    #ecolor=pointcolour, elinewidth=errLinewidth, capsize=capSize)

    plt.xlabel("MJD")
    plt.ylabel("Flux density (Janskys)")
    plt.title(graphTitle)
    ax.grid(True)


def calcRedX2():
    # calculates reduced chi squared
    chisqrd = 0
    dof = len(x) - 1
    for i in range(0, len(x)):
        E = slope * x[i] + intercept
        O = y[i]
        chisqrd += ((O - E) ** 2) / E

    redchisqrd = chisqrd / dof

    print('Reduced Chi squared:', redchisqrd)


######################################################


# INSERT LABEL STUFF HERE:
xlabel = str('xlabel')
ylabel = str('ylabel')
graphTitle = str('')
Label = str('')
capSize = 1
errLinewidth = 1
pointcolour = str('blue')

fileName = str('test1.csv')

################## Maincode ##########################


x = []
y = []
xerror = []
yerror = []

start = 1
try:
    # opens file to extract fields
    readFile = open(fileName, 'r')
except:
    start = 0

    # error if file not availible
    print('Error - File not found')

if start == 1:
    print('File Found')
    # reads each line, extracts fields and creates separate arrays

    data = pandas.read_csv(fileName)
    #print(data)
    dates = [Time(datetime.datetime.strptime(i, '%d-%b-%y')).mjd for i in data['date']]
    #print(dates)

    j=0
    for line in readFile:
        if j == 0:
            j+=1
            continue
        #print(line)
        # defines how to split each column
        splitUp = line.split(',')
        createArrays()
        x = dates
        #print(x)

    fig, ax = plt.subplots()
    # plots each point

    # prints parameters for fit

    print(x)
    print(y)

    try:
        model = np.polyfit(x, y, 1)
        # print(model)
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)

        print('Gradient =', slope)
        print('Intercept = ', intercept)

        plotGraph()
        plt.legend()
        plt.show()
        print('Graph plotted successfully\n')
    except:
        print('ERROR - Failed to plot graph')

    # calculates reduced chi squared
    calcRedX2()

#print(x)
# print(y)
# print(xerror)

#dt = datetime.strptime(i, '%d-%m-%y'))
#tm = Time(dt)
#print(tm.mjd)

