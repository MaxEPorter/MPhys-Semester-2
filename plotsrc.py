import pandas
import numpy as np
import matplotlib.pyplot as plt

df = pandas.read_excel('../J2217.xlsx', engine='openpyxl')
print(df)

vcols = df.shape[1]-2

dfnp = df.to_numpy()

values = dfnp[0][-vcols:]
print(values)

plt.hist(values)
plt.show()
