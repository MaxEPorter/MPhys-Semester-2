import pandas
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog


def get_sources(fname):

    df = pandas.read_excel(fname, engine='openpyxl')
    t1_end = df.shape[0]
    t2_start = 0

    for index, row in df.iterrows():

        if pandas.isna(row['NVSS id']):
            t1_end = index
            break

    for index, row in df.iterrows():

        if row['NVSS id'] == 'ERRORS':
            t2_start = index+2
            break

    f = df.iloc[:t1_end]
    e = df.iloc[t2_start:]

    return f, e


def get_fname():
    root = tk.Tk()
    root.withdraw()

    # get spreadsheets and find table
    file_path = filedialog.askopenfilename()
    return file_path


if __name__ == '__main__':

    fluxes, errors = get_sources(get_fname())

    print(fluxes)
    print(errors)

