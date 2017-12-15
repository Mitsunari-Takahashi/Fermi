#!/usr/bin/env python

import os
import sys
import pandas as pd


PATH_CATALOGUE = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/Swift_grb_table.csv'


def open_table(path=PATH_CATALOGUE):
    """Open your csv file and return its pandas dataframe.
"""
    if path in (None, ""):
        path=PATH_CATALOGUE
    dset = pd.read_csv(path, header=0)
    return dset


def read_one_row(dset, gcnname):
    """Returns a data series.
"""
    #print dset[dset.GRB==gcnname]
    #print dset
    #print [dset.GRB==gcnname]
    return dset[dset.GRB==gcnname].iloc[0]
