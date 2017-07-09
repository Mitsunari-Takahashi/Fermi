#!/usr/bin/env python
import pandas
import pMETandMJD


def ReadGBMCatalogueOneLine(nameGrb, path_catlogue_csv = "/disk/gamma/cta/store/takhsm/FermiData/catalogue/GBM-BusrtCatalogue_20170623.csv"):
    """Read one GRB's information from GBM catalogue in CSV format and return it in dictonary.
"""
    dct_info = {}
    num_lines = sum(1 for line in open(path_catlogue_csv))
    csv = pandas.read_csv(path_catlogue_csv)
    for iGrb in range(num_lines-1):
        if int(nameGrb) == int(csv.ix[iGrb,'name']):
            dct_info['ra'] = float(csv.ix[iGrb,'ra'])
            dct_info['dec'] = float(csv.ix[iGrb,'dec']) 
            dct_info['trigger_time'] = pMETandMJD.ConvertMjdToMet(float(csv.ix[iGrb,'trigger_time']))
            dct_info['error_radius'] = float(csv.ix[iGrb,'error_radius'])
    return dct_info
