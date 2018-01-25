#!/usr/bin/env python

import os
import sys
import numpy as np
from collections import OrderedDict
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
#import click
import ReadKonusWindRedshiftCatalogueInfo
import ReadLATCatalogueInfo
import ReadGBMCatalogueInfo


##### Matplotlib #####
#mpl.rcParams['text.usetex'] = True
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
plt.rcParams["font.size"] = 18
NMARKER_STYLE = 10


def read_catalogue():
    tb_lat = ReadLATCatalogueInfo.open_table()
    tb_gbm = ReadGBMCatalogueInfo.open_table()
    tb_kw = ReadKonusWindRedshiftCatalogueInfo.open_table()
    grbs_lat = ReadLATCatalogueInfo.read_all(tb_lat, tb_gbm, tb_kw)
    grbs_lat_gbm = ReadLATCatalogueInfo.select_gbm_exist(grbs_lat)
    grbs_lat_gbm_kw = ReadLATCatalogueInfo.select_konuswind_exist(grbs_lat_gbm)

    dict_all = {'efluence':{'GBM':OrderedDict(),
                            'KonusWind':OrderedDict()},
                'zenith':{'Fermi':OrderedDict()},
                'epeak':{'GBM':OrderedDict(),
                         'KonusWind':OrderedDict()},
                'alpha':{'GBM':OrderedDict(),
                         'KonusWind':OrderedDict()},
                'beta':{'GBM':OrderedDict(),
                         'KonusWind':OrderedDict()}}

    for grb in grbs_lat_gbm_kw:
        #print grb
        dict_all['efluence']['GBM'][grb['GRBNAME']] = grb['GBM']['FLUENCE']
        dict_all['efluence']['KonusWind'][grb['GRBNAME']] = grb['KonusWind']['S']*1e-6
        dict_all['zenith']['Fermi'][grb['GRBNAME']] = grb['ZENITH']
        dict_all['epeak']['GBM'][grb['GRBNAME']] = grb['GBM']['FLNC_BAND_EPEAK']
        dict_all['epeak']['KonusWind'][grb['GRBNAME']] = float(grb['KonusWind']['Ep'])
        dict_all['alpha']['GBM'][grb['GRBNAME']] = grb['GBM']['FLNC_BAND_ALPHA']
        dict_all['alpha']['KonusWind'][grb['GRBNAME']] = float(grb['KonusWind']['Alpha'])
        dict_all['beta']['GBM'][grb['GRBNAME']] = grb['GBM']['FLNC_BAND_BETA']
        dict_all['beta']['KonusWind'][grb['GRBNAME']] = float(grb['KonusWind']['Beta'])
    return dict_all
            

def plot(dict_info, ax):
    ax[0][0].scatter(dict_info['efluence']['GBM'].values(),dict_info['efluence']['KonusWind'].values(), marker='o')
    ax[0][0].set_xlim([1e-6, 0.01])
    ax[0][0].set_ylim([1e-6, 0.01])
    ax[0][0].set_xscale("log", nonposx='clip')
    ax[0][0].set_yscale("log", nonposy='clip')
    ax[0][0].grid()
    ax[0][0].set_xlabel('GBM [erg/cm^2]')
    ax[0][0].set_ylabel('Konus-Wind [erg/cm^2]')

    ax[0][1].scatter(np.array(dict_info['zenith']['Fermi'].values()), np.array(dict_info['efluence']['GBM'].values())/np.array(dict_info['efluence']['KonusWind'].values()), marker='o')
    #ax[0][1].set_ylim([1e-6, 0.01])
    ax[0][1].set_yscale("log", nonposy='clip')
    ax[0][1].grid()
    ax[0][1].set_xlabel('Zenith for Fermi [deg]')
    ax[0][1].set_ylabel('GBM/Konus-Wind [erg/cm^2]')

    ax[0][2].scatter(np.array(dict_info['efluence']['KonusWind'].values()), np.array(dict_info['efluence']['GBM'].values())/np.array(dict_info['efluence']['KonusWind'].values()), marker='o')
    ax[0][2].set_xlim([1e-6, 0.01])
    ax[0][2].set_xscale("log", nonposy='clip')
    ax[0][2].set_yscale("log", nonposy='clip')
    ax[0][2].grid()
    ax[0][2].set_xlabel('Konus-Wind [erg/cm^2]')
    ax[0][2].set_ylabel('GBM/Konus-Wind')

    ax[1][0].scatter(np.array(dict_info['epeak']['KonusWind'].values()), np.array(dict_info['efluence']['GBM'].values())/np.array(dict_info['efluence']['KonusWind'].values()), marker='o')
    #ax[1][0].set_xlim([1e-6, 0.01])
    ax[1][0].set_xscale("log", nonposy='clip')
    ax[1][0].set_yscale("log", nonposy='clip')
    ax[1][0].grid()
    ax[1][0].set_xlabel('Konus-Wind peak energy [keV]')
    ax[1][0].set_ylabel('GBM/Konus-Wind')

    ax[1][1].scatter(np.array(dict_info['epeak']['KonusWind'].values())-np.array(dict_info['epeak']['GBM'].values()), np.array(dict_info['efluence']['GBM'].values())/np.array(dict_info['efluence']['KonusWind'].values()), marker='o')
    #ax[1][1].set_xlim([1e-6, 0.01])
    #ax[1][1].set_xscale("log", nonposy='clip')
    ax[1][1].set_yscale("log", nonposy='clip')
    ax[1][1].grid()
    ax[1][1].set_xlabel('Difference in Epeak [keV]')
    ax[1][1].set_ylabel('GBM/Konus-Wind')

    ax[1][2].scatter(np.array(dict_info['beta']['KonusWind'].values()), np.array(dict_info['efluence']['GBM'].values())/np.array(dict_info['efluence']['KonusWind'].values()), marker='o')
    #ax[1][2].set_xlim([1e-6, 0.01])
    #ax[1][2].set_xscale("log", nonposy='clip')
    ax[1][2].set_yscale("log", nonposy='clip')
    ax[1][2].grid()
    ax[1][2].set_xlabel('Konus-Wind beta')
    ax[1][2].set_ylabel('GBM/Konus-Wind')


def savefig(fig, outdir, outname, figforms):
    for ff in figforms:
        outpath = '{0}/{1}.{2}'.format(outdir, outname, ff)
        fig.savefig(outpath)


###### Main #####
def main():
    fig, ax = plt.subplots(2,3, figsize=(30, 20))
    #fig = plt.figure(figsize=(5,5))
    #ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
    tbdata = read_catalogue()
    plot(tbdata, ax)
    fig.tight_layout()
    savefig(fig, '.', 'FIG_GBM_vs_KonusWind', ('png', 'pdf'))


if __name__ == '__main__':
    main()
