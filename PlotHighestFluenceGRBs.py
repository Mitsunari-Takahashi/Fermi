#!/usr/bin/env python

import sys
import os
import copy
#import ROOT
#from ROOT import TTree, TChain, TGraph, TH1
import numpy as np
#from array import array
#import math
#from math import pi, cos, sin, tan, acos, asin, atan, radians, degrees
#from pColor import *
import click
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from pLsList import ls_list


def PlotHighestFluenceGRBs(lst_path_input, str_suffix):
    NGRB = len(lst_path_input)
    if str_suffix!='':
        str_suffix = '_' + str_suffix
    NAX = 2
    NAY = 1
    fig, ax = plt.subplots(NAX, NAY, figsize=(10, 10))

    x = [[], []]
    y = [[], []]
    x0 = [0.5]
    y0 = [[], []]
    y_err_lo = [[], []]
    y_err_hi = [[], []]
    ax_title = ['Flux from 56.2 GeV to 316 GeV in the rest frame', 'Energy flux from 56.2 GeV to 316 GeV in the rest frame']
    ylabel = [r'$(F - F_{extrapolated})/\Delta F_{extrapolated}$', r'$(eF - eF_{extrapolated})/\Delta eF_{extrapolated}$']
    for (iin, pathin) in enumerate(lst_path_input):
        x[0].append(iin+1)
        x0.append(iin+1.5)
        x[1].append(os.path.basename(pathin)[3:11])
        df = pd.read_csv(pathin, names=('start', 'stop', 'emin_rest', 'emax_rest', 'emin_shifted', 'emax_shifted', 'ts', 'Integral', 'Integral_err', 'Integral_ul95', 'Integral_ll95', 'Integral_ul68', 'Integral_ll68', 'Index', 'Index_err', 'flux', 'flux_err', 'flux_ul95', 'flux_ll95', 'flux_ul68', 'flux_ll68', 'eflux', 'eflux_err', 'eflux_ul95', 'eflux_ll95', 'eflux_ul68', 'eflux_ll68', 'fluxhe', 'fluxhe_err', 'efluxhe', 'efluxhe_err'), comment='#') 
        fluxhe_extrapolated = df.ix[1, 'fluxhe']
        fluxhe_extrapolated_err = df.ix[1, 'fluxhe_err']
        efluxhe_extrapolated = df.ix[1, 'efluxhe']
        efluxhe_extrapolated_err = df.ix[1, 'efluxhe_err']
        fluxhe = df.ix[2, 'flux']
        fluxhe_sigma = (fluxhe - fluxhe_extrapolated) / fluxhe_extrapolated_err
        fluxhe0_sigma = (0 - fluxhe_extrapolated) / fluxhe_extrapolated_err
        fluxhe_ul95 = df.ix[2, 'flux_ul95']
        fluxhe_ul95_sigma = (fluxhe_ul95 - fluxhe_extrapolated) / fluxhe_extrapolated_err
        fluxhe_ll95 = df.ix[2, 'flux_ll95']
        fluxhe_ll95_sigma = (fluxhe_ll95 - fluxhe_extrapolated) / fluxhe_extrapolated_err
        y[0].append(fluxhe_sigma)
        y0[0].append(fluxhe0_sigma)
        if iin==0:
            y0[0].append(fluxhe0_sigma)
        #print fluxhe0_sigma
        y_err_lo[0].append(fluxhe_sigma-fluxhe_ll95_sigma)
        y_err_hi[0].append(fluxhe_ul95_sigma-fluxhe_sigma)
        efluxhe = df.ix[2, 'eflux']
        efluxhe_sigma = (efluxhe - efluxhe_extrapolated) / efluxhe_extrapolated_err
        efluxhe0_sigma = (0 - efluxhe_extrapolated) / efluxhe_extrapolated_err
        efluxhe_ul95 = df.ix[2, 'eflux_ul95']
        efluxhe_ll95 = df.ix[2, 'eflux_ll95']
        efluxhe_ul95_sigma = (efluxhe_ul95 - efluxhe_extrapolated) / efluxhe_extrapolated_err
        efluxhe_ll95_sigma = (efluxhe_ll95 - efluxhe_extrapolated) / efluxhe_extrapolated_err
        y[1].append(efluxhe_sigma)
        y0[1].append(efluxhe0_sigma)
        if iin==0:
            y0[1].append(efluxhe0_sigma)

        y_err_lo[1].append(efluxhe_sigma-efluxhe_ll95_sigma)
        y_err_hi[1].append(efluxhe_ul95_sigma-efluxhe_sigma)
    #x[0].append(NGRB+1)
    #x[1].append('')
#    for iax in range(NAX):
    x_ticks = copy.deepcopy(x)
    x_ticks[0].insert(0, 0)
    x_ticks[1].insert(0, '')
    x_ticks[0].insert(-1, NGRB+1)
    x_ticks[1].insert(-1, '')
    x_consistent = [0.5, NGRB+0.5, NGRB+0.5, 0.5] #np.array([0, NGRB+1])
    y_consistent = [-1, -1, 1, 1]
    y_consistent_lo = np.array([-1, -1])
    y_consistent_hi = np.array([1, 1])
    for iax in range(NAX):
        ax[iax].step(x0, y0[iax], color='black', ls='dashed', label='F=0')
        #ax[iax].fill_between(x_consistent, y_consistent_lo, y_consistent_hi, facecolor='y',alpha=0.5, label='Region consistent with PL')
        ax[iax].fill(x_consistent, y_consistent, facecolor='g',alpha=0.25, label=r'Consistent with PL ($\pm 1 \sigma$)', lw=0)
        ax[iax].errorbar(x[0], y[iax], yerr=[y_err_lo[iax], y_err_hi[iax]], fmt='o', label='Observed flux and \n 95% upper/lower limits', lw=1)
        ax[iax].set_title(ax_title[iax])
        ax[iax].set_xticks(x_ticks[0])
        ax[iax].set_xticklabels(x_ticks[1], rotation=15, fontsize='small')
        ax[iax].set_xlim([0, NGRB+1])
        ax[iax].set_ylim([-3, 12])
        ax[iax].set_ylabel(ylabel[iax])
        ax[iax].grid(axis='y')
        ax[iax].legend(loc=2, fontsize=8)
    fig.savefig('VHEflux_in_sigma{0}.png'.format(str_suffix))


@click.command()
@click.argument('inputs', type=str)
@click.option('--suffix', '-s', type=str, default='')
def main(inputs, suffix):
    lst_path_incsv = ls_list(inputs)
    PlotHighestFluenceGRBs(lst_path_incsv, suffix)


if __name__ == '__main__':
    main()
