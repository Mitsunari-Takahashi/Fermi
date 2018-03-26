#!/usr/bin/env python
"""Module for making light curves of multi-wavelength data.
The main class LightCurve is a chain of another module pLATLikelihoodConfig.py.
The authour: Mitsunari Takahashi
"""
import sys
import os
import os.path
path_upstairs = os.path.join(os.path.dirname(__file__), '../')
sys.path.append(path_upstairs)
import logging
import pickle
import datetime
import numpy as np
import math
from math import log10, log, sqrt, ceil, isnan, pi, factorial
from astropy.io import fits
import click
import csv
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import pickle_utilities
import pMatplot
import pMETandMJD

# mpl.rcParams['text.usetex'] = True
# mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
# plt.rcParams["font.size"] = 18
# pgf_with_rc_fonts = {"pgf.texsystem": "pdflatex"}
# mpl.rcParams.update(pgf_with_rc_fonts)
NMARKER_STYLE = 10

##### VERSION OF THIS MACRO #####
VERSION = 0.1


##### Conversion from MeV to erg ######
MEVtoERG = 1.6021766208E-6

##### Refered information #####
GRB_CATALOGUE = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LAT2CATALOG-v1-LTF.fits'


def convert_eflux_to_vFv(eflux, spindex, eref, emin, emax):
    if spindex!=-2:
        return eflux * (spindex+2)*eref**(spindex+2)/(emax**(spindex+2)-emin**(spindex+2))
    else:
        return eflux / (log(emax)-log(emin))


def predict_LAT_curve(eref, erange, fluenceGBM, tref, tmin, tmax, spindex):
    if erange=='highE':
        prefactor = 5.07E-04
        emax = 100000.
        emin = 10000.
        lcindex = -1.095
    elif erange=='midE':
        prefactor = 15.4881661891
        emax = 10000.
        emin = 1000.
        lcindex = -1.315
    if erange=='lowE':
        prefactor = 138.03842646
        emax = 1000.
        emin = 100.
        lcindex = -1.305
    return prefactor*fluenceGBM * pow(tref/10., lcindex)
#    return n_t95*fluxGBM*(spindex+1.)/(pow(emax, spindex+1) - pow(emin, spindex+1))*pow(eref, spindex+1) * pow(eref*MEVtoERG, 1) * pow(tref/tnorm, lcindex)


class SourceObject:
    def __init__(self, name, suffix):
        self.name = name
        self.datasets = {}
        self.suffix = '' if suffix=='' else '_'+suffix
        

    def load_csv(self, path_csv):
        with open(path_csv, 'r') as f:
            reader = csv.reader(f)
            header = next(reader)
            for row in reader:
                print row
                if not row[0] in self.datasets:
                    self.datasets[row[0]] = {}
                    for icol, col in enumerate(header[1:]):
                        self.datasets[row[0]][col] = []
                for icol, col in enumerate(header[1:]):
                    if row[0][:3]!='LAT' and col[:5]=='nuFnu':
                        self.datasets[row[0]][col].append(convert_eflux_to_vFv(float(row[icol+1]), float(row[9]), float(row[1]), float(row[2]), float(row[3])))
                    else:
                        self.datasets[row[0]][col].append(float(row[icol+1]))


    def load_pickle(self, path_pickle, nindex, lcindex):
        slots = pickle_utilities.load(path_pickle)
        for islot, slot in enumerate(slots):
            print '===== {0} - {1} s ====='.format(slot['period'][0], slot['period'][1])
            print 'Normalization: {0}'.format(slot['normalization'][nindex][lcindex])
            print 'E^2 dN/dE: {0}'.format(slot['normalization'][nindex][lcindex]*4.4723E-6*10000*MEVtoERG)
            print 'Energy flux: {0}'.format(slot['eflux'][nindex][lcindex]*MEVtoERG)


    def plot(self, pathout, figform):
        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_axes((0.07, 0.15, 0.9, 0.75))
        ax.set_title('GRB {name}'.format(name=self.name))
        ax.set_xscale("log", nonposx='clip')
        ax.set_yscale("log", nonposy='clip')
        #ax.set_xlabel('$T - T0 [s]}$')
        #ax.set_ylabel(r'$\nu F_{\nu} \rm{[erg/cm^2 \cdot s]}$')
        for keydata, dataset in self.datasets.items():
            ax.errorbar(dataset['Time'], dataset['nuFnu'], xerr=dataset['TimeErr'], yerr=(dataset['nuFnuErrNeg'], dataset['nuFnuErrPos']), label=keydata, fmt='x', lw=1) #marker=lst_markers[nmarker])
            # ax.set_xlim([max(0.1,tmin), tmax])
            # ax.set_ylim([emin/1000., emax/1000.])
        ax.grid()#axis='both', which='major',color='black',linestyle='--')
        #ax.grid(axis='y', which='minor',color='black',linestyle='-.', alpha=0.25)
        #ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
        #ax.yaxis.set_minor_formatter(FormatStrFormatter('%.0f'))

        #Prediction in 10-100GeV
        t_pred = [7.168, 100000]
        vFv_pred_highE = [predict_LAT_curve(10000., 'highE', fluenceGBM=4.4723e-06, tref=t, tmin=t_pred[0], tmax=t_pred[1], spindex=-2.1) for t in t_pred]
        ax.plot(t_pred, vFv_pred_highE, lw=1, c='b')
        #vFv_pred_lowE = [predict_LAT_curve(500., 'lowE', fluenceGBM=4.4723e-06, tref=t, tmin=t_pred[0], tmax=t_pred[1], spindex=-2.0) for t in t_pred]
        #ax.plot(t_pred, vFv_pred_lowE, lw=1)
        
        ax.legend(loc=0, fancybox=True, framealpha=0.5, fontsize=9) #, ncol=ngrb_plotted/10+1)
        for ff in figform:
            fig.savefig('{0}{1}.{2}'.format(pathout, self.suffix, ff))


@click.command()
@click.argument('name', type=str)
@click.argument('csv', type=str)
# @click.option('--emin', type=float, default=10000.)
# @click.option('--emax', type=float, default=10**5.5)
# @click.option('--tmin', type=float, default=0)
# @click.option('--tmax', type=float, default=100000)
# @click.option('--roi', type=float, default=1.)
# @click.option('--redshift', '-z', is_flag=True)
@click.option('--suffix', '-s', type=str, default='')
@click.option('--pathout', type=str, default='./MWLLightCurve')
@click.option('--figform', type=str, default=('png',), multiple=True)
def main(name, csv, suffix, pathout, figform):
    src = SourceObject(name, suffix)
    src.load_csv(csv)
    src.load_pickle('/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/140928437/E0010000-0100000MeV/r01deg/briefslots/LightCurve_140928437_ScaleFactor::PowerLaw2_IndexFree_400x180_count.pickle', 111, 79)
    src.plot(pathout, figform)
    

if __name__ == '__main__':
    main()
