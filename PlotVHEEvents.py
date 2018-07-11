#!/usr/bin/env python
"""Module for making light curves of LAT data.
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
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from FindGoodstatPeriods import find_goodstat_periods, get_entries, get_event_time_and_energy, get_event_time_energy_angsep
from FindCrossEarthlimb import find_cross_earthlimb
#import ReadLTFCatalogueInfo
import ReadLATCatalogueInfo
import ReadGBMCatalogueInfo
import pickle_utilities
import pMatplot
import pMETandMJD
import DownloadFermiData
from pMatplot import TPL_COLOR_WO_BLACK as TPL_COLOR

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
plt.rcParams["font.size"] = 15
NMARKER_STYLE = 10

##### VERSION OF THIS MACRO #####
VERSION = 0.1


##### Conversion from MeV to erg ######
MEVtoERG = 1.6021766208E-6

##### Refered information #####
GRB_CATALOGUE = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LAT2CATALOG-v1-LTF.fits'

##### Plot style #####
lst_markers = ['s', 'o', 'D', 'x']


##### CalOnly photons #####
DICT_ADDPHOTON={'090926A': {'time': (422.7,), 'energy':(50.49609375*1000.,)},
                '150902A': {'time': (2064.52053469,), 'energy':(83.6322578125*1000.,)},
                '160509A': {'time': (2035.85387415,5757.82151717), 'energy':(115.829539063*1000.,63.1624726562*1000.)}
                }


class SourceObject:
    def __init__(self, name, met0, ra, dec, gcnname=None, tmin=0, tmax=100000, emin=10000, emax=316228, deg_roi=1., zmax=100., suffix='', tprompt_start=0., redshift=0):
        self.config = {'name':name, 
                       'gcnname':gcnname if gcnname is not None else name,
                       'coord':{'ra':ra, 'dec':dec},
                       'time':{'min':tmin, 'max':tmax}, 
                       'met':{'o':met0, 'min':met0+tmin, 'max':met0+tmax}, 
                       'energy':{'min':emin, 'max':emax},
                       'roi':{'radius':deg_roi}, #, 'margin':rad_margin}, 
                       'zenith':{'max':zmax}, 
                       'suffix':suffix,
                       'path_evtfile':'/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{na}/GRB{na}_P8_P302_BASE_T00-999-101000_r030_ft1.fits'.format(na=name)
                       }
        self.gbm = {'t90_start':tprompt_start}
        self.target = {'coord':{'ra':ra, 'dec':dec},
                       'redshift':redshift}


    def get_events(self):
        events = get_event_time_and_energy(self.config['path_evtfile'], tStart=min(self.config['time']['min'], self.config['time']['min'] + self.gbm['t90_start']), tStop=self.config['time']['max'], rlim=self.config['roi']['radius'], ra=self.target['coord']['ra'], dec=self.target['coord']['dec'], torigin=self.config['met']['o'], emin=self.config['energy']['min'], emax=self.config['energy']['max'], zmax=self.config['zenith']['max'], z=self.target['redshift'])
        #if self.target['redshift']>0:
        #    self.energies = self.energies * (1.+self.target['redshift'])
        self.times = events[0]
        self.energies = events[1]

    def get_events_angsep(self):
        events = get_event_time_energy_angsep(self.config['path_evtfile'], tStart=min(self.config['time']['min'], self.config['time']['min'] + self.gbm['t90_start']), tStop=self.config['time']['max'], rlim=sqrt(10.), ra=self.target['coord']['ra'], dec=self.target['coord']['dec'], torigin=self.config['met']['o'], emin=self.config['energy']['min'], emax=self.config['energy']['max'], zmax=self.config['zenith']['max'], z=self.target['redshift'])
        self.angsep = events[2]


    def plot(self, ax, marker, color=None, alpha=1, addphoton=False):
        #nmarker = nplot / NMARKER_STYLE
        p = ax.scatter(self.times, self.energies/1000., alpha=alpha, label=self.config['gcnname'], marker=marker, facecolors='none', edgecolors=color)
        #p = ax.scatter(self.times, self.energies/1000., alpha=0.75, label=self.config['gcnname'], marker=lst_markers[nmarker])
        if addphoton==True and self.config['gcnname'] in DICT_ADDPHOTON.keys():
            print 'Additional photons:'
            t_added = np.array(DICT_ADDPHOTON[self.config['gcnname']]['time'])
            e_added = np.array(DICT_ADDPHOTON[self.config['gcnname']]['energy'])
            for t,e in zip(t_added, e_added):
                print '{e:.1f} MeV at {t:.1f}'.format(t=t, e=e)
            ax.scatter(t_added, e_added/1000., marker=marker, label='CalOnly', facecolors=p.get_edgecolors(), edgecolors=p.get_edgecolors(), s=p.get_sizes()[0]*1.5)


@click.command()
@click.option('--emin', type=float, default=10000.)
@click.option('--emax', type=float, default=10**5.5)
@click.option('--tmin', type=float, default=0)
@click.option('--tmax', type=float, default=100000)
@click.option('--roi', type=float, default=1.)
@click.option('--redshift', '-z', is_flag=True)
@click.option('--suffix', '-s', type=str, default='')
@click.option('--pathout', type=str, default='./Stacked_VHEevents')
@click.option('--figform', type=str, default=('png',), multiple=True)
@click.option('--longonly', '-l', is_flag=True)
@click.option('--exceptions', '-e', type=str, default=(), multiple=True)
@click.option('--tolerr', type=float, default=180)
@click.option('--addphoton', '-a', is_flag=True)
def main(emin, emax, tmin, tmax, roi, redshift, suffix, pathout, figform, longonly, exceptions, tolerr, addphoton):
    fig = plt.figure(figsize=(15,4.5))
    ax = fig.add_axes((0.07, 0.15, 0.9, 0.75))
    str_title = 'Events within {rad:2.1f} deg from LAT-catalogued {lo}GRBs'.format(rad=roi, lo='long-' if longonly==True else '')
    if tolerr<180:
        str_title = str_title +  ' with localization error smaller than {te} deg'.format(te=tolerr)
    ax.set_title(str_title)
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    ax.set_xlabel(r'$Time - \rm{ min( T_{0}, T_{0}+T_{05}) [s]}$')
    ax.set_ylabel(r'$Energy \rm{[GeV]}$')

    # fig_sep = plt.figure(figsize=(15, 4.5))
    # ax_sep = fig_sep.add_axes((0.07, 0.15, 0.9, 0.75))
    # angsep_stacked_sq = []

    tb_gbm = ReadGBMCatalogueInfo.open_table()
    tb_lat = ReadLATCatalogueInfo.open_table()
    lst_lat = ReadLATCatalogueInfo.read_all(tb_lat, tb_gbm)
    #lst_lat = ReadLATCatalogueInfo.select_gbm_exist(lst_lat)
    if longonly==True:
        lst_lat = ReadLATCatalogueInfo.select_long(lst_lat)
    lst_lat = ReadLATCatalogueInfo.select_small_error(lst_lat, tolerr)
    ngrb_plotted = 0
    for (irow , tb_one) in enumerate(lst_lat):
        #print tb_one
        name = tb_one['GRBNAME']
        print '##### No.{0} GRB{1} #####'.format(irow, name)
        if name in exceptions:
            print '{0} is skipped.'.format(name)
            continue
        #tb_one = ReadLTFCatalogueInfo.select_one_by_name(tb_lat, name)
        trigger_met = tb_one['TRIGGER_TIME'] #pMETandMJD.ConvertMjdToMet(tb_one['GBM']['TRIGGER_TIME'])
        if tb_one['GBM'] is not None:
            t05 = pMETandMJD.ConvertMjdToMet(tb_one['GBM']['T90_START'])
        else:
            t05 = 0.
        z = 0 #tb_one['REDSHIFT']
        if redshift==False or z>0:
            src = SourceObject(name=name, met0=trigger_met, ra=tb_one['RA'], dec=tb_one['DEC'], gcnname=tb_one['GCNNAME'], tmin=tmin, tmax=tmax, emin=emin, emax=emax, deg_roi=roi, tprompt_start=t05, redshift=z if redshift==True else 0)
            path_dir_data = '/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{na}'.format(na=name)
            if not os.path.exists(path_dir_data):
                os.mkdir(path_dir_data)
            src.get_events()
            #src.get_events_angsep()
            #angsep_stacked_sq.append(src.angsep**2)

            if len(src.times)>0:
                src.plot(ax, lst_markers[ngrb_plotted/NMARKER_STYLE], TPL_COLOR[ngrb_plotted%len(TPL_COLOR)], addphoton=addphoton)
                ngrb_plotted+=1
    ax.set_xlim([max(0.1,tmin), tmax])
    ax.set_ylim([emin/1000., emax/1000.])
    ax.grid(axis='both', which='major',color='black',linestyle='--')
    ax.grid(axis='y', which='minor',color='black',linestyle='-.', alpha=0.25)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
    ax.legend(loc=2, fancybox=True, framealpha=0.5, fontsize=10, ncol=ngrb_plotted/7+1, scatterpoints=2)

    #ax_sep.hist(angsep_stacked_sq, histtype='barstacked', bins=1000)

    for ff in figform:
        fig.savefig('{0}{1}{2}.{3}'.format(pathout, suffix if suffix=='' else '_'+suffix, '_redshift' if redshift==True else '', ff))
        #fig_sep.savefig('{0}_separation{1}{2}.{3}'.format(pathout, suffix if suffix=='' else '_'+suffix, '_redshift' if redshift==True else '', ff))




if __name__ == '__main__':
    main()
