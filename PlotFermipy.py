#!/usr/bin/env python

import sys
import os
import numpy as np
#from fermipy.gtanalysis import GTAnalysis
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi, log10, sqrt, ceil, isnan
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
ROOT.gROOT.SetBatch()
from pColor import *
from pMETandMJD import *


MEV2ERG = 1.6021766208e-6


class Datum:
    def __init__(self, time=None, enr_cfg=None, enr_sed=None, pl_index=None, flux=None, eflux=None, model_dnde=None):
        self.time = time # min, max
        self.enr_cfg = enr_cfg # min, max
        self.enr_sed = enr_sed # ref, min, max
        self.pl_index = pl_index # value, error
        self.flux =  flux # value, lower error, higher error
        self.eflux = eflux # value, lower error, higher error
        self.model_dnde = model_dnde
        

    def get_model_e2dnde_lowest_energy():
        NLOWEST = 0
        elowest = self.model_dnde['energies'][NLOWEST]
        e2dnde = self.model_dnde['dnde'][NLOWEST] * elowest * elowest
        e2dnde_lo = self.model_dnde['dnde_lo'][NLOWEST] * elowest * elowest
        e2dnde_hi = self.model_dnde['dnde_hi'][NLOWEST] * elowest * elowest
        return [elowest, e2dnde, e2dnde_lo, e2dnde_hi]


def PlotFermipy(name_src, lst_np_input_pairs, t0=0): # pair of output.npy and sed.npy
    # LC
    lst_enr = []
    lst_time = []
    lst_data = []

    for np_input_pair in lst_np_input_pairs:

        npi = np.load(np_input_pair[0]).flat[0]
        np_cfg = npi['config']
        np_src = npi['sources'][name_src]
        np_sed = np.load(np_input_pair[1]).flat[0]
        lst_time.append([np_cfg['selection']['tmin']-t0, np_cfg['selection']['tmax']-t0])
        lst_enr.append([np_cfg['selection']['emin'], np_cfg['selection']['emax']])

        lst_data.append(Datum(time=lst_time[-1], enr_cfg=lst_enr[-1], enr_sed=[np_sed['e_ref'], np_sed['e_min'], np_sed['e_max']], flux=[np_sed['flux'], np_sed['flux_err_lo'], np_sed['flux_err_hi'], np_sed['flux_ul95']], eflux=[np_sed['eflux'], np_sed['eflux_err_lo'], np_sed['eflux_err_hi'], np_sed['eflux_ul95']], model_dnde=np_sed['model_flux']))
        for (ipar, par) in enumerate(np_src['param_names']):
            if par=='Index':
                lst_data[-1].pl_index = [np_src['param_values'][ipar], np_src['param_errors'][ipar]]
    
    lst_enr_irr = []
    for e in lst_enr:
        if not e in lst_enr_irr:
            lst_enr_irr.append(e)
    lst_enr_irr.sort()
    nenr = len(lst_enr_irr)
    print nenr, 'energy bins.'
    lst_time_irr = []
    for t in lst_time:
        if not t in lst_time_irr:
            lst_time_irr.append(t)
    lst_time_irr.sort()
    ntime = len(lst_time_irr)
    print ntime, 'time bins.'
    
    # LC
    nrows_enr = int(sqrt(nenr))
    print '#row:', nrows_enr
    ncols_enr = int(ceil(nenr/nrows_enr))
    print '#cols:', ncols_enr
    fig_lc_flux, axes_lc_flux = plt.subplots(nrows=nrows_enr, ncols=ncols_enr, squeeze=False)
    fig_lc_eflux, axes_lc_eflux = plt.subplots(nrows=nrows_enr, ncols=ncols_enr, squeeze=False)
    fig_lc_index, axes_lc_index = plt.subplots(nrows=nrows_enr, ncols=ncols_enr, squeeze=False)
    lst_key_spec = ['Best', 'Softer', 'Harder']
    for (ienr, enr) in enumerate(lst_enr_irr):
        print 'Light curve for', enr, 'sec'
        xtime = []
        xtime_err = []
        yflux = []
        yflux_err_lo = []
        yflux_err_hi = []
        yflux_ul = []
        yeflux = []
        yeflux_err_lo = []
        yeflux_err_hi = []
        yeflux_ul = []
        yindex = []
        yindex_err = []
        ye2dnde = {}
        ye2dnde_lo = {}
        ye2dnde_hi = {}
        ye2dnde = []
        ye2dnde_lo = []
        ye2dnde_hi = []
        for datum in lst_data:
            if datum.enr_cfg==enr:
                # Time
                xtime.append((datum.time[1]+datum.time[0])/2.)
                xtime_err.append((datum.time[1]-datum.time[0])/2.)
                # flux
                if not isnan(datum.flux[2]):
                    yflux.append(datum.flux[0])
                    yflux_err_hi.append(datum.flux[1])
                    yflux_err_lo.append(datum.flux[2])
                    yflux_ul.append(datum.flux[3]) #-datum.flux[0])
                else:
                    yflux.append(datum.flux[0])
                    yflux_err_lo.append(0)
                    yflux_err_hi.append(0)
                    yflux_ul.append(datum.flux[3]) #-datum.flux[0])
                # Energy flux
                if not isnan(datum.eflux[2]):
                    yeflux.append(datum.eflux[0])
                    yeflux_err_hi.append(datum.eflux[1])
                    yeflux_err_lo.append(datum.eflux[2])
                    yeflux_ul.append(datum.eflux[3]) #-datum.eflux[0])
                else:
                    yeflux.append(datum.flux[0])
                    yeflux_err_lo.append(0)
                    yeflux_err_hi.append(0)
                    yeflux_ul.append(datum.flux[3]) #-datum.flux[0])
                # Power-law index
                yindex.append(datum.pl_index[0])
                yindex_err.append(datum.pl_index[1])
                # Model e2dnde 
                lst_e2dnde = datum.get_model_e2dnde_lowest_energy()
                ye2dnde.append(lst_e2dnde[1])
                ye2dnde_lo.append(lst_e2dnde[2])
                ye2dnde_hi.append(lst_e2dnde[3])
        print 'Time:', xtime
        print 'Time Error', xtime_err
        print 'Flux', yflux
        print 'Flux Error higher', yflux_err_hi
        print 'Flux Error lower', yflux_err_lo
        print 'Eflux', yeflux
        print 'Eflux Error higher', yeflux_err_hi
        print 'Eflux Error lower', yeflux_err_lo
        print 'Index', yindex
        print 'Index Error', yindex_err

        axes_lc_flux[ienr/ncols_enr,ienr%ncols_enr].set_xscale("log", nonposy='clip')
        axes_lc_flux[ienr/ncols_enr,ienr%ncols_enr].set_yscale("log", nonposy='clip')
        axes_lc_flux[ienr/ncols_enr,ienr%ncols_enr].errorbar(np.array(xtime), np.array(yflux), xerr=np.array(xtime_err), yerr=[np.array(yflux_err_lo), np.array(yflux_err_hi)], ls='')
        axes_lc_flux[ienr/ncols_enr,ienr%ncols_enr].errorbar(np.array(xtime), np.array(yflux_ul), xerr=np.array(xtime_err), yerr=np.array(yflux)*0.1, uplims=True, ls='')
      #  axes_lc_flux[ienr/ncols_enr,ienr%ncols_enr].errorbar(np.array(xtime), np.array(yflux_ul), xerr=np.array(xtime_err), ls='')
        axes_lc_flux[ienr/ncols_enr,ienr%ncols_enr].set_title('{0:d} - {1:d} MeV'.format(int(enr[0]), int(enr[1])))
        axes_lc_flux[ienr/ncols_enr,ienr%ncols_enr].set_ylim([5e-9, 5e-3])
        axes_lc_flux[ienr/ncols_enr,ienr%ncols_enr].set_xlabel('t - {0} [s]'.format(t0))
        axes_lc_flux[ienr/ncols_enr,ienr%ncols_enr].set_ylabel('Photon flux [cm^-2 s^-1]')

        axes_lc_eflux[ienr/ncols_enr,ienr%ncols_enr].set_xscale("log", nonposy='clip')
        axes_lc_eflux[ienr/ncols_enr,ienr%ncols_enr].set_yscale("log", nonposy='clip')
        axes_lc_eflux[ienr/ncols_enr,ienr%ncols_enr].errorbar(np.array(xtime), np.array(yeflux)*MEV2ERG, xerr=np.array(xtime_err), yerr=[np.array(yeflux_err_lo)*MEV2ERG, np.array(yeflux_err_hi)*MEV2ERG], ls='')
        axes_lc_eflux[ienr/ncols_enr,ienr%ncols_enr].errorbar(np.array(xtime), np.array(yeflux_ul)*MEV2ERG, xerr=np.array(xtime_err), yerr=np.array(yeflux)*MEV2ERG*0.1, uplims=True, ls='')
        axes_lc_eflux[ienr/ncols_enr,ienr%ncols_enr].set_title('{0:d} - {1:d} MeV'.format(int(enr[0]), int(enr[1])))
        axes_lc_eflux[ienr/ncols_enr,ienr%ncols_enr].set_ylim([5e-12, 1e-7])
        axes_lc_eflux[ienr/ncols_enr,ienr%ncols_enr].set_xlabel('t - {0} [s]'.format(t0))
        axes_lc_eflux[ienr/ncols_enr,ienr%ncols_enr].set_ylabel('Energy flux [erg cm^-2 s^-1]')

        axes_lc_index[ienr/ncols_enr,ienr%ncols_enr].set_xscale("log", nonposy='clip')
        axes_lc_index[ienr/ncols_enr,ienr%ncols_enr].errorbar(np.array(xtime), np.array(yindex), xerr=np.array(xtime_err), yerr=np.array(yindex_err), ls='')
        axes_lc_index[ienr/ncols_enr,ienr%ncols_enr].set_title('{0:d} - {1:d} MeV'.format(int(enr[0]), int(enr[1])))
        axes_lc_index[ienr/ncols_enr,ienr%ncols_enr].set_xlabel('t - {0} [s]'.format(t0))
        axes_lc_index[ienr/ncols_enr,ienr%ncols_enr].set_ylabel('Photon index')

        #plt.show()
    fig_lc_flux.tight_layout()
    fig_lc_eflux.tight_layout()
    fig_lc_index.tight_layout()
    fig_lc_flux.savefig('{0}_LightCurve_flux.png'.format(name_src))
    fig_lc_eflux.savefig('{0}_LightCurve_eflux.png'.format(name_src))
    fig_lc_index.savefig('{0}_LightCurve_index.png'.format(name_src))
        
    


# @click.command()
# @click.argument('srcname', type=str)
# @click.argument('inputfile', type=str)
# @click.option('--suffix', type=str, default='')
# def main(srcname, inputfile):
#     npi = np.load(inputfile).flat[0]
#     np_src = np['sources'][grbname]
#     PlotSED(np_src)


# if __name__ == '__main__':
#     main()
