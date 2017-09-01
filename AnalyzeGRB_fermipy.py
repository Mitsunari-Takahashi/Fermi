#!/usr/bin/env python

import os
import os.path
import subprocess
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
from sympy import *
from scipy import integrate
from fermipy.utils import get_parameter_limits
from fermipy.gtanalysis import GTAnalysis
import GtApp
import FluxDensity
from LikelihoodState import LikelihoodState
from fermipy.gtutils import BinnedAnalysis, SummedLikelihood
import BinnedAnalysis as ba
import pyLikelihood as pyLike
import ROOT

import numpy as np
from math import log10, log, sqrt, ceil, isnan, pi, factorial
import copy
from pLsList import ls_list
from pReadGBMCatalogueInfo import ReadGBMCatalogueOneLine
import click
from FindCrossEarthlimb import find_cross_earthlimb
from FindGoodstatPeriods import find_goodstat_periods, get_entries
import ReadLTFCatalogueInfo
import pMETandMJD
from DownloadFermiData import download_fermi_data_grb


def compute_GPoisson(x, m, s, n):
#    if n>0:
 #       return pow(x,n)*exp(-pow(x-m,2)/2./pow(s,2)-x)/sqrt(2.*pi)/s
    #  else:
    return pow(x,n)* ( exp(-pow(x-m,2)/2./pow(s,2)-x) + int(n==0)*exp(-pow(x+m,2)/2./pow(s,2)) ) /sqrt(2.*pi)/s

# def compute_GPoisson_zero(x, m, s):
#     return exp(-pow(x-m,2)/2./pow(s,2)-x)/sqrt(2.*pi)/s


def scale_limit(value, limit, norm):
    if norm>0:
        return value*limit/norm
    else:
        return 0.0

def zero_for_except(a):
    try:
        if isinstance(a, float) or isinstance(a, int):
            if isnan(a) is True:
                return 0.0
            else:
                return a
        else:
            return 0.0
    except (TypeError, NameError, IndexError):
        return 0.0


def profile(gta, name, parName, logemin=None, logemax=None, reoptimize=False, xvals=None, npts=None, savestate=True, **kwargs):
    """Profile the likelihood for the given source and parameter.
    Parameters
    ----------
        name : str
           Source name.
        parName : str
           Parameter name.
        reoptimize : bool
           Re-fit nuisance parameters at each step in the scan.  Note
           that enabling this option will only re-fit parameters that
           were free when the method was executed.
    Returns
    -------
        lnlprofile : dict
           Dictionary containing results of likelihood scan.
           """
        # Find the source
    name = gta.roi.get_source_by_name(name).name

    par = gta.like.normPar(name)
    idx = gta.like.par_index(name, parName)
    bounds = gta.like.model[idx].getBounds()
    value = gta.like.model[idx].getValue()
    loge_bounds = gta.loge_bounds
    
    optimizer = kwargs.get('optimizer', gta.config['optimizer'])

    if savestate:
        saved_state = gta._latch_state()
            
        # If parameter is fixed temporarily free it
    par.setFree(True)
    if optimizer['optimizer'] == 'NEWTON':
        gta._create_fitcache()

    if logemin is not None or logemax is not None:
        loge_bounds = gta.set_energy_range(logemin, logemax)
    else:
        loge_bounds = gta.loge_bounds

    loglike0 = -gta.like()

    if xvals is None:
        err = par.error()
        val = par.getValue()
        if err <= 0 or val <= 3 * err:
            xvals = 10 ** np.linspace(-2.0, 2.0, 51)
            if val < xvals[0]:
                xvals = np.insert(xvals, val, 0)
        else:
            xvals = np.linspace(0, 1, 25)
            xvals = np.concatenate((-1.0 * xvals[1:][::-1], xvals))
            xvals = val * 10 ** xvals

        # Update parameter bounds to encompass scan range
    gta.like[idx].setBounds(min(min(xvals), value, bounds[0]),
                                 max(max(xvals), value, bounds[1]))

    o = {'xvals': xvals,
         'npred': np.zeros(len(xvals)),
         'dnde': np.zeros(len(xvals)),
         'flux': np.zeros(len(xvals)),
         'eflux': np.zeros(len(xvals)),
         'dloglike': np.zeros(len(xvals)),
         'loglike': np.zeros(len(xvals))
         }

    if reoptimize and hasattr(gta.like.components[0].logLike,
                                  'setUpdateFixedWeights'):

        for c in gta.components:
            c.like.logLike.setUpdateFixedWeights(False)

    for i, x in enumerate(xvals):

        gta.like[idx] = x
        if gta.like.nFreeParams() > 1 and reoptimize:
                # Only reoptimize if not all frozen
            gta.like.freeze(idx)
            fit_output = gta._fit(errors=False, **optimizer)
            loglike1 = fit_output['loglike']
            gta.like.thaw(idx)
        else:
            loglike1 = -gta.like()

        flux = gta.like[name].flux(10 ** loge_bounds[0],
                                        10 ** loge_bounds[1])
        eflux = gta.like[name].energyFlux(10 ** loge_bounds[0],
                                               10 ** loge_bounds[1])
        prefactor = gta.like[idx]

        o['dloglike'][i] = loglike1 - loglike0
        o['loglike'][i] = loglike1
        o['dnde'][i] = prefactor.getTrueValue()
        o['flux'][i] = flux
        o['eflux'][i] = eflux

        cs = gta.model_counts_spectrum(name,
                                        loge_bounds[0],
                                        loge_bounds[1], summed=True)
        o['npred'][i] += np.sum(cs)

    gta.like[idx] = value

    if reoptimize and hasattr(gta.like.components[0].logLike,
                                  'setUpdateFixedWeights'):

        for c in gta.components:
            c.like.logLike.setUpdateFixedWeights(True)

        # Restore model parameters to original values
    if savestate:
        saved_state.restore()

    gta.like[idx].setBounds(*bounds)
    if logemin is not None or logemax is not None:
        gta.set_energy_range(*loge_bounds)

    return o


def AnalyzeGRB_fermipy(name, ft1_candidates, ft2_candidates, tmin, tmax, tbinedges, suffix, force, skipts, skipsed, skipresid, eranges, tb_masked, path_outdir='.', mode='unified', catalogues=['3FGL'], goodstat=16, shiftenergies=True, edisp=False, lst_spec_func=['PL'], index_fixed=None, ecut_fixed=None, radius_roi=12., sedadjusted=False) :
    NAME_TGT = name
    #dct_grb = ReadGBMCatalogueOneLine(NAME_TGT, '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/Highest-GBM-fluence-GRBs.csv')
    #print dct_grb
    path_first = os.getcwd()

    RA = tb_masked['RA'] #dct_grb['ra']
    DEC = tb_masked['DEC'] #dct_grb['dec']
    T0 = pMETandMJD.ConvertMjdToMet(tb_masked['TRIGGER_TIME']) #dct_grb['trigger_time']
    T90 = tb_masked['T90'] #dct_grb['t90']
    T90_START = tb_masked['T90_START'] #dct_grb['t90_start']
    if tmin is None:
        if mode in ('prompt', 'unified'):
            tmin = 0.
        elif mode in ('afterglow', 'earlyAG'):
            tmin = T90_START+T90
        elif mode in ('lateAG'):
            tmin = T90_START+T90*3.
    if tmax is None:
        if mode in ('prompt'):
            tmax =  T90_START+T90
        elif mode in ('earlyAG'):
            tmax = T90_START+T90*3.
        elif mode in ('lateAG', 'unified', 'afterglow'):
            tmax = 10000.
    eranges_shifted = []
    # erange_sed_lin = [100.0, 316.228, 1000.0, 3162.28, 10000.0, 31622.8, 100000.0]
    # if shiftenergies==True:
    #     erange_sed_lin = [316.228, 1778.28, 5623.41, 17782.8, 56234.1, 177828.0] #[316.228, 1778.28, 10000.0, 56234.1, 316228.] #[316.228, 1000.0, 3162.28, 10000.0, 56234.1, 316228.]
    # erange_sed_shifted_lin = []
    # erange_hiend = [10000., 100000.] #[17782.8, 100000.] #erange_sed_lin[-2:] #[56234., 316228.]
    # erange_hiend_shifted = []
    #if (dct_grb['z'] is not '') and (dct_grb['z'] is not '-') and (dct_grb['z'] is not '?'):
    if tb_masked['REDSHIFT']>0:
        REDSHIFT = tb_masked['REDSHIFT'] #dct_grb['z']
        if shiftenergies==True:
            for eset in eranges:
                eranges_shifted.append([e/(1+REDSHIFT) for e in eset])
            #erange_hiend_shifted = [e/(1+REDSHIFT) for e in erange_hiend]
            #erange_sed_shifted_lin = [e/(1+REDSHIFT) for e in erange_sed_lin]
        else:
            eranges_shifted = eranges
            #erange_hiend_shifted = erange_hiend
            #erange_sed_shifted_lin = erange_sed_lin
    else:
        REDSHIFT = 1
        eranges_shifted = eranges
        #erange_hiend_shifted = erange_hiend
     #   erange_sed_shifted_lin = erange_sed_lin
    #erange_sed_shifted = np.array([log10(x) for x in erange_sed_shifted_lin])
    print 'Shifted energy ranges:', eranges_shifted, 'MeV'
    SUFFIX = ''
    if suffix!='':
        SUFFIX = '_' + suffix

    # Model catalogues
    if len(catalogues)>0:
        str_catalogues = 'catalogs : ["'
        for cata in catalogues:
            str_catalogues = str_catalogues + cata + '", "'
        str_catalogues = str_catalogues[:-3] + ']'
    else:
        str_catalogues = ''

    ZCUT = 100.
    lst_tbin = []
    mstatag = 0
    str_suffix_csv = SUFFIX

    if mode in ('unified', 'prompt', 'afterglow', 'earlyAG', 'lateAG'):
        lst_tbin = [[tmin, tmax]]
        str_suffix_csv = mode + str_suffix_csv
    else:
        validtimes = find_cross_earthlimb(ft2_candidates[0], RA, DEC, T0+tmin, T0+tmax, ZCUT, T0)
        print validtimes
        for (ivt, vt) in enumerate(validtimes):
            if goodstat>0:
                if mode=='prompt':
                    periods_goodstat = find_goodstat_periods(ft1_candidates[0], T0+vt[0], T0+vt[1], goodstat)
                    print periods_goodstat
                    for igs in range(len(periods_goodstat)-1):
                        lst_tbin.append([periods_goodstat[igs]-T0, periods_goodstat[igs+1]-T0])
                elif mode=='afterglow':
                    mstatag += get_entries(ft1_candidates[0], vt[0]+T0, vt[1]+T0)
                    if ivt==0:
                        lst_tbin.append([vt[0]])
                    if mstatag>=goodstat:
                        lst_tbin[-1].append(vt[1])
                        mstatag = 0
                        if ivt<len(validtimes)-1:
                            lst_tbin.append([validtimes[ivt+1][0]])
                    if ivt==len(validtimes)-1:
                        if len(lst_tbin)>1:
                            lst_tbin = lst_tbin[:-1]
                        if len(lst_tbin[-1])<2:
                            lst_tbin[-1].append(vt[1])
                        else:
                            lst_tbin[-1][1] = vt[1]
            else:
                if ivt==0:
                    if tmin>=vt[0] and tmin<vt[1]:
                        lst_tbin = [[tmin]]
                    else:
                        lst_tbin = [[vt[0]]]
                    if tbinedges is not None and len(tbinedges)>0:
                        for tedge in tbinedges:
                            if vt[1]-tedge>tedge:
                                if tedge < vt[1]:
                                    lst_tbin[-1].append(tedge)
                                    lst_tbin.append([tedge])
                            else:
                                lst_tbin[-1].append(vt[1])
                    else:
                        lst_tbin[-1].append(vt[1])
                else: 
                    lst_tbin.append([vt[0]])
                    lst_tbin[-1].append(vt[1])
            
    print 'Time bin edges:', lst_tbin

    LST_RAN_TIME = lst_tbin
    NRAN_TIME = len(LST_RAN_TIME)
    LST_SED_ITEM_CSV = ['e_min', 'e_max', 'e_ref', 'index', 'ts', 'flux', 'flux_err_lo', 'flux_err_hi', 'eflux', 'eflux_err_lo', 'eflux_err_hi', 'ref_flux', 'ref_dnde', 'ref_dnde_e_min', 'ref_dnde_e_max']

    print 'Going to start standard analysis of', NAME_TGT
    if path_outdir=='.':
        path_outdir = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/' + NAME_TGT
    if not os.path.exists(path_outdir):
        os.makedirs(path_outdir)
    params_prev = [1e-10, 2.0, 1000.0, 10000.]
    fluxhe_frac_err_prev = 1.

    # Figure for comparison of extrapolated spectrum and observed counts
    fig_cspec = plt.figure()
    ax_cspec = fig_cspec.add_axes((0.1, 0.1, 0.8, 0.8))
    x_cspec_extrapolated_eref = []
    #x_cspec_extrapolated_ebin = []
    x_cspec_extrapolated_eref_errhi = []
    x_cspec_extrapolated_eref_errlo = []
    x_cspec_extrapolated_logemin = []
    x_cspec_extrapolated_logemax = []
    ie_half_decade = 0
    loge_append = log10(eranges_shifted[-1][0])
    #for ie_half_decade in range(12):#12):
    while loge_append<log10(eranges_shifted[-1][1]):
        #if loge_append<4.001:
        x_cspec_extrapolated_logemin.append(loge_append)
        x_cspec_extrapolated_logemax.append(loge_append+0.25)
        x_cspec_extrapolated_eref.append(10**(loge_append+0.125))
        x_cspec_extrapolated_eref_errhi.append(10**x_cspec_extrapolated_logemax[-1]-x_cspec_extrapolated_eref[-1])
        x_cspec_extrapolated_eref_errlo.append(x_cspec_extrapolated_eref[-1]-10**x_cspec_extrapolated_logemin[-1])
        #x_cspec_extrapolated_ebin.append(10**(loge_append))
        ie_half_decade+=1
        loge_append = log10(eranges_shifted[-1][0])+0.25*ie_half_decade
    x_cspec_extrapolated_logemin = np.array(x_cspec_extrapolated_logemin)
    x_cspec_extrapolated_logemax[-1] = log10(eranges_shifted[-1][1])
    x_cspec_extrapolated_logemax = np.array(x_cspec_extrapolated_logemax)
    x_cspec_extrapolated_eref[-1] = 10**((x_cspec_extrapolated_logemin[-1]+x_cspec_extrapolated_logemax[-1])/2.0)
    x_cspec_extrapolated_eref_errhi[-1] = 10**x_cspec_extrapolated_logemax[-1]-x_cspec_extrapolated_eref[-1]
    x_cspec_extrapolated_eref_errlo[-1] = x_cspec_extrapolated_eref[-1]-10**x_cspec_extrapolated_logemin[-1]
    x_cspec_extrapolated_eref = np.array(x_cspec_extrapolated_eref)
    x_cspec_extrapolated_eref_errhi = np.array(x_cspec_extrapolated_eref_errhi)
    x_cspec_extrapolated_eref_errlo = np.array(x_cspec_extrapolated_eref_errlo)
    #x_cspec_extrapolated_ebin.append(eranges_shifted[-1][1])
    #x_cspec_extrapolated_ebin = np.array(x_cspec_extrapolated_ebin)
    print 'Energy axis for count spectrum:'
    print x_cspec_extrapolated_logemin
    print x_cspec_extrapolated_logemax
    print x_cspec_extrapolated_eref
    print x_cspec_extrapolated_eref_errhi
    print x_cspec_extrapolated_eref_errlo
    
    cspec_extrapolated_prev = None
    cspec_extrapolated_tgt_prev = None
    cspec_extrapolated_other_prev = None
    cspec_extrapolated_err_prev = None
    flux_prev = None
    flux_err_prev = None
    flux_fracerr_prev = None
    fitresult_prev = None

    str_lc = """#name/C:function/C:start/F:stop:emin_rest:emax_rest:emin_shifted:emax_shifted:ts:Integral:Integral_err:Integral_ul95:Integral_ll95:Integral_ul68:Integral_ll68:Index1:Index1_err:Index2:Index2_err:BreakValue:BreakValue_err:flux:flux_err:flux_ul95:flux_ll95:flux_ul68:flux_ll68:eflux:eflux_err:eflux_ul95:eflux_ll95:eflux_ul68:eflux_ll68:fluxhe:fluxhe_err:efluxhe:efluxhe_err:flux_hiest_e:flux_hiest_e_err_hi:flux_hiest_e_err_lo:flux_hiest_e_ul:flux_hiest_e_ll:eflux_hiest_e:eflux_hiest_e_err_hi:eflux_hiest_e_err_lo:eflux_hiest_e_ul:eflux_hiest_e_ll:npred_hiest_e:npred_hiest_e_err:npred_all_hiest_e:nobs_hiest_e:deviation:sign_deviation:ecutoff:ecutoff_err_hi:ecutoff_err_lo:ecutoff_ul95:ecutoff_ll95:ecutoff_ul68:ecutoff_ll68:loglike
"""

    for (ieedges, eedges) in enumerate(eranges_shifted):
        #strenergies = 'E{0:0>6}-{1:0>6}MeV'.format(int(eedges[0]+0.5), int(eedges[1]+0.5))
        strenergies = 'E{0:0>6}-{1:0>6}MeV'.format(int(eranges[ieedges][0]+0.5), int(eranges[ieedges][1]+0.5))
        if shiftenergies==True:
            strenergies += '_shifted'
        print '%%%%%%%%%%%%%%%%%%'
        print strenergies
        print '%%%%%%%%%%%%%%%%%%'
        erange_sed_shifted = []
        erange_hiend_shifted = [0.1*eedges[1], eedges[1]]
        erange_extrapolate_shifted = [10000.0, 100000.0]
        ne_half_decade = int(2.*log10(eedges[1]/eedges[0]))
        for ie_half_decade in range(ne_half_decade):
            loge_append = log10(eedges[0])+0.5*ie_half_decade
            if loge_append<4.001:
                erange_sed_shifted.append(loge_append)
        erange_sed_shifted.append(log10(eedges[1]))
        erange_sed_shifted = np.array(erange_sed_shifted)
        print 'Initial parameters:', params_prev
        for itime in range(NRAN_TIME):
            strtime = 'T{0:0>6}-{1:0>6}s'.format(int(0.5+LST_RAN_TIME[itime][0]), int(0.5+LST_RAN_TIME[itime][1]))
            print '=====', strtime, '====='
            if mode!='special':
                strtime = mode
            dct_loglike = {}
            for (ispec, fspec) in enumerate(lst_spec_func):
                if fspec=='ExpCutoff' and eedges[1]<=10000.:
                    continue
                print '-----', fspec, '-----'
                str_index = 'IndexFree'
                if index_fixed is not None:
                    str_index = 'Index{0:0>3}'.format(int(100*index_fixed))
                path_subdir = '{0}/{1}/r{2:0>2}deg/{3}/{4}/{5}'.format(path_outdir, strenergies, int(radius_roi+0.5), strtime, fspec, str_index)
                if not os.path.exists(path_subdir):
                    os.makedirs(path_subdir)
                os.chdir(path_subdir)
                path_anlaysis = '{0}/fit_model{1}.npy'.format(path_subdir, SUFFIX)
                if os.path.exists(path_anlaysis) and force==False:
                    print 'Loading previous analysis...'
                    gta = GTAnalysis.create(path_anlaysis)
                else:
                    str_path_lt = path_subdir+'/ltcube_00.fits'
                    if not os.path.exists(str_path_lt):
                        str_path_lt = 'Null'
                    if fspec=='PL':
                        if not isinstance(index_fixed, float): #==None:
                            str_spectrum = """'SpectrumType' : 'PowerLaw', 'Prefactor' : {{ value : {0}, max : !!float 1e-3, min : !!float 1e-15, scale : 1.0, free : '1' }}, 'Index' : {{ value : {1}, min : -1.0, max : 8.0, scale : -1, free : '1' }}, 'Scale' : {{ value : {2}, max : 100000., min : 30, scale : 1, free : '0' }}""".format(params_prev[0], -params_prev[1], params_prev[2])
                        else:
                            str_spectrum = """'SpectrumType' : 'PowerLaw', 'Prefactor' : {{ value : !!float 1e-10, max : !!float 1e-3, min : !!float 1e-15, scale : 1, free : '1' }}, 'Index' : {{ value : {0}, min : -1.0, max : 8.0, scale : -1, free : '0' }}, 'Scale' : {{ value : 1000.0, max : 100000., min : 30, scale : 1, free : '0' }}""".format(index_fixed)
                    elif fspec=='ExpCutoff':
                        if not isinstance(ecut_fixed, float): #ecut_fixed==None:
                            str_spectrum = """'SpectrumType' : 'ExpCutoff', 'Prefactor' : {{ value : {0}, max : !!float 1e-3, min : !!float 1e-15, scale : 1.0, free : '1' }}, 'Index' : {{ value : {1}, min : -1.0, max : 8.0, scale : -1, free : '1' }}, 'Scale' : {{ value : {2}, max : 100000., min : 30, scale : 1, free : '0' }}, 'Ebreak' : {{ value : 10, min : 1, max : 300, scale : 1, free : '1' }}, 'P1' : {{ value : {3}, min : 10., max : 300000., scale : 1, free : '1' }}, 'P2' : {{ value : 0, min : 0.1, max : 10, scale : 1, free : '0' }}, 'P3' : {{ value : 0, min : 0.1, max : 10, scale : 1, free : '0' }}""".format(params_prev[0], -params_prev[1], params_prev[2], params_prev[3])
                        else:
                            str_spectrum = """'SpectrumType' : 'ExpCutoff', 'Prefactor' : {{ value : {0}, max : !!float 1e-3, min : !!float 1e-15, scale : 1.0, free : '1' }}, 'Index' : {{ value : {1}, min : -1.0, max : 8.0, scale : -1, free : '1' }}, 'Scale' : {{ value : {2}, max : 100000., min : 30, scale : 1, free : '0' }}, 'Ebreak' : {{ value : 0, min : 0, max : 50000, scale : 1, free : '0' }}, 'P1' : {{ value : {3}, min : 10., max : 100000., scale : 1, free : '0' }}, 'P2' : {{ value : 0, min : 0.1, max : 10, scale : 1, free : '0' }}, 'P3' : {{ value : 0, min : 0.1, max : 10, scale : 1, free : '0' }}""".format(params_prev[0], -params_prev[1], params_prev[2], ecut_fixed)
                    elif fspec=='BPL':
                        if index_fixed==None:
                            str_spectrum = """'SpectrumType' : 'BrokenPowerLaw', 'Prefactor' : { value : !!float 1e-10, max : !!float 1e-3, min : !!float 1e-15, scale : 1, free : '1' }, 'Index1' : { value : 2.3, min : -1.0, max : 8.0, scale : -1, free : '1' }, 'Index2' : { value : 1.6, min : -1.0, max : 8.0, scale : -1, free : '1' }, 'BreakValue' : { value : 5000, min : 500, max : 30000, scale : 1, free : '1' }"""
                        else:
                            str_spectrum = """'SpectrumType' : 'BrokenPowerLaw', 'Prefactor' : {{ value : !!float 1e-10, max : !!float 1e-3, min : !!float 1e-15, scale : 1, free : '1' }}, 'Index1' : {{ value : {0}, min : -1.0, max : 8.0, scale : -1, free : '0' }}, 'Index2' : {{ value : 1.6, min : -1.0, max : 8.0, scale : -1, free : '1' }}, 'BreakValue' : {{ value : 5000, min : 500, max : 30000, scale : 1, free : '1' }}""".format(index_fixed)
                    if fspec=='EblPL':
                        if index_fixed==None:
                            str_spectrum = """'SpectrumType' : 'EblAtten::PowerLaw2', 'Integral' : {{ value : 1.0, max : !!float 1e6, min : !!float 1e-6, scale : !!float 1e-6, free : '1' }}, 'Index' : {{ value : 2.0, min : -1.0, max : 8.0, scale : -1, free : '1' }}, 'LowerLimit' : {{ value : {0}, max : 100000., min : 30, scale : 1, free : '0' }}, 'UpperLimit' : {{ value : {1}, max : 1000000., min : 100., scale : 1, free : '0' }}, 'tau_norm' : {{ value : 1.0, max : 10, min : 0, scale : 1.0, free : '0' }}, 'redshift' : {{ value : {2}, max : 10, min : 0, scale : 1, free : '0' }}, 'ebl_model' : {{ value : 4, max : 8, min : 0, scale : 1.0, free : '0'}}""".format(eedges[0], eedges[1], REDSHIFT)
                        else:
                            str_spectrum = """'SpectrumType' : 'EblAtten::PowerLaw2', 'Integral' : {{ value : 1.0, max : !!float 1e6, min : !!float 1e-6, scale : !!float 1e-6, free : '1' }}, 'Index' : {{ value : {3}, min : -1.0, max : 8.0, scale : -1, free : '0' }}, 'LowerLimit' : {{ value : {0}, max : 100000., min : 30, scale : 1, free : '0' }}, 'UpperLimit' : {{ value : {1}, max : 1000000., min : 100., scale : 1, free : '0' }}, 'tau_norm' : {{ value : 1.0, max : 10, min : 0, scale : 1.0, free : '0' }}, 'redshift' : {{ value : {2}, max : 10, min : 0, scale : 1, free : '0' }}, 'ebl_model' : {{ value : 4, max : 8, min : 0, scale : 1.0, free : '0'}}""".format(eedges[0], eedges[1], REDSHIFT, index_fixed)
                    elif fspec=='EblBPL':
                        if index_fixed==None:
                            str_spectrum = """'SpectrumType' : 'EblAtten::BrokenPowerLaw2', 'Integral' : {{ value : 1.0, max : !!float 1e6, min : !!float 1e-6, scale : !!float 1e-6, free : '1' }}, 'Index1' : {{ value : 2.3, min : -1.0, max : 8.0, scale : -1, free : '1' }}, 'Index2' : {{ value : 1.6, min : -1.0, max : 8.0, scale : -1, free : '1' }}, 'BreakValue' : {{ value : 5000, min : 500, max : 30000, scale : 1, free : '1' }}, 'LowerLimit' : {{ value : {0}, max : 100000., min : 30, scale : 1, free : '0' }}, 'UpperLimit' : {{ value : {1}, max : 1000000., min : 100., scale : 1, free : '0' }}, 'tau_norm' : {{ value : 1.0, max : 10, min : 0, scale : 1.0, free : '0' }}, 'redshift' : {{ value : {2}, max : 10, min : 0, scale : 1, free : '0' }}, 'ebl_model' : {{ value : 4, max : 8, min : 0, scale : 1.0, free : '0'}}""".format(eedges[0], eedges[1], REDSHIFT)
                        else:
                            str_spectrum = """'SpectrumType' : 'EblAtten::BrokenPowerLaw2', 'Integral' : {{ value : 1.0, max : !!float 1e6, min : !!float 1e-6, scale : !!float 1e-6, free : '1' }}, 'Index1' : {{ value : {3}, min : -1.0, max : 8.0, scale : -1, free : '0' }}, 'Index2' : {{ value : 1.6, min : -1.0, max : 8.0, scale : -1, free : '1' }}, 'BreakValue' : {{ value : 5000, min : 500, max : 30000, scale : 1, free : '1' }}, 'LowerLimit' : {{ value : {0}, max : 100000., min : 30, scale : 1, free : '0' }}, 'UpperLimit' : {{ value : {1}, max : 1000000., min : 100., scale : 1, free : '0' }}, 'tau_norm' : {{ value : 1.0, max : 10, min : 0, scale : 1.0, free : '0' }}, 'redshift' : {{ value : {2}, max : 10, min : 0, scale : 1, free : '0' }}, 'ebl_model' : {{ value : 4, max : 8, min : 0, scale : 1.0, free : '0'}}""".format(eedges[0], eedges[1], REDSHIFT, index_fixed)

                    str_config = """fileio:
  outdir : {0}
data:
  evfile : {1}
  scfile : {2}
  ltcube : {3} #ltcube : ltcube_00.fits
binning:
  roiwidth   : {15} #21.
  binsz      : 0.1
  binsperdec : 4
selection :
  emin : {4}
  emax : {5}
  zmax    : {12}
  evclass : 128 # 8
  evtype  : 3
  tmin    : {6}
  tmax    : {7}
  filter  : null
  ra      : {8}
  dec     : {9}
gtlike:
  edisp : {14}
  irfs : 'P8R2_SOURCE_V6' #'P8R2_TRANSIENT020E_V6' #'P8R2_SOURCE_V6'
  edisp_disable : ['isodiff','galdiff']
model:
  src_radius  : {16} #25.0
  galdiff  : '/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v5r0/gll_iem_v06.fits'
  isodiff  : '/afs/slac.stanford.edu/g/glast/ground/GLAST_EXT/diffuseModels/v5r0/iso_P8R2_SOURCE_V6_v06.txt' #iso_P8R2_TRANSIENT020E_V6_v06.txt'
  {13}
  sources :
    - {{ 'name' : 'GRB{10}', 'ra' : {8}, 'dec' :{9}, {11}, 'SpatialModel': 'PointSource'}}
""".format(path_subdir, ft1_candidates[0], ft2_candidates[0], str_path_lt, eedges[0], eedges[1], int(T0+LST_RAN_TIME[itime][0]+0.5), int(T0+LST_RAN_TIME[itime][1]+0.5), RA, DEC, NAME_TGT, str_spectrum, ZCUT, str_catalogues, str(edisp), ceil(radius_roi*sqrt(2)), radius_roi+10.)

                    with open("{0}/config.yaml".format(path_subdir), 'w') as conf:
                        conf.write(str_config)
                    with open("{0}/config.yaml".format(path_subdir), 'r') as conf:
                        print conf

                    print 'Setting up...'
                    gta = GTAnalysis('{0}/config.yaml'.format(path_subdir),logging={'verbosity' : 3})
                    #try:
                    gta.setup()
                    #except RuntimeError:
                     #   print 'RuntimeError'
                    print 'Checking ft1 file...'
                    hdulist=fits.open('{0}/ft1_00.fits'.format(path_subdir))
                    print 'FT1 file has', len(hdulist['EVENTS'].data), 'events.'
                    #    return 1
                    if index_fixed is not None:
                        if fspec in ('PL', 'EblPL'):
                            gta.lock_parameter(NAME_TGT, 'Index', lock=True)
                        if fspec in ('BPL', 'EblBPL'):
                            gta.lock_parameter(NAME_TGT, 'Index1', lock=True)
                    if ieedges==0:
                        gta.optimize()
                    #else:
                    # Detive spectrum from previous results
                    npredhe_prev = sum(gta.model_counts_spectrum('GRB'+NAME_TGT, log10(erange_hiend_shifted[0]), log10(erange_hiend_shifted[1]))[0])
                    npredhe_prev_all = 0
                    for src in gta.get_sources():
                        npredhe_prev_all += sum(gta.model_counts_spectrum(src.name, log10(erange_hiend_shifted[0]), log10(erange_hiend_shifted[1]))[0])

                    if flux_fracerr_prev is not None:
                        cspec_extrapolated_prev = []
                        cspec_extrapolated_tgt_prev = []
                        cspec_extrapolated_other_prev = []
                        cspec_extrapolated_err_prev = []
                        for me in range(len(x_cspec_extrapolated_eref)):
                            nc_extra = 0
                            nc_other_extra = 0
                            for src in gta.get_sources():
                                nc_src = sum(gta.model_counts_spectrum(src.name, x_cspec_extrapolated_logemin[me], x_cspec_extrapolated_logemax[me])[0])
                                if src.name[3:]==NAME_TGT:
                                    cspec_extrapolated_tgt_prev.append(nc_src)
                                else:
                                    nc_other_extra += nc_src
                                if me==0:
                                    print src.name, nc_src
                                nc_extra += nc_src
                            cspec_extrapolated_prev.append(nc_extra)
                            cspec_extrapolated_other_prev.append(nc_other_extra)
                            cspec_extrapolated_err_prev.append(nc_extra*flux_fracerr_prev[me])
                        cspec_extrapolated_prev = np.array(cspec_extrapolated_prev)
                        print cspec_extrapolated_tgt_prev
                        cspec_extrapolated_tgt_prev = np.array(cspec_extrapolated_tgt_prev)
                        cspec_extrapolated_other_prev = np.array(cspec_extrapolated_other_prev)
                        cspec_extrapolated_err_prev = np.array(cspec_extrapolated_err_prev)

                    if ieedges>0:
                        gta.optimize()
                    gta.print_roi()
            
                    gta.free_sources(free=False)
                    #gta.free_source('galdiff', free=True, pars=['Prefactor', 'Index'])
                    if fspec=='PL':
                        if index_fixed==None:
                            gta.free_source('GRB'+NAME_TGT, free=True, pars=['Prefactor', 'Index'])
                        else:
                            gta.free_source('GRB'+NAME_TGT, free=True, pars=['Prefactor'])
                    elif fspec=='ExpCutoff':
                        if  not isinstance(ecut_fixed, float):
                            gta.free_source('GRB'+NAME_TGT, free=True, pars=['Prefactor', 'Index', 'Ebreak', 'P1'])
                        else:
                            gta.free_source('GRB'+NAME_TGT, free=True, pars=['Prefactor', 'Index', 'Ebreak'])
                    elif fspec=='BPL':
                        if index_fixed==None:
                            gta.free_source('GRB'+NAME_TGT, free=True, pars=['Prefactor', 'Index1', 'Index2', 'BreakValue'])
                        else:
                            gta.free_source('GRB'+NAME_TGT, free=True, pars=['Prefactor', 'Index2', 'BreakValue'])
                    elif fspec=='EblPL':
                        if index_fixed==None:
                            gta.free_source('GRB'+NAME_TGT, free=True, pars=['Integral', 'Index'])
                        else:
                            gta.free_source('GRB'+NAME_TGT, free=True, pars=['Integral'])
                    elif fspec=='EblBPL':
                        if index_fixed==None:
                            gta.free_source('GRB'+NAME_TGT, free=True, pars=['Integral', 'Index1', 'Index2', 'BreakValue'])
                        else:
                            gta.free_source('GRB'+NAME_TGT, free=True, pars=['Integral', 'Index2', 'BreakValue'])
                    else:
                        print fspec, 'is not available!!!'
                        return 1

                    gta.free_sources(distance=3.0,pars='norm')
                    gta.free_sources(minmax_npred=[0.1,None],pars='norm')
                    gta.free_source('galdiff',free=False)
                    #gta.free_source('isodiff',free=False)
                    print 'Fitting...'
                    fitresult = gta.fit()
                    gta.free_sources(minmax_npred=[None,0.5],free=False)
                    gta.free_sources(minmax_npred=[0.5,None],free=True)
                    #gta.free_source('galdiff',free=False)
                    gta.free_source('isodiff',free=False)
                    print 'Fitting...'
                    fitresult = gta.fit()
                    gta.write_roi('fit_model'+SUFFIX)
                    gta.print_roi()
                    print ' Fitting finished.'

                    #npout = np.load('{0}/fit_model.npy'.format(path_subdir)).flat[0]                
                    npout = np.load(path_anlaysis).flat[0]                
                    dct_loglike[fspec] = fitresult['loglike']
                    print '** log(likelihood):', dct_loglike[fspec]
                    src_model = gta.get_src_model('GRB'+NAME_TGT)
                    if not src_model['ts']==src_model['ts']:
                        print 'TS is NaN!!!'
                        print 'Checking ft1 file...'
                        hdulist=fits.open('{0}/ft1_00.fits'.format(path_subdir))
                        print 'FT1 file has', len(hdulist['EVENTS'].data), 'events.'
                        return 1
                    print '** TS:', src_model['ts']
                    norm_lims95 = get_parameter_limits(src_model['norm_scan'], src_model['loglike_scan'], 0.95)
                    norm_lims68 = get_parameter_limits(src_model['norm_scan'], src_model['loglike_scan'], 0.68)
                    if fspec=='PL' or fspec=='BPL':
                        print '** Prefactor:', src_model['param_values'][0], '+/-', src_model['param_errors'][0]
                    if fspec=='EblPL' or fspec=='EblBPL':
                        print '** Integral:', src_model['param_values'][0], '+/-', src_model['param_errors'][0]
                    print '  95% limits:', scale_limit(src_model['param_values'][0], norm_lims95['ll'], norm_lims95['x0']), '-', scale_limit(src_model['param_values'][0], norm_lims95['ul'], norm_lims95['x0'])
                    print '  68% limits:', scale_limit(src_model['param_values'][0], norm_lims68['ll'], norm_lims95['x0']), '-', scale_limit(src_model['param_values'][0], norm_lims68['ul'], norm_lims95['x0'])
                    print '** Flux:', src_model['flux'], '+/-', src_model['flux_err'], '(UL:', src_model['flux_ul95'], ')'
                    print '  95% limits:', scale_limit(src_model['flux'], norm_lims95['ll'], norm_lims95['x0']), '-', scale_limit(src_model['flux'], norm_lims95['ul'], norm_lims95['x0'])
                    print '  68% limits:', scale_limit(src_model['flux'], norm_lims68['ll'], norm_lims95['x0']), '-', scale_limit(src_model['flux'], norm_lims68['ul'], norm_lims95['x0'])
                    print '** Energy flux:', src_model['eflux'], '+/-', src_model['eflux_err'], '(UL:', src_model['eflux_ul95'], ')'
                    print '  95% limits:', scale_limit(src_model['eflux'], norm_lims95['ll'], norm_lims95['x0']), '-', scale_limit(src_model['eflux'], norm_lims95['ul'], norm_lims95['x0'])
                    print '  68% limits:', scale_limit(src_model['eflux'], norm_lims68['ll'], norm_lims95['x0']), '-', scale_limit(src_model['eflux'], norm_lims68['ul'], norm_lims95['x0'])
                    if fspec=='PL':
                        print '** Index:', src_model['param_values'][1], '+/-', src_model['param_errors'][1]
                    elif fspec=='BPL':
                        print '** Index1:', src_model['param_values'][1], '+/-', src_model['param_errors'][1]
                        print '** Index2:', src_model['param_values'][2], '+/-', src_model['param_errors'][2]
                        print '** Break energy:', src_model['param_values'][3], '+/-', src_model['param_errors'][3]
                    elif fspec=='EblPL':
                        print '** Index:', src_model['param_values'][1], '+/-', src_model['param_errors'][1]
                    elif fspec=='EblBPL':
                        print '** Index1:', src_model['param_values'][1], '+/-', src_model['param_errors'][1]
                        print '** Index2:', src_model['param_values'][2], '+/-', src_model['param_errors'][2]
                        print '** Break energy:', src_model['param_values'][3], '+/-', src_model['param_errors'][3]
                    fluxhe = gta.like.flux('GRB'+NAME_TGT, erange_extrapolate_shifted[0], erange_extrapolate_shifted[1])
                    fluxhe_err = gta.like.fluxError('GRB'+NAME_TGT, erange_extrapolate_shifted[0], erange_extrapolate_shifted[1])
                    print '** Extrapolated flux in', int(erange_extrapolate_shifted[0]+0.5), '-', int(erange_extrapolate_shifted[1]+0.5), 'MeV:', fluxhe, '+/-', fluxhe_err
                    
                    efluxhe = gta.like.energyFlux('GRB'+NAME_TGT, erange_extrapolate_shifted[0], erange_extrapolate_shifted[1])
                    efluxhe_err = gta.like.energyFluxError('GRB'+NAME_TGT, erange_extrapolate_shifted[0], erange_extrapolate_shifted[1])
                    print '** Extrapolated energy flux in', int(erange_extrapolate_shifted[0]+0.5), '-', int(erange_extrapolate_shifted[1]+0.5), 'MeV:', efluxhe, '+/-', efluxhe_err
                    print '** Extrapolated counts in', int(erange_hiend_shifted[0]+0.5), '-', int(erange_hiend_shifted[1]+0.5), 'MeV with Prefactor', params_prev[0], ', Index', params_prev[1], ':', npredhe_prev
                    npredhe_prev_err = npredhe_prev * fluxhe_frac_err_prev
                    print '** Extrapolated counts of all sources in', int(erange_hiend_shifted[0]+0.5), '-', int(erange_hiend_shifted[1]+0.5), 'MeV with Prefactor', params_prev[0], ', Index', params_prev[1], ':', npredhe_prev_all, '+/-', npredhe_prev_err
                    nobshe = 0 

                    for (le, loge_low) in enumerate(npout['roi']['log_energies'][:-1]):
                        if loge_low>=log10(erange_hiend_shifted[0]):
                            nobshe += npout['roi']['counts'][le]
                        
                    print '** Observed counts in', int(erange_hiend_shifted[0]+0.5), '-', int(erange_hiend_shifted[1]+0.5), 'MeV:', nobshe

                    # Calc deviation from Poisson fluctualtion, including Gaussian uncertainty of PL fitting
                    eref_hiend = sqrt(erange_hiend_shifted[0]*erange_hiend_shifted[1])
                    # Calc tentative uncertainty of observed count for TS evaluation
                    if fitresult_prev is not None:
                        npar_index_tgt = None
                        for (ipar, name_par) in enumerate(fitresult_prev['par_names']):
                            if fitresult_prev['src_names'][ipar][3:]==NAME_TGT and name_par=='Index':
                                npar_index_tgt = ipar
                                continue
                        nobshe_sigma = npredhe_prev * sqrt(pow(fluxhe_frac_err_prev,2) + pow(nobshe/npredhe_prev*(log10(eref_hiend)-log10(params_prev[2])) ,2) * fitresult_prev['covariance'][npar_index_tgt][npar_index_tgt])
                        print 'Tentative uncertainty of observed count ({0}): {1}'.format(nobshe, nobshe_sigma)
                        npredhe_sigma = npredhe_prev * sqrt(pow(fluxhe_frac_err_prev,2) + pow(log10(eref_hiend)-log10(params_prev[2] ), 2) * fitresult_prev['covariance'][npar_index_tgt][npar_index_tgt])
                        print 'Check of uncertainty of predicted count ({0}): {1}'.format(npredhe_prev_err, npredhe_sigma)
                    else:
                        nobshe_sigma = nobshe*fluxhe_frac_err_prev

                    # Use simpy
                    (integral_PGauss_mu, integral_PGauss_mu_err)  = integrate.quad(compute_GPoisson, 0, npredhe_prev_all+5.*npredhe_prev_err, args=(npredhe_prev_all, npredhe_prev_err, nobshe))
                    if nobshe>0:
                        (integral_PGauss_n, integral_PGauss_n_err)  = integrate.quad(compute_GPoisson, 0, nobshe+5.*nobshe*fluxhe_frac_err_prev, args=(nobshe, nobshe_sigma, nobshe))
                    else:
                        (integral_PGauss_n, integral_PGauss_n_err)  = (1, 0)
                    print 'Integral of modified Gaussian:', integral_PGauss_mu, '+/-', integral_PGauss_mu_err
                    print 'Integral of modified Gaussian for mu=n:', integral_PGauss_n, '+/-', integral_PGauss_n_err
                    if integral_PGauss_mu<=0:
                        print 'Not positive!!!'
                        if ieedges>0:
                            return 1
                        else:
                            integral_PGauss_mu = 0.001
                    elif integral_PGauss_mu_err>integral_PGauss_mu/100. and ieedges>0:
                        print 'Uncertainty is too large!!!'
                        if ieedges>0:
                            return 1
                        else:
                            integral_PGauss_mu = 0.001
                    deviation = 2. * ( log(integral_PGauss_n / integral_PGauss_mu) )
                    # if nobshe==0:
                    #     deviation = 2. * ( log(integral_PGauss_n) -log(integral_PGauss_mu) )
                    # else:
                    #     deviation = 2. * ( nobshe*log(nobshe) - nobshe + log(integral_PGauss_n) -log(integral_PGauss_mu) )
                    sign_deviation = int(nobshe>=npredhe_prev_all)*2-1
                    print '** TS of Deviation in', int(erange_hiend_shifted[0]+0.5), '-', int(erange_hiend_shifted[1]+0.5), 'MeV from power-law with Prefactor', params_prev[0], ', Index', params_prev[1], ':', deviation
                    print ''

                    for (ipar, par) in enumerate(src_model['param_names']):
                        print par, ':', src_model['param_values'][ipar], '+/-', src_model['param_errors'][ipar]
                        if ipar<len(params_prev):
                            if src_model['param_values'][ipar] == src_model['param_values'][ipar]:
                                params_prev[ipar] = src_model['param_values'][ipar]
                    if src_model['ts']<9 and ieedges==0:
                        print 'TS is not enough!!!'
                        return 1
                    if flux_fracerr_prev is not None:
                        #ax_cspec.plot(x_cspec_extrapolated_eref, npout['roi']['counts'], 'D', label='Observed')
                        ax_cspec.errorbar(x_cspec_extrapolated_eref, npout['roi']['counts'], xerr=(x_cspec_extrapolated_eref_errlo, x_cspec_extrapolated_eref_errhi), fmt='D', label='Observed')
                        ax_cspec.plot(x_cspec_extrapolated_eref, cspec_extrapolated_prev, 'g', label='Sum of all models')
                        ax_cspec.plot(x_cspec_extrapolated_eref, cspec_extrapolated_tgt_prev, 'c', ls='--', label='Model of GRB{0}'.format(NAME_TGT))
                        ax_cspec.plot(x_cspec_extrapolated_eref, cspec_extrapolated_other_prev, 'm', ls=':', label='Model of other sources'.format(NAME_TGT))
                        #sigma1_hi = []
                        #sigma1_lo = []
                        #cspec_extrapolated_err_total = []
                        # for (im, m) in enumerate(cspec_extrapolated_prev):
                        #     n = Symbol('n')
                        #     sols = solve(m-n+n*ln(n/m)-0.5, n)
                        #     if len(sols)==1 and sols[0].is_Float==True:
                        #         sigma1 = sols[0] - m
                        #         if sigma1>0:
                        #             sigma1_hi.append(sols[0])
                        #             sigma1_lo.append(2.*m-sols[0])
                        #             cspec_extrapolated_err_total.append(sqrt(pow(sigma1,2) + pow(cspec_extrapolated_err_prev[im],2)))
                        #         else:
                        #             print 'Standard deviation is not positive!!!'
                        #             return 1
                        #     else:
                        #         print 'Solution:', sols, '!!!!'
                        #         return 1
                        # sigma1_hi = np.array(sigma1_hi)
                        # sigma1_lo = np.array(sigma1_lo)
                        # cspec_extrapolated_err_total = np.array(cspec_extrapolated_err_total)
                        # ax_cspec.fill_between(x_cspec_extrapolated_eref, sigma1_lo, sigma1_hi, alpha=0.75, color='g')
                        ax_cspec.fill_between(x_cspec_extrapolated_eref, cspec_extrapolated_prev+cspec_extrapolated_err_prev, cspec_extrapolated_prev-cspec_extrapolated_err_prev, alpha=0.2, color='g', label='Fitting uncertainty')
                        ax_cspec.grid(c='black', ls='--', lw=0.5, alpha=0.5)
                        ax_cspec.set_xlim(100, 100000)
                        #ax_cspec.fill_between(x=[100,5623], y1=0, y2=100, color='c', alpha=0.2, label='for count deviation')
                        #ax_cspec.fill_between(x=[10000,100000], y1=0, y2=100, color='r', alpha=0.2, label='for fitting')
                        ax_cspec.legend(loc=1, fontsize=12, fancybox=True, framealpha=0.5)
                        ax_cspec.set_xscale('log')
                        ax_cspec.set_yscale('log')
                        ax_cspec.set_xlabel('Energy [MeV]')
                        ax_cspec.set_ylabel('[counts]')
                        ax_cspec.set_title('RoI of GRB'+NAME_TGT)
                        fig_cspec.savefig("{0}/Count_spectrum_{1}{2}.png".format(path_subdir, NAME_TGT, SUFFIX))

                    flux_prev = []
                    flux_err_prev = []
                    flux_fracerr_prev = []
                    for me in range(len(x_cspec_extrapolated_eref)):
                        flux_prev.append(gta.like.flux('GRB'+NAME_TGT, x_cspec_extrapolated_logemin[me], x_cspec_extrapolated_logemax[me]))
                        flux_err_prev.append(gta.like.fluxError('GRB'+NAME_TGT, x_cspec_extrapolated_logemin[me], x_cspec_extrapolated_logemax[me]))
                        flux_fracerr_prev.append(flux_err_prev[-1]/flux_prev[-1])
                    fitresult_prev = copy.deepcopy(fitresult)

                    # Cutoff
                    ecutoff = 0
                    ecutoff_err_lo = 0
                    ecutoff_err_hi = 0
                    ecutoff_ul95 = 0
                    ecutoff_ll95 = 0
                    ecutoff_ul68 = 0
                    ecutoff_ll68 = 0
                    if fspec=='ExpCutoff':
                        gta.free_sources(free=False)
                        gta.free_source('GRB'+NAME_TGT, free=True, pars=['Prefactor', 'Index', 'P1'])
                        val = gta.like.model[gta.like.par_index('GRB'+NAME_TGT, 'P1')].getValue()
                        xvals = 10 ** np.linspace(2.0, 6.0, 101)
                        if val < xvals[0]:
                            xvals = np.insert(xvals, val, 0)
                        prof_ecutoff = profile(gta, 'GRB'+NAME_TGT, 'P1', reoptimize=True, xvals=xvals)
                        dloglike_ecutoff = prof_ecutoff['loglike'] - dct_loglike[fspec]
                        ecutoff_limits95 = get_parameter_limits(xvals, dloglike_ecutoff, 0.95)
                        ecutoff = ecutoff_limits95['x0']
                        ecutoff_err_lo = ecutoff_limits95['err_lo']
                        ecutoff_err_hi = ecutoff_limits95['err_hi']
                        ecutoff_ul95 = ecutoff_limits95['ul']
                        ecutoff_ll95 = ecutoff_limits95['ll']

                        ecutoff_limits68 = get_parameter_limits(xvals, dloglike_ecutoff, 0.68)
                        ecutoff = ecutoff_limits68['x0']
                        ecutoff_err_lo = ecutoff_limits68['err_lo']
                        ecutoff_err_hi = ecutoff_limits68['err_hi']
                        ecutoff_ul68 = ecutoff_limits68['ul']
                        ecutoff_ll68 = ecutoff_limits68['ll']

                        print '** Likelihood for cutoff energy values:'
                        print '  P1:', val
                        for (ix, x) in enumerate(prof_ecutoff['xvals']):
                            print '  Ecutoff:', x, 'Delta loglikelihood:', dloglike_ecutoff[ix]
                        print '  95% Limits', ecutoff_limits95
                        print '  68% Limits', ecutoff_limits68
                    print 'SED with adjusted energy bins.'
                    print erange_sed_shifted
                        #gta_cloned = gta.clone(gta.config)
                    if sedadjusted==True:
                        if index_fixed == None:
                            sed_ad  = gta.sed('GRB'+NAME_TGT, prefix='AD', use_local_index=True, make_plots=True, outfile='sed_GRB{0}{1}_ad.fits'.format(NAME_TGT, SUFFIX), loge_bins=erange_sed_shifted) #[2., 2.5, 3., 3.5, 4., 4.5, 5.]) #erange_sed_shifted)
                        else:
                            sed_ad  = gta.sed('GRB'+NAME_TGT, prefix='AD', bin_index=index_fixed, make_plots=True, outfile='sed_GRB{0}{1}_ad.fits'.format(NAME_TGT, SUFFIX), loge_bins=erange_sed_shifted)
                    else:
                        if index_fixed == None:
                            sed_ad  = gta.sed('GRB'+NAME_TGT, prefix='AD', use_local_index=True, make_plots=True, outfile='sed_GRB{0}{1}_ad.fits'.format(NAME_TGT, SUFFIX))
                        else:
                            sed_ad  = gta.sed('GRB'+NAME_TGT, prefix='AD', bin_index=index_fixed, make_plots=True, outfile='sed_GRB{0}{1}_ad.fits'.format(NAME_TGT, SUFFIX))

                    if skipresid==False:
                        model_resid = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
                        maps = gta.residmap(strenergies, model=model_resid, make_plots=True)

                    for (je, e_ref) in enumerate(sed_ad['e_ref']):
                        print '** SED', int(0.5+sed_ad['e_min'][je]), '-', int(0.5+sed_ad['e_max'][je]), 'MeV'
                        print '* TS:', sed_ad['ts'][je]
                        norm_sed_lims95 = get_parameter_limits(sed_ad['norm_scan'][je], sed_ad['dloglike_scan'][je], 0.95)
                        norm_sed_lims68 = get_parameter_limits(sed_ad['norm_scan'][je], sed_ad['dloglike_scan'][je], 0.68)
                        print '* Flux:', sed_ad['flux'][je], '+', sed_ad['flux_err_hi'][je], '-', sed_ad['flux_err_lo'][je], '(UL:', sed_ad['flux_ul95'][je], ')'
                        print '  Max likelihood:', sed_ad['ref_flux'][je]*norm_sed_lims95['x0']
                        print '  95% limits:', sed_ad['ref_flux'][je]*norm_sed_lims95['ll'], '-', sed_ad['ref_flux'][je]*norm_sed_lims95['ul'] 
                        print '  68% limits:', sed_ad['ref_flux'][je]*norm_sed_lims68['ll'], '-', sed_ad['ref_flux'][je]*norm_sed_lims68['ul'] 
                        print '* Energy flux:', sed_ad['eflux'][je], '+', sed_ad['eflux_err_hi'][je], '-', sed_ad['eflux_err_lo'][je], '(UL:', sed_ad['eflux_ul95'][je], ')'
                        print '  Max likelihood:', sed_ad['ref_eflux'][je]*norm_sed_lims95['x0']
                        print '  95% limits:', sed_ad['ref_eflux'][je]*norm_sed_lims95['ll'], '-', sed_ad['ref_eflux'][je]*norm_sed_lims95['ul']
                        print '  68% limits:', sed_ad['ref_eflux'][je]*norm_sed_lims68['ll'], '-', sed_ad['ref_eflux'][je]*norm_sed_lims68['ul']
                    #print 'SED with equivalent energy bins.'
                    #sed_eq  = gta.sed('GRB'+NAME_TGT, prefix='EQ', use_local_index=True, make_plots=True, outfile='sed_GRB{0}{1}_eq.fits'.format(NAME_TGT, SUFFIX))
                    #if fspec=='PL':
                    str_lc += """{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47},{48},{49},{50},{51},{52},{53},{54},{55},{56},{57},{58},{59},{60}
""".format(NAME_TGT, fspec, LST_RAN_TIME[itime][0], LST_RAN_TIME[itime][1], eranges[ieedges][0], eranges[ieedges][1], eedges[0], eedges[1], src_model['ts'], src_model['param_values'][0], src_model['param_errors'][0], src_model['param_values'][0]*norm_lims95['ul']/norm_lims95['x0'], src_model['param_values'][0]*norm_lims95['ll']/norm_lims95['x0'], src_model['param_values'][0]*norm_lims68['ul']/norm_lims68['x0'], src_model['param_values'][0]*norm_lims68['ll']/norm_lims68['x0'], zero_for_except(src_model['param_values'][1]), zero_for_except(src_model['param_errors'][1]), zero_for_except(src_model['param_values'][2]), zero_for_except(src_model['param_errors'][2]), zero_for_except(src_model['param_values'][3]), zero_for_except(src_model['param_errors'][3]), src_model['flux'], src_model['flux_err'], src_model['flux_ul95'], scale_limit(src_model['flux'], norm_lims95['ll'], norm_lims95['x0']), scale_limit(src_model['flux'], norm_lims68['ul'], norm_lims68['x0']), scale_limit(src_model['flux'], norm_lims68['ll'], norm_lims68['x0']), src_model['eflux'], src_model['eflux_err'], src_model['eflux_ul95'], scale_limit(src_model['eflux'], norm_lims95['ll'], norm_lims95['x0']), scale_limit(src_model['eflux'], norm_lims68['ul'], norm_lims68['x0']), scale_limit(src_model['eflux'], norm_lims68['ll'], norm_lims68['x0']), fluxhe, fluxhe_err, efluxhe, efluxhe_err, sed_ad['ref_flux'][je]*norm_sed_lims68['x0'], sed_ad['ref_flux'][je]*(norm_sed_lims68['ul']-norm_sed_lims68['x0']), sed_ad['ref_flux'][je]*(norm_sed_lims68['x0']-norm_sed_lims68['ll']), sed_ad['ref_flux'][je]*norm_sed_lims95['ul'], sed_ad['ref_flux'][je]*norm_sed_lims95['ll'], sed_ad['ref_eflux'][je]*norm_sed_lims68['x0'], sed_ad['ref_eflux'][je]*(norm_sed_lims68['ul']-norm_sed_lims68['x0']), sed_ad['ref_eflux'][je]*(norm_sed_lims68['x0']-norm_sed_lims68['ll']), sed_ad['ref_eflux'][je]*norm_sed_lims95['ul'], sed_ad['ref_eflux'][je]*norm_sed_lims95['ll'], npredhe_prev, npredhe_prev_err, npredhe_prev_all, nobshe, deviation, sign_deviation, ecutoff, ecutoff_err_lo, ecutoff_err_hi, ecutoff_ul95, ecutoff_ll95, ecutoff_ul68, ecutoff_ll68, dct_loglike[fspec]) 
            
                fluxhe_frac_err_prev = fluxhe_err / fluxhe
                    #continue

            if 'PL' in lst_spec_func:
                try:
                    if 'BPL' in lst_spec_func:
                        ts_ebreak = 2*(dct_loglike['BPL'] - dct_loglike['PL'])
                        print 'TS of PL fit respect to BPL:', ts_ebreak
                        p_ebreak = ROOT.TMath.Prob(ts_ebreak, 2)
                        print 'p-value of PL fit respect to BPL:', p_ebreak
                    if 'ExpCutoff' in lst_spec_func:
                        ts_ecut = 2*(dct_loglike['ExpCutoff'] - dct_loglike['PL'])
                        print 'TS of PL fit respect to ExpCutoff:', ts_ecut
                        p_ecut = ROOT.TMath.Prob(ts_ecut, 1)
                        print 'p-value of PL fit respect to ExpCutoff:', p_ecut
                except KeyError:
                    print 'Comparable analysis is not done.'

    #withopt = 'a'
    #if ieedges==0 and ispec==0:
    #    withopt = 'w'
    with open("{0}/Summary_GRB{1}_{2}.csv".format(path_outdir, NAME_TGT, str_suffix_csv), 'w') as text:
        print str_lc
        text.write(str_lc)


@click.command()
@click.argument('name', type=str)
@click.option('--tmin', type=float, default=None)
@click.option('--tmax', type=float, default=None)
@click.option('-s', '--suffix', type=str, default='')
@click.option('-f', '--force', is_flag=True)
@click.option('--skipts', is_flag=True)
@click.option('--skipsed', is_flag=True)
@click.option('--skipresid', is_flag=True)
@click.option('--emin', default=0, type=float)
@click.option('--emax', default=0, type=float)
@click.option('--nebindecade', default=0)
@click.option('--tbinedges', '-t', multiple=True, default=None, type=float)
@click.option('--outpath', '-o', default=None)
@click.option('--mode', '-m', type=click.Choice(['prompt', 'afterglow', 'unified', 'earlyAG', 'lateAG']))
@click.option('--catalogues', '-c', multiple=True, default=None, type=str)
@click.option('--goodstat', '-g', type=int, default=0)
@click.option('--shiftenergies', is_flag=True)
@click.option('--edisp', is_flag=True)
#@click.option('--bpl', '-b', is_flag=True)
@click.option('--func', multiple=True, default=None, type=str)
@click.option('--fixindex', default=None, type=float, help='Fix the spectral index of power-law to the assigned value.')
@click.option('--fixecut', default=None, type=float, help='Fix the break energy of power-law with exponential cutoff to the assigned value.')
#@click.option('--eshift', '-e', type=click.Choice(['fixed', 'shifted', 'both']))
@click.option('--roi', '-r', type=float, default=12)
@click.option('--reftable', type=str, default='/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LAT2CATALOG-v1-LTF.fits')
@click.option('--download', is_flag=True)
@click.option('--sedadjusted', is_flag=True)
def main(name, tmin, tmax, tbinedges, suffix, force, skipts, skipsed, skipresid, emin, emax, nebindecade, outpath, mode, catalogues, goodstat, edisp, shiftenergies, func, fixindex, fixecut, roi, reftable, download, sedadjusted):
    lst_ebin = []
    #if bpl==True:
    #lst_ebin = [[316.228, 316228.0]]
    if emax==0:
        lst_ebin = [[177.828, 5623.41], [177.828, 100000.0]] #[[177.828, 5623.41], [177.828, 100000.0]] #[[100.0, 5623.41], [100.0, 100000.0]] #[[316.228, 5623.41], [316.228, 100000.0]] #[[100.0, 5623.41], [100.0, 100000.0]] #[[316.228, 177828.0]]
    #    lst_ebin = [[316.228, 316228.0], [316.228, 10000.0], [56234.1, 316228.0]] #[[562.34, 316228.0], [562.34, 10000.0], [56234.0, 316228.0]]
    elif nebindecade>0:
        nebin = int((log10(emax)-log10(emin))*nebindecade+0.5)
        webin = (log10(emax)-log10(emin))/float(nebin)
        for iebin in range(1, 1+nebin):
            lst_ebin.append([pow(10, log10(emin)), pow(10, log10(emin)+iebin*webin)])
    else:
        lst_ebin.append([emin, emax])
    print 'Energy bin edges:', lst_ebin
    ft1_candidates = [None]
    ft2_candidates = [None]
    if outpath == None:
        outpath = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/' + name
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    tbfits = ReadLTFCatalogueInfo.open_table(1, reftable)
    tb_masked = ReadLTFCatalogueInfo.select_one_by_name(tbfits, name)

    path_ft1_exist = ls_list(outpath+'/*_ft1*.fits')[0]
    path_ft2_exist = ls_list(outpath+'/*_ft2*.fits')[0]
    ft1_exist = path_ft1_exist[0]=='/'
    ft2_exist = path_ft2_exist[0]=='/'

    if ft1_exist==True and download==True:
        subprocess.call(['mv', path_ft1_exist, path_ft1_exist.replace('.fits', '_old.fits')])
    if ft2_exist==True and download==True:
        subprocess.call(['mv', path_ft2_exist, path_ft2_exist.replace('.fits', '_old.fits')])
    if ft1_exist==False or download==True:
        print 'Downloading FT1 data...'
        os.chdir(outpath)
        download_fermi_data_grb(name, lst_ft=[1], path_catalogue=reftable, path_outdir=outpath)
    if ft2_exist==False or download==True:
        print 'Downloading FT2 data...'
        os.chdir(outpath)
        download_fermi_data_grb(name, lst_ft=[2], path_catalogue=reftable, path_outdir=outpath)
    else:
        print 'Downloading data is skipped.'

    if mode=='prompt':
        ft1_candidates = ls_list(outpath+'/*_ft1*.fits')
        ft2_candidates = ls_list(outpath+'/*_ft2*.fits')
    elif mode in ('afterglow', 'unified', 'earlyAG', 'lateAG', 'special'):
        ft1_candidates = ls_list(outpath+'/*_ft1*.fits')
        ft2_candidates = ls_list(outpath+'/*_ft2*.fits')
        #ft1_candidates = ls_list(outpath+'/*_PH??.fits')
        #ft2_candidates = ls_list(outpath+'/*_SC??.fits')
    lst_assum_spec = ['PL']
    if not any(func):
        func = ['PL', 'BPL']
    AnalyzeGRB_fermipy(name, ft1_candidates, ft2_candidates, tmin, tmax, tbinedges, suffix, force, skipts, skipsed, skipresid, lst_ebin, tb_masked, outpath, mode, catalogues, goodstat, shiftenergies, edisp, func, fixindex, fixecut, roi, sedadjusted)


if __name__ == '__main__':
    main()
