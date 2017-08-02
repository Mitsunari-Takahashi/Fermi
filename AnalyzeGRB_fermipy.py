#!/usr/bin/env python

import os
import os.path
import subprocess
import matplotlib as mpl
import matplotlib.pyplot as plt
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
from math import log10, sqrt, ceil, isnan
from pLsList import ls_list
from pReadGBMCatalogueInfo import ReadGBMCatalogueOneLine
import click
from FindCrossEarthlimb import find_cross_earthlimb
from FindGoodstatPeriods import find_goodstat_periods, get_entries
import ReadLTFCatalogueInfo
import pMETandMJD
from DownloadFermiData import download_fermi_data_grb


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


def AnalyzeGRB_fermipy(name, ft1_candidates, ft2_candidates, tmin, tmax, tbinedges, suffix, force, skipts, skipsed, skipresid, eranges, tb_masked, path_outdir='.', mode='unified', catalogues=['3FGL'], goodstat=16, shiftenergies=True, edisp=False, lst_spec_func=['PL'], radius_roi=12.) :
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
        tmin = T90_START+T90
    if tmax is None:
        tmax = 10000.
    eranges_shifted = []
    erange_sed_lin = [100.0, 316.228, 1000.0, 3162.28, 10000.0, 31622.8, 100000.0]
    if shiftenergies==True:
        erange_sed_lin = [316.228, 1778.28, 5623.41, 17782.8, 56234.1, 177828.0] #[316.228, 1778.28, 10000.0, 56234.1, 316228.] #[316.228, 1000.0, 3162.28, 10000.0, 56234.1, 316228.]
    erange_sed_shifted_lin = []
    erange_hiend = erange_sed_lin[-2:] #[56234., 316228.]
    erange_hiend_shifted = []
    #if (dct_grb['z'] is not '') and (dct_grb['z'] is not '-') and (dct_grb['z'] is not '?'):
    if tb_masked['REDSHIFT']>0:
        REDSHIFT = tb_masked['REDSHIFT'] #dct_grb['z']
        if shiftenergies==True:
            for eset in eranges:
                eranges_shifted.append([e/(1+REDSHIFT) for e in eset])
            erange_hiend_shifted = [e/(1+REDSHIFT) for e in erange_hiend]
            erange_sed_shifted_lin = [e/(1+REDSHIFT) for e in erange_sed_lin]
        else:
            eranges_shifted = eranges
            erange_hiend_shifted = erange_hiend
            erange_sed_shifted_lin = erange_sed_lin
    else:
        REDSHIFT = 1
        eranges_shifted = eranges
        erange_hiend_shifted = erange_hiend
        erange_sed_shifted_lin = erange_sed_lin
    erange_sed_shifted = np.array([log10(x) for x in erange_sed_shifted_lin])
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
    if mode == 'unified' or mode == 'afterglow':
        lst_tbin = [[tmin, tmax]]
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

    for (ieedges, eedges) in enumerate(eranges_shifted):
        #strenergies = 'E{0:0>6}-{1:0>6}MeV'.format(int(eedges[0]+0.5), int(eedges[1]+0.5))
        strenergies = 'E{0:0>6}-{1:0>6}MeV'.format(int(eranges[ieedges][0]+0.5), int(eranges[ieedges][1]+0.5))
        if shiftenergies==True:
            strenergies += '_shifted'
        print '%%%%%%%%%%%%%%%%%%'
        print strenergies
        print '%%%%%%%%%%%%%%%%%%'
        for itime in range(NRAN_TIME):
            strtime = 'T{0:0>6}-{1:0>6}s'.format(int(0.5+LST_RAN_TIME[itime][0]), int(0.5+LST_RAN_TIME[itime][1]))
            print '=====', strtime, '====='
            dct_loglike = {}
            for (ispec, fspec) in enumerate(lst_spec_func):
                print '-----', fspec, '-----'
                path_subdir = '{0}/{1}/{2}/{3}'.format(path_outdir, strenergies, strtime, fspec)
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
                        str_spectrum = """'SpectrumType' : 'PowerLaw', 'Prefactor' : { value : 1.0, max : !!float 1e6, min : !!float 1e-3, scale : !!float 1e-6, free : '1' }, 'Index' : { value : 2.0, min : -1.0, max : 8.0, scale : -1, free : '1' }, 'Scale' : { value : 100.0, max : 100000., min : 30, scale : 1, free : '0' }"""
                    elif fspec=='BPL':
                        str_spectrum = """'SpectrumType' : 'BrokenPowerLaw', 'Prefactor' : { value : 1.0, max : !!float 1e6, min : !!float 1e-3, scale : !!float 1e-6, free : '1' }, 'Index1' : { value : 2.0, min : -1.0, max : 8.0, scale : -1, free : '1' }, 'Index2' : { value : 1.6, min : -1.0, max : 8.0, scale : -1, free : '1' }, 'BreakValue' : { value : 5000, min : 500, max : 50000, scale : 1, free : '1' }"""
                    if fspec=='EblPL':
                        str_spectrum = """'SpectrumType' : 'EblAtten::PowerLaw2', 'Integral' : {{ value : 1.0, max : !!float 1e6, min : !!float 1e-6, scale : !!float 1e-6, free : '1' }}, 'Index' : {{ value : 2.0, min : -1.0, max : 8.0, scale : -1, free : '1' }}, 'LowerLimit' : {{ value : {0}, max : 100000., min : 30, scale : 1, free : '0' }}, 'UpperLimit' : {{ value : {1}, max : 1000000., min : 100., scale : 1, free : '0' }}, 'tau_norm' : {{ value : 1.0, max : 10, min : 0, scale : 1.0, free : '0' }}, 'redshift' : {{ value : {2}, max : 10, min : 0, scale : 1, free : '0' }}, 'ebl_model' : {{ value : 4, max : 8, min : 0, scale : 1.0, free : '0'}}""".format(eedges[0], eedges[1], REDSHIFT)
                    elif fspec=='EblBPL':
                        str_spectrum = """'SpectrumType' : 'EblAtten::BrokenPowerLaw2', 'Integral' : {{ value : 1.0, max : !!float 1e6, min : !!float 1e-6, scale : !!float 1e-6, free : '1' }}, 'Index1' : {{ value : 2.0, min : -1.0, max : 8.0, scale : -1, free : '1' }}, 'Index2' : {{ value : 1.6, min : -1.0, max : 8.0, scale : -1, free : '1' }}, 'BreakValue' : {{ value : 5000, min : 500, max : 50000, scale : 1, free : '1' }}, 'LowerLimit' : {{ value : {0}, max : 100000., min : 30, scale : 1, free : '0' }}, 'UpperLimit' : {{ value : {1}, max : 1000000., min : 100., scale : 1, free : '0' }}, 'tau_norm' : {{ value : 1.0, max : 10, min : 0, scale : 1.0, free : '0' }}, 'redshift' : {{ value : {2}, max : 10, min : 0, scale : 1, free : '0' }}, 'ebl_model' : {{ value : 4, max : 8, min : 0, scale : 1.0, free : '0'}}""".format(eedges[0], eedges[1], REDSHIFT)

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

                    gta.setup()
                    gta.optimize()
                    gta.print_roi()
            
                    gta.free_sources(free=False)
                    if fspec=='PL':
                        gta.free_source('GRB'+NAME_TGT, free=True, pars=['Prefactor', 'Index'])
                    elif fspec=='BPL':
                        gta.free_source('GRB'+NAME_TGT, free=True, pars=['Prefactor', 'Index1', 'Index2', 'BreakValue'])
                    elif fspec=='EblPL':
                        gta.free_source('GRB'+NAME_TGT, free=True, pars=['Integral', 'Index'])
                    elif fspec=='EblBPL':
                        gta.free_source('GRB'+NAME_TGT, free=True, pars=['Integral', 'Index1', 'Index2', 'BreakValue'])

                    print 'Fitting...'
                    fitresult = gta.fit()
                    gta.write_roi('fit_model'+SUFFIX)
                    gta.print_roi()
                    print ' Fitting finished.'
                
                    dct_loglike[fspec] = fitresult['loglike']
                    print '** log(likelihood):', dct_loglike[fspec]
                    src_model = gta.get_src_model('GRB'+NAME_TGT)
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
                    fluxhe = gta.like.flux('GRB'+NAME_TGT, erange_hiend_shifted[0], erange_hiend_shifted[1])
                    fluxhe_err = gta.like.fluxError('GRB'+NAME_TGT, erange_hiend_shifted[0], erange_hiend_shifted[1])
                    print '** Extrapolated flux in', int(erange_hiend_shifted[0]+0.5), '-', int(erange_hiend_shifted[1]+0.5), 'MeV:', fluxhe, '+/-', fluxhe_err
                    efluxhe = gta.like.energyFlux('GRB'+NAME_TGT, erange_hiend_shifted[0], erange_hiend_shifted[1])
                    efluxhe_err = gta.like.energyFluxError('GRB'+NAME_TGT, erange_hiend_shifted[0], erange_hiend_shifted[1])
                    print '** Extrapolated energy flux in', int(erange_hiend_shifted[0]+0.5), '-', int(erange_hiend_shifted[1]+0.5), 'MeV:', efluxhe, '+/-', efluxhe_err
                    print ''
                    for (ipar, par) in enumerate(src_model['param_names']):
                        print par, ':', src_model['param_values'][ipar], '+/-', src_model['param_errors'][ipar]
                    #if ieedges==0:
                    print 'SED with adjusted energy bins.'
                    print erange_sed_shifted
                        #gta_cloned = gta.clone(gta.config)
                    sed_ad  = gta.sed('GRB'+NAME_TGT, prefix='AD', use_local_index=True, make_plots=True, outfile='sed_GRB{0}{1}_ad.fits'.format(NAME_TGT, SUFFIX), loge_bins=erange_sed_shifted) #[2., 2.5, 3., 3.5, 4., 4.5, 5.]) #erange_sed_shifted)
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
                    str_lc = """#name:function:start:stop:emin_rest:emax_rest:emin_shifted:emax_shifted:ts:Integral:Integral_err:Integral_ul95:Integral_ll95:Integral_ul68:Integral_ll68:Index1:Index1_err:Index2:Index2_err:BreakValue:BreakValue_err:flux:flux_err:flux_ul95:flux_ll95:flux_ul68:flux_ll68:eflux:eflux_err:eflux_ul95:eflux_ll95:eflux_ul68:eflux_ll68:fluxhe:fluxhe_err:efluxhe:efluxhe_err:flux_hiest_e:flux_hiest_e_err_hi:flux_hiest_e_err_lo:flux_hiest_e_ul:flux_hiest_e_ll:eflux_hiest_e:eflux_hiest_e_err_hi:eflux_hiest_e_err_lo:eflux_hiest_e_ul:eflux_hiest_e_ll:loglike
{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16},{17},{18},{19},{20},{21},{22},{23},{24},{25},{26},{27},{28},{29},{30},{31},{32},{33},{34},{35},{36},{37},{38},{39},{40},{41},{42},{43},{44},{45},{46},{47}
""".format(NAME_TGT, fspec, LST_RAN_TIME[itime][0], LST_RAN_TIME[itime][1], eranges[ieedges][0], eranges[ieedges][1], eedges[0], eedges[1], src_model['ts'], src_model['param_values'][0], src_model['param_errors'][0], src_model['param_values'][0]*norm_lims95['ul']/norm_lims95['x0'], src_model['param_values'][0]*norm_lims95['ll']/norm_lims95['x0'], src_model['param_values'][0]*norm_lims68['ul']/norm_lims68['x0'], src_model['param_values'][0]*norm_lims68['ll']/norm_lims68['x0'], zero_for_except(src_model['param_values'][1]), zero_for_except(src_model['param_errors'][1]), zero_for_except(src_model['param_values'][2]), zero_for_except(src_model['param_errors'][2]), zero_for_except(src_model['param_values'][3]), zero_for_except(src_model['param_errors'][3]), src_model['flux'], src_model['flux_err'], src_model['flux_ul95'], scale_limit(src_model['flux'], norm_lims95['ll'], norm_lims95['x0']), scale_limit(src_model['flux'], norm_lims68['ul'], norm_lims68['x0']), scale_limit(src_model['flux'], norm_lims68['ll'], norm_lims68['x0']), src_model['eflux'], src_model['eflux_err'], src_model['eflux_ul95'], scale_limit(src_model['eflux'], norm_lims95['ll'], norm_lims95['x0']), scale_limit(src_model['eflux'], norm_lims68['ul'], norm_lims68['x0']), scale_limit(src_model['eflux'], norm_lims68['ll'], norm_lims68['x0']), fluxhe, fluxhe_err, efluxhe, efluxhe_err, sed_ad['ref_flux'][je]*norm_sed_lims95['x0'], sed_ad['ref_flux'][je]*norm_sed_lims68['ul'], sed_ad['ref_flux'][je]*norm_sed_lims68['ll'], sed_ad['ref_flux'][je]*norm_sed_lims95['ul'], sed_ad['ref_flux'][je]*norm_sed_lims95['ll'], sed_ad['ref_eflux'][je]*norm_sed_lims95['x0'], sed_ad['ref_eflux'][je]*norm_sed_lims68['ul'], sed_ad['ref_eflux'][je]*norm_sed_lims68['ll'], sed_ad['ref_eflux'][je]*norm_sed_lims95['ul'], sed_ad['ref_eflux'][je]*norm_sed_lims95['ll'], dct_loglike[fspec]) #.format(NAME_TGT, fspec, LST_RAN_TIME[itime][0], LST_RAN_TIME[itime][1], eranges[ieedges][0], eranges[ieedges][1], eedges[0], eedges[1], src_model['ts'], src_model['param_values'][0], src_model['param_errors'][0], src_model['param_values'][0]*norm_lims95['ul']/norm_lims95['x0'], src_model['param_values'][0]*norm_lims95['ll']/norm_lims95['x0'], src_model['param_values'][0]*norm_lims68['ul']/norm_lims68['x0'], src_model['param_values'][0]*norm_lims68['ll']/norm_lims68['x0'], zero_for_except(src_model['param_values'][1]), zero_for_except(src_model['param_errors'][1]), zero_for_except(src_model['param_values'][2]), zero_for_except(src_model['param_errors'][2]), zero_for_except(src_model['param_values'][3]), zero_for_except(src_model['param_errors'][3]), src_model['flux'], src_model['flux_err'], src_model['flux_ul95'], scale_limit(src_model['flux'], norm_lims95['ll'], norm_lims95['x0']), scale_limit(src_model['flux'], norm_lims68['ul'], norm_lims68['x0']), scale_limit(src_model['flux'], norm_lims68['ll'], norm_lims68['x0']), src_model['eflux'], src_model['eflux_err'], src_model['eflux_ul95'], scale_limit(src_model['eflux'], norm_lims95['ll'], norm_lims95['x0']), scale_limit(src_model['eflux'], norm_lims68['ul'], norm_lims68['x0']), scale_limit(src_model['eflux'], norm_lims68['ll'], norm_lims68['x0']), fluxhe, fluxhe_err, efluxhe, efluxhe_err, sed_ad['flux'][-1], sed_ad['flux_err_hi'][-1], sed_ad['flux_err_lo'][-1], sed_ad['flux_ul95'][-1], sed_ad['ref_flux'][-1]*norm_sed_lims68['ul'], sed_ad['eflux'][-1], sed_ad['eflux_err_hi'][-1], sed_ad['eflux_err_lo'][-1], sed_ad['eflux_ul95'][-1], sed_ad['ref_eflux'][-1]*norm_lims68['ul'], dct_loglike[fspec])
                
                    withopt = 'a'
                    if ieedges==0 and ispec==0:
                        withopt = 'w'
                    with open("{0}/GRB{1}_{2}_lc_summary.csv".format(path_outdir, NAME_TGT, strtime), withopt) as text:
                        print str_lc
                        text.write(str_lc)
                    #continue
            try:
                ts_ebreak = 2*(dct_loglike['BPL'] - dct_loglike['PL'])
                print 'TS of PL fit respect to BPL:', ts_ebreak
                p_ebreak = ROOT.TMath.Prob(ts_ebreak, 2)
                print 'p-value of PL fit respect to BPL:', p_ebreak
            except KeyError:
                print 'Comparable analysis is not done.'
    
    #             print '===== Fitting parameters ====='
    #             c = np.load('{0}/fit_model{1}.npy'.format(path_subdir, SUFFIX)).flat[0]
    #             np_src = c['sources']['GRB'+NAME_TGT]
    #             norm_lims95 = get_parameter_limits(np_src['norm_scan'], np_src['loglike_scan'], 0.95)
    #             norm_lims68 = get_parameter_limits(np_src['norm_scan'], np_src['loglike_scan'], 0.68)
    #             print '* Integral:', np_src['param_values'][0], '+/-', np_src['param_errors'][0]
    #             print '  95% limits:', np_src['param_values'][0]*norm_lims95['ll']/norm_lims95['x0'], '-', np_src['param_values'][0]*norm_lims95['ul']/norm_lims95['x0']
    #             print '  68% limits:', np_src['param_values'][0]*norm_lims68['ll']/norm_lims95['x0'], '-', np_src['param_values'][0]*norm_lims68['ul']/norm_lims95['x0']
    #             print '* Flux:', np_src['flux'], '+/-', np_src['flux_err']
    #             print '  95% limits:', np_src['flux']*norm_lims95['ll']/norm_lims95['x0'], '-', np_src['flux']*norm_lims95['ul']/norm_lims95['x0']
    #             print '  68% limits:', np_src['flux']*norm_lims68['ll']/norm_lims95['x0'], '-', np_src['flux']*norm_lims68['ul']/norm_lims95['x0']
    #             print '* Energy flux:', np_src['eflux'], '+/-', np_src['eflux_err']
    #             print '  95% limits:', np_src['eflux']*norm_lims95['ll']/norm_lims95['x0'], '-', np_src['eflux']*norm_lims95['ul']/norm_lims95['x0']
    #             print '  68% limits:', np_src['eflux']*norm_lims68['ll']/norm_lims95['x0'], '-', np_src['eflux']*norm_lims68['ul']/norm_lims95['x0']
    #             print '* Index:', np_src['param_values'][1], '+/-', np_src['param_errors'][1]
    #             print ''
    #             for (ipar, par) in enumerate(np_src['param_names']):
    #                 print par, ':', np_src['param_values'][ipar], '+/-', np_src['param_errors'][ipar]
    # #TS map
    #             if skipts==False:
    #                 model_ts = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
    #                 maps_ts = gta.tsmap('fit_ts',model=model_ts, make_plots=True)
    #                 gta.plotter.make_tsmap_plots(maps_ts, roi=gta.roi)                
    # #SED
    #             if skipsed==False:
    #                 if c['sources']['GRB'+NAME_TGT]['npred']>0:
    #                     print 'Going to SED analysis...'
    #                     print 'SED with equivalent energy bins.'
    #                     sed_eq  = gta.sed('GRB'+NAME_TGT, prefix='EQ', use_local_index=True, make_plots=True, outfile='sed_GRB{0}{1}_eq.fits'.format(NAME_TGT, SUFFIX))

    # # Residual map
    #             if skipresid==False:
    #                 model_resid = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
    #                 maps_resid = gta.residmap('fit_resid',model=model_resid, make_plots=True)
    #                 gta.plotter.make_residmap_plots(maps_resid, roi=gta.roi)


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
@click.option('--mode', '-m', type=click.Choice(['prompt', 'afterglow', 'unified']))
@click.option('--catalogues', '-c', multiple=True, default=None, type=str)
@click.option('--goodstat', '-g', type=int, default=0)
@click.option('--shiftenergies', is_flag=True)
@click.option('--edisp', is_flag=True)
#@click.option('--bpl', '-b', is_flag=True)
@click.option('--func', multiple=True, default=None, type=str)
#@click.option('--eshift', '-e', type=click.Choice(['fixed', 'shifted', 'both']))
@click.option('--roi', '-r', type=float, default=12)
@click.option('--reftable', type=str, default='/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LAT2CATALOG-v1-LTF.fits')
@click.option('--download', is_flag=True)
def main(name, tmin, tmax, tbinedges, suffix, force, skipts, skipsed, skipresid, emin, emax, nebindecade, outpath, mode, catalogues, goodstat, edisp, shiftenergies, func, roi, reftable, download):
    lst_ebin = []
    #if bpl==True:
    #lst_ebin = [[316.228, 316228.0]]
    if emax==0:
        lst_ebin = [[316.228, 177828.0]]
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
    elif mode=='afterglow' or mode=='unified':
        ft1_candidates = ls_list(outpath+'/*_ft1*.fits')
        ft2_candidates = ls_list(outpath+'/*_ft2*.fits')
        #ft1_candidates = ls_list(outpath+'/*_PH??.fits')
        #ft2_candidates = ls_list(outpath+'/*_SC??.fits')
    lst_assum_spec = ['PL']
    if not any(func):
        func = ['PL', 'BPL']
    AnalyzeGRB_fermipy(name, ft1_candidates, ft2_candidates, tmin, tmax, tbinedges, suffix, force, skipts, skipsed, skipresid, lst_ebin, tb_masked, outpath, mode, catalogues, goodstat, shiftenergies, edisp, func, roi)


if __name__ == '__main__':
    main()
