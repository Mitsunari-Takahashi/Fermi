#!/usr/bin/env python

import sys
import os
path_upstairs = os.path.join(os.path.dirname(__file__), '../')
sys.path.append(path_upstairs)
import logging
import numpy as np
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi, sqrt
import click
import matplotlib as mpl
import matplotlib.pyplot as plt
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TTree, TH1, TGraphErrors, kWhite, kBlack, kGray, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink
ROOT.gROOT.SetBatch()
from scipy.optimize import leastsq
import pickle_utilities
import pMatplot
from pLsList import ls_list
import ReadLTFCatalogueInfo
from STLikelihoodAnalysis import get_module_logger


##### Logger #####
logger = get_module_logger(__name__)


##### GBM ######
NCAT_GBM = 3


###### Matplotlib setting ######
mpl.rcParams['font.size'] = 15


def plot_grbs_deviation(lst_path_pickle, path_table, path_out):
    # GBM data
    tb = ReadLTFCatalogueInfo.open_table(1, path_table)
    FLUENCE_CUT = [1.09e-04, 3.06e-05] #[1.45E-4, 3.70E-5]

    # Common numbers
    d0 = pickle_utilities.load(lst_path_pickle[0])
    PHASE = d0['phase']
    EMIN_FIT = d0['lower_energies']['emin']
    EMAX_FIT = d0['lower_energies']['emax']
    EMIN_EXTRA = d0['highest_energies']['emin']
    EMAX_EXTRA = d0['highest_energies']['emax']
    EMIN_WHOLE = d0['whole_energies']['emin']
    EMAX_WHOLE = d0['whole_energies']['emax']

    # string for x/y-axis label
    #str_allgrbs = '{phs} of {nall} long-GRBs'.format(phs=PHASE.capitalize(), nall=lst_ngrb_valid[0])
    str_fluence_gbm = r'GBM fluence $[\mathrm{{erg/cm^2}}]$'
    str_index_lower_energies = 'Index in {emin:3.3f} - {emax:3.3f} GeV'.format(emin=EMIN_FIT/1000, emax=EMAX_FIT/1000)
    str_t90_gbm = 'GBM T90 [s]'
    str_index_whole_energies = 'Index in {emin:3.3f} - {emax:3.3f} GeV'.format(emin=EMIN_FIT/1000, emax=EMAX_EXTRA/1000)
    str_deviation = r'Deviation in {emin:3.0f} - {emax:3.0f} GeV $\mathrm{{[\sigma]}}$'.format(emin=EMIN_EXTRA/1000, emax=EMAX_EXTRA/1000)
    str_flux_whole_energies = r'Flux in {emin:3.3f} - {emax:3.0f} GeV $\mathrm{{[/cm^2 s]}}$'.format(emin=EMIN_FIT/1000, emax=EMAX_EXTRA/1000)
    str_flux_highest_energies = r'Flux in {emin:3.0f} - {emax:3.0f} GeV $\mathrm{{[/cm^2 s]}}$'.format(emin=EMIN_EXTRA/1000, emax=EMAX_EXTRA/1000)
    str_nobs_highest_energies = 'Observed photon number in {emin:3.0f} - {emax:3.0f} GeV [counts]'.format(emin=EMIN_EXTRA/1000, emax=EMAX_EXTRA/1000)
    str_npred_all_highest_energies = 'Predicted photon number in {emin:3.0f} - {emax:3.0f} GeV [counts]'.format(emin=EMIN_EXTRA/1000, emax=EMAX_EXTRA/1000)
    str_flux_highest_energies_per_fluence_gbm = 'Flux in {emin:3.0f} - {emax:3.0f} GeV / GBM fluence'.format(emin=EMIN_EXTRA/1000, emax=EMAX_EXTRA/1000)
    str_ts_lower_energies = 'TS of fitting in {emin:3.3f} - {emax:3.3f} GeV'.format(emin=EMIN_FIT/1000, emax=EMAX_FIT/1000)
    str_ts_highest_energies = 'TS of extrapolated power-law in {emin:3.0f} - {emax:3.0f} GeV'.format(emin=EMIN_EXTRA/1000, emax=EMAX_EXTRA/1000)
    str_lightcurve_index = 'Decay index'

    # Lists for results
    lst_names = []

    lst_fluence_gbm = []
    lst_fluence_gbm_err = []

    lst_t90_gbm = []
    lst_t90_gbm_err = []

    lst_index_lower_energies = []
    lst_index_err_lower_energies = []

    lst_index_whole_energies = []
    lst_index_err_whole_energies = []

    lst_deviation_ts = []

    lst_flux_whole_energies = []
    lst_flux_err_lo_whole_energies = []
    lst_flux_err_hi_whole_energies = []

    lst_flux_highest_energies = []
    lst_flux_ul_highest_energies = []
    lst_flux_err_lo_highest_energies = []
    lst_flux_err_hi_highest_energies = []

    lst_nobs_highest_energies = []
    lst_npred_all_highest_energies = []

    lst_bool_flux_finite_highest_energies = []

    lst_bool_lightcurve = []
    lst_lightcurve_index = []
    lst_lightcurve_index_err = []

    lst_index_scanned = []
    lst_dev_scanned = []
    lst_ts_fit_scanned = []

    for icat in range(NCAT_GBM+1):
        lst_names.append([])

        lst_fluence_gbm.append([])
        lst_fluence_gbm_err.append([])

        lst_t90_gbm.append([])
        lst_t90_gbm_err.append([])

        lst_index_lower_energies.append([])
        lst_index_err_lower_energies.append([])

        lst_index_whole_energies.append([])
        lst_index_err_whole_energies.append([])

        lst_deviation_ts.append([])

        lst_flux_whole_energies.append([])
        lst_flux_err_lo_whole_energies.append([])
        lst_flux_err_hi_whole_energies.append([])

        lst_flux_highest_energies.append([])
        lst_flux_ul_highest_energies.append([])
        lst_flux_err_lo_highest_energies.append([])
        lst_flux_err_hi_highest_energies.append([])

        lst_nobs_highest_energies.append([])
        lst_npred_all_highest_energies.append([])

        lst_bool_flux_finite_highest_energies.append([])

        lst_bool_lightcurve.append([])
        lst_lightcurve_index.append([])
        lst_lightcurve_index_err.append([])

        lst_index_scanned.append([])
        lst_dev_scanned.append([])
        lst_ts_fit_scanned.append([])

    dp_ts_scanned = []

    for igrb, path_pickle in enumerate(lst_path_pickle):
        d = pickle_utilities.load(path_pickle)
        name = d['target']
        logger.info('===== GRB{name} ====='.format(name=name))
        ts = d['lower_energies']['TS']
        if ts<25:
            logger.warning('TS={ts} is NOT enough!!'.format(ts=ts))
            continue

        #GBM category
        tb1 = ReadLTFCatalogueInfo.select_one_by_name(tb, name)
        ncategory = ReadLTFCatalogueInfo.judge_category_fluence(tb, name, FLUENCE_CUT)+1
        logger.info('GBM fluence category: {0}'.format(ncategory))

        # Check outstanding values
        if d['deviation_ts']>2.25 or d['deviation_ts']<-4.0:
            logger.info('TS of deviation: {ts} !!'.format(ts=d['deviation_ts']))
        if d['highest_energies']['flux']['value']/tb1['FLUENCE']>1e-3:
            logger.info('{st}: {flu} +/- {fluer} !!'.format(st=str_flux_highest_energies, flu=d['highest_energies']['flux']['value'], fluer=d['highest_energies']['flux']['error']))
            logger.info(' cf. GBM fluence: {flu}'.format(flu=tb1['FLUENCE']))

        #Light curve data
        path_lc_pickle = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/LongGRBs/{name}/E0000100-0100000MeV/r12deg/lightcurve/LightCurve_{name}_indexfree_fit.pickle'.format(name=name)
        lc = pickle_utilities.load(path_lc_pickle)

        # Fill data
        for icat in range(NCAT_GBM+1):
            if icat in (0, ncategory):
                lst_names[icat].append(name)

                lst_fluence_gbm[icat].append(tb1['FLUENCE'])
                lst_fluence_gbm_err[icat].append(tb1['FLUENCE_ERROR'])

                lst_t90_gbm[icat].append(tb1['T90'])
                lst_t90_gbm_err[icat].append(tb1['T90_ERROR'])

                lst_index_lower_energies[icat].append(d['lower_energies']['Index']['value'])
                lst_index_err_lower_energies[icat].append(d['lower_energies']['Index']['error'])

                lst_index_whole_energies[icat].append(d['whole_energies']['Index']['value'])
                lst_index_err_whole_energies[icat].append(d['whole_energies']['Index']['error'])

                lst_deviation_ts[icat].append(sqrt(abs(d['deviation_ts'])) * (1 if d['deviation_ts']>=0 else -1))

                lst_flux_whole_energies[icat].append(d['whole_energies']['flux']['value'])
                lst_flux_err_lo_whole_energies[icat].append(d['whole_energies']['limits']['best']['flux']['err_lo']-d['whole_energies']['limits']['best']['flux']['x0']+d['whole_energies']['flux']['value'])
                lst_flux_err_hi_whole_energies[icat].append(d['whole_energies']['limits']['best']['flux']['err_hi']-d['whole_energies']['limits']['best']['flux']['x0']+d['whole_energies']['flux']['value'])

                lst_nobs_highest_energies[icat].append(d['highest_energies']['nobs'])
                lst_npred_all_highest_energies[icat].append(d['highest_energies']['npred_all']['value'])

                lst_ts_fit_scanned[icat].append(d['scan']['lower_energies']['TS'])
                xscan = d['scan']['lower_energies']['index']
                yscan = d['scan']['lower_energies']['norm']
                xmesh, ymesh = np.meshgrid(xscan, yscan)
                lst_index_scanned[icat].append(xmesh)
                lst_dev_scanned[icat].append(d['scan']['highest_energies']['TS']*d['scan']['highest_energies']['sign'])
                logger.debug("""
TS of deviation from power-law:""")
                logger.debug(lst_dev_scanned[icat][-1])

                # Load light curve
                if 'index' in lc['fit']['flux']['afterglow']:
                    lst_bool_lightcurve[icat].append(1)
                    lst_lightcurve_index[icat].append(lc['fit']['flux']['afterglow']['index']['value'])
                    lst_lightcurve_index_err[icat].append(lc['fit']['flux']['afterglow']['index']['error'])
                else:
                    lst_bool_lightcurve[icat].append(0)

                if d['highest_energies']['flux']['value']==d['highest_energies']['flux']['value'] and d['highest_energies']['limits']['best']['flux']['err_lo']==d['highest_energies']['limits']['best']['flux']['err_lo'] and d['highest_energies']['limits']['best']['flux']['err_hi']==d['highest_energies']['limits']['best']['flux']['err_hi'] and d['highest_energies']['limits']['best']['flux']['err_hi']+d['highest_energies']['limits']['best']['flux']['x0']<d['highest_energies']['limits']['best']['flux']['ul']:
                    lst_flux_highest_energies[icat].append(d['highest_energies']['flux']['value'])
                    lst_flux_err_lo_highest_energies[icat].append(d['highest_energies']['limits']['best']['flux']['err_lo']-d['highest_energies']['limits']['best']['flux']['x0']+d['highest_energies']['flux']['value'])
                    lst_flux_err_hi_highest_energies[icat].append(d['highest_energies']['limits']['best']['flux']['err_hi']-d['highest_energies']['limits']['best']['flux']['x0']+d['highest_energies']['flux']['value'])
                    lst_bool_flux_finite_highest_energies[icat].append(1)
                elif d['highest_energies']['limits']['best']['flux']['ul']==d['highest_energies']['limits']['best']['flux']['ul']:
                    lst_flux_ul_highest_energies[icat].append(d['highest_energies']['limits']['best']['flux']['ul'])
                    lst_bool_flux_finite_highest_energies[icat].append(0)
                else:
                    logger.warning('GRB{name} does NOT have valid flux value or UL!!'.format(name=name))

    #print lst_flux_highest_energies
    #print lst_flux_err_lo_highest_energies
    #print lst_flux_err_hi_highest_energies
    #print lst_flux_ul_highest_energies
    #print lst_bool_flux_finite_highest_energies

    # Number of valid GRBs
    lst_ngrb_valid = [len(lst_names[i]) for i in range(NCAT_GBM+1)]
    logger.info('Numbers of valid GRBs (all, category 1, 2, 3) : {0}'.format(lst_ngrb_valid))

    # List of ndarrays
    names = []

    indices_lower_energies = []
    indices_err_lower_energies = []

    indices_whole_energies = []
    indices_err_whole_energies = []

    deviations_ts = []

    fluences_gbm = []
    fluences_gbm_err = []
    #fluences_gbm_hiEv = []
    #fluences_gbm_err_hiEv = []
    #fluences_gbm_hiEu = []
    #fluences_gbm_err_hiEu = []

    t90_gbm = []
    t90_gbm_err = []

    fluxes_whole_energies = []
    fluxes_err_lo_whole_energies = []
    fluxes_err_hi_whole_energies = []

    fluxes_highest_energies = []
    fluxes_err_lo_highest_energies = []
    fluxes_err_hi_highest_energies = []
    fluxes_ul_highest_energies = []

    nobs_highest_energies = []
    npred_all_highest_energies = []

    bools_flux_finite_highest_energies = []
    bools_non_flux_finite_highest_energies = []

    bools_lightcurve = []
    bools_non_lightcurve = []
    lightcurve_indices = []
    lightcurve_indices_err = []

    # List of Data_plotted
    dp_deviation_vs_fluence_gbm = []
    dp_deviation_vs_indices_lower_energies = []
    dp_deviation_vs_t90_gbm = []
    dp_deviation = []
    dp_flux_highest_energies_vs_fluence_gbm = []
    dp_flux_ul_highest_energies_vs_fluence_gbm = []
    dp_flux_highest_energies_per_fluence_gbm_vs_fluence_gbm = []
    dp_flux_ul_highest_energies_per_fluence_gbm_vs_fluence_gbm = []
    dp_index_whole_energies = []
    dp_index_whole_energies_vs_t90_gbm = []
    dp_deviation_vs_fluence_gbm_vs_index_lower_energies = []
    dp_index_whole_energies_vs_nobs_highest_energies_vs_flux_whole_energies_vs_fluence_gbm = []
    dp_index_whole_energies_vs_nobs_highest_energies_vs_npred_all_highest_energies = []
    dp_index_whole_energies_vs_fluence_gbm_vs_nobs_highest_energies_vs_npred_all_highest_energies = []
    dp_index_whole_energies_vs_fluence_gbm = []
    dp_flux_whole_energies_vs_fluence_gbm = []
    dp_lightcurve_index = []
    dp_lightcurve_index_vs_spectral_index = []
    dp_flux_highest_energies_vs_lightcurve_index = []
    dp_flux_ul_highest_energies_vs_lightcurve_index = []
#    dp_nobs_highest_energies_vs_npred_all_highest_energies_vs_
    #dp_ts_scanned = []

    for icat in range(NCAT_GBM+1):
        str_category = 'Category {0}'.format(icat) if icat>0 else 'All categories'
        indices_lower_energies.append(np.array(lst_index_lower_energies[icat]))
        indices_err_lower_energies.append(np.array(lst_index_err_lower_energies[icat]))

        indices_whole_energies.append(np.array(lst_index_whole_energies[icat]))
        indices_err_whole_energies.append(np.array(lst_index_err_whole_energies[icat]))

        deviations_ts.append(np.array(lst_deviation_ts[icat]))

        fluences_gbm.append(np.array(lst_fluence_gbm[icat]))
        fluences_gbm_err.append(np.array(lst_fluence_gbm_err[icat]))

        t90_gbm.append(np.array(lst_t90_gbm[icat]))
        t90_gbm_err.append(np.array(lst_t90_gbm_err[icat]))

        fluxes_whole_energies.append(np.array(lst_flux_whole_energies[icat]))
        fluxes_err_lo_whole_energies.append(np.array(lst_flux_err_lo_whole_energies[icat]))
        fluxes_err_hi_whole_energies.append(np.array(lst_flux_err_hi_whole_energies[icat]))

        fluxes_highest_energies.append(np.array(lst_flux_highest_energies[icat]))
        fluxes_err_lo_highest_energies.append(np.array(lst_flux_err_lo_highest_energies[icat]))
        fluxes_err_hi_highest_energies.append(np.array(lst_flux_err_hi_highest_energies[icat]))
        fluxes_ul_highest_energies.append(np.array(lst_flux_ul_highest_energies[icat]))

        nobs_highest_energies.append(np.array(lst_nobs_highest_energies[icat]))
        npred_all_highest_energies.append(np.array(lst_npred_all_highest_energies[icat]))

        bools_flux_finite_highest_energies.append(np.array(lst_bool_flux_finite_highest_energies[icat]))
        bools_non_flux_finite_highest_energies.append(1.-bools_flux_finite_highest_energies[icat])

        bools_lightcurve.append(np.array(lst_bool_lightcurve[icat]))
        bools_non_lightcurve.append(1.-bools_lightcurve[icat])
        lightcurve_indices.append(np.array(lst_lightcurve_index[icat]))
        lightcurve_indices_err.append(np.array(lst_lightcurve_index_err[icat]))

        # Data_plotted
        dp_deviation_vs_fluence_gbm.append(pMatplot.Data_plotted(label=str_category, gr_type='errorbar', xdata=fluences_gbm[icat], ydata=deviations_ts[icat], xdata_err=fluences_gbm_err[icat]))

        dp_flux_whole_energies_vs_fluence_gbm.append(pMatplot.Data_plotted(label=str_category, gr_type='errorbar', xdata=fluences_gbm[icat], ydata=fluxes_whole_energies[icat], xdata_err=fluences_gbm_err[icat], ydata_err=[fluxes_err_lo_whole_energies[icat], fluxes_err_hi_whole_energies[icat]]))

        dp_deviation_vs_fluence_gbm_vs_index_lower_energies.append(pMatplot.Data_plotted(label=str_category, gr_type='scatter', xdata=indices_lower_energies[icat], ydata=fluences_gbm[icat], zdata=deviations_ts[icat]))

        dp_deviation_vs_t90_gbm.append(pMatplot.Data_plotted(label=str_category, gr_type='errorbar', xdata=t90_gbm[icat], xdata_err=t90_gbm_err[icat], ydata=deviations_ts[icat]))

        dp_index_whole_energies_vs_t90_gbm.append(pMatplot.Data_plotted(label=str_category, gr_type='errorbar', xdata=t90_gbm[icat], xdata_err=t90_gbm_err[icat], ydata=indices_whole_energies[icat], ydata_err=indices_err_whole_energies[icat]))

        dp_index_whole_energies_vs_nobs_highest_energies_vs_flux_whole_energies_vs_fluence_gbm.append(pMatplot.Data_plotted(label=str_nobs_highest_energies, gr_type='scatter', xdata=fluences_gbm[icat], ydata=fluxes_whole_energies[icat], zdata=nobs_highest_energies[icat], wdata=indices_whole_energies[icat]))

        dp_index_whole_energies_vs_fluence_gbm_vs_nobs_highest_energies_vs_npred_all_highest_energies.append(pMatplot.Data_plotted(label=str_index_whole_energies, gr_type='scatter', xdata=npred_all_highest_energies[icat], ydata=nobs_highest_energies[icat], wdata=indices_whole_energies[icat], zdata=np.log(fluences_gbm[icat]/(min(fluences_gbm[icat])*np.ones_like(fluences_gbm[icat])))))
        #print dp_index_whole_energies_vs_fluence_gbm_vs_nobs_highest_energies_vs_npred_all_highest_energies[-1].xdata

        dp_index_whole_energies_vs_fluence_gbm.append(pMatplot.Data_plotted(label=str_index_whole_energies, gr_type='errorbar', xdata=fluences_gbm[icat], ydata=indices_whole_energies[icat], xdata_err=fluences_gbm_err[icat], ydata_err=indices_err_whole_energies[icat]))

        dp_deviation_vs_indices_lower_energies.append(pMatplot.Data_plotted(label=str_category, gr_type='errorbar', xdata=indices_lower_energies[icat], ydata=deviations_ts[icat], xdata_err=indices_err_lower_energies[icat]))

        dp_deviation.append(pMatplot.Data_plotted(label=str_category, gr_type='hist', xdata=deviations_ts[icat]))

        dp_flux_highest_energies_vs_fluence_gbm.append(pMatplot.Data_plotted(label=str_category, 
                                                                             gr_type='errorbar', 
                                                                             xdata=np.array([f for f in bools_flux_finite_highest_energies[icat]*fluences_gbm[icat] if f>0]), ydata=fluxes_highest_energies[icat], xdata_err=np.array([f for f in bools_flux_finite_highest_energies[icat]*fluences_gbm_err[icat] if f>0]), 
                                                                             ydata_err=[fluxes_err_lo_highest_energies[icat], fluxes_err_hi_highest_energies[icat]]))
        dp_flux_ul_highest_energies_vs_fluence_gbm.append(pMatplot.Data_plotted(label=str_category+' UL', 
                                                                                gr_type='errorbar', 
                                                                                xdata=np.array([f for f in (np.ones_like(bools_flux_finite_highest_energies[icat])-bools_flux_finite_highest_energies[icat])*fluences_gbm[icat] if f>0]), 
                                                                                ydata=fluxes_ul_highest_energies[icat], xdata_err=np.array([f for f in (np.ones_like(bools_flux_finite_highest_energies[icat])-bools_flux_finite_highest_energies[icat])*fluences_gbm_err[icat] if f>0]), 
                                                                                ul=True))
        dp_index_whole_energies.append(pMatplot.Data_plotted(label=str_category, gr_type='hist', xdata=indices_whole_energies[icat]))

        fluences_gbm_finite_flux_highest_energies = np.array([f for f in bools_flux_finite_highest_energies[icat]*fluences_gbm[icat] if f>0])
        fluences_gbm_err_finite_flux_highest_energies = np.array([f for f in bools_flux_finite_highest_energies[icat]*fluences_gbm_err[icat] if f>0])
        fluences_gbm_ul_flux_highest_energies = np.array([f for f in (np.ones_like(bools_flux_finite_highest_energies[icat])-bools_flux_finite_highest_energies[icat])*fluences_gbm[icat] if f>0])
        fluences_gbm_err_ul_flux_highest_energies = np.array([f for f in (np.ones_like(bools_flux_finite_highest_energies[icat])-bools_flux_finite_highest_energies[icat])*fluences_gbm_err[icat] if f>0])
        dp_flux_highest_energies_per_fluence_gbm_vs_fluence_gbm.append(pMatplot.Data_plotted(label=str_category, 
                                                                                             gr_type='errorbar', 
                                                                                             xdata=fluences_gbm_finite_flux_highest_energies, 
                                                                                             ydata=fluxes_highest_energies[icat]/fluences_gbm_finite_flux_highest_energies, 
                                                                                             xdata_err=fluences_gbm_err_finite_flux_highest_energies, 
                                                                                             ydata_err=[fluxes_highest_energies[icat]/fluences_gbm_finite_flux_highest_energies*np.sqrt(fluxes_err_lo_highest_energies[icat]*fluxes_err_lo_highest_energies[icat]/fluxes_highest_energies[icat]/fluxes_highest_energies[icat]+fluences_gbm_err_finite_flux_highest_energies*fluences_gbm_err_finite_flux_highest_energies/fluences_gbm_finite_flux_highest_energies/fluences_gbm_finite_flux_highest_energies), 
                                                                                                        fluxes_highest_energies[icat]/fluences_gbm_finite_flux_highest_energies*np.sqrt(fluxes_err_hi_highest_energies[icat]*fluxes_err_hi_highest_energies[icat]/fluxes_highest_energies[icat]/fluxes_highest_energies[icat] + fluences_gbm_err_finite_flux_highest_energies*fluences_gbm_err_finite_flux_highest_energies/fluences_gbm_finite_flux_highest_energies/fluences_gbm_finite_flux_highest_energies)]))

        dp_flux_ul_highest_energies_per_fluence_gbm_vs_fluence_gbm.append(pMatplot.Data_plotted(label='95% UL', 
                                                                             gr_type='errorbar', 
                                                                             xdata=fluences_gbm_ul_flux_highest_energies, ydata=fluxes_ul_highest_energies[icat]/fluences_gbm_ul_flux_highest_energies, xdata_err=fluences_gbm_err_ul_flux_highest_energies, ul=True))

        dp_lightcurve_index.append(pMatplot.Data_plotted(label=str_category, gr_type='errorbar', xdata=np.array([f for f in bools_lightcurve[icat]*fluences_gbm[icat] if f>0]), ydata=lightcurve_indices[icat], xdata_err=np.array([f for f in bools_lightcurve[icat]*fluences_gbm_err[icat] if f>0]), ydata_err=lightcurve_indices_err[icat]))

        dp_lightcurve_index_vs_spectral_index.append(pMatplot.Data_plotted(label=str_category, gr_type='errorbar', xdata=np.array([f for f in bools_lightcurve[icat]*indices_whole_energies[icat] if f!=0]), ydata=lightcurve_indices[icat], xdata_err=np.array([f for f in bools_lightcurve[icat]*indices_err_whole_energies[icat] if f!=0]), ydata_err=lightcurve_indices_err[icat]))

        print 'X'
        print np.array([h for h in np.array([g for f, g in zip(bools_lightcurve[icat], bools_flux_finite_highest_energies[icat]) if f==1])*lightcurve_indices[icat] if h!=0])
        print 'Xerr'
        print np.array([h for h in np.array([g for f, g in zip(bools_lightcurve[icat], bools_flux_finite_highest_energies[icat]) if f==1])*lightcurve_indices_err[icat] if h!=0])
        print 'Y'
        print np.array([h for h in np.array([f for f, g in zip(bools_lightcurve[icat], bools_flux_finite_highest_energies[icat]) if g==1])*fluxes_highest_energies[icat] if h!=0])
        print 'Yerr'
        print [np.array([h for h in np.array([f for f, g in zip(bools_lightcurve[icat], bools_flux_finite_highest_energies[icat]) if g==1])*fluxes_err_lo_highest_energies[icat] if h!=0]), np.array([h for h in np.array([f for f, g in zip(bools_lightcurve[icat], bools_flux_finite_highest_energies[icat]) if g==1])*fluxes_err_hi_highest_energies[icat] if h!=0])]


        dp_flux_highest_energies_vs_lightcurve_index.append(pMatplot.Data_plotted(label=str_category, gr_type='errorbar', ydata=np.array([h for h in np.array([f for f, g in zip(bools_lightcurve[icat], bools_flux_finite_highest_energies[icat]) if g==1])*fluxes_highest_energies[icat] if h!=0]), xdata=np.array([h for h in np.array([g for f, g in zip(bools_lightcurve[icat], bools_flux_finite_highest_energies[icat]) if f==1])*lightcurve_indices[icat] if h!=0]), ydata_err=[np.array([h for h in np.array([f for f, g in zip(bools_lightcurve[icat], bools_flux_finite_highest_energies[icat]) if g==1])*fluxes_err_lo_highest_energies[icat] if h!=0]), np.array([h for h in np.array([f for f, g in zip(bools_lightcurve[icat], bools_flux_finite_highest_energies[icat]) if g==1])*fluxes_err_hi_highest_energies[icat] if h!=0])], xdata_err=np.array([h for h in np.array([g for f, g in zip(bools_lightcurve[icat], bools_flux_finite_highest_energies[icat]) if f==1])*lightcurve_indices_err[icat] if h!=0])))

        dp_flux_ul_highest_energies_vs_lightcurve_index.append(pMatplot.Data_plotted(label='Flux UL', gr_type='errorbar', ydata=np.array([h for h in np.array([f for f, g in zip(bools_lightcurve[icat], bools_non_flux_finite_highest_energies[icat]) if g==1])*fluxes_ul_highest_energies[icat] if h!=0]), xdata=np.array([h for h in np.array([g for f, g in zip(bools_lightcurve[icat], bools_non_flux_finite_highest_energies[icat]) if f==1])*lightcurve_indices[icat] if h!=0]), xdata_err=np.array([h for h in np.array([g for f, g in zip(bools_lightcurve[icat], bools_non_flux_finite_highest_energies[icat]) if f==1])*lightcurve_indices_err[icat] if h!=0]), ul=True))

#        dp_flux_ul_highest_energies_vs_lightcurve_index.append(pMatplot.Data_plotted(label=str_category, gr_type='errorbar', ydata=np.array([f for f in bools_lightcurve[icat]*bools_non_flux_finite_highest_energies[icat]*fluxes_ul_highest_energies[icat] if f!=0]), xdata=np.array([f for f in lightcurve_indices[icat]*bools_non_flux_finite_highest_energies[icat] if f!=0]), xdata_err=np.array([f for f in bools_non_flux_finite_highest_energies[icat]*lightcurve_indices_err[icat] if f!=0]), ul=True))

    # Fit index
    def residuals_weighted(const, y, x, yerr):
        err = (y - const) / yerr
        return(err)

    index_fitted = leastsq(residuals_weighted, -2.,
                                          args=(indices_whole_energies[0], fluences_gbm[0], indices_err_whole_energies[0]), 
                                          full_output=True)
    chisq_index = 0

    for jgrb in range(len(lst_dev_scanned[0])):
        dp_ts_scanned.append(pMatplot.Data_plotted(label=name, gr_type='contour', xdata=lst_index_scanned[0][jgrb], ydata=lst_dev_scanned[0][jgrb], zdata=lst_ts_fit_scanned[0][jgrb]))
        chisq_index += pow(residuals_weighted(index_fitted[0][0], indices_whole_energies[0][jgrb], fluences_gbm[0][jgrb], indices_err_whole_energies[0][jgrb]), 2)
    chisq_index = chisq_index/len(lst_dev_scanned[0])
    logger.info('Averaged index: {0} +/- {1}'.format(index_fitted[0][0], sqrt(index_fitted[1][0][0])))
    logger.info('Chi^2/NDF: {0}'.format(chisq_index))
    
    # plot
    ps_deviation_vs_fluence_gbm = pMatplot.Plot_single((dp_deviation_vs_fluence_gbm[0],), xtitle=str_fluence_gbm, ytitle=str_deviation, xlog=True)
    ps_deviation_vs_fluence_gbm.plot(path_save='{dire}/fig_deviation_vs_fluence_gbm'.format(dire=path_out))

    ps_flux_whole_energies_vs_fluence_gbm = pMatplot.Plot_single(dp_flux_whole_energies_vs_fluence_gbm[1:], xtitle=str_fluence_gbm, ytitle=str_flux_whole_energies, xlog=True, ylog=True)
    ps_flux_whole_energies_vs_fluence_gbm.plot(path_save='{dire}/fig_flux_whole_energies_vs_fluence_gbm'.format(dire=path_out))

    ps_deviation_vs_fluence_gbm_vs_index_lower_energies = pMatplot.Plot_single((dp_deviation_vs_fluence_gbm_vs_index_lower_energies[0],), title=str_deviation, xtitle=str_index_lower_energies, ytitle=str_fluence_gbm, ylog=True)
    ps_deviation_vs_fluence_gbm_vs_index_lower_energies.plot(path_save='{dire}/fig_deviation_vs_fluence_gbm_vs_index_lower_energies'.format(dire=path_out))

    ps_deviation_vs_t90_gbm = pMatplot.Plot_single(dp_deviation_vs_t90_gbm[1:], xtitle=str_t90_gbm, ytitle=str_deviation, xlog=True)
    ps_deviation_vs_t90_gbm.plot(path_save='{dire}/fig_deviation_vs_t90_gbm'.format(dire=path_out))

    ps_index_whole_energies_vs_t90_gbm = pMatplot.Plot_single(dp_index_whole_energies_vs_t90_gbm[1:], xtitle=str_t90_gbm, ytitle=str_index_whole_energies, xlog=True)
    ps_index_whole_energies_vs_t90_gbm.plot(path_save='{dire}/fig_index_whole_energies_vs_t90_gbm'.format(dire=path_out))

    ps_index_whole_energies_vs_nobs_highest_energies_vs_flux_whole_energies_vs_fluence_gbm = pMatplot.Plot_single((dp_index_whole_energies_vs_nobs_highest_energies_vs_flux_whole_energies_vs_fluence_gbm[0],), title=str_nobs_highest_energies, xtitle=str_fluence_gbm, ytitle=str_flux_whole_energies, xlog=True, ylog=True)
    ps_index_whole_energies_vs_nobs_highest_energies_vs_flux_whole_energies_vs_fluence_gbm.plot(path_save='{dire}/fig_index_whole_energies_vs__nobs_highest_energies_vs_flux_whole_energies_vs_fluence_gbm'.format(dire=path_out))

    ps_index_whole_energies_vs_fluence_gbm_vs_nobs_highest_energies_vs_npred_all_highest_energies = pMatplot.Plot_single((dp_index_whole_energies_vs_fluence_gbm_vs_nobs_highest_energies_vs_npred_all_highest_energies[0],), xtitle=str_npred_all_highest_energies, ytitle=str_nobs_highest_energies)
    ps_index_whole_energies_vs_fluence_gbm_vs_nobs_highest_energies_vs_npred_all_highest_energies.plot(path_save='{dire}/fig_index_whole_energies_vs_fluence_gbm_vs_nobs_highest_energies_vs_npred_all_highest_energies'.format(dire=path_out))

    ps_index_whole_energies_vs_fluence_gbm = pMatplot.Plot_single((dp_index_whole_energies_vs_fluence_gbm[0],), xtitle=str_fluence_gbm, ytitle=str_index_whole_energies, xlog=True)
    ps_index_whole_energies_vs_fluence_gbm.plot(path_save='{dire}/fig_index_whole_energies_vs_fluence_gbm'.format(dire=path_out))

    ps_deviation_vs_indices_lower_energies = pMatplot.Plot_single(dp_deviation_vs_indices_lower_energies[1:], xtitle=str_index_lower_energies, ytitle=str_deviation)
    ps_deviation_vs_indices_lower_energies.plot(path_save='{dire}/fig_deviation_vs_indices_lower_energies'.format(dire=path_out))

    ps_deviation = pMatplot.Hist_single(dp_deviation[1:], xtitle=str_deviation)
    ps_deviation.plot(bins=12, range=(-3, 3), path_save='{dire}/fig_deviation'.format(dire=path_out))

    ps_flux_highest_energies_vs_fluence_gbm = pMatplot.Plot_single((dp_flux_highest_energies_vs_fluence_gbm[0], dp_flux_ul_highest_energies_vs_fluence_gbm[0]), xtitle=str_fluence_gbm, ytitle=str_flux_highest_energies, xlog=True, ylog=True)
#    ps_flux_highest_energies_vs_fluence_gbm.plot(ylim=(1E-9, 2E-6), path_save='{dire}/fig_flux_highest_energies_vs_fluence_gbm'.format(dire=path_out))
    ps_flux_highest_energies_vs_fluence_gbm.plot(ylim=(0.1*min(fluxes_highest_energies[0]), min(1e-3, 2.*max(fluxes_ul_highest_energies[0]))), path_save='{dire}/fig_flux_highest_energies_vs_fluence_gbm'.format(dire=path_out))

    ps_index_whole_energies = pMatplot.Hist_single(dp_index_whole_energies[1:], xtitle=str_index_whole_energies)
    ps_index_whole_energies.plot(bins=25, range=(-3.5, -1), path_save='{dire}/fig_index_whole_energies'.format(dire=path_out))

    ps_flux_highest_energies_per_fluence_gbm_vs_fluence_gbm = pMatplot.Plot_single((dp_flux_highest_energies_per_fluence_gbm_vs_fluence_gbm[0], dp_flux_ul_highest_energies_per_fluence_gbm_vs_fluence_gbm[0]), xtitle=str_fluence_gbm, ytitle=str_flux_highest_energies_per_fluence_gbm, xlog=True, ylog=True)
    ps_flux_highest_energies_per_fluence_gbm_vs_fluence_gbm.plot(ylim=(1E-5, 10), path_save='{dire}/fig_flux_highest_energies_per_fluence_gbm_vs_fluence_gbm'.format(dire=path_out))

    ps_lightcurve_index_vs_fluence_gbm = pMatplot.Plot_single((dp_lightcurve_index[0],), xtitle=str_fluence_gbm, ytitle=str_lightcurve_index, xlog=True, ylog=False)
    ps_lightcurve_index_vs_fluence_gbm.plot(path_save='{dire}/fig_lightcurve_index_vs_fluence_gbm'.format(dire=path_out))

    ps_lightcurve_index_vs_spectral_index = pMatplot.Plot_single(dp_lightcurve_index_vs_spectral_index[1:], xtitle=str_index_whole_energies, ytitle=str_lightcurve_index, xlog=False, ylog=False)
    ps_lightcurve_index_vs_spectral_index.plot(path_save='{dire}/fig_lightcurve_index_vs_spectral_index'.format(dire=path_out))

    ps_flux_highest_energies_vs_lightcurve_index = pMatplot.Plot_single((dp_flux_highest_energies_vs_lightcurve_index[0], dp_flux_ul_highest_energies_vs_lightcurve_index[0]), xtitle=str_lightcurve_index, ytitle=str_flux_highest_energies, xlog=False, ylog=True)
    ps_flux_highest_energies_vs_lightcurve_index.plot(ylim=(0.1*min(fluxes_highest_energies[0]), min(1e-3, 2.*max(fluxes_ul_highest_energies[0]))), path_save='{dire}/fig_flux_highest_energies_vs_lightcurve_index'.format(dire=path_out))

    ps_ts_scanned1 = pMatplot.Plot_single(dp_ts_scanned, xtitle=str_index_lower_energies, ytitle=str_ts_highest_energies, title=str_ts_lower_energies, xfigsize=20, yfigsize=20, zlevels=[2.30])
    ps_ts_scanned1.plot(xlim=(-3.5, -1.0), ylim=(-10, 25), path_save='{dire}/fig_ts_scanned_1sigma'.format(dire=path_out)) #ylim=(-100, 100), 
    ps_ts_scanned1.ax.axhline(y=2.30, label=r'$1\sigma$', ls='--', c='k')
    ps_ts_scanned1.ax.axhline(y=4.61, label=r'$2\sigma$', ls='--', c='k')
    ps_ts_scanned1.ax.axhline(y=-2.30, label=r'$-1\sigma$', ls='--', c='k')
    ps_ts_scanned1.ax.axhline(y=-4.61, label=r'$-2\sigma$', ls='--', c='k')
    #ps_ts_scanned1.savefig('{dire}/fig_ts_scanned_2sigma'.format(dire=path_out))
    ps_ts_scanned1.fig.savefig('{dire}/fig_ts_scanned_1sigma.png'.format(dire=path_out))
    ps_ts_scanned1.fig.savefig('{dire}/fig_ts_scanned_1sigma.pdf'.format(dire=path_out))

    ps_ts_scanned4 = pMatplot.Plot_single(dp_ts_scanned, xtitle=str_index_lower_energies, ytitle=str_ts_highest_energies, title=str_ts_lower_energies, xfigsize=20, yfigsize=20, zlevels=[6.18])
    ps_ts_scanned4.plot(xlim=(-3.5, -1.0), ylim=(-10, 25), path_save='{dire}/fig_ts_scanned_2sigma'.format(dire=path_out)) #ylim=(-100, 100), 
    ps_ts_scanned4.ax.axhline(y=2.30, label=r'$1\sigma$', ls='--', c='k', lw=2)
    ps_ts_scanned4.ax.axhline(y=4.61, label=r'$2\sigma$', ls='--', c='k', lw=2)
    ps_ts_scanned4.ax.axhline(y=-2.30, label=r'$-1\sigma$', ls='--', c='k', lw=2)
    ps_ts_scanned4.ax.axhline(y=-4.61, label=r'$-2\sigma$', ls='--', c='k', lw=2)
    #ps_ts_scanned4.savefig('{dire}/fig_ts_scanned_2sigma'.format(dire=path_out))
    ps_ts_scanned4.fig.savefig('{dire}/fig_ts_scanned_2sigma.png'.format(dire=path_out))
    ps_ts_scanned4.fig.savefig('{dire}/fig_ts_scanned_2sigma.pdf'.format(dire=path_out))
                                                     

@click.command()
@click.argument('inputs', type=str)
@click.argument('output', type=str)
@click.option('--table', type=str, default='/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LAT2CATALOG-v1-LTF.fits')
def main(inputs, output, table):
    lst_path_inputs = ls_list(inputs)
    plot_grbs_deviation(lst_path_pickle=lst_path_inputs, path_table=table, path_out=output)


if __name__ == '__main__':
    main()
