#!/usr/bin/env python
"""Module for making light curves of LAT data.
The main class LightCurve is a chain of another module pLATLikelihoodConfig.py.
The authour: Mitsunari Takahashi
* Version: 2.1 (2017.12.10) 
  Automatic change of energy for e2dnde
* Version: 2.0 (2017.12.09) 
  Can set start time
* Version: 1.8 (2017.10.02) 
  Parameter limits with ndf=2 for free index
* Version: 1.7 (2017.09.27) 
  Fit with weighting by yerr/ydata
* Version: 1.6 (2017.09.26) 
  More homogeneous time binning
* Version: 1.5 (2017.09.22) 
  Energy vs. time plot.
* Version: 1.4 (2017.09.22) 
  Fitting light curves.
* Version: 1.3 (2017.09.21) 
  Fixed bug of freeing index
* Version: 1.2 (2017.09.19) 
  Save results as CSV file.
* Version: 1.1 (2017.09.19) 
  Allowed time bin shorter than 1s.
* Version: 1.0 (2017.09.19) 
  Added plot function.
* Version: 0.0 (2017.09.16)
"""
import sys
import os
import os.path
path_upstairs = os.path.join(os.path.dirname(__file__), '../')
sys.path.append(path_upstairs)
import logging
import pickle
import subprocess
import datetime
import numpy as np
import itertools
from collections import OrderedDict
import math
from math import log10, log, sqrt, ceil, isnan, pi, factorial, cos, sin, degrees, radians
from sympy import *
from scipy import integrate
from scipy.interpolate import interp2d, interp1d
from scipy.optimize import minimize_scalar
from astropy.io import fits
import click
import ROOT
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import pyLikelihood
from UnbinnedAnalysis import *
from BinnedAnalysis import *
import STLikelihoodAnalysis.pLATLikelihoodConfig as pLATLikelihoodConfig
from STLikelihoodAnalysis.pLATLikelihoodConfig import TABLE_DLOGLIKE_SIGNIFICANCE
from FindGoodstatPeriods import find_goodstat_periods, get_entries, get_event_time_and_energy, get_event_time_energy_angsep, find_goodstat_periods_back
from FindCrossEarthlimb import find_cross_earthlimb#, get_good_intervals
#from STLikelihoodAnalysis import get_module_logger
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
import ReadLATCatalogueInfo
import ReadGBMCatalogueInfo
import GetGTI
import pickle_utilities
import pMatplot
import pMETandMJD
import PlotGRBClosureRelations


mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

##### VERSION OF THIS MACRO #####
VERSION = 2.1 # 2017.12.10


##### Logger #####
#logger = get_module_logger(__name__)
#logger = getLogger('__main__')
logger = getLogger(__name__)
handler = StreamHandler()
loglevel = 'INFO'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)


##### Conversion from MeV to erg ######
MEVtoERG = 1.6021766208E-6


##### Scanned range of parameters #####
DICT_SCANNED_RANGES_PARAMETERS = {'flux': (1E-12, 1E-2),
                                  'eflux': (1E-6, 1.),
                                  'e2dnde': (1E-6, 1.),
                                  'Index': (-5., 3.)}


class LightCurve:
    def __init__(self, name, met0, tmin, tmax, emin, emax, eref, evclass=128, evtype=3, ft2interval='30s', deg_roi=12., rad_margin=10., zmax=100., index_fixed=None, suffix='', grbcatalogue=pLATLikelihoodConfig.GRB_CATALOGUE_LAT, psForce=False, refit=True, force=False, outdir=None, phase='lightcurve'):

        self.config = {'name':name, 
                       'time':{'min':tmin, 'max':tmax}, 
                       'met':{'o':met0, 'min':met0+tmin, 'max':met0+tmax}, 
                       'energy':{'min':emin, 'max':emax, 'ref':eref},
                       'evclass':evclass,
                       'evtype':evtype, 
                       'ft2interval':ft2interval, 
                       'roi':{'radius':deg_roi, 'margin':rad_margin}, 
                       'zenith':{'max':zmax}, 
                       'index_fixed':index_fixed, 
                       'suffix':suffix, 
                       'model':{'psForce':psForce},
                       'refit':refit,
                       'force':force,
                       'phase':phase}
        self.outdir = outdir if outdir is not None else os.getcwd()


class LightCurveGRB(LightCurve):
    def __init__(self, name, wholephase='special', tmin=0.0, tmax=10000.0, emin=100.0, emax=100000.0, eref=1000., evclass=128, evtype=3, ft2interval='30s', deg_roi=12., rad_margin=10., zmax=100., index_fixed=None, redshift=0.0, suffix='', grbcatalogue=pLATLikelihoodConfig.GRB_CATALOGUE_LAT, psForce=False, rlim_goodstat=2.0, ngoodstat=10, ntbinperdecade=4, refit=True, force=False, outdir=None, phase='lightcurve', spectraltype='PowerLaw', spectralpars={'Prefactor':1e-10, 'Index':-2, 'Scale':1000}, synccutoff=None):
        # Target GRB
        self.grb = pLATLikelihoodConfig.GRBTarget(name, grbcatalogue, spectraltype=spectraltype, spectralpars=spectralpars, redshift=redshift, synccutoff=synccutoff)
        #self.grb = pLATLikelihoodConfig.GRBTarget(name, grbcatalogue, spectraltype='ScaleFactor::PowerLaw2', spectralpars={'Integral':norm, 'Index':index, 'LowerLimit':emin, 'UpperLimit':emax, 'ScaleFactor':scalefactor})

        # Light curve instance
        LightCurve.__init__(self, name=name, met0=self.grb.met0, tmin=tmin, tmax=tmax, emin=emin, emax=emax, eref=eref, evclass=evclass, evtype=evtype, ft2interval=ft2interval, deg_roi=deg_roi, rad_margin=rad_margin, zmax=zmax, index_fixed=index_fixed, suffix=suffix, psForce=psForce, refit=refit, force=force, phase=phase)
        self.config['goodstat'] = {'n':ngoodstat, 'r':rlim_goodstat}
        self.config['ntbinperdecade'] = ntbinperdecade
        self.config['wholephase'] = wholephase
        self.periods_goodstat = []
        self.periods_goodstat_zabove90 = []
        self.analyses = []
        self.counts_time = None
        self.counts_energy = None
        self.counts_angsep = None
        self.calonly_counts_time = None
        self.calonly_counts_energy = None
        self.calonly_counts_angsep = None
        self.synccutoff = synccutoff


    def setup(self):
        # Analysis instance which covers the whole time range
        self.analysis_whole = pLATLikelihoodConfig.GRBConfig(target=self.grb, phase=self.config['wholephase'], tstop=self.config['time']['max'], emin=self.config['energy']['min'], emax=self.config['energy']['max'], deg_roi=self.config['roi']['radius'], zmax=self.config['zenith']['max'] if self.config['energy']['min']>=300. else min(90., self.config['zenith']['max']), suffix=self.config['suffix'], tmin_special=self.config['time']['min'], tmax_special=self.config['time']['max'], ft2interval='30s')
        self.analysis_whole.setup(force={'download':False, 'filter':self.config['force'], 'maketime':True, 'evtbin':self.config['force'], 'livetime':self.config['force'], 'exposure':self.config['force'], 'model_3FGL_sources':self.config['force'], 'diffuse_responses':self.config['force'], 'srcmaps':self.config['force']})
        if self.config['force']==True or not os.path.exists(self.analysis_whole.path_model_xml.replace('.xml', '_new.xml')):
            self.analysis_whole.fit(bredo=self.config['refit'])
        else:
            self.analysis_whole.path_model_xml_new = self.analysis_whole.path_model_xml.replace('.xml', '_new.xml')

        self.outdir = '{base}/{target}/{energy}/{roi}/{phase}'.format(base=self.analysis_whole.dir_base, target=self.analysis_whole.target.name, energy=self.analysis_whole.str_energy, roi=self.analysis_whole.str_roi, phase=self.config['phase'])
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
        self.outbasename = 'LightCurve_{target}_{spectype}_{index}{suffix}'.format(target=self.analysis_whole.target.name, spectype=self.analysis_whole.target.spectraltype, index=self.analysis_whole.str_index, suffix=self.analysis_whole.suffix)
        self.path_periods_save = '{0}/{1}_TimeBins_T{2:0>6.0f}-{3:0>6.0f}sec_{4}.pickle'.format(self.outdir, self.outbasename, self.analysis_whole.tmin, self.analysis_whole.tmax, self.analysis_whole.ft2interval)

        if os.path.exists(self.path_periods_save):
            periods_loaded = pickle_utilities.load(self.path_periods_save)
            logger.info('Loading analysis periods from {0}...'.format(self.path_periods_save))
            if periods_loaded['config']['time']['min']==self.config['time']['min'] and periods_loaded['config']['time']['max']==self.config['time']['max'] and periods_loaded['config']['roi']['radius']==self.config['roi']['radius'] and periods_loaded['config']['goodstat']['n']==self.config['goodstat']['n'] and periods_loaded['config']['goodstat']['r']==self.config['goodstat']['r'] and periods_loaded['config']['ntbinperdecade']==self.config['ntbinperdecade']:
                self.periods_goodstat = periods_loaded['periods']
                self.periods_goodstat_zabove90 = periods_loaded['periods_zabove90']
        if len(self.periods_goodstat)<1:
            validtimes_extra_zmax = []
            if self.config['zenith']['max']<=90 or self.config['energy']['min']>=300:
                validtimes = find_cross_earthlimb(self.analysis_whole.path_ft2, self.analysis_whole.target.ra, self.analysis_whole.target.dec, self.analysis_whole.tmin, self.analysis_whole.tmax, self.analysis_whole.zmax, thcut=self.analysis_whole.thetamax, torigin=self.analysis_whole.target.met0, )
            else:
                validtimes = find_cross_earthlimb(self.analysis_whole.path_ft2, self.analysis_whole.target.ra, self.analysis_whole.target.dec, self.analysis_whole.tmin, self.analysis_whole.tmax, min(90., self.config['zenith']['max']), thcut=self.analysis_whole.thetamax, torigin=self.analysis_whole.target.met0, )
                validtimes_zmax = find_cross_earthlimb(self.analysis_whole.path_ft2, self.analysis_whole.target.ra, self.analysis_whole.target.dec, self.analysis_whole.tmin, self.analysis_whole.tmax, self.config['zenith']['max'], thcut=self.analysis_whole.thetamax, torigin=self.analysis_whole.target.met0, )
                for vt_zmax in validtimes_zmax:
                    for vt in validtimes:
                        if vt_zmax[0]<vt[0] and vt[1]<=vt_zmax[1]:
                            validtimes_extra_zmax.append([vt_zmax[0], vt[0]])
                        if vt[1]<vt_zmax[1] and vt[0]>=vt_zmax[0]:
                            validtimes_extra_zmax.append([vt[1], vt_zmax[1]])
            if self.config['goodstat']['n']>0:
                logger.info("""Good time intervals (zenith < {zen} deg)
{vt}""".format(zen=min(90., self.config['zenith']['max']), vt=validtimes))
           # Find good-statistics interval
                for vt in validtimes:
                    logger.info(vt)
                    if vt[0]<vt[1]:
                        logger.debug('Event file:', self.analysis_whole.path_filtered)
                        self.periods_goodstat += find_goodstat_periods(self.analysis_whole.path_filtered, vt[0], vt[1], nthreshold=self.config['goodstat']['n'], rlim=self.config['goodstat']['r'], ra=self.grb.ra, dec=self.grb.dec, torigin=self.analysis_whole.target.met0)
                    else:
                        logger.warning('The end of the interval {end} is earlier than the start {sta}! This interval is skipped.'.format(end=vt[1], sta=vt[0]))

                logger.info("""Good-statistics peroids (>{nth} events within {rlim} deg)
{vt}""".format(nth=self.config['goodstat']['n'], rlim=self.config['goodstat']['r'], vt=self.periods_goodstat))
                # z>90
                self.periods_goodstat_zabove90 = []
                for vtabove90 in validtimes_extra_zmax:
                    logger.info(vtabove90)
                    if vtabove90[0]<vtabove90[1]:
                        logger.debug('Event file:', self.analysis_whole.path_filtered)
                        self.periods_goodstat_zabove90 += find_goodstat_periods(self.analysis_whole.path_filtered, vt[0], vt[1], nthreshold=self.config['goodstat']['n'], rlim=self.config['goodstat']['r'], ra=self.grb.ra, dec=self.grb.dec, torigin=self.analysis_whole.target.met0)
                    else:
                        logger.warning('The end of the interval {end} is earlier than the start {sta}! This interval is skipped.'.format(end=vt[1], sta=vt[0]))

            else: # Logarithmically equivalent time binning
                # Prompt phase
                tpromptstart = 0. #self.analysis_whole.target.t05
                tpromptend = (self.analysis_whole.target.t25 + self.analysis_whole.target.t50) #*3./2.)
                if validtimes[0][0]<tpromptend and tpromptend<validtimes[0][1]:
                    if tpromptstart>validtimes[0][0] and tpromptstart<tpromptend:
                        #self.periods_goodstat.append([validtimes[0][0], tpromptstart])
                        tp0, tp1 = tpromptstart, tpromptend
                    else:
                        tp0, tp1 = validtimes[0][0], tpromptend
                    self.periods_goodstat += find_goodstat_periods_back(self.analysis_whole.path_filtered, tp0, tp1, nthreshold=8, rlim=12., ra=self.grb.ra, dec=self.grb.dec, torigin=self.analysis_whole.target.met0)
                    #self.periods_goodstat += find_goodstat_periods(self.analysis_whole.path_filtered, tp0, tp1, nthreshold=10, rlim=12., ra=self.grb.ra, dec=self.grb.dec, torigin=self.analysis_whole.target.met0)
                        
                    validtimes[0][0] = tpromptend # Forget about the prompt phase
                # Afterglow phase
                tpow_index = 0.5*(np.log10(2))
                tpow_unit = (pow(10000, -tpow_index) - pow(100000, -tpow_index))/2.# 2 bins in 10^4-10^5 s for index -0.5*(np.log10(2))
                ntpow_bin = 1./tpow_unit
                tbins = []
                validtimes_1 = []
                validtimes_not1 = []
                tbin_trans1 = 0.
                for vt in validtimes:
                    #if np.log10(vt[1]/vt[0])>=1./self.config['ntbinperdecade']: #flag_1==True:
                    if pow(vt[0], -tpow_index) - pow(vt[1], -tpow_index) >= tpow_unit:
                        validtimes_1.append(vt)
                    else:
                        validtimes_not1.append(vt)

                logger.debug('validtimes_1: {0}'.format(validtimes_1))
                logger.debug('validtimes_not1: {0}'.format(validtimes_not1))

                for vt in validtimes_1:
                    logtmin_1 = pow(max(1, vt[0]), -tpow_index) #np.log10(max(1, vt[0]))
                    logtmax_1 = pow(vt[1], -tpow_index) #np.log10(vt[1])
                    tedges_1 = pow(np.linspace(logtmin_1, logtmax_1, int((-logtmax_1+logtmin_1)/tpow_unit+1.5)), -1./tpow_index)
                    #tedges_1 = 10**np.linspace(logtmin_1, logtmax_1, int(self.config['ntbinperdecade']*(logtmax_1-logtmin_1)+1.5))
                    for t0, t1 in zip(tedges_1[:-1],tedges_1[1:]):
                        tbins.append([t0, t1])
                    tbin_trans1 = max(tbin_trans1, vt[1])

                validtimes_3 = []
                validtimes_not3 = []
                print 'tpow_unit: {0}'.format(tpow_unit)
                for vt in validtimes_not1:
                    logger.debug('vt[0]: {0}, vt[1]: {1}'.format(vt[0], vt[1]))
                    logger.debug((pow(vt[0], -tpow_index) - pow(vt[1], -tpow_index) < tpow_unit))
                    if pow(vt[0], -tpow_index) - pow(vt[1], -tpow_index) < tpow_unit:
                        validtimes_3.append(vt)
                    else:
                        validtimes_not3.append(vt)

                logger.debug('validtimes_3: {0}'.format(validtimes_3))
                logger.debug('validtimes_not3: {0}'.format(validtimes_not3))
                tbins += validtimes_not3

                logtmin_3 = pow(max(1, validtimes_3[0][0]), -tpow_index)
                logtmax_3 = pow(validtimes_3[-1][1], -tpow_index) 
                tedges_3 = pow(np.linspace(logtmin_3, logtmax_3, int((-logtmax_3+logtmin_3)/tpow_unit+1.5)), -1./tpow_index)
                for t0, t1 in zip(tedges_3[:-1], tedges_3[1:]):
                    tbins.append([t0, t1])

                self.periods_goodstat += tbins

                #z>90
                self.periods_goodstat_zabove90 = []
                if len(validtimes_extra_zmax)>0:
                    for vt_ex in validtimes_extra_zmax:
                        if pow(vt_ex[0], -tpow_index) - pow(vt_ex[1], -tpow_index) >= tpow_unit/2.:
                            self.periods_goodstat_zabove90.append(vt_ex)

        print """Time slots:
{0}""".format(self.periods_goodstat)
        for ip, pds in enumerate(self.periods_goodstat):
            logger.info('Time slot No.{0}: {1} - {2} s'.format(ip, pds[0], pds[1]))
        print """Time slots (z>90 deg):
{0}""".format(self.periods_goodstat)
        for ip, pds in enumerate(self.periods_goodstat_zabove90):
            logger.info('Time slot No.{0}: {1} - {2} s'.format(ip, pds[0], pds[1]))
        if self.path_periods_save is not None:
            pickle_utilities.dump(self.path_periods_save, {'periods':self.periods_goodstat, 'config':self.config, 'periods_zabove90':self.periods_goodstat_zabove90})
        self.summary_results = []

        # Analysis instances
        for ip, pds in enumerate(self.periods_goodstat):
            logger.debug(pds)
            self.analyses.append(pLATLikelihoodConfig.GRBConfig(target=self.grb, phase=self.config['phase'], emin=self.config['energy']['min'], emax=self.config['energy']['max'], deg_roi=self.config['roi']['radius'], zmax=self.config['zenith']['max'] if self.config['energy']['min']>=300. else min(90., self.config['zenith']['max']), suffix=self.config['suffix'], tmin_special=pds[0], tmax_special=pds[1]))
            self.summary_results.append({})
            self.summary_results[-1]['time'] = {'min':self.analyses[-1].tmin, 'max':self.analyses[-1].tmax}
        #z>90
        for ip, pds in enumerate(self.periods_goodstat_zabove90):
            logger.debug(pds)
            self.analyses.append(pLATLikelihoodConfig.GRBConfig(target=self.grb, phase=self.config['phase'], emin=max(10**2.5, self.config['energy']['min']), emax=self.config['energy']['max'], deg_roi=min(6., self.config['roi']['radius']), zmax=self.config['zenith']['max'], suffix=self.config['suffix'], tmin_special=pds[0], tmax_special=pds[1]))
            self.summary_results.append({})
            self.summary_results[-1]['time'] = {'min':self.analyses[-1].tmin, 'max':self.analyses[-1].tmax}

        for ip, pds in enumerate(self.periods_goodstat+self.periods_goodstat_zabove90):
            self.analyses[ip].setup(force={'download':False, 'filter':self.config['force'], 'maketime':True, 'livetime':self.config['force'], 'exposure':self.config['force'], 'model_3FGL_sources':True, 'diffuse_responses':self.config['force']}, skip_zero_data=True)
            #self.analyses[-1].use_external_model(self.analysis_whole.path_model_xml_new)


    def add_calonly(self, calonly_info):
        path_onevt = calonly_info[0] 
        path_onexp = calonly_info[1] 
        path_bkgevt = calonly_info[2] 
        rclass = calonly_info[3]
        file_onevt = ROOT.TFile(path_onevt, "READ")
        tr_onevt = file_onevt.Get('EVENTS_GRB{0}'.format(self.analysis_whole.target.name))
        file_onexp = ROOT.TFile(path_onexp, "READ")
        file_bkgevt = ROOT.TFile(path_bkgevt, "READ")
        calonly_counts_time = []
        calonly_counts_energy = []
        calonly_counts_angsep = []
        for evt in tr_onevt:
            if evt.FLAG_PSF68 and evt.EVENT_CLASS>=4096 and evt.ZENITH_ANGLE<self.config['zenith']['max']:
                calonly_counts_time.append(evt.TIME_GRB)
                calonly_counts_energy.append(10**evt.ENERGY)
                calonly_counts_angsep.append(evt.ANG_SEP)
        self.calonly_counts_time = np.array(calonly_counts_time)
        self.calonly_counts_energy = np.array(calonly_counts_energy)
        self.calonly_counts_angsep = np.array(calonly_counts_angsep)

        for iana, ana in enumerate(self.analyses):
            dir_bkgevt = file_bkgevt.Get('TimeBin{0}'.format(iana))
            ana.calonly = pLATLikelihoodConfig.CalOnlyData(tr_onevt=tr_onevt, htg_onexp=file_onexp.Get('htgExp_{0}_projYX'.format(iana)), htg_bkgevt=dir_bkgevt.Get('htgExBKG_CalOnly_{rcla}_PSF68_projE'.format(rcla=rclass)), on_energy=(ana.emin, ana.emax), on_time=(ana.target.met0+ana.tmin, ana.target.met0+ana.tmax), on_classes=pLATLikelihoodConfig.DICT_EVCLASSES_BIT_CALONLY[rclass], on_zenith=(0., ana.zmax), on_theta=(0., ana.thetamax))
        #     calonly_counts_time = np.hstack((calonly_counts_time, ana.calonly.onevt_unbinned['time']))
        #     calonly_counts_energy = np.hstack((calonly_counts_time, ana.calonly.onevt_unbinned['energy']))


    def count_energy(self, rlim=6.0):
        events = get_event_time_energy_angsep(self.analysis_whole.path_filtered, self.analysis_whole.tmin, self.analysis_whole.tmax, rlim, ra=self.grb.ra, dec=self.grb.dec, torigin=self.analysis_whole.target.met0, zmax=self.analysis_whole.zmax)
        self.counts_time = events[0]
        self.counts_energy = events[1]
        self.counts_angsep = events[2]


    def run_analysis(self, use_calonly=False):
        for ip, pds in enumerate(self.periods_goodstat+self.periods_goodstat_zabove90):
            logger.info('Analyzing {0}...'.format(pds))
            nevt_rough = self.analyses[ip].nevt_rough
            self.summary_results[ip]['Nobs'] = nevt_rough
            gti_table = GetGTI.get_gti_table(self.analyses[ip].path_filtered_gti)
            if len(gti_table)<1:
                logger.warning('No GTIs in {0} - {1} s'.format(pds[0], pds[1]))
                logger.warning('Skipped...')
                self.summary_results[ip]['TS'] = np.nan
                self.summary_results[ip].update(self.analyses[ip].dct_summary_results)
                continue
            self.analyses[ip].set_likelihood_external_model(self.analysis_whole.path_model_xml_new)
            #self.analyses[ip].diffuse_responses()
            logger.info('Fixing normalization of other sources.')
            for source in self.analyses[ip].like.sourceNames():
                if source not in (self.analyses[ip].target.name):
                    self.analyses[ip].like.normPar(source).setFree(False)
            self.analyses[ip].likeobj = pyLike.NewMinuit(self.analyses[ip].like.logLike)
            try:
                self.analyses[ip].dofit()
                self.analyses[ip].summarize_fit_results()
            except RuntimeError:
                logger.warning('Analysis for {0} did not converge.'.format(pds))
                logger.warning('Relaxing the tolerance...')
                self.analyses[ip].reset_target_norm()
                try:
                    self.analyses[ip].dofit(tol=self.analyses[ip].like.tol*10.)
                    self.analyses[ip].summarize_fit_results()
                except RuntimeError:
                    logger.error('Analysis for {0} did not converge finally!'.format(pds))
                    self.analyses[ip].dct_summary_results['TS'] = np.nan
                    self.analyses[ip].reset_target_norm()
            if sum(self.analyses[ip].like._Nobs())>0:
                self.analyses[ip].plot_countspectra_fitted()

            self.analyses[ip].scan_norm_and_index(eref=self.config['energy']['ref'], use_calonly=use_calonly)
            self.summary_results[ip].update(self.analyses[ip].dct_summary_results)


    def scan_parameters(self, lcintegrals=None, lcindices=None, tnorm=10., rescaler=1.):
        lst_periods = []
        for ip, pds in enumerate(self.periods_goodstat+self.periods_goodstat_zabove90):
            logger.info('==========')
            logger.info(pds)
            lst_periods.append({'period':pds,
                                'duration':self.analyses[ip].duration,
                                'tref':np.zeros(shape=(len(lcintegrals), len(lcindices))),
                                'loglike':np.zeros(shape=(len(lcintegrals), len(lcindices))),
                                'normalization':np.zeros(shape=(len(lcintegrals), len(lcindices))),
                                'flux':np.zeros(shape=(len(lcintegrals), len(lcindices))),
                                'eflux':np.zeros(shape=(len(lcintegrals), len(lcindices))),
                                'fluence':np.zeros(shape=(len(lcintegrals), len(lcindices))),
                                'efluence':np.zeros(shape=(len(lcintegrals), len(lcindices))),
                                'conversion':{'flux':np.zeros(shape=(len(lcintegrals), len(lcindices))), 
                                              'eflux':np.zeros(shape=(len(lcintegrals), len(lcindices)))},
                                'nobs':0,
                                'npred':{'target':np.zeros(shape=(len(lcintegrals), len(lcindices))),
                                         'others':np.zeros(shape=(len(lcintegrals), len(lcindices)))}
                                })
            tb_gti = GetGTI.get_gti_table(self.analyses[ip].path_filtered_gti)
            for ilcintegral, jlcindex in itertools.product(range(len(lcintegrals)), range(len(lcindices))):
                lst_periods[ip]['normalization'][ilcintegral][jlcindex], lst_periods[ip]['tref'][ilcintegral][jlcindex] = lightcurve_prefactor(tmin=pds[0], tmax=pds[1], integral=lcintegrals[ilcintegral], index=lcindices[jlcindex], tb_gti=tb_gti, torigin=self.analyses[ip].target.met0, tnorm=tnorm, rescaler=rescaler) #, self.grb.t95, 100000
            logger.debug("""Normalization:
{0}""".format(lst_periods[ip]['normalization'][ilcintegral][jlcindex]))
            logger.debug("""Tref:
{0}""".format(lst_periods[ip]['tref'][ilcintegral][jlcindex]))
            if self.analyses[ip].duration<=0:
                continue

            self.analyses[ip].set_likelihood()
            if lcintegrals is None:
                lcintegrals = [self.analyses[ip].like.normPar(self.analyses[ip].target.name).getValue()]
            if lcindices is None:
                lcindex_idx = self.analyses[ip].like.par_index(self.analyses[ip].target.name, 'Index')
                lcindices = [np.array([-1.0, -1.3])]

            logger.debug('  Nobs: {0}'.format(sum(self.analyses[ip].like._Nobs())))
            lst_periods[ip]['nobs'] = sum(self.analyses[ip].like._Nobs())

            for ilcintegral, jlcindex in itertools.product(range(len(lcintegrals)), range(len(lcindices))):
                norm_factor = lst_periods[ip]['normalization'][ilcintegral][jlcindex]
                if ilcintegral==0 and jlcindex==0:
                    logger.debug('  Normalization factor: {0}'.format(norm_factor))
                try:
                    self.analyses[ip].like.normPar(self.analyses[ip].target.name).setValue(norm_factor)
                except RuntimeError:
                    logger.critical('Setting normalization factor {0} ({1},{2}) failed!!!'.format(norm_factor, ilcintegral, jlcindex))
                    sys.exit(1)

                if self.analyses[ip].target.spectraltype[-8:]=='PowerLaw':
                    norm_idx = self.analyses[ip].like.par_index(self.analyses[ip].target.name, 'Prefactor')
                elif self.analyses[ip].target.spectraltype[-9:]=='PowerLaw2':
                    norm_idx = self.analyses[ip].like.par_index(self.analyses[ip].target.name, 'Integral')
                else:
                    logger.critical('Spectral function must be PowerLaw or PowerLaw2!!!')
                    sys.exit(1)
                self.analyses[ip].like.freeze(norm_idx)
                index_idx = self.analyses[ip].like.par_index(self.analyses[ip].target.name, 'Index')
                self.analyses[ip].like.freeze(index_idx)

               # Current loglike value
                lst_periods[ip]['loglike'][ilcintegral][jlcindex] = -self.analyses[ip].like()
                if ilcintegral==0 and jlcindex==0:
                    logger.debug('  log Likelihood: {0}'.format(lst_periods[ip]['loglike'][ilcintegral][jlcindex]))
               # Current values
                lst_periods[ip]['flux'][ilcintegral][jlcindex] = self.analyses[ip].like[self.analyses[ip].target.name].flux(self.config['energy']['min'], self.config['energy']['max'])
                lst_periods[ip]['conversion']['flux'][ilcintegral][jlcindex] = lst_periods[ip]['flux'][ilcintegral][jlcindex] / norm_factor
                lst_periods[ip]['eflux'][ilcintegral][jlcindex] = self.analyses[ip].like[self.analyses[ip].target.name].energyFlux(self.config['energy']['min'], self.config['energy']['max'])
                lst_periods[ip]['conversion']['eflux'][ilcintegral][jlcindex] = lst_periods[ip]['eflux'][ilcintegral][jlcindex] / norm_factor

                lst_periods[ip]['fluence'][ilcintegral][jlcindex] = lst_periods[ip]['flux'][ilcintegral][jlcindex] * self.analyses[ip].duration
                lst_periods[ip]['efluence'][ilcintegral][jlcindex] = lst_periods[ip]['eflux'][ilcintegral][jlcindex] * self.analyses[ip].duration
                if ilcintegral==0 and jlcindex==0:
                    logger.debug('  Npred: {0}'.format(sum(self.analyses[ip].like._srcCnts(str(self.analyses[ip].target.name)))))
                for src in self.analyses[ip].like.sourceNames():
                    if src==str(self.analyses[ip].target.name):
                        lst_periods[ip]['npred']['target'][ilcintegral][jlcindex] += sum(self.analyses[ip].like._srcCnts(src))
                    else:
                        lst_periods[ip]['npred']['others'][ilcintegral][jlcindex] += sum(self.analyses[ip].like._srcCnts(src))
        return lst_periods



    def pickle(self, stuff=None):
        self.dct_stored = {'config':self.config, 'results':self.summary_results, 'counts':{'time':self.counts_time, 'energy':self.counts_energy, 'angsep':self.counts_angsep}, 'counts_CalOnly':{'time':self.calonly_counts_time, 'energy':self.calonly_counts_energy, 'angsep':self.calonly_counts_angsep}} if stuff is None else stuff
        path_pickle = '{0}/{1}.pickle'.format(self.outdir, self.outbasename)
        with open(path_pickle, mode='wb') as f:
            pickle.dump(self.dct_stored, f)
        logger.info('Result summary has been serialized as {0}'.format(path_pickle))


def scan_norm_beta_alpha(dict_summary, torigin, outdir=None, suffix='', norms=None, betas=None, alphas=None, tnorm=10., tmin=0., tmax=sys.maxint, mwlindices=None):

    if not 'scan3D' in dict_summary:
        dict_summary['scan3D'] = {}
        dict_summary['scan3D']['alphas'] = alphas
        dict_summary['scan3D']['betas'] = betas
        dict_summary['scan3D']['norms'] = norms
        alpha_mesh, beta_mesh, norm_mesh = np.meshgrid(alphas, betas, norms)
        loglike_mesh = np.full_like(norm_mesh, sys.maxint)
        logger.info('3D mesh shape: {0}'.format(loglike_mesh.shape))

        list_periods_used = []
        t0_min = sys.maxint
        t1_max = -sys.maxint
        for period in dict_summary['results']:
            t0 = period['time']['min']
            t1 = period['time']['max']
            if t0>=tmin and t1<=tmax:
                list_periods_used.append(period)
                list_periods_used[-1]['tb_gti'] = GetGTI.get_gti_table(period['data_path'])
                if t0<t0_min:
                    t0_min = t0
                if t1>t1_max:
                    t1_max = t1

        for inorm, jbeta, kalpha  in itertools.product(range(len(norms)), range(len(betas)), range(len(alphas))):
            if inorm%10+jbeta%10+kalpha%10<2:
                logger.debug('No.{0},{1},{2}'.format(inorm, jbeta, kalpha))
            dll_sum = 0.
            flag_fail = False

            for period in list_periods_used: 
                t0 = period['time']['min']
                t1 = period['time']['max']
                if t0>=tmin and t1<=tmax and flag_fail == False:
                    logger.debug('Time interval {0} - {1} s'.format(t0, t1))
                    logger.debug(period.keys())
                    norm = lightcurve_prefactor(tmin=t0, tmax=t1, integral=norms[inorm], index=alphas[kalpha], tb_gti=period['tb_gti'], torigin=torigin, tnorm=tnorm, rescaler=1.)[0]
                    normalizations = period['dloglike']['normalization']
                    indices = period['dloglike']['Index']
                    dloglikes = period['dloglike']['dloglike']
                    dloglikes[np.isinf(dloglikes)] = sys.maxint
                    func_interp = interp2d(x=indices, y=normalizations, z=dloglikes, kind='linear', bounds_error=True)
                    try:
                        dll = func_interp(betas[jbeta], norm)
                        dll_sum += dll
                    except ValueError:
                        logger.error('ValueError!!')
                        logger.error('Time interval {0} - {1} s'.format(t0, t1))
                        logger.error('Norm: {0} Should be within {1} - {2}'.format(norm, normalizations[0], normalizations[-1]))
                        logger.error('Beta: {0} Should be within {1} - {2}'.format(betas[jbeta], indices[0], indices[-1]))
                        logger.error('Norm at t_scale: {0}, alpha: {1}'.format(norms[inorm], alphas[kalpha]))
                        flag_fail = True
                        break
                else:
                    logger.debug('Time interval {0} - {1} s is skipped.'.format(t0, t1))
            logger.debug('Grid norm: {0}, beta: {1}, alpha:{2}'.format(inorm, jbeta, kalpha))
            logger.debug('Summed d-loglike: {0}'.format(dll_sum))
            if flag_fail==False: 
                loglike_mesh[jbeta][kalpha][inorm] = dll_sum
            else:
                loglike_mesh[jbeta][kalpha][inorm] = sys.maxint

        loglike_min = np.nanmin(loglike_mesh)
        dloglike_mesh = loglike_mesh - loglike_min
        dict_summary['scan3D']['dloglike_mesh'] = dloglike_mesh
    else:        
        dloglike_mesh = dict_summary['scan3D']['dloglike_mesh']
        alphas = dict_summary['scan3D']['alphas']
        betas = dict_summary['scan3D']['betas']
        norms = dict_summary['scan3D']['norms']
        alpha_mesh, beta_mesh, norm_mesh = np.meshgrid(alphas, betas, norms)
        t0_min = sys.maxint
        t1_max = -sys.maxint
        for period in dict_summary['results']:
            t0 = period['time']['min']
            t1 = period['time']['max']
            if t0>=tmin and t1<=tmax:
                if t0<t0_min:
                    t0_min = t0
                if t1>t1_max:
                    t1_max = t1

    print 'd-loglike mesh:'
    print dloglike_mesh

    # Projection to alpha-beta
    dloglike_alpha_beta = np.full_like(dloglike_mesh[:,:,0], float(sys.maxint))
    for jb, ka  in itertools.product(range(len(betas)), range(len(alphas))):
        dll_local = dloglike_mesh[jb,ka,:]
        dloglike_alpha_beta[jb, ka] = np.nanmin(dll_local)
    dloglike_alpha_beta = np.ma.masked_invalid(dloglike_alpha_beta)

    # Projection to norm-alpha
    dloglike_norm_alpha = np.full_like(dloglike_mesh[0,:,:], float(sys.maxint))
    for ka, ino in itertools.product(range(len(alphas)), range(len(norms))):
        dll_local = dloglike_mesh[:,ka,ino]
        dloglike_norm_alpha[ka,ino] = np.nanmin(dll_local)
    dloglike_norm_alpha = np.ma.masked_invalid(dloglike_norm_alpha)

    # Projection to beta-norm
    dloglike_beta_norm = np.full_like(dloglike_mesh[:,0,:], float(sys.maxint))
    for ino, jb  in itertools.product(range(len(norms)), range(len(betas))):
        dll_local = dloglike_mesh[jb,:,ino]
        dloglike_beta_norm[jb,ino] = np.nanmin(dll_local)
    dloglike_beta_norm = np.ma.masked_invalid(dloglike_beta_norm)

    # Plot for only result
    fig_result = plt.figure(figsize=(5, 5))
    ax_result = fig_result.add_axes((0.15, 0.1, 0.8, 0.8))
    cont_levels_NDF3 = sorted(TABLE_DLOGLIKE_SIGNIFICANCE[2].as_void())
    cont_alpha_beta = ax_result.contour(-alpha_mesh[:,:,0], -1.-beta_mesh[:,:,0], dloglike_alpha_beta*2., levels=cont_levels_NDF3, colors='k', linestyles=pMatplot.TPL_LINE)
    if mwlindices is not None:
        mwlobs = PlotGRBClosureRelations.ObservedIndices()
        mwlobs.read(mwlindices, instruments=['XRT'])
        mwlobs.draw(ax_result)
    ax_result.set_title(dict_summary['config']['name'])
    ax_result.set_xlabel(r'Temporal index $\alpha$')
    ax_result.set_ylabel(r'Spectral index $\beta$')
    ax_result.set_xlim((0, 2.5))
    ax_result.set_ylim((-0.5, 1.5))
    ax_result.grid()
    for ff in ('png','pdf'):
        path_save = "{dire}/scanned3D_result_{targ}{suff}.{form}".format(dire=outdir, targ=dict_summary['config']['name'], suff=suffix, form=ff)
        fig_result.savefig(path_save)
        logger.info('{0} has been saved.'.format(path_save))
    fig_result.clf()


    # Plot for comparison with closure relations
    fig = plt.figure(figsize=(12, 12))
    ax = [[fig.add_axes((0.075, 0.6, 0.4, 0.35)),
          fig.add_axes((0.55, 0.6, 0.4, 0.35))],
          [fig.add_axes((0.075, 0.15, 0.4, 0.35)),
          fig.add_axes((0.55, 0.15, 0.4, 0.35))]]
    cbaxes = fig.add_axes([0.1, 0.05, 0.8, 0.02]) 
    print 'alpha:', alpha_mesh[:,:,0].shape
    print alpha_mesh[:,:,0]
    print 'beta:', beta_mesh[:,:,0].shape
    print beta_mesh[:,:,0]
    print 'd-loglike:', dloglike_alpha_beta.shape
    print dloglike_alpha_beta

    # closure_relations = {'Synchrotron':PlotGRBClosureRelations.ClosureRelation(alpha=PlotGRBClosureRelations.DICT_ALPHA['Synchrotron']['ISM']['Slow']['Highest-E'], beta=PlotGRBClosureRelations.DICT_BETA['Synchrotron']['ISM']['Slow']['Highest-E'], name='Synchrotron Highest-E'),
    #                      'SSC':PlotGRBClosureRelations.ClosureRelation(alpha=PlotGRBClosureRelations.DICT_ALPHA['SSC']['ISM']['Slow']['2nd highest-E'], beta=PlotGRBClosureRelations.DICT_BETA['Synchrotron']['ISM']['Slow']['2nd highest-E'], name='SSC 2nd highest-E')}
    # for clrel in closure_relations.values():
    #     clrel.draw(ax)
    flag_cbar = False
    cbar = None
    for iax,cb in enumerate(('ISM', 'Wind')):
        for jax,coo in enumerate(('Fast', 'Slow')):
            str_title = '{0}, {1}-cooling'.format(cb, coo)
            for em in ('Synchrotron', 'SSC'):
                for eseg, formula in PlotGRBClosureRelations.DICT_ALPHA[em][cb][coo].items():
                    if eseg in ('1st HE', '2nd HE', '1st HE (IC-dom)'):
                        str_name = """{em}
{eseg}""".format(em=em[:3], eseg=eseg)
                        clrel = PlotGRBClosureRelations.ClosureRelation(alpha=formula, beta=PlotGRBClosureRelations.DICT_BETA[em][cb][coo][eseg], emission=em, cbmprof=cb, cooling=coo, esegment=eseg, name=str_name)
                        clrel.draw(ax[iax][jax])
                        if em is 'Synchrotron' and coo is 'Fast' and eseg in ('1st HE', '2nd HE'):
                            str_name = """{em}
{eseg}""".format(em=em[:3], eseg=eseg)
                            clrel = PlotGRBClosureRelations.ClosureRelation(alpha=PlotGRBClosureRelations.DICT_ALPHA[em][cb]['Radiative'][eseg], beta=PlotGRBClosureRelations.DICT_BETA[em][cb]['Radiative'][eseg], emission=em, cbmprof=cb, cooling='Radiative', esegment=eseg, name=str_name)
                            clrel.draw(ax[iax][jax])
                        if flag_cbar==False and (coo!='Fast' and eseg!='2nd HE'):
                            cbar = plt.colorbar(clrel.im, cax = cbaxes, orientation="horizontal")
                            cbar.set_label('p of the electrons')
                            #cbar = fig.colorbar(clrel.im)
                            #cbar.set_label('p of the electrons')
                            flag_cbar = True
            cont_alpha_beta = ax[iax][jax].contour(-alpha_mesh[:,:,0], -1.-beta_mesh[:,:,0], dloglike_alpha_beta*2., levels=cont_levels_NDF3, colors='k', linestyles=pMatplot.TPL_LINE)
            if mwlindices is not None: #and coo is 'Slow':
                mwlobs = PlotGRBClosureRelations.ObservedIndices()
                mwlobs.read(mwlindices, instruments=['XRT'])
                mwlobs.draw(ax[iax][jax])
            ax[iax][jax].set_title(str_title)
            ax[iax][jax].set_xlabel(r'Temporal index $\alpha$')
            ax[iax][jax].set_ylabel(r'Spectral index $\beta$')
            ax[iax][jax].set_xlim((-0.25, 3.25))
            ax[iax][jax].set_ylim((-0.5, 1.75))
            ax[iax][jax].grid()

    for ff in ('png','pdf'):
        path_save = "{dire}/scanned3D_{targ}{suff}.{form}".format(dire=outdir, targ=dict_summary['config']['name'], suff=suffix, form=ff)
        fig.savefig(path_save)
        logger.info('{0} has been saved.'.format(path_save))
    fig.clf()

    dict_norm_alpha_allowd = {}
    for sgnf in [1.0]:
        dict_norm_alpha_allowd[sgnf] = dloglike_norm_alpha*2. <= TABLE_DLOGLIKE_SIGNIFICANCE[str(sgnf)][2] # NDF=3
    dict_summary['scan3D']['bowtie'] = get_lightcurve_bowtie(norm_mesh[0,:,:], alpha_mesh[0,:,:], dict_norm_alpha_allowd, (t0_min, t1_max), tnorm)

    dict_summary['scan3D']['allowed_range'] = {}
    dict_summary['scan3D']['allowed_range']['norm'] = {}
    dict_summary['scan3D']['allowed_range']['beta'] = {}
    dict_summary['scan3D']['allowed_range']['alpha'] = {}
    for sgnf in [1.0, 2.0]:
        norm_mesh_allowed = norm_mesh[dloglike_mesh*2.<=TABLE_DLOGLIKE_SIGNIFICANCE[str(sgnf)][2]] # NDF=3
        beta_mesh_allowed = beta_mesh[dloglike_mesh*2.<=TABLE_DLOGLIKE_SIGNIFICANCE[str(sgnf)][2]] # NDF=3
        alpha_mesh_allowed = alpha_mesh[dloglike_mesh*2.<=TABLE_DLOGLIKE_SIGNIFICANCE[str(sgnf)][2]] # NDF=3
        dict_summary['scan3D']['allowed_range']['norm'][sgnf] = (np.nanmin(norm_mesh_allowed), np.nanmax(norm_mesh_allowed))
        dict_summary['scan3D']['allowed_range']['beta'][sgnf] = (np.nanmin(beta_mesh_allowed), np.nanmax(beta_mesh_allowed))
        dict_summary['scan3D']['allowed_range']['alpha'][sgnf] = (np.nanmin(alpha_mesh_allowed), np.nanmax(alpha_mesh_allowed))

    return dict_summary

        
def lightcurve_prefactor(tmin, tmax, integral, index, tb_gti, torigin, tnorm=10., rescaler=1.): #, torigin=0): #, tend=100000
    logger.debug('Integral:{0}'.format(integral))
    logger.debug('LC index:{0}'.format(index))

    powerlaw = lambda t, alpha: pow(t/tnorm, alpha)
    itgl = 0.
    weight = 0.
    
    ti0, ti1 = tb_gti['START']-torigin, tb_gti['STOP']-torigin
    tdiff = tb_gti['STOP']-tb_gti['START']
    if index!=-1:
        #tiref = (index+1) / (index+2) * (pow(ti1, index+2)-pow(ti0, index+2)) / (pow(ti1, index+1)-pow(ti0, index+1))
        tiref = pow((pow(ti1,index+1) + pow(ti0,index+1))/2., 1./(index+1))
    else:
        #tiref = (ti1-ti0) / (np.log(ti1/ti0))
        tiref = np.exp((np.log(ti1)+np.log(ti0))/2.)
    amp = powerlaw(tiref, index) * tdiff
    itgl = np.sum(tiref * amp)
    weight = np.sum(amp)
    tref = itgl / weight
    return (integral * pow(tref/tnorm, index) * rescaler), tref


def make_lightcurves(name, wholephase, emin, emax, eref, roi, ngoodstat, rgoodstat, ntbinperdecade, suffix, grbcatalogue, refit, force, outdir, index, redshift=False, tmin=0, tmax=10000, addphoton=None, calonly=tuple([None]*4), synccutoff=None):
    if synccutoff is not None and synccutoff is not False and synccutoff==synccutoff:
        if redshift==0 or redshift!=redshift:
            sptype = 'ExpCutoff'
            sppars = {'Prefactor':1E-10, 'Index':-2., 'Scale':1e3, 'Ebreak':5e4, 'P1':5e4, 'P2':0, 'P3':0}
        else:
            sptype = 'EblAtten::ExpCutoff'
            sppars = {'Prefactor':1E-10, 'Index':-2., 'Scale':1e3, 'Ebreak':5e4, 'P1':5e4, 'P2':0, 'P3':0, 'tau_norm':1., 'redshift':redshift, 'ebl_model':4}
    else:
        if redshift==0 or redshift!=redshift:
            sptype = 'PowerLaw2'
            sppars = {'Integral':1E-7, 'Index':-2., 'LowerLimit':emin, 'UpperLimit':emax}
        else:
            sptype = 'EblAtten::PowerLaw2'
            sppars = {'Integral':1E-7, 'Index':-2., 'LowerLimit':emin, 'UpperLimit':emax, 'tau_norm':1., 'redshift':redshift, 'ebl_model':4}
        
    lc = LightCurveGRB(name=name, wholephase=wholephase, tmin=tmin, tmax=tmax, emin=emin, emax=emax, eref=eref, deg_roi=roi, ngoodstat=ngoodstat, rlim_goodstat=rgoodstat, ntbinperdecade=ntbinperdecade, suffix=suffix, grbcatalogue=grbcatalogue, refit=refit, force=force, outdir=None, spectraltype=sptype, spectralpars=sppars, zmax=100., synccutoff=synccutoff)
    lc.setup()
    if tuple(calonly) != tuple([None]*4):
        logger.info('CalOnly data:')
        logger.info(calonly)
        lc.add_calonly(calonly)
    else:
        logger.warning('No CalOnly data.')
    lc.count_energy()
    lc.run_analysis(use_calonly=True if calonly!=tuple([None]*4) else False)
    lc.pickle()
    plot_lightcurves(lc.dct_stored, outdir=outdir, index=index, addphoton=addphoton, redshift=redshift)


def dataframe(dct_summary, outdir=None, ts_threshold=4.0, index='free'):
    if isinstance(dct_summary, basestring) and dct_summary[-7:]=='.pickle':
        dct_summary = pickle_utilities.load(dct_summary)
    elif not isinstance(dct_summary, dict):
        logger.critical('The input {0} was NOT a dictionary or path of pickle file!!!'.format(dct_summary))
        sys.exit(1)

    results = dct_summary['results']
    dct_frame = {'time_min': [],
                 'time_max': [],
                 'TS': [], 
                 'flux[cm^-2 s^-1]': [],
                 'flux_err_lo': [],
                 'flux_err_hi': [],
                 'eflux[erg cm^-2 s^-1]': [],
                 'eflux_err_lo': [],
                 'eflux_err_hi': [],
                 'e2dnde[erg cm^-2 s^-1]': [],
                 'e2dnde_err_lo': [],
                 'e2dnde_err_hi': [],
                 'Index': [],
                 'Index_err': []}
    for period in results:
        if period['TS']>=ts_threshold:
            dct_frame['time_min'].append(period['time']['min'])
            dct_frame['time_max'].append(period['time']['max'])
            dct_frame['TS'].append(period['TS'])
            dct_frame['flux[cm^-2 s^-1]'].append(period['limits'][index]['flux']['x0'])
            dct_frame['flux_err_lo'].append(period['limits'][index]['flux']['err_lo'])
            dct_frame['flux_err_hi'].append(period['limits'][index]['flux']['err_hi'])
            dct_frame['eflux[erg cm^-2 s^-1]'].append(period['limits'][index]['eflux']['x0']*MEVtoERG)
            dct_frame['eflux_err_lo'].append(period['limits'][index]['eflux']['err_lo']*MEVtoERG)
            dct_frame['eflux_err_hi'].append(period['limits'][index]['eflux']['err_hi']*MEVtoERG)
            dct_frame['e2dnde[erg cm^-2 s^-1]'].append(period['limits'][index]['e2dnde']['x0']*MEVtoERG)
            dct_frame['e2dnde_err_lo'].append(period['limits'][index]['e2dnde']['err_lo']*MEVtoERG)
            dct_frame['e2dnde_err_hi'].append(period['limits'][index]['e2dnde']['err_hi']*MEVtoERG)
            dct_frame['Index'].append(period['Index']['value'])
            dct_frame['Index_err'].append(period['Index']['error'])

    df = pd.DataFrame(dct_frame)
    outdir = outdir if outdir is not None else '{base}/{target}/E{emin:0>7.0f}-{emax:0>7.0f}MeV/r{roi:0>2.0f}deg/{phase}'.format(base=pLATLikelihoodConfig.PATH_BASEDIR, target=str(dct_summary['config']['name']), emin=dct_summary['config']['energy']['min'], emax=dct_summary['config']['energy']['max'], roi=dct_summary['config']['roi']['radius'], phase='lightcurve')
    outbasename = 'LightCurve_{target}_index{idx}{suffix}'.format(target=str(dct_summary['config']['name']), idx=index, suffix=str(dct_summary['config']['suffix']))
    df.to_csv('{dire}/{name}.csv'.format(dire=outdir, name=outbasename))


def plot_lightcurves(dct_summary, outdir=None, ts_threshold=4.0, index='free', grbcatalogue=pLATLikelihoodConfig.GRB_CATALOGUE_LAT, addphoton=None, fitlc=False, fit_phase='afterglow', mwlindices=None, redshift=False):

    # Plotting light curve composite
    dct_curves = None
    path_dct_summary = None
    if isinstance(dct_summary, basestring) and dct_summary[-7:]=='.pickle':
        if dct_summary[-18:]=='plotmediate.pickle':
            dct_mediate = pickle_utilities.load(dct_summary)
            logger.info(dct_mediate.items())
            dct_curves = dct_mediate['curves']
            dct_summary = dct_mediate['summary']
        else:
            path_dct_summary = dct_summary
            dct_summary = pickle_utilities.load(dct_summary)
    elif not isinstance(dct_summary, dict):
        logger.critical('The input {0} was NOT a dictionary or path of pickle file!!!'.format(dct_summary[-7:]))
        sys.exit(1)

    tb_lat = ReadLATCatalogueInfo.open_table(grbcatalogue)
    tb_gbm = ReadGBMCatalogueInfo.open_table()
    tb_one = ReadLATCatalogueInfo.select_one_by_name(tb_lat, str(dct_summary['config']['name']), tb_gbm)

    # Time-evolution
    gbm_fluence = 1E-5
    if 'GBM'in tb_one:
        if tb_one['GBM'] is not None:
            if 'FLUENCE' in tb_one['GBM']:
                gbm_fluence = tb_one['GBM']['FLUENCE']
    if fit_phase=='afterglow':
        if dct_summary['config']['name'] == '130907904':
            tmin = 360
        else:
            tmin = tb_one['GBM']['T90'] + tb_one['GBM']['T90_START']
        tmax = 100000.
    else:
        tmin = 0.
        tmax = 100000.
    if dct_summary['config']['name'] == '130907904':
        t_triggered = 400282876.000
    else:
        t_triggered = pMETandMJD.ConvertMjdToMet(float(tb_one['GBM']['TRIGGER_TIME']))
    dct_summary = scan_norm_beta_alpha(dict_summary=dct_summary, torigin=t_triggered, outdir=outdir, norms=gbm_fluence*(10**np.linspace(-3, 2, 51)), betas=np.linspace(-2.75, -0.5, 63), alphas=np.linspace(-2.5, 0.0, 76), tnorm=100., tmin=tmin, tmax=tmax, suffix=str(dct_summary['config']['suffix']), mwlindices=mwlindices) #GRB 160509374
    #dct_summary = scan_norm_beta_alpha(dict_summary=dct_summary, torigin=t_triggered, outdir=outdir, norms=gbm_fluence*(10**np.linspace(-3, 2, 51)), betas=np.linspace(-2.5, -1.25, 63), alphas=np.linspace(-2.0, -1.0, 76), tnorm=100., tmin=tmin, tmax=tmax, suffix=str(dct_summary['config']['suffix']), mwlindices=mwlindices) #GRB 090926181
    bowtie_times, bowtie_curves = dct_summary['scan3D']['bowtie'][0], dct_summary['scan3D']['bowtie'][1]

    #print path_dct_summary
    #pickle_utilities.dump(path_dct_summary, dct_summary)

    outbasename = 'LightCurve_{target}_index{idx}{suffix}'.format(target=str(dct_summary['config']['name']), idx=index, suffix=str(dct_summary['config']['suffix']))
    pickle_utilities.dump('{0}/{1}.pickle'.format(outdir, outbasename), dct_summary)

    if dct_curves is None:
    # Config
        str_energies = '{emin:3.3f} - {emax:3.0f}'.format(emin=dct_summary['config']['energy']['min'], emax=dct_summary['config']['energy']['max'])
    #dct_summary['fit'] = {'flux':{}}

    # Starting point of fitting

    # Characteristices
        dct_curves = OrderedDict() #{}

        #dct_curves['TS'] = pMatplot.Curve('TS', xlabel='Time - T0 [s]', ylabel=r'$\sqrt{\rm{max}(TS, 0)}$', xerr_asym=True, yerr_asym=False)

        dct_curves['flux'] = pMatplot.Curve('flux', xlabel='Time - T0 [s]', ylabel=r'Photon flux $\mathrm{[/cm^2 s]}$', xerr_asym=True, yerr_asym=True)
        dct_curves['flux_ul'] = pMatplot.Curve('flux', xlabel='Time - T0 [s]', ylabel=r'Photon flux $\mathrm{[/cm^2 s]}$', xerr_asym=True, yerr_asym=False, ul=True)

        #dct_curves['eflux'] = pMatplot.Curve('eflux', xlabel='Time - T0 [s]', ylabel=r'Energy flux $\mathrm{[erg/cm^2 s]}$', xerr_asym=True, yerr_asym=True)
        #dct_curves['eflux_ul'] = pMatplot.Curve('eflux', xlabel='Time - T0 [s]', ylabel=r'Energy flux $\mathrm{[erg/cm^2 s]}$', xerr_asym=True, yerr_asym=False, ul=True)

        #dct_curves['e2dnde'] = pMatplot.Curve('e2dnde', xlabel='Time - T0 [s]', ylabel=r'$E^2 dN/dE \, \rm{{at}} \, {ene:3.1f} \rm{{GeV}} \, \mathrm{{[erg/cm^2 s]}}$'.format(ene=dct_summary['config']['energy']['ref']/1000.), xerr_asym=True, yerr_asym=True)
        #dct_curves['e2dnde_ul'] = pMatplot.Curve('e2dnde', xerr_asym=True, yerr_asym=False, ul=True, xlabel='Time - T0 [s]', ylabel=r'$E^2 dN/dE \, \rm{{at}} \, {ene:3.1f} \rm{{GeV}} \, \mathrm{{[erg/cm^2 s]}}$'.format(ene=dct_summary['config']['energy']['ref']/1000.)) 

        dct_curves['Index'] = pMatplot.Curve('Index', xlabel='Time - T0 [s]', ylabel='Spectral index', xerr_asym=True, yerr_asym=True)
        dct_curves['Index_ul'] = pMatplot.Curve('Index', xlabel='Time - T0 [s]', ylabel='Spectral index', xerr_asym=True, yerr_asym=False, ul=True)
        dct_curves['Index_ll'] = pMatplot.Curve('Index', xlabel='Time - T0 [s]', ylabel='Spectral index', xerr_asym=True, yerr_asym=False, ll=True)
    
        for ic, curve in enumerate(dct_curves.values()):
            logger.info('===== {0} ====='.format(curve.quantity))
            for period in dct_summary['results']:
                t0 = max(1, period['time']['min'])
                t1 = period['time']['max']
                tref = 10**((np.log10(t0)+np.log10(t1))/2.0) #(t0+t1)/2.0 #sqrt(t0*t1)
                logger.info('----- {tmin:.1f} - {tmax:.1f} s -----'.format(tmin=t0, tmax=t1))            
                dloglike = period['dloglike']['dloglike']
                dict_meshes_allowed = {}
                for sgnf in [1.0, 2.0, 3.0]:
                    dict_meshes_allowed[sgnf] = dloglike*2. <= TABLE_DLOGLIKE_SIGNIFICANCE[str(sgnf)][1] # NDF=2
                indices = period['dloglike']['Index']
                norms = period['dloglike']['normalization']
                indices_mesh, norms_mesh = np.meshgrid(indices, norms)
                dloglikes_mesh = period['dloglike']['dloglike']
                flag_detect = True if dloglikes_mesh[0][0]*2>TABLE_DLOGLIKE_SIGNIFICANCE['2.0'][1] else False

                if period['TS']!=period['TS']:
                    logger.warning('Ananlysis result is NaN! Skipping...')
                if curve.quantity == 'TS':
                    y = sqrt(max(0, period[curve.quantity]))
                    yerr = 0
                    logger.info('{v:.2}'.format(v=y))
                    curve.set_point(tref, y, {'hi':t1-tref, 'lo':tref-t0}, yerr)
                elif curve.quantity in ('Index', 'flux', 'eflux', 'e2dnde'):
                    y_values = {}

                    if curve.quantity in ('flux', 'eflux', 'e2dnde'):
                        for sgnf in [1.0, 2.0, 3.0]:
                            y_values[sgnf] = set(period['dloglike'][curve.quantity][dict_meshes_allowed[sgnf]])
                    elif curve.quantity in ('Index',):
                        for sgnf in [1.0, 2.0, 3.0]:
                            y_values[sgnf] = set(indices_mesh[dict_meshes_allowed[sgnf]])
                            
                    if curve.ul==False and curve.ll==False and flag_detect==True:
                        logger.info('Pointed')
                        curve.set_point(tref, period['dloglike']['best'][curve.quantity]*(MEVtoERG if curve.quantity in ('eflux', 'e2dnde') else 1), {'hi':t1-tref, 'lo':tref-t0}, {'hi': (max(y_values[1.0])-period['dloglike']['best'][curve.quantity])*(MEVtoERG if curve.quantity in ('eflux', 'e2dnde') else 1), 'lo': (period['dloglike']['best'][curve.quantity]-min(y_values[1.0]))*(MEVtoERG if curve.quantity in ('eflux', 'e2dnde') else 1)})
                    elif curve.ul==True and curve.ll==False and flag_detect==False:
                        logger.info('Upper limit')
                        curve.set_point(tref, max(y_values[2.0])*(MEVtoERG if curve.quantity in ('eflux', 'e2dnde') else 1), {'hi':t1-tref, 'lo':tref-t0})
                    elif curve.ul==False and curve.ll==True and flag_detect==False:
                        logger.info('Lower limit')
                        curve.set_point(tref, min(y_values[2.0])*(MEVtoERG if curve.quantity in ('eflux', 'e2dnde') else 1), {'hi':t1-tref, 'lo':tref-t0})
            #    curve.make_meshes()
        dct_plotmediate = {'summary':dct_summary, 'curves':dct_curves}
        pickle_utilities.dump('{0}/{1}_plotmediate.pickle'.format(outdir, outbasename), dct_plotmediate)

    NPLOTS = 1+sum([int(v.ul==False and v.ll==False) for v in dct_curves.values()])
    fig, ax = plt.subplots(NPLOTS, 1, figsize=(10, 16./6.*NPLOTS), sharex=True)
    iax = -1
    ax[0].set_title('GRB '+str(dct_summary['config']['name']))
    ax[0].set_xlim((1.0, 100000.0))

    if 'TS' in dct_curves:
        iax+=1
        ax[iax].errorbar(dct_curves['TS'].get_xdata(), dct_curves['TS'].get_ydata(), xerr=dct_curves['TS'].get_xerr(), yerr=dct_curves['TS'].get_yerr(), fmt=dct_curves['TS'].fmt)
        ax[iax].set_ylabel(dct_curves['TS'].ylabel)
        ax[iax].set_xscale("log", nonposx='clip')
        ax[iax].grid(ls='-', lw=0.5, alpha=0.5)

    if 'flux' in dct_curves:
        iax+=1
        ax[iax].errorbar(dct_curves['flux'].get_xdata(), dct_curves['flux'].get_ydata(), xerr=dct_curves['flux'].get_xerr(), yerr=dct_curves['flux'].get_yerr(), fmt=dct_curves['flux'].fmt, ms=2)
        ax[iax].errorbar(dct_curves['flux_ul'].get_xdata(), dct_curves['flux_ul'].get_ydata(), xerr=dct_curves['flux_ul'].get_xerr(), fmt=dct_curves['flux_ul'].fmt)

        for name_shown in bowtie_times.keys():
            times = bowtie_times[name_shown]
            curve_lo = bowtie_curves[0][name_shown]
            curve_hi = bowtie_curves[1][name_shown]
            if times is not None:
                ax[iax].fill_between(times, curve_lo, curve_hi, alpha=0.2, facecolor='r' if name_shown=='1.0' else 'b') #, label=name_shown)

        ax[iax].set_ylabel(dct_curves['flux'].ylabel)
        ax[iax].set_xscale("log", nonposx='clip')
        ax[iax].set_yscale("log", nonposx='clip')
        ax[iax].grid(ls='-', lw=0.5, alpha=0.5)
        logger.debug('Axis label:')
        logger.debug(ax[iax].get_yticks())
        ax[iax].set_yticks([y for y in ax[iax].get_yticks() if y<0.5*ax[iax].get_ylim()[1]])
    # Fit
#     if fitlc==True:
#         dct_fit = {'fit':{'flux':{'lightcurve':{}}}}
#         flux_max = dct_curves['flux'].get_maximum()
#         if flux_max!=0:
#             t_flux_max = flux_max[1]
#             fit_result = dct_curves['flux'].fit_lin(t_flux_max, t_flux_max)
#             if isinstance(fit_result, tuple) and len(fit_result)==2:
#                 params, butterfly = fit_result[0], fit_result[1]
#                 for iparam, param, in enumerate(params[0]):
#                     logger.info("""Parameter {0}: {1} +/- {2}""".format(iparam, params[0][iparam], params[1][iparam]))
#                 str_powerlaw_fitted = r"""$ \displaystyle F(t) = \rm{{ F_{{0}} }} \left( \frac{{ \it{{ t }} }}{{ {t_scale:.1E} }} \right)^{{- \rm{{ \alpha }} }}$
# $\rm{{ F_{{0}} = {f0:.2E} \pm {f0err:.1E} }}$
# $\rm{{ \alpha = {alpha:.2E} \pm {alphaerr:.2E} }}$""".format(f0=params[0][0], f0err=params[1][0], alpha=params[0][1], alphaerr=params[1][1], t_scale=t_flux_max)
#                 ax[iax].plot(butterfly[0], butterfly[1], c='g')
#                 ax[iax].text(butterfly[0][0]*10, butterfly[1][0]/5., str_powerlaw_fitted)
#                 dct_fit['fit']['flux']['lightcurve']['amplitude'] = {'value':params[0][0], 'error':params[1][0]}
#                 dct_fit['fit']['flux']['lightcurve']['index'] = {'value':params[0][1], 'error':params[1][1]}
#                 dct_fit['fit']['flux']['lightcurve']['t_scale'] = t_flux_max #t95
#                 dct_fit['fit']['flux']['lightcurve']['tmin'] = t_flux_max #t95
#                 dct_fit['fit']['flux']['lightcurve']['tmax'] = None
#                 pickle_utilities.dump('{0}/{1}_fit.pickle'.format(outdir, outbasename), dct_fit)

    if 'eflux' in dct_curves:
        iax+=1
        ax[iax].errorbar(dct_curves['eflux'].get_xdata(), dct_curves['eflux'].get_ydata(), xerr=dct_curves['eflux'].get_xerr(), yerr=dct_curves['eflux'].get_yerr(), fmt=dct_curves['eflux'].fmt)
        ax[iax].errorbar(dct_curves['eflux_ul'].get_xdata(), dct_curves['eflux_ul'].get_ydata(), xerr=dct_curves['eflux_ul'].get_xerr(), fmt=dct_curves['eflux_ul'].fmt)
        ax[iax].set_ylabel(dct_curves['eflux'].ylabel)
        ax[iax].set_xscale("log", nonposx='clip')
        ax[iax].set_yscale("log", nonposx='clip')
        ax[iax].grid(ls='-', lw=0.5, alpha=0.5)
        logger.debug('Axis label:')
        logger.debug(ax[iax].get_yticks())
        ax[iax].set_yticks([y for y in ax[iax].get_yticks() if y<0.5*ax[iax].get_ylim()[1]])

    if 'e2dnde' in dct_curves:
        iax+=1
        ax[iax].errorbar(dct_curves['e2dnde'].get_xdata(), dct_curves['e2dnde'].get_ydata(), xerr=dct_curves['e2dnde'].get_xerr(), yerr=dct_curves['e2dnde'].get_yerr(), fmt=dct_curves['e2dnde'].fmt)
        ax[iax].errorbar(dct_curves['e2dnde_ul'].get_xdata(), dct_curves['e2dnde_ul'].get_ydata(), xerr=dct_curves['e2dnde_ul'].get_xerr(), fmt=dct_curves['e2dnde_ul'].fmt)
        ax[iax].set_ylabel(dct_curves['e2dnde'].ylabel)
        ax[iax].set_xscale("log", nonposx='clip')
        ax[iax].set_yscale("log", nonposx='clip')
        ax[iax].grid(ls='-', lw=0.5, alpha=0.5)
        logger.debug('Axis label:')
        logger.debug(ax[iax].get_yticks())
        ax[iax].set_yticks([y for y in ax[iax].get_yticks() if y<0.5*ax[iax].get_ylim()[1]])

    if 'Index' in dct_curves:
        iax+=1
        ax[iax].errorbar(dct_curves['Index'].get_xdata(), dct_curves['Index'].get_ydata(), xerr=dct_curves['Index'].get_xerr(), yerr=dct_curves['Index'].get_yerr(), fmt=dct_curves['Index'].fmt)
        #ax[iax].errorbar(dct_curves['Index_ul'].get_xdata(), dct_curves['Index_ul'].get_ydata(), xerr=dct_curves['Index_ul'].get_xerr(), yerr=dct_curves['Index_ul'].get_yerr(), fmt=dct_curves['Index_ul'].fmt)
        #ax[iax].errorbar(dct_curves['Index_ll'].get_xdata(), dct_curves['Index_ll'].get_ydata(), xerr=dct_curves['Index_ll'].get_xerr(), yerr=dct_curves['Index_ll'].get_yerr(), fmt=dct_curves['Index_ll'].fmt)
        for sgnf, beta_range in dct_summary['scan3D']['allowed_range']['beta'].items():
            if sgnf==2.0:
                ax[iax].fill_between(x=[bowtie_times.values()[0][0], bowtie_times.values()[0][-1]], y1=[beta_range[0]]*2, y2=[beta_range[1]]*2, alpha=0.2, facecolor='c')
        ax[iax].set_ylabel(dct_curves['Index'].ylabel)
        ax[iax].set_xlabel(dct_curves['Index'].xlabel)
        ax[iax].set_xscale("log", nonposx='clip')
        ax[iax].grid(ls='-', lw=0.5, alpha=0.5)
        logger.debug('Axis label:')
        logger.debug(ax[iax].get_yticks())
        ax[iax].set_yticks([y for y in ax[iax].get_yticks() if y<ax[iax].get_ylim()[1]-0.1])

    iax+=1
    im_angsep = ax[iax].scatter(dct_summary['counts']['time'], dct_summary['counts']['energy'], alpha=0.5, c=np.clip(dct_summary['counts']['angsep'], 0.0, 6.0), cmap=cm.Greens_r, vmin=0, vmax=6)
    if 'counts_CalOnly' in dct_summary:
        if 'time' in dct_summary['counts_CalOnly'] and 'energy' in dct_summary['counts_CalOnly'] and 'angsep' in dct_summary['counts_CalOnly']:
            if isinstance(dct_summary['counts_CalOnly']['time'], np.ndarray) and isinstance(dct_summary['counts_CalOnly']['energy'], np.ndarray) and isinstance(dct_summary['counts_CalOnly']['angsep'], np.ndarray):
                if len(dct_summary['counts_CalOnly']['time'])*len(dct_summary['counts_CalOnly']['energy'])*len(dct_summary['counts_CalOnly']['angsep'])>0:
                    ax[iax].scatter(dct_summary['counts_CalOnly']['time'], dct_summary['counts_CalOnly']['energy'], alpha=0.5, c=np.clip(dct_summary['counts_CalOnly']['angsep'], 0.0, 6.0), cmap=cm.Greens_r, marker='D', s=40, vmin=0, vmax=6, edgecolors='red', linewidth=1)
    ax[iax].set_xlabel('Time - T0 [s]')
    ax[iax].set_ylabel('Energy [MeV]')
    if len(dct_summary['counts']['time'])>0 or addphoton is not None:
        ax[iax].set_yscale("log", nonposx='clip')
    ax[iax].grid(ls='-', lw=0.5, alpha=0.5)
    ax[iax].set_yticks([y for y in ax[iax].get_yticks() if y<0.5*ax[iax].get_ylim()[1]])
    ax[iax].set_ylim(bottom=max(300,dct_summary['config']['energy']['min']), top=1.1*dct_summary['config']['energy']['max'])
    if (not redshift in (None, False)) and redshift==redshift:
        ax_intrinsic_energy = ax[iax].twinx()
        ax_intrinsic_energy.set_yscale("log", nonposx='clip')
        ax_intrinsic_energy.set_ylim(np.array(ax[iax].get_ylim()) * (1.+redshift))
        ax_intrinsic_energy.set_ylabel('Energy * (1+{0:1.2f}) [MeV]'.format(redshift), color='magenta')
        ax_intrinsic_energy.grid(ls='-', lw=0.5, alpha=0.5, axis='y', c='magenta')
        ax_intrinsic_energy.tick_params(axis='y', colors='magenta')
        # ax_intrinsic_time = ax[0].twiny()
        # ax_intrinsic_time.set_yscale("log", nonposx='clip')
        # ax_intrinsic_time.set_xlim(np.array(ax[0].get_xlim()) / (1.+redshift))        
    ax_cbar_angsep = fig.add_axes([0.1, 1./(1.+iax), 0.18, 0.02]) 
    cbar_angsep = plt.colorbar(im_angsep, cax=ax_cbar_angsep, orientation="horizontal", ticks=[6.0, 4.0, 2.0, 0.0])
    cbar_angsep.set_label('Angular separation [deg]')

    fig.tight_layout() #subplots_adjust(hspace=0)
    fig.subplots_adjust(hspace=0)
    outdir = outdir if outdir is not None else '{base}/{target}/E{emin:0>7.0f}-{emax:0>7.0f}MeV/r{roi:0>2.0f}deg/{phase}'.format(base=pLATLikelihoodConfig.PATH_BASEDIR, target=str(dct_summary['config']['name']), emin=dct_summary['config']['energy']['min'], emax=dct_summary['config']['energy']['max'], roi=dct_summary['config']['roi']['radius'], phase='lightcurve')
    for ff in ['png', 'pdf']:
        fig.savefig('{0}/{1}.{2}'.format(outdir, outbasename, ff))
    fig.clf()


def get_lightcurve_bowtie(norm_mesh, indices_mesh, dict_meshes_shown, trange, tscale, ntperdec=20):
    logt_min = np.log10(trange[0])
    logt_max = np.log10(trange[1])
    tvals = 10 ** np.linspace(logt_min, logt_max, int((logt_max-logt_min)*ntperdec+1))
    def flux(t, norm, index):
        return norm*pow(t/tscale, index)

    odict_curves_lo = OrderedDict()
    odict_curves_hi = OrderedDict()
    odict_times = OrderedDict()
    list_flux_maps = [flux(t, norm_mesh, indices_mesh) for t in tvals]
            
    for name_shown,shown in dict_meshes_shown.items():
        bound_lo = []
        bound_hi = []
        times = []
        for jt,t in enumerate(tvals):
            flux_shown = list_flux_maps[jt][shown]
            if len(flux_shown)>0:
                bound_lo.append(min(flux_shown.flatten())) 
                bound_hi.append(max(flux_shown.flatten())) 
                times.append(t)
            if len(times)>0:
                odict_curves_lo[name_shown] = np.array(bound_lo)
                odict_curves_hi[name_shown] = np.array(bound_hi)
                odict_times[name_shown] = np.array(times)
            else:
                odict_curves_lo[name_shown] = None
                odict_curves_hi[name_shown] = None
                odict_energies[name_shown] = None
        return (odict_times, (odict_curves_lo, odict_curves_hi))


@click.command()
@click.option('--namemin', type=str, default='000000000')
@click.option('--namemax', type=str, default='999999999')
@click.option('--mode', '-m', type=click.Choice(['unified', 'prompt', 'primary', 'intermittent', 'afterglow', 'earlyAG', 'lateAG', 'farAG', 'T95to01ks', 'T95to03ks', '01ksto10ks', '01ksto100ks', 'T95to03ks', '03ksto10ks', '03ksto100ks', 'lightcurve', 'special']))
@click.option('--grbcatalogue', '-c', type=str, default=pLATLikelihoodConfig.GRB_CATALOGUE_LAT)
@click.option('--emin', type=float, default=100.)
@click.option('--emax', type=float, default=100000.)
@click.option('--eref', type=float, default=1000.)
@click.option('--tmin', type=float, default=0.)
@click.option('--tmax', type=float, default=10000.)
@click.option('--roi', type=float, default=12.)
@click.option('--ngoodstat', type=int, default=10)
@click.option('--rgoodstat', type=int, default=2)
@click.option('--ntbinperdecade', type=float, default=5.0, help='Number of time bins in onde decade. This is active only if --ngoodstat 0.')
@click.option('--index', type=click.Choice(['free', 'best']), default='free')
@click.option('--redshift', '-z', type=float, default=0.)
@click.option('--calonly', type=(str, str, str, str), default=[None]*4, help='path_onevt, path_onexp, path_offevt, path_offexp, rclass')
@click.option('--suffix', '-s', type=str, default='')
@click.option('--force', '-f', is_flag=True)
@click.option('--refit', '-r', is_flag=True)
@click.option('--outdir', '-o', type=str, default='.')
@click.option('--plotonly', '-p', type=str, default=None, help='Path of result pickle file if you skip analyses.')
@click.option('--mwlindices', type=str, default=None, help='Path of a CSV file of the observed temporal/spectral indices in other wavelengh.')
@click.option('--synccutoff', type=float, default=None, help='Cutoff at synchrotron maximum energy at 1000 s in MeV')
@click.option('--bsub', '-b', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(namemin, namemax, mode, emin, emax, eref, tmin, tmax, roi, ngoodstat, rgoodstat, ntbinperdecade, refit, force, suffix, grbcatalogue, outdir, plotonly, index, redshift, calonly, mwlindices, synccutoff, bsub, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    tb_lat = ReadLATCatalogueInfo.open_table(grbcatalogue)
    tb_gbm = ReadGBMCatalogueInfo.open_table()
    if bsub==True:
        tb_lat = ReadLATCatalogueInfo.select_by_name(tb_lat, namemin, namemax, tb_gbm)
        #tb_lat = ReadLATCatalogueInfo.select_gbm_exist(tb_lat)
        if len(tb_lat)<1:
            print 'No GRBs.'
            return 1
        for (irow , row) in enumerate(tb_lat):
            name = row['GRBNAME']
            print '##### No.{0} GRB{1} #####'.format(irow, name)
            if not os.path.exists(name):
                os.mkdir(name)
            acmd = ['bsub', '-o','{0}/GRB{0}_lightcurve{1}.log'.format(name, suffix if suffix=='' else '_'+suffix), '-J','lc{0}'.format(name[:-3]), '-W','1200', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/STLikelihoodAnalysis/LightCurve.py', '--mode', mode, '--emin', str(emin), '--emax', str(emax), '--eref', str(eref), '--tmin', str(tmin), '--tmax', str(tmax), '-s', suffix, '--index', 'free', '--redshift', str(redshift), '--roi', str(roi), '--ngoodstat', str(ngoodstat), '--rgoodstat', str(rgoodstat), '--ntbinperdecade', str(ntbinperdecade), '--grbcatalogue', grbcatalogue, '--namemin', name, '--outdir', outdir]
            if force==True:
                acmd.append('--force')
            if refit==True:
                acmd.append('--refit')
            if plotonly is not None:
                acmd.append('--plotonly')
                acmd.append(plotonly)
            if synccutoff is not None:
                acmd.append('--synccutoff')
                acmd.append(str(synccutoff))
            if calonly != tuple([None]*4):
                acmd.append('--calonly')
                for i in range(4):
                    acmd.append(str(calonly[i]))
            if mwlindices is not None:
                acmd.append('--mwlindices')
                acmd.append(mwlindices)
            print acmd
            subprocess.call(acmd)
        return 0
    else:
        tb_lat = ReadLATCatalogueInfo.select_one_by_name(tb_lat, namemin, tb_gbm)

        # CalOnly photons
        addphoton = None

        if plotonly==None:
            make_lightcurves(name=namemin, wholephase=mode, emin=emin, emax=emax, eref=eref, tmin=tmin, tmax=tmax, roi=roi, ngoodstat=ngoodstat, rgoodstat=rgoodstat, ntbinperdecade=ntbinperdecade, suffix=suffix, grbcatalogue=grbcatalogue, refit=refit, force=force, outdir=outdir, index=index, redshift=redshift, addphoton=addphoton, calonly=calonly, synccutoff=synccutoff)#, tmin=tb_lat['LAT_TRIGGER_TIME']-tb_lat['TRIGGER_TIME'])
        else:
            plot_lightcurves(plotonly, outdir=outdir, index=index, addphoton=addphoton, mwlindices=mwlindices, redshift=redshift)


if __name__ == '__main__':
    main()
