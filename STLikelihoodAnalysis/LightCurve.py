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
import math
from math import log10, log, sqrt, ceil, isnan, pi, factorial
from sympy import *
from scipy import integrate
from astropy.io import fits
import click
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import pyLikelihood
from UnbinnedAnalysis import *
from BinnedAnalysis import *
import STLikelihoodAnalysis.pLATLikelihoodConfig as pLATLikelihoodConfig
from FindGoodstatPeriods import find_goodstat_periods, get_entries, get_event_time_and_energy, get_event_time_energy_angsep
from FindCrossEarthlimb import find_cross_earthlimb#, get_good_intervals
#from STLikelihoodAnalysis import get_module_logger
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
import ReadLATCatalogueInfo
import ReadGBMCatalogueInfo
import GetGTI
import pickle_utilities
import pMatplot

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
    def __init__(self, name, wholephase='special', tmin=0.0, tmax=10000.0, emin=100.0, emax=100000.0, eref=1000., evclass=128, evtype=3, ft2interval='30s', deg_roi=12., rad_margin=10., zmax=100., index_fixed=None, redshift=0.0, suffix='', grbcatalogue=pLATLikelihoodConfig.GRB_CATALOGUE_LAT, psForce=False, rlim_goodstat=2.0, ngoodstat=10, ntbinperdecade=4, refit=True, force=False, outdir=None, phase='lightcurve', spectraltype='PowerLaw', spectralpars={'Prefactor':1e-10, 'Index':-2, 'Scale':1000}):
        # Target GRB
        self.grb = pLATLikelihoodConfig.GRBTarget(name, grbcatalogue, spectraltype=spectraltype, spectralpars=spectralpars, redshift=redshift)
        #self.grb = pLATLikelihoodConfig.GRBTarget(name, grbcatalogue, spectraltype='ScaleFactor::PowerLaw2', spectralpars={'Integral':norm, 'Index':index, 'LowerLimit':emin, 'UpperLimit':emax, 'ScaleFactor':scalefactor})

        # Light curve instance
        LightCurve.__init__(self, name=name, met0=self.grb.met0, tmin=tmin, tmax=tmax, emin=emin, emax=emax, eref=eref, evclass=evclass, evtype=evtype, ft2interval=ft2interval, deg_roi=deg_roi, rad_margin=rad_margin, zmax=zmax, index_fixed=index_fixed, suffix=suffix, psForce=psForce, refit=refit, force=force, phase=phase)
        self.config['goodstat'] = {'n':ngoodstat, 'r':rlim_goodstat}
        self.config['ntbinperdecade'] = ntbinperdecade
        self.config['wholephase'] = wholephase
        self.periods_goodstat = []
        self.analyses = []
        self.counts_time = None
        self.counts_energy = None
        self.counts_angsep = None


    def setup(self):
        # Analysis instance which covers the whole time range
        self.analysis_whole = pLATLikelihoodConfig.GRBConfig(target=self.grb, phase=self.config['wholephase'], tstop=self.config['time']['max'], emin=self.config['energy']['min'], emax=self.config['energy']['max'], deg_roi=self.config['roi']['radius'], zmax=self.config['zenith']['max'], suffix=self.config['suffix'], tmin_special=self.config['time']['min'], tmax_special=self.config['time']['max'], ft2interval='1s')
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
        if len(self.periods_goodstat)<1:
            validtimes = find_cross_earthlimb(self.analysis_whole.path_ft2, self.analysis_whole.target.ra, self.analysis_whole.target.dec, self.analysis_whole.tmin, self.analysis_whole.tmax, self.analysis_whole.zmax, thcut=self.analysis_whole.thetamax, torigin=self.analysis_whole.target.met0, )
            if self.config['goodstat']['n']>0:
                logger.info("""Good time intervals (zenith < {zen} deg)
{vt}""".format(zen=self.config['zenith']['max'], vt=validtimes))
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

            else: # Logarithmically equivalent time binning
                # Prompt phase
                tpromptstart = self.analysis_whole.target.t05
                tpromptend = self.analysis_whole.target.t25 + self.analysis_whole.target.t50 # *3./2.
                if validtimes[0][0]<tpromptend and tpromptend<validtimes[0][1]:
                    if tpromptstart>validtimes[0][0] and tpromptstart<tpromptend:
                        self.periods_goodstat.append([validtimes[0][0], tpromptstart])
                        tp0, tp1 = tpromptstart, tpromptend
                    else:
                        tp0, tp1 = validtimes[0][0], tpromptend
                    self.periods_goodstat += find_goodstat_periods(self.analysis_whole.path_filtered, tp0, tp1, nthreshold=10, rlim=5., ra=self.grb.ra, dec=self.grb.dec, torigin=self.analysis_whole.target.met0)
                        
                    validtimes[0][0] = tpromptend # Forget about the prompt phase
                # Afterglow phase
                tbins = []
                validtimes_1 = []
                validtimes_not1 = []
                tbin_trans1 = 0.
                for vt in validtimes:
                    if np.log10(vt[1]/vt[0])>=1./self.config['ntbinperdecade']: #flag_1==True:
                        validtimes_1.append(vt)
                    else:
                        validtimes_not1.append(vt)

                logger.debug('validtimes_1: {0}'.format(validtimes_1))
                logger.debug('validtimes_not1: {0}'.format(validtimes_not1))

                for vt in validtimes_1:
                    logtmin_1 = np.log10(max(1, vt[0]))
                    logtmax_1 = np.log10(vt[1])
                    tedges_1 = 10**np.linspace(logtmin_1, logtmax_1, int(self.config['ntbinperdecade']*(logtmax_1-logtmin_1)+1.5))
                    for t0, t1 in zip(tedges_1[:-1],tedges_1[1:]):
                        tbins.append([t0, t1])
                    tbin_trans1 = max(tbin_trans1, vt[1])

                validtimes_3 = []
                validtimes_not3 = []
                for vt0, vt1 in zip(validtimes_not1[:-1], validtimes_not1[1:]):
                    if np.log10(vt1[1]/vt0[0])<1./self.config['ntbinperdecade']:
                        validtimes_3.append(vt0)
                        if vt1==validtimes_not1[-1]:
                            validtimes_3.append(vt1)
                    else:
                        validtimes_not3.append(vt0)                    
                        if vt1==validtimes_not1[-1]:
                            validtimes_not3.append(vt1)

                logger.debug('validtimes_3: {0}'.format(validtimes_3))
                logger.debug('validtimes_not3: {0}'.format(validtimes_not3))
                tbins += validtimes_not3

                logtmin_3 = np.log10(max(1, validtimes_3[0][0]))
                logtmax_3 = np.log10(validtimes_3[-1][1])
                tedges_3 = 10**np.linspace(logtmin_3, logtmax_3, int(self.config['ntbinperdecade']*(logtmax_3-logtmin_3)+1.5))
                for t0, t1 in zip(tedges_3[:-1], tedges_3[1:]):
                    tbins.append([t0, t1])

                self.periods_goodstat += tbins

        print """Time slots:
{0}""".format(self.periods_goodstat)
        for ip, pds in enumerate(self.periods_goodstat):
            logger.info('Time slot No.{0}: {1} - {2} s'.format(ip, pds[0], pds[1]))
            #for sta, sto in pds:
            #    logger.debug('  {0} - {1} s'.format(sta, sto))
        if self.path_periods_save is not None:
            pickle_utilities.dump(self.path_periods_save, {'periods':self.periods_goodstat, 'config':self.config})
        self.summary_results = []
        for ip, pds in enumerate(self.periods_goodstat):
            # Analysis instances
            logger.debug(pds)
            self.analyses.append(pLATLikelihoodConfig.GRBConfig(target=self.grb, phase=self.config['phase'], emin=self.config['energy']['min'], emax=self.config['energy']['max'], deg_roi=self.config['roi']['radius'], zmax=self.config['zenith']['max'], suffix=self.config['suffix'], tmin_special=pds[0], tmax_special=pds[1]))
            self.summary_results.append({})
            self.summary_results[-1]['time'] = {'min':self.analyses[-1].tmin, 'max':self.analyses[-1].tmax}

            #nevt_rough = 
            self.analyses[-1].setup(force={'download':False, 'filter':self.config['force'], 'maketime':True, 'livetime':self.config['force'], 'exposure':self.config['force'], 'model_3FGL_sources':True, 'diffuse_responses':self.config['force']}, skip_zero_data=True)
            #self.analyses[-1].use_external_model(self.analysis_whole.path_model_xml_new)


    def add_calonly(self, path_onevt, path_onexp, path_offevt, path_offexp, rclass='R100'):
        file_onevt = ROOT.TFile(path_onevt, "READ")
        tr_onevt = file_onevt.Get('EVENTS_GRB{0}'.format(self.analysis_whole.target.name))
        file_onexp = ROOT.TFile(path_onexp, "READ")
        file_offevt = ROOT.TFile(path_offevt, "READ")
        file_offexp = ROOT.TFile(path_offexp, "READ")
        for iana, ana in enumerate(self.analyses):
            ana.calonly = pLATLikelihoodConfig.CalOnlyData(energy_bins=ana.like.energies, tr_onevt=tr_onevt, htg_onexp=file_onexp.Get('htgExp_{0}_scaled'.format(iana)), htg_offevt=file_offevt.Get('htgEvt_GalOffCalOnly_{rcla}'.format(rcla=rclass)), htg_offexp=file_offexp.Get('htgExp_GalacticOFF_yx_CalOnly{rcla}'.format(rcla=rclass)), on_time=(ana.target.met0+ana.tmin, ana.target.met0+ana.tmax), on_classes=pLATLikelihoodConfig.DICT_EVCLASSES_BIT_CALONLY[rclass], on_zenith=(0., ana.zmax), on_theta=(0., ana.thetamax))


    def count_energy(self, rlim=6.0):
        events = get_event_time_energy_angsep(self.analysis_whole.path_filtered, self.analysis_whole.tmin, self.analysis_whole.tmax, rlim, ra=self.grb.ra, dec=self.grb.dec, torigin=self.analysis_whole.target.met0, zmax=self.analysis_whole.zmax)
        self.counts_time = events[0]
        self.counts_energy = events[1]
        self.counts_angsep = events[2]


    def run_analysis(self):
        for ip, pds in enumerate(self.periods_goodstat):
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
            #if self.analyses[ip].dct_summary_results['TS']>=4:
            #    self.analyses[ip].eval_flux_and_error()
            if self.analyses[ip].dct_summary_results['TS']>0: #4:
                self.analyses[ip].eval_limits_powerlaw(str_index_fixed=['best', 'free'], emin=self.config['energy']['min'], emax=self.config['energy']['max'], eref=self.config['energy']['ref'])
                self.analyses[ip].reset_target_norm()
                self.analyses[ip].eval_limits_powerlaw_index(emin=self.config['energy']['min'], emax=self.config['energy']['max'], eref=self.config['energy']['ref'])
            else:
                self.analyses[ip].eval_limits_powerlaw(str_index_fixed=['best'], emin=self.config['energy']['min'], emax=self.config['energy']['max'], eref=self.config['energy']['ref'])
            self.analyses[ip].scan_norm_and_index(eref=self.config['energy']['ref'], use_calonly=False)

            dict_allowed_intervals = {}
            dict_allowed_intervals['1sigma'] = pLATLikelihoodConfig.get_allowed_intervals(ana.dct_summary_results['dloglike'], ana.dct_summary_results['dloglike']['dloglike']<=2.30)
            dict_allowed_intervals['2sigma'] = pLATLikelihoodConfig.get_allowed_intervals(ana.dct_summary_results['dloglike'], ana.dct_summary_results['dloglike']['dloglike']<=6.18)
            self.analyses[ip].dct_summary_results['allowed_intervals'] = dict_allowed_intervals

            self.summary_results[ip].update(self.analyses[ip].dct_summary_results)


    def scan_parameters(self, lcintegrals=None, lcindices=None, tnorm=10., rescaler=1.):
        lst_periods = []
        for ip, pds in enumerate(self.periods_goodstat):
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
        self.dct_stored = {'config':self.config, 'results':self.summary_results, 'counts':{'time':self.counts_time, 'energy':self.counts_energy, 'angsep':self.counts_angsep}} if stuff is None else stuff
        #path_pickle = '{base}/{target}/{energy}/{roi}/{phase}/LightCurve_{target}_{spectype}_{index}{suffix}.pickle'.format(base=self.analysis_whole.dir_base, target=self.analysis_whole.target.name, energy=self.analysis_whole.str_energy, roi=self.analysis_whole.str_roi, phase='lightcurve', spectype=self.analysis_whole.target.spectraltype, index=self.analysis_whole.str_index, suffix=self.analysis_whole.suffix)
        path_pickle = '{0}/{1}.pickle'.format(self.outdir, self.outbasename)
        #logger.info("""Object contents: 
#{0}""".format(self.dct_stored))
        with open(path_pickle, mode='wb') as f:
            pickle.dump(self.dct_stored, f)
        logger.info('Result summary has been serialized as {0}'.format(path_pickle))

        
def lightcurve_prefactor(tmin, tmax, integral, index, tb_gti, torigin, tnorm=10., rescaler=1.): #, torigin=0): #, tend=100000
    logger.debug('Integral:{0}'.format(integral))
    logger.debug('LC index:{0}'.format(index))

    powerlaw = lambda t, alpha: pow(t/tnorm, alpha)
    itgl = 0.
    weight = 0.
    
    ti0, ti1 = tb_gti['START']-torigin, tb_gti['STOP']-torigin
    tdiff = tb_gti['STOP']-tb_gti['START']
    if index!=-1:
        tiref = (index+1) / (index+2) * (pow(ti1, index+2)-pow(ti0, index+2)) / (pow(ti1, index+1)-pow(ti0, index+1))
    else:
        tiref = (ti1-ti0) / (np.log(ti1/ti0))
    amp = powerlaw(tiref, index) * tdiff
    itgl = np.sum(tiref * amp)
    weight = np.sum(amp)
    tref = itgl / weight
    return (integral * pow(tref/tnorm, index) * rescaler), tref


def make_lightcurves(name, wholephase, emin, emax, eref, roi, ngoodstat, rgoodstat, ntbinperdecade, suffix, grbcatalogue, refit, force, outdir, index, redshift=False, tmin=0, tmax=10000, addphoton=None):
    if redshift==0 or redshift!=redshift:
        sptype = 'PowerLaw2'
        sppars = {'Integral':1E-7, 'Index':-2., 'LowerLimit':emin, 'UpperLimit':emax}
    else:
        sptype = 'EblAtten::PowerLaw2'
        sppars = {'Integral':1E-7, 'Index':-2., 'LowerLimit':emin, 'UpperLimit':emax, 'tau_norm':1., 'redshift':redshift, 'ebl_model':4}
    lc = LightCurveGRB(name=name, wholephase=wholephase, tmin=tmin, tmax=tmax, emin=emin, emax=emax, eref=eref, deg_roi=roi, ngoodstat=ngoodstat, rlim_goodstat=rgoodstat, ntbinperdecade=ntbinperdecade, suffix=suffix, grbcatalogue=grbcatalogue, refit=refit, force=force, outdir=None, spectraltype=sptype, spectralpars=sppars)
    lc.setup()
    lc.count_energy()
    lc.run_analysis()
    lc.pickle()
    plot_lightcurves(lc.dct_stored, index=index, addphoton=addphoton)
    dataframe(lc.dct_stored, index=index)


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


def plot_lightcurves(dct_summary, outdir=None, ts_threshold=4.0, index='free', grbcatalogue=pLATLikelihoodConfig.GRB_CATALOGUE_LAT, addphoton=None, fitlc=False):
    if isinstance(dct_summary, basestring) and dct_summary[-7:]=='.pickle':
        dct_summary = pickle_utilities.load(dct_summary)
    elif not isinstance(dct_summary, dict):
        logger.critical('The input {0} was NOT a dictionary or path of pickle file!!!'.format(dct_summary[-7:]))
        sys.exit(1)
    outbasename = 'LightCurve_{target}_index{idx}{suffix}'.format(target=str(dct_summary['config']['name']), idx=index, suffix=str(dct_summary['config']['suffix']))

    # Config
    str_energies = '{emin:3.3f} - {emax:3.0f}'.format(emin=dct_summary['config']['energy']['min'], emax=dct_summary['config']['energy']['max'])
    #dct_summary['fit'] = {'flux':{}}

    # Starting point of fitting
    tb_lat = ReadLATCatalogueInfo.open_table(grbcatalogue)
    tb_gbm = ReadGBMCatalogueInfo.open_table()
    tb_one = ReadLATCatalogueInfo.select_one_by_name(tb_lat, dct_summary['config']['name'], tb_gbm)
    #t95 = tb_one['GBM']['T90'] + tb_one['GBM']['T90_START']

    # Characteristices
    dct_curves = {}
    dct_curves['TS'] = pMatplot.Curve('TS', xlabel='Time - T0 [s]', ylabel=r'$\sqrt{\rm{max}(TS, 0)}$', xerr_asym=True, yerr_asym=False)

    dct_curves['flux'] = pMatplot.Curve('flux', xlabel='Time - T0 [s]', ylabel=r'Photon flux $\mathrm{[/cm^2 s]}$', xerr_asym=True, yerr_asym=True)
    dct_curves['flux_ul'] = pMatplot.Curve('flux', xlabel='Time - T0 [s]', ylabel=r'Photon flux $\mathrm{[/cm^2 s]}$', xerr_asym=True, yerr_asym=False, ul=True)

    dct_curves['eflux'] = pMatplot.Curve('eflux', xlabel='Time - T0 [s]', ylabel=r'Energy flux $\mathrm{[erg/cm^2 s]}$', xerr_asym=True, yerr_asym=True)
    dct_curves['eflux_ul'] = pMatplot.Curve('eflux', xlabel='Time - T0 [s]', ylabel=r'Energy flux $\mathrm{[erg/cm^2 s]}$', xerr_asym=True, yerr_asym=False, ul=True)

    dct_curves['e2dnde'] = pMatplot.Curve('e2dnde', xerr_asym=True, yerr_asym=True, xlabel='Time - T0 [s]', ylabel=r'$E^2 dN/dE \, \rm{{at}} \, {ene:3.1f} \rm{{GeV}} \, \mathrm{{[erg/cm^2 s]}}$'.format(ene=dct_summary['config']['energy']['ref']/1000.))
    dct_curves['e2dnde_ul'] = pMatplot.Curve('e2dnde', xerr_asym=True, yerr_asym=False, ul=True, xlabel='Time - T0 [s]', ylabel=r'$E^2 dN/dE \, \rm{{at}} \, {ene:3.1f} \rm{{GeV}} \, \mathrm{{[erg/cm^2 s]}}$'.format(ene=dct_summary['config']['energy']['ref']/1000.)) #ene=dct_summary['results'][0]['Scale']['value']/1000.), )

    dct_curves['Index'] = pMatplot.Curve('Index', xlabel='Time - T0 [s]', ylabel='Spectral index', xerr_asym=True, yerr_asym=False)
    #dct_curves['Index_ul'] = pMatplot.Curve('Index', xlabel='Time - T0 [s]', ylabel='Spectral index', xerr_asym=True, yerr_asym=False, ul=True)
    #dct_curves['Index_ll'] = pMatplot.Curve('Index', xlabel='Time - T0 [s]', ylabel='Spectral index', xerr_asym=True, yerr_asym=False, ll=True)

#    dct_curves['Energy'] = pMatplot.Curve('Energy', xlabel='Time - T0 [s]', ylabel=r'$\log Energy \, \rm{{[GeV]}}$')

    fig, ax = plt.subplots(1+sum([int(v.ul==False and v.ll==False) for v in dct_curves.values()]), 1, figsize=(10, 16), sharex=True)
    
    for ic, curve in enumerate(dct_curves.values()):
        logger.info('===== {0} ====='.format(curve.quantity))
        for period in dct_summary['results']:
            t0 = max(1, period['time']['min'])
            t1 = period['time']['max']
            tref = 10**((np.log10(t0)+np.log10(t1))/2.0) #(t0+t1)/2.0 #sqrt(t0*t1)
            logger.info('----- {tmin:.1f} - {tmax:.1f} s -----'.format(tmin=t0, tmax=t1))            
            #logger.debug(tref)
            if period['TS']!=period['TS']:
                logger.warning('Ananlysis result is NaN! Skipping...')
            if curve.quantity == 'TS':
                y = sqrt(max(0, period[curve.quantity]))
                yerr = 0
                logger.info('{v:.2}'.format(v=y))
            if curve.ul==False and period['TS']>=ts_threshold:
                if curve.quantity in ('Index'):
                    y = period[curve.quantity]['value']
                    yerr = period[curve.quantity]['error']
                    logger.info('{v:.2} +/- {e:.2}'.format(v=y, e=yerr))
                elif curve.quantity in ('flux', 'eflux', 'e2dnde'):
                    logger.debug('Index: '+index)
                    logger.debug(period['limits'])
                    logger.debug(period['limits'][index])
                    logger.debug(period['limits'][index][curve.quantity])
                    y = period['limits'][index][curve.quantity]['x0']
                    yerr_hi = period['limits'][index][curve.quantity]['err_hi']
                    yerr_lo = period['limits'][index][curve.quantity]['err_lo']
                    if curve.quantity in ('eflux', 'e2dnde'):
                        y = y*MEVtoERG
                        yerr_hi = yerr_hi*MEVtoERG
                        yerr_lo = yerr_lo*MEVtoERG
                    logger.info('{v:.2} + {eh:.2} - {el:.2}'.format(v=y, eh=yerr_hi, el=yerr_lo))

                if curve.yerr_asym:
                    curve.set_point(tref, y, {'hi':t1-tref, 'lo':tref-t0}, {'lo':yerr_lo, 'hi':yerr_hi})
                else:
                    curve.set_point(tref, y, {'hi':t1-tref, 'lo':tref-t0}, yerr)

            elif (curve.ul==True or curve.quantity in ('TS')) and period['TS']<ts_threshold: #and period['Nobs']>0:
                if curve.quantity in ('flux', 'eflux', 'e2dnde'):
                    y = period['limits']['best'][curve.quantity]['ul'] #y = period['limits'][index][curve.quantity]['ul']
                    yerr = 0
                    if curve.quantity in ('eflux', 'e2dnde'):
                        y = y*MEVtoERG
                if curve.quantity in ('Index'):
                    y = period['index_limit']['free']['index']['ul']
                    yerr = 0
                logger.info('UL: {ul:.2}'.format(ul=y))
                curve.set_point(tref, y, {'hi':t1-tref, 'lo':tref-t0}, yerr)
            elif (curve.ll==True) and period['TS']<ts_threshold:
                if curve.quantity in ('Index'):
                    y = period['index_limit']['free']['index']['ll']
                    yerr = 0
                logger.info('LL: {ll:.2}'.format(ll=y))
                curve.set_point(tref, y, {'hi':t1-tref, 'lo':tref-t0}, yerr)
    
    ax[0].set_title('GRB '+str(dct_summary['config']['name']))
    ax[0].errorbar(dct_curves['TS'].get_xdata(), dct_curves['TS'].get_ydata(), xerr=dct_curves['TS'].get_xerr(), yerr=dct_curves['TS'].get_yerr(), fmt=dct_curves['TS'].fmt)
    ax[0].set_ylabel(dct_curves['TS'].ylabel)
    ax[0].set_xscale("log", nonposx='clip')
    ax[0].grid(ls='-', lw=0.5, alpha=0.5)
    ax[0].set_xlim((1.0, 100000.0))

    ax[1].errorbar(dct_curves['flux'].get_xdata(), dct_curves['flux'].get_ydata(), xerr=dct_curves['flux'].get_xerr(), yerr=dct_curves['flux'].get_yerr(), fmt=dct_curves['flux'].fmt, ms=2)
    ax[1].errorbar(dct_curves['flux_ul'].get_xdata(), dct_curves['flux_ul'].get_ydata(), xerr=dct_curves['flux_ul'].get_xerr(), fmt=dct_curves['flux_ul'].fmt)
    ax[1].set_ylabel(dct_curves['flux'].ylabel)
    ax[1].set_xscale("log", nonposx='clip')
    ax[1].set_yscale("log", nonposx='clip')
    ax[1].grid(ls='-', lw=0.5, alpha=0.5)
    logger.debug('Axis label:')
    logger.debug(ax[1].get_yticks())
    ax[1].set_yticks([y for y in ax[1].get_yticks() if y<0.5*ax[1].get_ylim()[1]])
    # Fit
    if fitlc==True:
        dct_fit = {'fit':{'flux':{'lightcurve':{}}}}
        flux_max = dct_curves['flux'].get_maximum()
        if flux_max!=0:
            t_flux_max = flux_max[1]
            fit_result = dct_curves['flux'].fit_lin(t_flux_max, t_flux_max)
            if isinstance(fit_result, tuple) and len(fit_result)==2:
                params, butterfly = fit_result[0], fit_result[1]
                for iparam, param, in enumerate(params[0]):
                    logger.info("""Parameter {0}: {1} +/- {2}""".format(iparam, params[0][iparam], params[1][iparam]))
                str_powerlaw_fitted = r"""$ \displaystyle F(t) = \rm{{ F_{{0}} }} \left( \frac{{ \it{{ t }} }}{{ {t_scale:.1E} }} \right)^{{- \rm{{ \alpha }} }}$
$\rm{{ F_{{0}} = {f0:.2E} \pm {f0err:.1E} }}$
$\rm{{ \alpha = {alpha:.2E} \pm {alphaerr:.2E} }}$""".format(f0=params[0][0], f0err=params[1][0], alpha=params[0][1], alphaerr=params[1][1], t_scale=t_flux_max)
                ax[1].plot(butterfly[0], butterfly[1], c='g')
       #ax[1].fill_between(butterfly[0],butterfly[1]+butterfly[2],butterfly[1]-butterfly[2] ,alpha=0.2, facecolor='g')
                ax[1].text(butterfly[0][0]*10, butterfly[1][0]/5., str_powerlaw_fitted)
                dct_fit['fit']['flux']['lightcurve']['amplitude'] = {'value':params[0][0], 'error':params[1][0]}
                dct_fit['fit']['flux']['lightcurve']['index'] = {'value':params[0][1], 'error':params[1][1]}
                dct_fit['fit']['flux']['lightcurve']['t_scale'] = t_flux_max #t95
                dct_fit['fit']['flux']['lightcurve']['tmin'] = t_flux_max #t95
                dct_fit['fit']['flux']['lightcurve']['tmax'] = None
                pickle_utilities.dump('{0}/{1}_fit.pickle'.format(outdir, outbasename), dct_fit)

    ax[2].errorbar(dct_curves['eflux'].get_xdata(), dct_curves['eflux'].get_ydata(), xerr=dct_curves['eflux'].get_xerr(), yerr=dct_curves['eflux'].get_yerr(), fmt=dct_curves['eflux'].fmt)
    ax[2].errorbar(dct_curves['eflux_ul'].get_xdata(), dct_curves['eflux_ul'].get_ydata(), xerr=dct_curves['eflux_ul'].get_xerr(), fmt=dct_curves['eflux_ul'].fmt)
    ax[2].set_ylabel(dct_curves['eflux'].ylabel)
    ax[2].set_xscale("log", nonposx='clip')
    ax[2].set_yscale("log", nonposx='clip')
    ax[2].grid(ls='-', lw=0.5, alpha=0.5)
    logger.debug('Axis label:')
    logger.debug(ax[2].get_yticks())
    ax[2].set_yticks([y for y in ax[2].get_yticks() if y<0.5*ax[2].get_ylim()[1]])

    ax[3].errorbar(dct_curves['e2dnde'].get_xdata(), dct_curves['e2dnde'].get_ydata(), xerr=dct_curves['e2dnde'].get_xerr(), yerr=dct_curves['e2dnde'].get_yerr(), fmt=dct_curves['e2dnde'].fmt)
    ax[3].errorbar(dct_curves['e2dnde_ul'].get_xdata(), dct_curves['e2dnde_ul'].get_ydata(), xerr=dct_curves['e2dnde_ul'].get_xerr(), fmt=dct_curves['e2dnde_ul'].fmt)
    ax[3].set_ylabel(dct_curves['e2dnde'].ylabel)
    ax[3].set_xscale("log", nonposx='clip')
    ax[3].set_yscale("log", nonposx='clip')
    ax[3].grid(ls='-', lw=0.5, alpha=0.5)
    logger.debug('Axis label:')
    logger.debug(ax[3].get_yticks())
    ax[3].set_yticks([y for y in ax[3].get_yticks() if y<0.5*ax[3].get_ylim()[1]])

    ax[4].errorbar(dct_curves['Index'].get_xdata(), dct_curves['Index'].get_ydata(), xerr=dct_curves['Index'].get_xerr(), yerr=dct_curves['Index'].get_yerr(), fmt=dct_curves['Index'].fmt)
    #ax[4].errorbar(dct_curves['Index_ul'].get_xdata(), dct_curves['Index_ul'].get_ydata(), xerr=dct_curves['Index_ul'].get_xerr(), yerr=dct_curves['Index_ul'].get_yerr(), fmt=dct_curves['Index_ul'].fmt)
    #ax[4].errorbar(dct_curves['Index_ll'].get_xdata(), dct_curves['Index_ll'].get_ydata(), xerr=dct_curves['Index_ll'].get_xerr(), yerr=dct_curves['Index_ll'].get_yerr(), fmt=dct_curves['Index_ll'].fmt)
    ax[4].set_ylabel(dct_curves['Index'].ylabel)
    ax[4].set_xlabel(dct_curves['Index'].xlabel)
    ax[4].set_xscale("log", nonposx='clip')
    ax[4].grid(ls='-', lw=0.5, alpha=0.5)
    logger.debug('Axis label:')
    logger.debug(ax[4].get_yticks())
    ax[4].set_yticks([y for y in ax[4].get_yticks() if y<ax[4].get_ylim()[1]-0.1])

    ax[5].scatter(dct_summary['counts']['time'], dct_summary['counts']['energy'], alpha=0.5, c=np.clip(1./dct_summary['counts']['angsep'], 0.0, 1.0), cmap=cm.Purples)
    if addphoton is not None:
        logger.info('Additional photons:')
        t_added = np.array(addphoton['time'])
        e_added = np.array(addphoton['energy'])
        for t,e in zip(t_added, e_added):
            logger.info('{e:.1f} MeV at {t:.1f}'.format(t=t, e=e))
        ax[5].scatter(t_added, e_added, marker='D', c='r', edgecolors='face')#, s=np.log10(e_added))
    ax[5].set_xlabel('Time - T0 [s]')
    ax[5].set_ylabel('Energy [MeV]')
    #ax[5].set_xscale("log", nonposx='clip')
    if len(dct_summary['counts']['time'])>0 or addphoton is not None:
        ax[5].set_yscale("log", nonposx='clip')
    ax[5].grid(ls='-', lw=0.5, alpha=0.5)
    ax[5].set_yticks([y for y in ax[5].get_yticks() if y<0.5*ax[5].get_ylim()[1]])
    ax[5].set_ylim(bottom=max(300,dct_summary['config']['energy']['min']), top=1.1*dct_summary['config']['energy']['max'])

    fig.tight_layout() #subplots_adjust(hspace=0)
    fig.subplots_adjust(hspace=0)
    outdir = outdir if outdir is not None else '{base}/{target}/E{emin:0>7.0f}-{emax:0>7.0f}MeV/r{roi:0>2.0f}deg/{phase}'.format(base=pLATLikelihoodConfig.PATH_BASEDIR, target=str(dct_summary['config']['name']), emin=dct_summary['config']['energy']['min'], emax=dct_summary['config']['energy']['max'], roi=dct_summary['config']['roi']['radius'], phase='lightcurve')
    for ff in ['png', 'pdf']:
        fig.savefig('{0}/{1}.{2}'.format(outdir, outbasename, ff))


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
@click.option('--suffix', '-s', type=str, default='')
@click.option('--force', '-f', is_flag=True)
@click.option('--refit', '-r', is_flag=True)
@click.option('--outdir', '-o', type=str, default='')
@click.option('--plotonly', '-p', type=str, default=None, help='Path of result pickle file if you skip analyses.')
@click.option('--bsub', '-b', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(namemin, namemax, mode, emin, emax, eref, tmin, tmax, roi, ngoodstat, rgoodstat, ntbinperdecade, refit, force, suffix, grbcatalogue, outdir, plotonly, index, redshift, bsub, loglevel):
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
            acmd = ['bsub', '-o','{0}/GRB{0}_lightcurve{1}.log'.format(name, suffix if suffix=='' else '_'+suffix), '-J','lc{0}'.format(name[:-3]), '-W','400', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/STLikelihoodAnalysis/LightCurve.py', '--mode', mode, '--emin', str(emin), '--emax', str(emax), '--eref', str(eref), '--tmin', str(tmin), '--tmax', str(tmax), '-s', suffix, '--index', 'free', '--redshift', str(redshift), '--roi', str(roi), '--ngoodstat', str(ngoodstat), '--rgoodstat', str(rgoodstat), '--ntbinperdecade', str(ntbinperdecade), '--grbcatalogue', grbcatalogue, '--namemin', name]
            if force==True:
                acmd.append('--force')
            if refit==True:
                acmd.append('--refit')
            print acmd
            subprocess.call(acmd)
        return 0
    else:
        tb_lat = ReadLATCatalogueInfo.select_one_by_name(tb_lat, namemin, tb_gbm)

        # CalOnly photons
        addphoton = None
        if namemin=='090926181':
            addphoton={'time': (422.7,), 'energy':(50.49609375*1000.,)}
        elif namemin=='150902733':
            addphoton={'time': (2064.52053469,), 'energy':(83.6322578125*1000.,)}
        elif namemin=='160509374':
            addphoton={'time': (2035.85387415,5757.82151717), 'energy':(115.829539063*1000.,63.1624726562*1000.)}

        if plotonly==None:
            make_lightcurves(name=namemin, wholephase=mode, emin=emin, emax=emax, eref=eref, tmin=tmin, tmax=tmax, roi=roi, ngoodstat=ngoodstat, rgoodstat=rgoodstat, ntbinperdecade=ntbinperdecade, suffix=suffix, grbcatalogue=grbcatalogue, refit=refit, force=force, outdir=outdir, index=index, redshift=redshift, addphoton=addphoton)#, tmin=tb_lat['LAT_TRIGGER_TIME']-tb_lat['TRIGGER_TIME'])
        else:
            plot_lightcurves(plotonly, index=index, addphoton=addphoton) #, addphoton={'time': (422.7,), 'energy':(50.49609375*1000.,)})
            dataframe(plotonly, index=index)


if __name__ == '__main__':
    main()
