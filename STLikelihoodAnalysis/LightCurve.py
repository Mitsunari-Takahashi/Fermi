#!/usr/bin/env python
"""Module for making light curves of LAT data.
The main class LightCurve is a chain of another module pLATLikelihoodConfig.py.
The authour: Mitsunari Takahashi
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
import datetime
import numpy as np
import math
from math import log10, log, sqrt, ceil, isnan, pi, factorial
from sympy import *
from scipy import integrate
import click
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import STLikelihoodAnalysis.pLATLikelihoodConfig as pLATLikelihoodConfig
from FindGoodstatPeriods import find_goodstat_periods
from FindCrossEarthlimb import find_cross_earthlimb
from STLikelihoodAnalysis import get_module_logger
import pickle_utilities
import pMatplot


##### VERSION OF THIS MACRO #####
VERSION = 1.3 # 2017.09.21


##### Logger #####
logger = get_module_logger(__name__)


##### Conversion from MeV to erg ######
MEVtoERG = 1.6021766208E-6

class LightCurve:
    def __init__(self, name, met0, tmin, tmax, emin, emax, evclass=128, evtype=3, ft2interval='30s', deg_roi=12., rad_margin=10., zmax=100., index_fixed=None, suffix='', grbcatalogue=pLATLikelihoodConfig.GRB_CATALOGUE_LTF, psForce=False, refit=True, force=False, outdir=None):
        self.config = {'name':name, 
                       'time':{'min':tmin, 'max':tmax}, 
                       'met':{'o':met0, 'min':met0+tmin, 'max':met0+tmax}, 
                       'energy':{'min':emin, 'max':emax},
                       'evclass':evclass,
                       'evtype':evtype, 
                       'ft2interval':ft2interval, 
                       'roi':{'radius':deg_roi, 'margin':rad_margin}, 
                       'zenith':{'max':zmax}, 
                       'index_fixed':index_fixed, 
                       'suffix':suffix, 
                       'model':{'psForce':psForce},
                       'refit':refit,
                       'force':force}
        self.outdir = outdir if outdir is not None else os.getcwd()


class LightCurveGRB(LightCurve):
    def __init__(self, name, tmin=0.0, tmax=10000.0, emin=100.0, emax=100000.0, evclass=128, evtype=3, ft2interval='30s', deg_roi=12., rad_margin=10., zmax=100., index_fixed=None, suffix='', grbcatalogue=pLATLikelihoodConfig.GRB_CATALOGUE_LTF, psForce=False, rlim_goodstat=1.0, ngoodstat=10., refit=True, force=False, outdir=None):
        # Target GRB
        self.grb = pLATLikelihoodConfig.GRBTarget(name, grbcatalogue)

        # Light curve instance
        LightCurve.__init__(self, name=name, met0=self.grb.met0, tmin=tmin, tmax=tmax, emin=emin, emax=emax, evclass=evclass, evtype=evtype, ft2interval=ft2interval, deg_roi=deg_roi, rad_margin=rad_margin, zmax=zmax, index_fixed=index_fixed, suffix=suffix, psForce=psForce, refit=True, force=False)
        self.config['goodstat'] = {'n':ngoodstat, 'r':rlim_goodstat}
        self.periods_goodstat = []
        self.analyses = []


    def setup(self):
        # Analysis instance which covers the whole time range
        self.analysis_whole = pLATLikelihoodConfig.GRBConfig(target=self.grb, phase='unified', tstop=self.config['time']['max'], emin=self.config['energy']['min'], emax=self.config['energy']['max'], deg_roi=self.config['roi']['radius'], zmax=self.config['zenith']['max'], suffix=self.config['suffix'])
        self.analysis_whole.download()
        self.analysis_whole.set_directories()
        self.analysis_whole.filter(False)
        self.analysis_whole.maketime(False)

        self.outdir = '{base}/{target}/{energy}/{roi}/{phase}'.format(base=self.analysis_whole.dir_base, target=self.analysis_whole.target.name, energy=self.analysis_whole.str_energy, roi=self.analysis_whole.str_roi, phase='lightcurve')
        self.outbasename = 'LightCurve_{target}_{spectype}_{index}{suffix}'.format(target=self.analysis_whole.target.name, spectype=self.analysis_whole.target.spectraltype, index=self.analysis_whole.str_index, suffix=self.analysis_whole.suffix)

        # Find good-statistics interval
        validtimes = find_cross_earthlimb(self.analysis_whole.path_ft2, self.analysis_whole.target.ra, self.analysis_whole.target.dec, self.analysis_whole.tmin, self.analysis_whole.tmax, self.analysis_whole.zmax, self.analysis_whole.target.met0)
        logger.info("""Good time intervals (zenith < {zen} deg)
{vt}""".format(zen=self.config['zenith']['max'], vt=validtimes))

        for vt in validtimes:
            logger.info(vt)
            logger.debug('Event file:', self.analysis_whole.path_filtered_gti)
            self.periods_goodstat += find_goodstat_periods(self.analysis_whole.path_filtered_gti, vt[0], vt[1], nthreshold=self.config['goodstat']['n'], rlim=self.config['goodstat']['r'], ra=self.grb.ra, dec=self.grb.dec, torigin=self.analysis_whole.target.met0)

        logger.info("""Good-statistics peroids (>{nth} events within {rlim} deg)
{vt}""".format(nth=self.config['goodstat']['n'], rlim=self.config['goodstat']['r'], vt=self.periods_goodstat))


    def run_analysis(self):
        self.summary_results = []
        for ip, pds in enumerate(self.periods_goodstat):
            # Analysis instances
            self.analyses.append(pLATLikelihoodConfig.GRBConfig(target=self.grb, phase='lightcurve', emin=self.config['energy']['min'], emax=self.config['energy']['max'], deg_roi=self.config['roi']['radius'], zmax=self.config['zenith']['max'], suffix=self.config['suffix'], tmin_special=pds[0], tmax_special=pds[1]))
            self.summary_results.append({})
            self.summary_results[-1]['time'] = {'min':self.analyses[-1].tmin, 'max':self.analyses[-1].tmax}

            self.analyses[-1].setup(force={'download':self.config['force'], 'filter':self.config['force'], 'maketime':self.config['force'], 'livetime':self.config['force'], 'exposure':self.config['force'], 'model_3FGL_sources':True, 'diffuse_responses':self.config['force']})

            self.analyses[-1].fit(bredo=self.config['refit'])
            self.analyses[-1].summarize_fit_results()
            self.analyses[-1].plot_countspectra_fitted()
            self.analyses[-1].eval_flux_and_error()
            self.analyses[-1].eval_limits_powerlaw(str_index_fixed=['free'])
            self.summary_results[-1].update(self.analyses[-1].dct_summary_results)


    def pickle(self):
        self.dct_stored = {'config':self.config, 'results':self.summary_results}
        #path_pickle = '{base}/{target}/{energy}/{roi}/{phase}/LightCurve_{target}_{spectype}_{index}{suffix}.pickle'.format(base=self.analysis_whole.dir_base, target=self.analysis_whole.target.name, energy=self.analysis_whole.str_energy, roi=self.analysis_whole.str_roi, phase='lightcurve', spectype=self.analysis_whole.target.spectraltype, index=self.analysis_whole.str_index, suffix=self.analysis_whole.suffix)
        path_pickle = '{0}/{1}.pickle'.format(self.outdir, self.outbasename)
        logger.info("""Object contents: 
{0}""".format(self.dct_stored))
        with open(path_pickle, mode='wb') as f:
            pickle.dump(self.dct_stored, f)
        logger.info('Result summary has been serialized as {0}'.format(path_pickle))



def make_lightcurves(name, emin, emax, roi, ngoodstat, suffix, grbcatalogue, refit, force, outdir, index):
    lc = LightCurveGRB(name=name, emin=emin, emax=emax, deg_roi=roi, ngoodstat=ngoodstat, suffix=suffix, grbcatalogue=grbcatalogue, refit=refit, force=force, outdir=None)
    lc.setup()
    lc.run_analysis()
    lc.pickle()
    plot_lightcurves(lc.dct_stored, index=index)
    dataframe(lc.dct_stored, index=index)


def dataframe(dct_summary, outdir=None, ts_threshold=4.0, index='free'):
    if isinstance(dct_summary, basestring) and dct_summary[-7:]=='.pickle':
        dct_summary = pickle_utilities.load(dct_summary)
    elif not isinstance(dct_summary, dict):
        logger.critical('The input {0} was NOT a dictionary or path of pickle file!!!'.format(dct_summary[-7:]))
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
                 'e2dnde[erg cm^-2 s^-1] @1GeV': [],
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
            dct_frame['e2dnde[erg cm^-2 s^-1] @1GeV'].append(period['limits'][index]['e2dnde']['x0']*MEVtoERG)
            dct_frame['e2dnde_err_lo'].append(period['limits'][index]['e2dnde']['err_lo']*MEVtoERG)
            dct_frame['e2dnde_err_hi'].append(period['limits'][index]['e2dnde']['err_hi']*MEVtoERG)
            dct_frame['Index'].append(period['Index']['value'])
            dct_frame['Index_err'].append(period['Index']['error'])

    df = pd.DataFrame(dct_frame)
    outdir = outdir if outdir is not None else '{base}/{target}/E{emin:0>7.0f}-{emax:0>7.0f}MeV/r{roi:0>2.0f}deg/{phase}'.format(base=pLATLikelihoodConfig.PATH_BASEDIR, target=str(dct_summary['config']['name']), emin=dct_summary['config']['energy']['min'], emax=dct_summary['config']['energy']['max'], roi=dct_summary['config']['roi']['radius'], phase='lightcurve')
    outbasename = 'LightCurve_{target}_index{idx}{suffix}'.format(target=str(dct_summary['config']['name']), idx=index, suffix=str(dct_summary['config']['suffix']))
    df.to_csv('{dire}/{name}.csv'.format(dire=outdir, name=outbasename))


def plot_lightcurves(dct_summary, outdir=None, ts_threshold=4.0, index='free'):
    if isinstance(dct_summary, basestring) and dct_summary[-7:]=='.pickle':
        dct_summary = pickle_utilities.load(dct_summary)
    elif not isinstance(dct_summary, dict):
        logger.critical('The input {0} was NOT a dictionary or path of pickle file!!!'.format(dct_summary[-7:]))
        sys.exit(1)

    # Config
    str_energies = '{emin:3.3f} - {emax:3.0f}'.format(emin=dct_summary['config']['energy']['min'], emax=dct_summary['config']['energy']['max'])

    # Characteristices
    dct_curves = {}
    dct_curves['TS'] = pMatplot.Curve('TS', xlabel='Time - T0 [s]', ylabel=r'$\sqrt{\rm{max}(TS, 0)}$', xerr_asym=True, yerr_asym=False)

    dct_curves['flux'] = pMatplot.Curve('flux', xlabel='Time - T0 [s]', ylabel=r'Photon flux $\mathrm{[/cm^2 s]}$', xerr_asym=True, yerr_asym=True)
    dct_curves['flux_ul'] = pMatplot.Curve('flux', xlabel='Time - T0 [s]', ylabel=r'Photon flux $\mathrm{[/cm^2 s]}$', xerr_asym=True, yerr_asym=False, ul=True)

    dct_curves['eflux'] = pMatplot.Curve('eflux', xlabel='Time - T0 [s]', ylabel=r'Energy flux $\mathrm{[erg/cm^2 s]}$', xerr_asym=True, yerr_asym=True)
    dct_curves['eflux_ul'] = pMatplot.Curve('eflux', xlabel='Time - T0 [s]', ylabel=r'Energy flux $\mathrm{[erg/cm^2 s]}$', xerr_asym=True, yerr_asym=False, ul=True)

    dct_curves['e2dnde'] = pMatplot.Curve('e2dnde', xlabel='Time - T0 [s]', ylabel=r'$E^2 dN/dE \, \rm{{at}} \, {ene:3.1f} \rm{{GeV}} \, \mathrm{{[erg/cm^2 s]}}$'.format(ene=dct_summary['results'][0]['Scale']['value']/1000.), xerr_asym=True, yerr_asym=True)
    dct_curves['e2dnde_ul'] = pMatplot.Curve('e2dnde', xlabel='Time - T0 [s]', ylabel=r'$E^2 dN/dE \, \rm{{at}} \, {ene:3.1f} \rm{{GeV}} \, \mathrm{{[erg/cm^2 s]}}$'.format(ene=dct_summary['results'][0]['Scale']['value']/1000.), xerr_asym=True, yerr_asym=False, ul=True)

    dct_curves['Index'] = pMatplot.Curve('Index', xlabel='Time - T0 [s]', ylabel='Spectral index', xerr_asym=True, yerr_asym=False)

    fig, ax = plt.subplots(sum([int(v.ul==False) for v in dct_curves.values()]), 1, figsize=(10, 16), sharex=True)
    
    for ic, curve in enumerate(dct_curves.values()):
        logger.info('===== {0} ====='.format(curve.quantity))
        for period in dct_summary['results']:
            t0 = max(1, period['time']['min'])
            t1 = period['time']['max']
            tref = (t0+t1)/2.0 #sqrt(t0*t1)
            logger.info('----- {tmin:.1f} - {tmax:.1f} s -----'.format(tmin=t0, tmax=t1))            
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
                    curve.set_point(tref, y, {'lo':t1-tref, 'hi':tref-t0}, {'lo':yerr_lo, 'hi':yerr_hi})
                else:
                    curve.set_point(tref, y, {'lo':t1-tref, 'hi':tref-t0}, yerr)
            elif (curve.ul==True or curve.quantity in ('TS')) and period['TS']<ts_threshold:
                if curve.quantity in ('flux', 'eflux', 'e2dnde'):
                    y = period['limits'][index][curve.quantity]['ul']
                    yerr = 0
                    if curve.quantity in ('eflux', 'e2dnde'):
                        y = y*MEVtoERG
                curve.set_point(tref, y, {'lo':t1-tref, 'hi':tref-t0}, yerr)

    
    ax[0].set_title('GRB '+str(dct_summary['config']['name']))
    ax[0].errorbar(dct_curves['TS'].get_xdata(), dct_curves['TS'].get_ydata(), xerr=dct_curves['TS'].get_xerr(), yerr=dct_curves['TS'].get_yerr(), fmt=dct_curves['TS'].fmt)
    ax[0].set_ylabel(dct_curves['TS'].ylabel)
    ax[0].set_xscale("log", nonposx='clip')
    ax[0].grid(ls='-', lw=0.5, alpha=0.5)
    ax[0].set_xlim((1.0, 10000.0))

    ax[1].errorbar(dct_curves['flux'].get_xdata(), dct_curves['flux'].get_ydata(), xerr=dct_curves['flux'].get_xerr(), yerr=dct_curves['flux'].get_yerr(), fmt=dct_curves['flux'].fmt)
    ax[1].errorbar(dct_curves['flux_ul'].get_xdata(), dct_curves['flux_ul'].get_ydata(), xerr=dct_curves['flux_ul'].get_xerr(), fmt=dct_curves['flux_ul'].fmt)
    ax[1].set_ylabel(dct_curves['flux'].ylabel)
    ax[1].set_xscale("log", nonposx='clip')
    ax[1].set_yscale("log", nonposx='clip')
    ax[1].grid(ls='-', lw=0.5, alpha=0.5)
    logger.debug('Axis label:')
    logger.debug(ax[1].get_yticks())
    ax[1].set_yticks([y for y in ax[1].get_yticks() if y<ax[1].get_ylim()[1]])

    ax[2].errorbar(dct_curves['eflux'].get_xdata(), dct_curves['eflux'].get_ydata(), xerr=dct_curves['eflux'].get_xerr(), yerr=dct_curves['eflux'].get_yerr(), fmt=dct_curves['eflux'].fmt)
    ax[2].errorbar(dct_curves['eflux_ul'].get_xdata(), dct_curves['eflux_ul'].get_ydata(), xerr=dct_curves['eflux_ul'].get_xerr(), fmt=dct_curves['eflux_ul'].fmt)
    ax[2].set_ylabel(dct_curves['eflux'].ylabel)
    ax[2].set_xscale("log", nonposx='clip')
    ax[2].set_yscale("log", nonposx='clip')
    ax[2].grid(ls='-', lw=0.5, alpha=0.5)
    logger.debug('Axis label:')
    logger.debug(ax[2].get_yticks())
    ax[2].set_yticks([y for y in ax[2].get_yticks() if y<ax[2].get_ylim()[1]])

    ax[3].errorbar(dct_curves['e2dnde'].get_xdata(), dct_curves['e2dnde'].get_ydata(), xerr=dct_curves['e2dnde'].get_xerr(), yerr=dct_curves['e2dnde'].get_yerr(), fmt=dct_curves['e2dnde'].fmt)
    ax[3].errorbar(dct_curves['e2dnde_ul'].get_xdata(), dct_curves['e2dnde_ul'].get_ydata(), xerr=dct_curves['e2dnde_ul'].get_xerr(), fmt=dct_curves['e2dnde_ul'].fmt)
    ax[3].set_ylabel(dct_curves['e2dnde'].ylabel)
    ax[3].set_xscale("log", nonposx='clip')
    ax[3].set_yscale("log", nonposx='clip')
    ax[3].grid(ls='-', lw=0.5, alpha=0.5)
    logger.debug('Axis label:')
    logger.debug(ax[3].get_yticks())
    ax[3].set_yticks([y for y in ax[3].get_yticks() if y<ax[3].get_ylim()[1]])

    ax[4].errorbar(dct_curves['Index'].get_xdata(), dct_curves['Index'].get_ydata(), xerr=dct_curves['Index'].get_xerr(), yerr=dct_curves['Index'].get_yerr(), fmt=dct_curves['Index'].fmt)
    ax[4].set_ylabel(dct_curves['Index'].ylabel)
    ax[4].set_xlabel(dct_curves['Index'].xlabel)
    ax[4].set_xscale("log", nonposx='clip')
    ax[4].grid(ls='-', lw=0.5, alpha=0.5)
    logger.debug('Axis label:')
    logger.debug(ax[4].get_yticks())
    ax[4].set_yticks([y for y in ax[4].get_yticks() if y<ax[4].get_ylim()[1]])

    fig.tight_layout() #subplots_adjust(hspace=0)
    fig.subplots_adjust(hspace=0)
    outdir = outdir if outdir is not None else '{base}/{target}/E{emin:0>7.0f}-{emax:0>7.0f}MeV/r{roi:0>2.0f}deg/{phase}'.format(base=pLATLikelihoodConfig.PATH_BASEDIR, target=str(dct_summary['config']['name']), emin=dct_summary['config']['energy']['min'], emax=dct_summary['config']['energy']['max'], roi=dct_summary['config']['roi']['radius'], phase='lightcurve')
    outbasename = 'LightCurve_{target}_index{idx}{suffix}'.format(target=str(dct_summary['config']['name']), idx=index, suffix=str(dct_summary['config']['suffix']))
    for ff in ['png', 'pdf']:
        fig.savefig('{0}/{1}.{2}'.format(outdir, outbasename, ff))


@click.command()
@click.argument('name', type=str)
@click.option('--grbcatalogue', '-c', type=str, default=pLATLikelihoodConfig.GRB_CATALOGUE_LTF)
@click.option('--emin', type=float, default=100.)
@click.option('--emax', type=float, default=100000.)
@click.option('--roi', type=float, default=12.)
@click.option('--ngoodstat', type=int, default=10)
@click.option('--index', type=click.Choice(['free', 'best']), default='free')
@click.option('--suffix', '-s', type=str, default='')
@click.option('--force', '-f', is_flag=True)
@click.option('--refit', '-r', is_flag=True)
@click.option('--outdir', '-o', type=str, default='')
@click.option('--plotonly', '-p', type=str, default=None, help='Path of result pickle file if you skip analyses.')
def main(name, emin, emax, roi, ngoodstat, refit, force, suffix, grbcatalogue, outdir, plotonly, index):
    if plotonly==None:
        make_lightcurves(name, emin, emax, roi, ngoodstat, suffix, grbcatalogue, refit, force, outdir, index)
    else:
        plot_lightcurves(plotonly, index=index)
        dataframe(plotonly, index=index)


if __name__ == '__main__':
    main()
