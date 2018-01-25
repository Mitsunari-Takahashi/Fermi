#!/usr/bin/env python
"""Module for making spectral index light curve collection of GRB LAT data.
The main class LightCurve is a chain of another module pLATLikelihoodConfig.py.
The authour: Mitsunari Takahashi
"""
import sys
import os
import os.path
path_upstairs = os.path.join(os.path.dirname(__file__), '../')
sys.path.append(path_upstairs)
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
import pickle
import csv
import numpy as np
import click
import itertools
import gc
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import ReadLATCatalogueInfo
import ReadGBMCatalogueInfo
import ReadSwiftCatalogueInfo
import ReadKonusWindRedshiftCatalogueInfo
import pConvertFluxValues
from pConvertFluxValues import MEVtoERG, KEVtoERG
import pickle_utilities
import pMatplot
from pMatplot import find_range_shown, find_indexrange, TPL_MARKER, TPL_COLOR, TPL_LINE
import PlotVHEEvents
import pMETandMJD


##### Logging #####
logger = getLogger(__name__)
handler = StreamHandler()


##### Matplotlib #####
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
plt.rcParams["font.size"] = 12
NMARKER_STYLE = 10

##### VERSION OF THIS MACRO #####
VERSION = 0.0


##### Catalogues #####
GRB_CATALOGUE_LAT = "/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LATBurstCatalogue.xml"
GRB_CATALOGUE_GBM = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/GBMcatalogue20171005.fits'

TABLE_LAT = ReadLATCatalogueInfo.open_table() #GRB_CATALOGUE_LAT)
TABLE_GBM = ReadGBMCatalogueInfo.open_table() #GRB_CATALOGUE_GBM)
TABLE_KW = ReadKonusWindRedshiftCatalogueInfo.open_table()

TIME_INTERVALS = {'prompt':'T0 to T95', 'T95to03ks':'T95 to T95+3ks', '03ksto100ks':'T95+3ks to 100ks'}
ENERGY_RANGES = {'whole': (100, 100000), 'low': (100, 1000), 'mid': (1000, 10000), 'high':(10000, 100000)}

powerlaws = {}
powerlaws[0.1] = lambda t: 1.60E-03*pow(t/10., -1.345)
powerlaws[1.0] = lambda t: 1.68E-03*pow(t/10., -1.325)
powerlaws[10.] = lambda t: 4.62E-04*pow(t/10., -1.085)

lst_markers = ['s', 'o', 'D', 'x', 'd', 'P', '*', 'h']

class ICCollection:
    def __init__(self, index, erange, brighten=[], interdirs='E0000100-0100000MeV/r12deg/lightcurve', suffix='LAT_N10in2deg', suffix2='', table_catalogue=None, figforms=['pdf', 'png'], ts_threshold=4):
        #self.name = name
        self.index = index
        self.skip=('130427324', '150314205')
        self.erange = erange
        self.interdirs = interdirs
        if table_catalogue is None:
            tb_lat = ReadLATCatalogueInfo.read_all(TABLE_LAT, TABLE_GBM, TABLE_KW)
            tb_lat = ReadLATCatalogueInfo.select_gbm_exist(tb_lat)
            tb_lat = ReadLATCatalogueInfo.select_long(tb_lat)
            #tb_lat = ReadLATCatalogueInfo.select_zenith_angle(tb_lat, 100.)
            table_catalogue = ReadLATCatalogueInfo.select_small_error(tb_lat, 0.3)
        self.table_grb = table_catalogue
        self.ts_threshold = ts_threshold
        self.suffix = suffix if suffix=='' else '_'+suffix
        self.suffix2 = suffix2 if suffix2=='' else '_'+suffix2
        self.figforms = figforms
        self.grbs_brighten = []
        for g in table_catalogue:
            if g['GRBNAME'] in self.skip:
                continue
            if g['GRBNAME'] in brighten:
                self.grbs_brighten.append(g['GRBNAME'])
        self.dct_grbs = {}
        self.dct_nmarker = {}
        self.dct_ncolor = {}


    def load_grbs(self, name, dirs='E0000100-0100000MeV/r12deg/lightcurve'):
        self.dataset = {}
        for grb in self.table_grb:
            if grb['GRBNAME'] in skip:
                continue
            logger.debug('Loading likelihood valus of {0}'.format(grb['GRBNAME']))
            path_likelihood = '/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{name}/{dirs}/LightCurve_{name}_PowerLaw_IndexFree_pm3dec_180x090.pickle'.format(name=grb['GRBNAME'], dirs=dirs)
            self.dataset[grb['GRBNAME']] = pickle_utilities.load(path_likelihood)


    def set_curves(self, grbname, dataset):
        for ic, curve in enumerate(self.dct_grbs[grbname].values()):
            logger.info('===== {0} ====='.format(curve.quantity))
            for period in dataset:
                t0 = max(1, period['time']['min'])
                t1 = period['time']['max']
                tref = 10**((np.log10(t0)+np.log10(t1))/2.0)
                logger.debug('----- {tmin:.1f} - {tmax:.1f} s -----'.format(tmin=t0, tmax=t1))
                if not 'limits' in period:
                    logger.warning('GRB{0} {1}-{2} s: No limits! Skipped.'.format(grbname, period['time']['min'], period['time']['max']))
                    #continue

                if curve.quantity == 'TS':
                    y = sqrt(max(0, period[curve.quantity]))
                    yerr = 0
                    logger.debug('{v:.2}'.format(v=y))
                if curve.ul==False and period['TS']>=self.ts_threshold:
                    if curve.quantity in ('Index'):
                        y = period[curve.quantity]['value']
                        yerr = period[curve.quantity]['error']
                        logger.debug('{v:.2} +/- {e:.2}'.format(v=y, e=yerr))
                    elif curve.quantity in ('flux', 'eflux', 'e2dnde'):
                        y = period['limits'][self.index][curve.quantity]['x0']
                        yerr_hi = period['limits'][self.index][curve.quantity]['err_hi']
                        yerr_lo = period['limits'][self.index][curve.quantity]['err_lo']
                        if curve.quantity in ('eflux', 'e2dnde'):
                            y = y*MEVtoERG
                            yerr_hi = yerr_hi*MEVtoERG
                            yerr_lo = yerr_lo*MEVtoERG
                        logger.debug('{v:.2} + {eh:.2} - {el:.2}'.format(v=y, eh=yerr_hi, el=yerr_lo))

                    if curve.yerr_asym:
                        curve.set_point(tref, y, {'hi':t1-tref, 'lo':tref-t0}, {'lo':yerr_lo, 'hi':yerr_hi})
                    else:
                        curve.set_point(tref, y, {'hi':t1-tref, 'lo':tref-t0}, yerr)

                elif (curve.ul==True or curve.quantity in ('TS')) and period['TS']<self.ts_threshold:
                    if curve.quantity in ('flux', 'eflux', 'e2dnde'):
                        y = period['limits']['best'][curve.quantity]['ul']
                        yerr = 0
                        if curve.quantity in ('eflux', 'e2dnde'):
                            y = y*MEVtoERG
                    logger.debug('UL: {ul:.2}'.format(ul=y))
                    curve.set_point(tref, y, {'lo':t1-tref, 'hi':tref-t0}, yerr)


    def plot(self, outdir='.'):
        fig, ax = plt.subplots(3, 1, figsize=(15, 15), sharex=True, sharey=False)
        ngrb_vhe = 0
        #ngrb_brighten  = 0
        for igrb,grb in enumerate(self.table_grb):
            logger.info('###### {0} #####'.format(grb['GRBNAME']))
            if grb['GRBNAME'] in self.skip:
                logger.warning('{0} is skipped!'.format(grb['GRBNAME']))
                continue

            logger.debug('Loading likelihood valus of {0}'.format(grb['GRBNAME']))
            path_likelihood = '/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{name}/{dirs}/LightCurve_{name}_PowerLaw_IndexFree{suf}.pickle'.format(name=grb['GRBNAME'], dirs=self.interdirs, suf=self.suffix)
            grbdata = pickle_utilities.load(path_likelihood)['results']

            self.dct_grbs[grb['GRBNAME']] = {}
            # Index
            self.dct_grbs[grb['GRBNAME']]['flux'] = pMatplot.Curve('flux', xerr_asym=True, yerr_asym=True, xlabel='Time - T0 [s]', ylabel=r'Photon flux \mathrm{[cm^2 \cdot s]}$')
            self.dct_grbs[grb['GRBNAME']]['Index'] = pMatplot.Curve('Index', xerr_asym=True, yerr_asym=False, xlabel='Time - T0 [s]', ylabel='Photon index')
            self.set_curves(grb['GRBNAME'], grbdata)
            logger.debug('X: {0}'.format(self.dct_grbs[grb['GRBNAME']]['Index'].get_xdata()))
            logger.debug('Y: {0}'.format(self.dct_grbs[grb['GRBNAME']]['Index'].get_ydata()))
            self.dct_nmarker[grb['GRBNAME']] = lst_markers[igrb / NMARKER_STYLE]
            self.dct_ncolor[grb['GRBNAME']] = TPL_COLOR[1:][igrb % len(TPL_COLOR[1:])]
            logger.info('Marker style: {0}'.format(self.dct_nmarker[grb['GRBNAME']]))
            logger.info('Marker color: {0}'.format(self.dct_ncolor[grb['GRBNAME']]))
            if grb['GRBNAME'] in self.grbs_brighten:
                alpha = 1
                lwidth = 1
                grblabel = grb['GCNNAME'] 
                #ngrb_brighten+=1
            else:
                alpha = 0.25
                lwidth = 0.25
                grblabel = None
            ax[0].errorbar(self.dct_grbs[grb['GRBNAME']]['flux'].get_xdata(), self.dct_grbs[grb['GRBNAME']]['flux'].get_ydata(), xerr=self.dct_grbs[grb['GRBNAME']]['flux'].get_xerr(), yerr=self.dct_grbs[grb['GRBNAME']]['flux'].get_yerr(), lw=lwidth, ms=5, c=self.dct_ncolor[grb['GRBNAME']], fmt=self.dct_nmarker[grb['GRBNAME']], ls='-', alpha=alpha, label=grblabel)
            ax[1].errorbar(self.dct_grbs[grb['GRBNAME']]['Index'].get_xdata(), self.dct_grbs[grb['GRBNAME']]['Index'].get_ydata(), xerr=self.dct_grbs[grb['GRBNAME']]['Index'].get_xerr(), yerr=self.dct_grbs[grb['GRBNAME']]['Index'].get_yerr(), lw=lwidth, ms=5, c=self.dct_ncolor[grb['GRBNAME']], fmt=self.dct_nmarker[grb['GRBNAME']], ls='-', alpha=alpha, label=grblabel)

           # >10GeV events
            vhe = PlotVHEEvents.SourceObject(name=grb['GRBNAME'], met0=pMETandMJD.ConvertMjdToMet(float(grb['TRIGGER_TIME'])), ra=grb['RA'], dec=grb['DEC'], gcnname=grb['GCNNAME'], tmin=0, tmax=100000, emin=10000, emax=120000, deg_roi=1.0, tprompt_start=grb['GBM']['T90_START'], redshift=0)
            #vhe = PlotVHEEvents.SourceObject(name=grb['GRBNAME'], met0=grb['TRIGGER_TIME'], ra=grb['RA'], dec=grb['DEC'], gcnname=grb['GCNNAME'], tmin=0, tmax=100000, emin=10000, emax=120000, deg_roi=6.0, tprompt_start=grb['GBM']['T90_START'], redshift=0)
            vhe.get_events()
            #vhe.get_events_angsep()
            logger.debug('>10GeV events:')
            logger.debug(vhe.times)
            logger.debug(vhe.energies)
            if len(vhe.times)>0:
                vhe.plot(ax[2], self.dct_nmarker[grb['GRBNAME']], self.dct_ncolor[grb['GRBNAME']], alpha=alpha) #, addphoton)
                ngrb_vhe+=1

        ax[0].set_xlabel(r'$T - \mathrm{T_0 [s]}$')
        ax[0].set_ylabel(r'Photon flux $\mathrm{[cm^{2} / s]}$')
        ax[0].set_xlim(0.1, 120000)
        #ax[0].set_ylim(-4, 0)
        ax[0].set_xscale("log", nonposx='clip')
        ax[0].set_yscale("log", nonposx='clip')
        ax[0].grid(ls='-', lw=0.5, alpha=0.5)
        #ax[0].set_yticks([y for y in ax[0].get_yticks() if y<0.5*ax[0].get_ylim()[1]])
        ax[0].legend(loc=0, fancybox=True, framealpha=0.5, fontsize=9, ncol=max(1,len(self.grbs_brighten)/3))

        ax[1].set_xlabel(r'$T - \mathrm{T_0 [s]}$')
        ax[1].set_ylabel('Photon index')
        ax[1].set_xlim(0.1, 120000)
        ax[1].set_ylim(-4, 0)
        ax[1].set_xscale("log", nonposx='clip')
        ax[1].grid(ls='-', lw=0.5, alpha=0.5)
        #ax[1].set_yticks([y for y in ax[1].get_yticks() if y<0.5*ax[1].get_ylim()[1]])
        ax[1].legend(loc=0, fancybox=True, framealpha=0.5, fontsize=9, ncol=max(1,len(self.grbs_brighten)/5))

        ax[2].set_xlabel(r'$T - \mathrm{T_0 [s]}$')
        ax[2].set_ylabel('Energy [GeV]')
        ax[2].set_xlim(0.1, 120000)
        ax[2].set_ylim(10., 120.)
        ax[2].set_xscale("log", nonposx='clip')
        ax[2].set_yscale("log", nonposx='clip')
        ax[2].grid(ls='-', lw=0.5, alpha=0.5)
        #ax[2].set_yticks([y for y in ax[2].get_yticks() if y<0.5*ax[2].get_ylim()[1]])
        ax[2].legend(loc=0, fancybox=True, framealpha=0.5, fontsize=9, ncol=ngrb_vhe/10)
        ax[2].set_ylim(10., 120.)

        fig.tight_layout()
        #fig.subplots_adjust(hspace=0)
        outbasename = 'IndexCurveCollection{suffix}{suffix2}'.format(suffix=self.suffix, suffix2=self.suffix2)
        for ff in self.figforms:
            fig.savefig('{0}/{1}.{2}'.format(outdir, outbasename, ff))


    def plot_swift(self, outdir='.'):
        fig, ax = plt.subplots(2, 2, figsize=(10, 10), sharex=False, sharey=False)

        curve_flux = {}
        curve_flux_kw = {}
        curve_flux_ratio = {}
        curve_flux_kw_ratio = {}
        dframe_xrt = ReadSwiftCatalogueInfo.open_table()
        tb_grb_swift = ReadLATCatalogueInfo.select_by_swift(self.table_grb, swift_dframe=dframe_xrt, columns=['XRT11HourFlux'])
        for igrb,grb in enumerate(tb_grb_swift):
            logger.info('###### {0} #####'.format(grb['GRBNAME']))
            if grb['GRBNAME'] in self.skip:
                logger.warning('{0} is skipped!'.format(grb['GRBNAME']))
                continue
            ##### XRT flux ####
            curve_flux[grb['GRBNAME']] = pMatplot.Curve('xrt_vs_gbm', xerr_asym=False, yerr_asym=False, xlabel=r'GBM Fluence (prompt) $\mathrm{[erg / cm^2]}$', ylabel=r'XRT Flux (@11/24hrs) $\mathrm{[10^{-11} erg/cm^{2}/s]$')
            curve_flux_ratio[grb['GRBNAME']] = pMatplot.Curve('xrt_vs_gbm', xerr_asym=False, yerr_asym=False, xlabel=r'GBM Fluence (prompt) $\mathrm{[erg/cm^2]}$', ylabel=r'XRT Flux / GBM Fluence $\mathrm{[10^-11 /s]}$')

            dseries = ReadSwiftCatalogueInfo.read_one_row(dframe_xrt, grb['GCNNAME'])
            fluen_prompt = grb['GBM']['FLUENCE'] if grb['ZENITH']<100 else None
            curve_flux[grb['GRBNAME']].set_point(grb['GBM']['FLUENCE'], float(dseries['XRT11HourFlux']))
            curve_flux[grb['GRBNAME']].set_point(grb['GBM']['FLUENCE'], float(dseries['XRT24HourFlux']))
            curve_flux_ratio[grb['GRBNAME']].set_point(grb['GBM']['FLUENCE'], float(dseries['XRT11HourFlux'])/grb['GBM']['FLUENCE'])
            curve_flux_ratio[grb['GRBNAME']].set_point(grb['GBM']['FLUENCE'], float(dseries['XRT24HourFlux'])/grb['GBM']['FLUENCE'])
            alpha = 1 if grb['GRBNAME'] in self.grbs_brighten else 0.25
            lwidth = 1 if grb['GRBNAME'] in self.grbs_brighten else 0.25
            grblabel = grb['GCNNAME'] if grb['GRBNAME'] in self.grbs_brighten else None
            ax[0][0].plot(curve_flux[grb['GRBNAME']].get_xdata(), curve_flux[grb['GRBNAME']].get_ydata(), marker=self.dct_nmarker[grb['GRBNAME']], c=self.dct_ncolor[grb['GRBNAME']], label=grblabel, alpha=alpha, lw=lwidth, ls='-')
            ax[1][0].plot(curve_flux_ratio[grb['GRBNAME']].get_xdata(), curve_flux_ratio[grb['GRBNAME']].get_ydata(), marker=self.dct_nmarker[grb['GRBNAME']], c=self.dct_ncolor[grb['GRBNAME']], label=grblabel, alpha=alpha, lw=lwidth, ls='-')
            logger.debug(grb)

            flag_KW_efluence = False
            if 'KonusWind' in grb:
                if isinstance(grb['KonusWind'],dict):
                    logger.debug(grb['KonusWind'])
                    if 'S' in grb['KonusWind']:
                        logger.debug('{0} has Konus-Wind data.'.format(grb['GRBNAME']))
                        kw_efluence = grb['KonusWind']['S']*1e-6
                        flag_KW_efluence = True
            if flag_KW_efluence==False and grb['GCNNAME'] in ReadKonusWindRedshiftCatalogueInfo.DCT_FLUENCE_GCN:
                kw_efluence = ReadKonusWindRedshiftCatalogueInfo.DCT_FLUENCE_GCN[grb['GCNNAME']]['value']
                flag_KW_efluence = True
            if flag_KW_efluence==True:
                curve_flux_kw[grb['GRBNAME']] = pMatplot.Curve('xrt_vs_kw', xerr_asym=False, yerr_asym=False, xlabel=r'Wind Fluence (prompt) $\mathrm{[erg/cm^2]}$', ylabel=r'XRT Flux (@11/24hrs) $\mathrm{[10^{-11} erg/cm^2/s]}$')
                curve_flux_kw_ratio[grb['GRBNAME']] = pMatplot.Curve('xrt_vs_kw', xerr_asym=False, yerr_asym=False, xlabel=r'Wind Fluence (prompt) $\mathrm{[erg/cm^2]}$', ylabel=r'XRT Flux / Wind Fluence $\mathrm{[10^{-11} /s]$')
                curve_flux_kw[grb['GRBNAME']].set_point(kw_efluence, float(dseries['XRT11HourFlux'])*1e-11)
                curve_flux_kw[grb['GRBNAME']].set_point(kw_efluence, float(dseries['XRT24HourFlux'])*1e-11)
                curve_flux_kw_ratio[grb['GRBNAME']].set_point(kw_efluence, float(dseries['XRT11HourFlux'])/kw_efluence*1e-11)
                curve_flux_kw_ratio[grb['GRBNAME']].set_point(kw_efluence, float(dseries['XRT24HourFlux'])/kw_efluence*1e-11)
                ax[0][1].plot(curve_flux_kw[grb['GRBNAME']].get_xdata(), curve_flux_kw[grb['GRBNAME']].get_ydata(), marker=self.dct_nmarker[grb['GRBNAME']], c=self.dct_ncolor[grb['GRBNAME']], label=grblabel, alpha=alpha, lw=lwidth, ls='-')
                ax[1][1].plot(curve_flux_kw_ratio[grb['GRBNAME']].get_xdata(), curve_flux_kw_ratio[grb['GRBNAME']].get_ydata(), marker=self.dct_nmarker[grb['GRBNAME']], c=self.dct_ncolor[grb['GRBNAME']], label=grblabel, alpha=alpha, lw=lwidth, ls='-')
            else:
                logger.warning('No Konus-Wind fluence data of {0}'.format(grb['GCNNAME']))
            
                    
        ax[0][0].set_xlabel('GBM Fluence (prompt) [erg/cm2]') #'XRT Flux (@11/24hrs) $\mathrm{[10^{-11} erg/cm^{2}/s]$')#curve_flux.values()[0].xlabel)#'GBM Fluence (prompt)')
        ax[0][0].set_ylabel('XRT Flux (@11/24hrs) [erg/cm2/s]') #curve_flux.values()[0].ylabel)#'XRT Flux @11/24hrs')
        ax[0][0].set_xlim(1e-6, 1e-2)
        ax[0][0].set_ylim(1e-15, 1e-10)
        ax[0][0].set_xscale("log", nonposx='clip')
        ax[0][0].set_yscale("log", nonposx='clip')
        ax[0][0].grid(ls='-', lw=0.5, alpha=0.5)
        ax[0][0].legend(loc=0, fancybox=True, framealpha=0.5, fontsize=9, ncol=max(1,len(self.grbs_brighten)/5))

        ax[0][1].set_xlabel('Wind Fluence (prompt) [erg/cm2]') #curve_flux_kw.values()[0].xlabel)#'Wind Fluence (prompt)')
        ax[0][1].set_ylabel('XRT Flux (@11/24 hrs) [erg/cm2/s]')#'XRT Flux @11/24hrs')
        ax[0][1].set_xlim(1e-6, 1e-2)
        ax[0][1].set_ylim(1e-15, 1e-10)
        ax[0][1].set_xscale("log", nonposx='clip')
        ax[0][1].set_yscale("log", nonposx='clip')
        ax[0][1].grid(ls='-', lw=0.5, alpha=0.5)
        ax[0][1].legend(loc=0, fancybox=True, framealpha=0.5, fontsize=9, ncol=max(1,len(self.grbs_brighten)/5))

        ax[1][0].set_xlabel('GBM Fluence (prompt) [erg/cm2]') #curve_flux_ratio.values()[0].xlabel)
        ax[1][0].set_ylabel('XRT Flux / GBM Fluence [erg/cm2/s]') #curve_flux_ratio.values()[0].ylabel)
        ax[1][0].set_xlim(1e-6, 1e-2)
        ax[1][0].set_ylim(0.3e-11, 3e-5)
        ax[1][0].set_xscale("log", nonposx='clip')
        ax[1][0].set_yscale("log", nonposx='clip')
        ax[1][0].grid(ls='-', lw=0.5, alpha=0.5)
        ax[1][0].legend(loc=0, fancybox=True, framealpha=0.5, fontsize=9, ncol=max(1,len(self.grbs_brighten)/5))

        ax[1][1].set_xlabel('Wind Fluence (prompt) [erg/cm2]') #curve_flux_kw_ratio.values()[0].xlabel)
        ax[1][1].set_ylabel('XRT Flux / Wind Fluence [/s]') #curve_flux_kw_ratio.values()[0].xlabel)
        ax[1][1].set_xlim(1e-6, 1e-2)
        ax[1][1].set_ylim(0.3e-11, 3e-5)
        ax[1][1].set_xscale("log", nonposx='clip')
        ax[1][1].set_yscale("log", nonposx='clip')
        ax[1][1].grid(ls='-', lw=0.5, alpha=0.5)
        ax[1][1].legend(loc=0, fancybox=True, framealpha=0.5, fontsize=9, ncol=max(1,len(self.grbs_brighten)/5))

        fig.tight_layout()
        outbasename = 'XRTCollection{suffix}{suffix2}'.format(suffix=self.suffix, suffix2=self.suffix2)
        for ff in self.figforms:
            fig.savefig('{0}/{1}.{2}'.format(outdir, outbasename, ff))


@click.command()
@click.option('--index', '-i', type=click.Choice(['free', 'best']), default='free')
@click.option('--erange', '-e', type=click.Choice(['whole', 'low', 'mid', 'high', 'lowmid', 'midhigh']), default='whole')
@click.option('--brighten', '-b', type=str, default=None, multiple=True)
@click.option('--suffix', '-s', type=str, default='LAT_N10in2deg')
@click.option('--suffix2', type=str, default='')
@click.option('--outdir', '-o', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/LongGRBs/Stacking/lightcurve')
@click.option('--figform', type=str, default=('png','pdf'), multiple=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(index, erange, brighten, outdir, suffix, suffix2, figform, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    DICT_DIRS = {'whole':'E0000100-0100000MeV/r12deg/lightcurve',
                 'low':'E0000100-0001000MeV/r12deg/lightcurve',
                 'mid':'E0001000-0010000MeV/r03deg/lightcurve',
                 'high':'E0010000-0100000MeV/r01deg/lightcurve',
                 'lowmid':'E0000100-0010000MeV/r12deg/lightcurve',
                 'midhigh':'E0001000-0100000MeV/r03deg/lightcurve'}

    ##### Catalogue information #####
    tb_lat = ReadLATCatalogueInfo.read_all(TABLE_LAT, TABLE_GBM, TABLE_KW)
    tb_lat = ReadLATCatalogueInfo.select_gbm_exist(tb_lat)
    tb_lat = ReadLATCatalogueInfo.select_long(tb_lat)
    tb_lat = ReadLATCatalogueInfo.select_small_error(tb_lat, 0.3)
    #tb_lat = ReadLATCatalogueInfo.select_zenith_angle(tb_lat, 100.)

    ##### Index curve Collection #####
    lcc = ICCollection(index=index, erange=erange, brighten=brighten, interdirs=DICT_DIRS[erange], suffix=suffix, suffix2=suffix2, table_catalogue=tb_lat, figforms=figform)
    lcc.plot(outdir=outdir)
    lcc.plot_swift(outdir=outdir)

if __name__ == '__main__':
    main()
