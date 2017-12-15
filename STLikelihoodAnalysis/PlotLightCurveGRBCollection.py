#!/usr/bin/env python
"""Module for making light curve collection of GRB LAT data.
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
import pConvertFluxValues
from pConvertFluxValues import MEVtoERG, KEVtoERG
import pickle_utilities
import pMatplot
from pMatplot import find_range_shown, find_indexrange, TPL_MARKER, TPL_COLOR, TPL_LINE


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

TIME_INTERVALS = {'prompt':'T0 to T95', 'T95to03ks':'T95 to T95+3ks', '03ksto100ks':'T95+3ks to 100ks'}
ENERGY_RANGES = {'whole': (100, 100000), 'low': (100, 1000), 'mid': (1000, 10000), 'high':(10000, 100000)}

powerlaws = {}
powerlaws[0.1] = lambda t: 1.60E-03*pow(t/10., -1.345)
powerlaws[1.0] = lambda t: 1.68E-03*pow(t/10., -1.325)
powerlaws[10.] = lambda t: 4.62E-04*pow(t/10., -1.085)


class LCCollection:
    def __init__(self, index, erange, dframe_xrt, interdirs='E0000100-0100000MeV/r12deg/lightcurve', suffix='LAT_N10in2deg', table_catalogue=None, figforms=['pdf', 'png']):
        #self.name = name
        self.index = index
        self.skip=('130427324', '150314205')
        self.erange = erange
        self.interdirs = interdirs
        if table_catalogue is None:
            tb_lat = ReadLATCatalogueInfo.read_all(TABLE_LAT, TABLE_GBM)
            tb_lat = ReadLATCatalogueInfo.select_gbm_exist(tb_lat)
            tb_lat = ReadLATCatalogueInfo.select_long(tb_lat)
            tb_lat = ReadLATCatalogueInfo.select_zenith_angle(tb_lat, 100.)
            table_catalogue = ReadLATCatalogueInfo.select_small_error(tb_lat, 0.3)
        self.table_grb = table_catalogue
        self.dframe_xrt = dframe_xrt
        self.suffix = suffix if suffix=='' else '_'+suffix
        self.str_scalefactor='GBMfluence'
        self.figforms = figforms


    def load_grbs(self, name, dirs='E0000100-0100000MeV/r12deg/lightcurve'):
        self.dataset = {}
        for grb in self.table_grb:
            if grb['GRBNAME'] in skip:
                continue
            logger.debug('Loading likelihood valus of {0}'.format(grb['GRBNAME']))
            path_likelihood = '/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{name}/{dirs}/LightCurve_{name}_PowerLaw_IndexFree_LAT_N10in2deg.pickle'.format(name=grb['GRBNAME'], dirs=dirs)
            self.dataset[grb['GRBNAME']] = pickle_utilities.load(path_likelihood)


    # def load_swift_table(self, remove_nan_thr=(11, 24)):
    #     self.dframe_xrt = ReadSwiftCatalogueInfo.open_table()
    #     # for thr in remove_nan_thr:
    #     #     self.dframe_xrt = self.dframe_xrt.query("XRT{t}HourFlux != 'n/a'".format(t=thr))
    #     # self.dframe_xrt = self.dframe_xrt.query("XRTSpectralIndex != 'n/a'")
    #     self.grbs_xrt = self.dframe_xrt['GRB'].values
    #     logger.info('There exist {0} GRBs with Swift data.'.format(len(self.grbs_xrt)))


    def load_model(self, path_model='/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/LongGRBs/Stacking/lightcurve/AsanoModel_example20171207.txt'):
        lst_time = []
        lst_vFv_100MeV = []
        lst_vFv_001GeV = []
        lst_vFv_010GeV = []
        with open(path_model, 'r') as m:
            reader = csv.reader(m, delimiter=' ')
            header = next(reader)
            for row in reader:
                lst_time.append(float(row[0]))
                lst_vFv_100MeV.append(float(row[1]))
                lst_vFv_001GeV.append(float(row[2]))
                lst_vFv_010GeV.append(float(row[3]))
        npa_time = np.array(lst_time)
        dct_npa_vFv = {0.1:np.array(lst_vFv_100MeV), 1.0:np.array(lst_vFv_001GeV), 10.0:np.array(lst_vFv_010GeV)}
        return (npa_time, dct_npa_vFv)


    def plot(self, outdir='.', ts_threshold=4, xrtlimit=True):
        lst_markers = ['s', 'o', 'D', 'x', 'd', 'P', '*', 'h']

        fig, ax = plt.subplots(1, 3, figsize=(30, 10), sharex=True)

        dct_grbs = {}
        for igrb,grb in enumerate(self.table_grb):
            logger.info('###### {0} #####'.format(grb['GRBNAME']))
            if grb['GRBNAME'] in self.skip:
                logger.warning('Skipped.')
                continue
            #if xrtlimit==True and not grb['GCNNAME'] in self.grbs_xrt:
            #    logger.warning('No Swift data. Skipped.')
            #    continue
            
            logger.debug('Loading likelihood valus of {0}'.format(grb['GRBNAME']))
            path_likelihood = '/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{name}/{dirs}/LightCurve_{name}_PowerLaw_IndexFree{suf}.pickle'.format(name=grb['GRBNAME'], dirs=self.interdirs, suf=self.suffix)
            grbdata = pickle_utilities.load(path_likelihood)['results']

            dct_grbs[grb['GRBNAME']] = {}
            # e2dnde and UL
            dct_grbs[grb['GRBNAME']]['e2dnde'] = pMatplot.Curve('e2dnde', xerr_asym=True, yerr_asym=True, xlabel='Time - T0 [s]', ylabel=r'$E^2 dN/dE \, \rm{{at}} \, {ene:3.1f} \rm{{GeV}} \, \mathrm{{[erg/cm^2 s]}}$'.format(ene=1.0))
            dct_grbs[grb['GRBNAME']]['e2dnde_ul'] = pMatplot.Curve('e2dnde', xerr_asym=True, yerr_asym=False, ul=True, xlabel='Time - T0 [s]', ylabel=r'$E^2 dN/dE \, \rm{{at}} \, {ene:3.1f} \rm{{GeV}} \, \mathrm{{[erg/cm^2 s]}}$'.format(ene=1.0))
            for ic, curve in enumerate(dct_grbs[grb['GRBNAME']].values()):
                logger.info('===== {0} ====='.format(curve.quantity))
                for period in grbdata:
                    t0 = max(1, period['time']['min'])
                    t1 = period['time']['max']
                    tref = 10**((np.log10(t0)+np.log10(t1))/2.0)
                    logger.debug('----- {tmin:.1f} - {tmax:.1f} s -----'.format(tmin=t0, tmax=t1))
                    if not 'limits' in period:
                        logger.warning('No limits. Skipped.')
                        continue

                    if curve.quantity == 'TS':
                        y = sqrt(max(0, period[curve.quantity]))
                        yerr = 0
                        logger.debug('{v:.2}'.format(v=y))
                    if curve.ul==False and period['TS']>=ts_threshold:
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

                    elif (curve.ul==True or curve.quantity in ('TS')) and period['TS']<ts_threshold: #and period['Nobs']>0:
                        if curve.quantity in ('flux', 'eflux', 'e2dnde'):
                            y = period['limits']['best'][curve.quantity]['ul']
                            yerr = 0
                            if curve.quantity in ('eflux', 'e2dnde'):
                                y = y*MEVtoERG
                            logger.debug('UL: {ul:.2}'.format(ul=y))
                        curve.set_point(tref, y, {'lo':t1-tref, 'hi':tref-t0}, yerr)
            ##### XRT flux ####
            dct_grbs[grb['GRBNAME']]['xrt'] = pMatplot.Curve('xrt', xerr_asym=False, yerr_asym=False, xlabel='Time - T0 [s]', ylabel=r'$E^2 dN/dE \, \rm{{at}} \, {ene:3.1f} \rm{{keV}} \, \mathrm{{[erg/cm^2 s]}}$'.format(ene=1.0))
            dseries = ReadSwiftCatalogueInfo.read_one_row(self.dframe_xrt, grb['GCNNAME'])
            t_xrt = np.array([11, 24])*3600.
            #logger.info(dseries['XRT11HourFlux'])
            flux_xrt = np.array([float(dseries['XRT11HourFlux']), 
                                 float(dseries['XRT24HourFlux'])])
            spindex_xrt = -float(dseries['XRTSpectralIndex'])

            yxrt = pConvertFluxValues.convert_eflux_to_vFv(flux_xrt, spindex=spindex_xrt, eref=1.0, emin=0.3, emax=10)*1E-11
            dct_grbs[grb['GRBNAME']]['xrt'].set_point(t_xrt, yxrt)

            nmarker = igrb / NMARKER_STYLE
            ncolor = igrb % len(TPL_COLOR[1:])
            logger.info('Marker style: {0}'.format(lst_markers[nmarker]))
            logger.info('Marker color: {0}'.format(TPL_COLOR[1:][ncolor]))
            ax[0].errorbar(dct_grbs[grb['GRBNAME']]['e2dnde'].get_xdata(), dct_grbs[grb['GRBNAME']]['e2dnde'].get_ydata(), xerr=dct_grbs[grb['GRBNAME']]['e2dnde'].get_xerr(), yerr=dct_grbs[grb['GRBNAME']]['e2dnde'].get_yerr(), lw=0.25, ms=5, c=TPL_COLOR[1:][ncolor], fmt=lst_markers[nmarker])#, fmt=dct_grbs[grb['GRBNAME']]['e2dnde'].fmt
            ax[0].errorbar(dct_grbs[grb['GRBNAME']]['e2dnde_ul'].get_xdata(), dct_grbs[grb['GRBNAME']]['e2dnde_ul'].get_ydata(), ms=5, fmt=dct_grbs[grb['GRBNAME']]['e2dnde_ul'].fmt, markeredgecolor=TPL_COLOR[1:][ncolor], markerfacecolor='none', alpha=0.5)
            ax[0].plot(dct_grbs[grb['GRBNAME']]['xrt'].get_xdata(), dct_grbs[grb['GRBNAME']]['xrt'].get_ydata(), lw=0.25, ms=5, markeredgecolor=TPL_COLOR[1:][ncolor], markerfacecolor='none', marker=lst_markers[nmarker])

            ax[1].errorbar(dct_grbs[grb['GRBNAME']]['e2dnde'].get_xdata(), dct_grbs[grb['GRBNAME']]['e2dnde'].get_ydata()/grb['GBM']['FLUENCE'], xerr=dct_grbs[grb['GRBNAME']]['e2dnde'].get_xerr(), yerr=dct_grbs[grb['GRBNAME']]['e2dnde'].get_yerr()/grb['GBM']['FLUENCE'], ms=5, lw=0.25, c=TPL_COLOR[1:][ncolor], fmt=lst_markers[nmarker]) #, fmt=dct_grbs[grb['GRBNAME']]['e2dnde'].fmt
            ax[1].errorbar(dct_grbs[grb['GRBNAME']]['e2dnde_ul'].get_xdata(), dct_grbs[grb['GRBNAME']]['e2dnde_ul'].get_ydata()/grb['GBM']['FLUENCE'], ms=5, fmt=dct_grbs[grb['GRBNAME']]['e2dnde_ul'].fmt, markeredgecolor=TPL_COLOR[1:][ncolor], markerfacecolor='none', alpha=0.5)#, xerr=dct_grbs[grb['GRBNAME']]['e2dnde_ul'].get_xerr())   
            ax[1].plot(dct_grbs[grb['GRBNAME']]['xrt'].get_xdata(), dct_grbs[grb['GRBNAME']]['xrt'].get_ydata()/grb['GBM']['FLUENCE'], lw=0.25, ms=5, markeredgecolor=TPL_COLOR[1:][ncolor], markerfacecolor='none', marker=lst_markers[nmarker])

            ax[2].errorbar(dct_grbs[grb['GRBNAME']]['e2dnde'].get_xdata(), dct_grbs[grb['GRBNAME']]['e2dnde'].get_ydata()/grb['GBM']['FLUENCE']/powerlaws[1.0](dct_grbs[grb['GRBNAME']]['e2dnde'].get_xdata()), xerr=dct_grbs[grb['GRBNAME']]['e2dnde'].get_xerr(), yerr=dct_grbs[grb['GRBNAME']]['e2dnde'].get_yerr()/grb['GBM']['FLUENCE']/powerlaws[1.0](dct_grbs[grb['GRBNAME']]['e2dnde'].get_xdata()), ms=5, lw=0.25, c=TPL_COLOR[1:][ncolor], fmt=lst_markers[nmarker])
            ax[2].errorbar(dct_grbs[grb['GRBNAME']]['e2dnde_ul'].get_xdata(), dct_grbs[grb['GRBNAME']]['e2dnde_ul'].get_ydata()/grb['GBM']['FLUENCE']/powerlaws[1.0](dct_grbs[grb['GRBNAME']]['e2dnde_ul'].get_xdata()), ms=5, fmt=dct_grbs[grb['GRBNAME']]['e2dnde_ul'].fmt, markeredgecolor=TPL_COLOR[1:][ncolor], markerfacecolor='none', alpha=0.5)

        #ax[0].set_ylabel(dct_grbs[grb['GRBNAME']]['e2dnde'].ylabel)
        ax[0].set_ylabel(r'$\nu F_{\nu} \mathrm{[erg/cm^2 \cdot s]}$')
        ax[0].set_xlim(1, 120000)
        ax[0].set_xscale("log", nonposx='clip')
        ax[0].set_yscale("log", nonposx='clip')
        ax[0].grid(ls='-', lw=0.5, alpha=0.5)
        ax[0].set_yticks([y for y in ax[0].get_yticks() if y<0.5*ax[0].get_ylim()[1]])
        #ax[0].legend(loc=0, fancybox=True, framealpha=0.5, fontsize=9, ncol=10)
        ax[1].set_ylabel(r'$\nu F_{\nu} / F_{GBM} \mathrm{[/s]}$')
        ax[1].set_xscale("log", nonposx='clip')
        ax[1].set_yscale("log", nonposx='clip')
        ax[1].grid(ls='-', lw=0.5, alpha=0.5)
        ax[1].set_yticks([y for y in ax[1].get_yticks() if y<0.5*ax[1].get_ylim()[1]])
        tpl = np.array([10,10000])
        vFv_pl_100MeV = 1.60E-03*pow(tpl/10., -1.345)
        vFv_pl_001GeV = 1.68E-03*pow(tpl/10., -1.325)
        vFv_pl_010GeV = 4.62E-04*pow(tpl/10., -1.085)
        ax[1].plot(tpl, vFv_pl_100MeV, '-', label='PL 0.1 GeV', lw=2, c='k', alpha=0.5)
        ax[1].plot(tpl, vFv_pl_001GeV, '--', label='PL 1.0 GeV', lw=2, c='k', alpha=0.5)
        ax[1].plot(tpl, vFv_pl_010GeV, ':', label='PL 10.0 GeV', lw=2, c='k', alpha=0.5)
        tssc, vFv_ssc = self.load_model()
        for ie, e in enumerate((0.1, 1.0, 10.0)):
            ax[1].plot(tssc, vFv_ssc[e], TPL_LINE[ie], label='Model {0:1.1f} GeV'.format(e), lw=2, c='r', alpha=0.5)
        ax[1].legend()

        ax[2].set_ylabel(r'$\nu F_{\nu} / F_{GBM} \mathrm{[/s]} / (1.68 \times 10^{-3}(T/10 \mathrm{[s]})^{-1.325})$')
        ax[2].set_xscale("log", nonposx='clip')
        ax[2].set_yscale("log", nonposx='clip')
        ax[2].grid(ls='-', lw=0.5, alpha=0.5)
        ax[2].set_yticks([y for y in ax[2].get_yticks() if y<0.5*ax[2].get_ylim()[1]])
        tpl = np.array([10,10000])
        tssc, vFv_ssc = self.load_model()
        for ie, e in enumerate((0.1, 1.0, 10.0)):
            ax[2].plot(tssc, vFv_ssc[e]/powerlaws[e](tssc), TPL_LINE[ie], label='Model {0:1.1f} GeV'.format(e), lw=2, c='r', alpha=0.5)

        fig.tight_layout() #subplots_adjust(hspace=0)
        fig.subplots_adjust(hspace=0)
    
        outbasename = 'LightCurveCollection{suffix}'.format(suffix=self.suffix)
        for ff in self.figforms:
            fig.savefig('{0}/{1}.{2}'.format(outdir, outbasename, ff))


@click.command()
@click.option('--index', '-i', type=click.Choice(['free', 'best']), default='free')
@click.option('--erange', '-e', type=click.Choice(['whole', 'low', 'mid', 'high']), default='whole')
@click.option('--suffix', '-s', type=str, default='LAT_N10in2deg')
@click.option('--outdir', '-o', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/LongGRBs/Stacking/lightcurve')
@click.option('--figform', type=str, default=('png','pdf'), multiple=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(index, erange, outdir, suffix, figform, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    DICT_DIRS = {'whole':'E0000100-0100000MeV/r12deg/lightcurve',
                 'low':'E0000100-0001000MeV/r12deg/lightcurve',
                 'mid':'E0001000-0010000MeV/r03deg/lightcurve',
                 'high':'E0010000-0100000MeV/r01deg/lightcurve'}

    ##### Catalogue information #####
    tb_lat = ReadLATCatalogueInfo.read_all(TABLE_LAT, TABLE_GBM)
    tb_lat = ReadLATCatalogueInfo.select_gbm_exist(tb_lat)
    tb_lat = ReadLATCatalogueInfo.select_long(tb_lat)
    tb_lat = ReadLATCatalogueInfo.select_small_error(tb_lat, 0.3)
    tb_lat = ReadLATCatalogueInfo.select_zenith_angle(tb_lat, 100.)
    dframe_xrt = ReadSwiftCatalogueInfo.open_table()
    tb_lat = ReadLATCatalogueInfo.select_by_swift(tb_lat, dframe_xrt, ['XRT11HourFlux', 'XRT24HourFlux', 'XRTSpectralIndex'])

    ##### Lightcurve Collection #####
    lcc = LCCollection(index=index, erange=erange, dframe_xrt=dframe_xrt, interdirs=DICT_DIRS[erange], suffix=suffix, table_catalogue=tb_lat, figforms=figform)
    #lcc.load_swift_table()    
    lcc.plot(outdir=outdir)

if __name__ == '__main__':
    main()
