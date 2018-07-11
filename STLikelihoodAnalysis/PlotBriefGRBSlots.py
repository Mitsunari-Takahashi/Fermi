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
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
import pickle
import numpy as np
from scipy.interpolate import interp1d
from scipy import optimize
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


##### Conversion from MeV to erg ######
MEVtoERG = 1.6021766208E-6


##### Catalogues #####
GRB_CATALOGUE_LAT = "/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/LATBurstCatalogue.xml"
GRB_CATALOGUE_GBM = '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/Regualr/catalogue/GBMcatalogue20171005.fits'

TABLE_LAT = ReadLATCatalogueInfo.open_table() #GRB_CATALOGUE_LAT)
TABLE_GBM = ReadGBMCatalogueInfo.open_table() #GRB_CATALOGUE_GBM)

TIME_INTERVALS = {'prompt':'T0 to T95', 'afterglow':'T0 to T95', 'T95to03ks':'T95 to T95+3ks', '03ksto100ks':'T95+3ks to 100ks'}
ENERGY_RANGES = {'whole': (100, 100000), 'low': (100, 1000), 'mid': (1000, 10000), 'high':(10000, 100000), 'lowmid':(100, 10000), 'midhigh':(1000, 100000)}


class ModelEnsemble:
    def __init__(self, name, erange, twindow, indir, norms, lcindices, specindices, suffix='', str_scalefactor='GBMfluence', table_catalogue=None):
        self.name = name
        self.erange = erange
        self.str_input_dir = indir
        self.norms = norms
        self.lcindices = lcindices
        self.lcindices_mesh, self.norms_mesh = np.meshgrid(self.lcindices, self.norms)
        logger.debug('{0}'.format(self.lcindices_mesh))
        self.specindices = specindices
        self.str_time_window = twindow
        if table_catalogue is None:
            tb_lat = ReadLATCatalogueInfo.read_all(TABLE_LAT, TABLE_GBM)
            tb_lat = ReadLATCatalogueInfo.select_gbm_exist(tb_lat)
            tb_lat = ReadLATCatalogueInfo.select_long(tb_lat)
            table_catalogue = ReadLATCatalogueInfo.select_small_error(tb_lat, 0.3)
        self.table_grb = table_catalogue
        self.suffix = suffix if suffix=='' else '_'+suffix
        self.cont_levels = (2.30, 5.99, 11.83) #[2.31, 4.61, 9.21]
        self.table_xrt = ReadSwiftCatalogueInfo.open_table()


    def load_values(self, target, quantity):
        path_likelihood = '/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{grb}/{dirs}/briefslots/LightCurve_{grb}_PowerLaw_IndexFree_pm3dec_180x090.pickle'.format(grb=target, dirs=self.str_input_dir)
#        path_likelihood = '/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/{grb}/{dirs}/briefslots/LightCurve_{grb}_PowerLaw_IndexFree_E100MeV-010GeV_normx300_100x090.pickle'.format(grb=target, dirs=self.str_input_dir)
        lst_slots = pickle_utilities.load(path_likelihood)
        lst_loglike = []
        flag = False
        for slot in lst_slots:
            if self.str_time_window=='T95to03ks':
                if slot['period'][1]<=3000:
                    lst_loglike.append(slot[quantity])
                    logger.debug('  Loaded {0} - {1} s'.format(slot['period'][0], slot['period'][1]))
            elif self.str_time_window=='03ksto100ks':
                if slot['period'][1]>3000 and slot['period'][1]<=100000:
                    lst_loglike.append(slot[quantity])
                    logger.debug('  Loaded {0} - {1} s'.format(slot['period'][0], slot['period'][1]))
            elif self.str_time_window=='T95to01ks':
                if slot['period'][1]<=1000:
                    lst_loglike.append(slot[quantity])
                    logger.debug('  Loaded {0} - {1} s'.format(slot['period'][0], slot['period'][1]))
            elif self.str_time_window=='01ksto100ks':
                if slot['period'][1]>1000 and slot['period'][1]<=100000:
                    lst_loglike.append(slot[quantity])
                    logger.debug('  Loaded {0} - {1} s'.format(slot['period'][0], slot['period'][1]))
            else:
                lst_loglike.append(slot[quantity])
        return lst_loglike


    def joint_likelihood(self):
        loglike_sum_all = np.zeros(shape=(len(self.norms), len(self.lcindices)))
        for grb in self.table_grb:
            logger.debug('Loading likelihood valus of {0}'.format(grb['GRBNAME']))
            loglikes = self.load_values(grb['GRBNAME'], 'loglike')
            loglike_sum_all += np.sum(loglikes, axis=0)
        self.loglike_sum_all = loglike_sum_all
        return loglike_sum_all


    def delta_loglikelihood(self):
        self.loglike_max = self.loglike_sum_all.max()
        self.dloglike = -self.loglike_sum_all + self.loglike_max
        self.dloglike_doubled = self.dloglike * 2.
        self.args_max = zip(np.where(self.loglike_sum_all==self.loglike_sum_all.max())[0], np.where(self.loglike_sum_all==self.loglike_sum_all.max())[1])[0]
        logger.debug('Best indexes: ({0}, {1}'.format(self.args_max[0], self.args_max[1]))
        logger.info('Maximum likelihood: {v} at norm={n}, index={i}'.format(v=self.loglike_max, n=self.norms[self.args_max[0]], i=self.lcindices[self.args_max[1]]))
        return (self.args_max, self.dloglike)


    # def eval_dispersion(self):
    #     dct_disp = {}
    #     for grb in self.table_grb:
    #         loglikes = self.load_values(grb['GRBNAME'], 'loglike')
    #         dct_disp[grb['GRBNAME']] = np.sum(loglikes, axis=0)[self.args_max[0]][self.args_max[1]]
    #     self.dispersion


    def get_e2dnde(self, teval, eeval):
        #val = self.norms[self.args_max[0]]*eeval*eeval*MEVtoERG * pow(teval/10., self.lcindices[self.args_max[1]])
        vals = self.norms_mesh*eeval*eeval*MEVtoERG * pow(teval/10., self.lcindices_mesh)
        logger.info('E^2dN/dE(T={t:.1f}s,E={e:1.1E}MeV) = {e2dnde:1.2E} erg/cm^2/s'.format(t=teval, e=eeval, e2dnde=vals[self.args_max[0], self.args_max[1]]))
        return vals


    def get_reference_time(self, t0, t1):
        tref = (self.lcindices_mesh!=-1) * pow((pow(t0, self.lcindices_mesh+1) + pow(t1, self.lcindices_mesh+1))/2., 1./((self.lcindices_mesh!=-1)*self.lcindices_mesh+1)) + (self.lcindices_mesh==-1) * np.exp((np.log(t1) + np.log(t0))/2. )
        logger.debug('Period: {0} - {1}'.format(t0, t1))
        logger.info('Reference time: {0}'.format(tref))
        return tref


    def model_courve(self):
        tb_lat = ReadLATCatalogueInfo.read_all(TABLE_LAT, TABLE_GBM)
        for grb in self.table_grb:
            logger.debug(grb['GRBNAME'])
            tb_one = ReadLATCatalogueInfo.select_one_by_name(tb_lat, grb['GRBNAME'])
            t95 = tb_one['GBM']['T90_START'] + tb_one['GBM']['T90']
            loglikes = self.load_values(grb['GRBNAME'], 'loglike')
            fluences = self.load_values(grb['GRBNAME'], 'efluence')
            eflux = self.load_values(grb['GRBNAME'], 'eflux')
            prefactors = self.load_values(grb['GRBNAME'], 'normalization')
            periods = self.load_values(grb['GRBNAME'], 'period')

            for iperiod, period in enumerate(periods):
                tref = self.get_reference_time(periods[0][0], periods[0][1])
                prefactors_tref = prefactors * pow(tref/t95, self.lcindices_mesh)


    def renormalize_prefactor(self, tnorm, enorm, emin, emax, specindex, path_save, name_save, figforms=('png', 'pdf')):
        #tb_lat = ReadLATCatalogueInfo.read_all(TABLE_LAT, TABLE_GBM)
        prefactor_norm_best = {}
        gbm_fluence = {}
        fig, ax = plt.subplots(1, 2, sharex=True, sharey=False, figsize=(10, 7))
        for grb in self.table_grb:
            logger.debug(grb['GRBNAME'])
            t95 = grb['GBM']['T90_START'] + grb['GBM']['T90']
            prefactors = self.load_values(grb['GRBNAME'], 'normalization')[0]
            firstslot = self.load_values(grb['GRBNAME'], 'period')[0]
            logger.debug('First slot: {0}'.format(firstslot))
            gbm_fluence[grb['GRBNAME']] = grb['GBM']['FLUENCE']

            logger.debug("Prefactors' shape: {0}".format(prefactors.shape))
            logger.debug("LC indices' shape: {0}".format(self.lcindices_mesh.shape))
            #prefactors_tnorm = prefactors * pow(tnorm/t95, self.lcindices_mesh)
            prefactors_tnorm = prefactors * (self.lcindices_mesh+1) * pow(tnorm, self.lcindices_mesh) / (pow(firstslot[1], self.lcindices_mesh+1) - pow(firstslot[0], self.lcindices_mesh+1))
            prefactors_tnorm_enorm = prefactors_tnorm * (specindex+1) / (pow(emax, specindex+1) - pow(emin, specindex+1)) * pow(enorm, specindex)
            logger.debug("Normalized prefactors' shape: {0}".format(prefactors_tnorm_enorm.shape))
            logger.debug('Best indexes: {0}'.format(self.args_max))
            prefactor_norm_best[grb['GRBNAME']] = prefactors_tnorm_enorm[self.args_max[0]][self.args_max[1]] * MEVtoERG

        prefactor_norm_best_sorted = np.array([ v for k, v in sorted(prefactor_norm_best.items()) ])
        gbm_fluence_sorted = np.array([ v for k, v in sorted(gbm_fluence.items()) ])
        ax[0].hist(np.log10(prefactor_norm_best_sorted))
        ax[0].set_title('Distrubution of the normalization factor at {0} s'.format(tnorm))
        ax[0].set_xlabel(r'$ \log(\nu F_{\nu} / E_{GBM}) \rm{[/s]} $')
        ax[0].set_ylabel('[bursts]')
        ax[0].grid()
        #ax[1].hist(prefactor_norm_best_sorted, weights=gbm_fluence_sorted, normed=True)
        ax[1].hist(np.log10(prefactor_norm_best_sorted), weights=gbm_fluence_sorted, normed=True)
        ax[1].set_title('Weighted by GBM fluence')
        ax[1].set_xlabel(r'$ \log(\nu F_{\nu} / E_{GBM}) \rm{[/s]} $')
        ax[1].set_ylabel('normalized')
        ax[1].grid()
        for ff in figforms:
            fig.savefig('{dire}/{name}_{tn:0>5.0f}sec_{en:0>7.0f}MeV.{form}'.format(dire=path_save, name=name_save, tn=tnorm, en=enorm, form=ff))


    def plot_count(self, path_save, name_save, figforms=('png', 'pdf')):
        fig, ax = plt.subplots(1, 4, sharex=False, sharey=False, figsize=(20, 5))
        resid_best = {}
        ts_best = {}
        ts_sqrt_best = {}
        fluxes_xrt = {}
        lst_fluxesxrt = []
        lst_resid_xrt = []
        nresids = {}
        trefs = {}
        for grb in self.table_grb:
            logger.debug(grb['GRBNAME'])
            nobss = self.load_values(grb['GRBNAME'], 'nobs')
            nobs = sum(nobss)
            npred = self.load_values(grb['GRBNAME'], 'npred')
            npred_target = 0
            npred_others = 0

            for ip, npp in enumerate(npred):
                npred_target = npp['target'][self.args_max[0]][self.args_max[1]] if ip==0 else npred_target+npp['target'][self.args_max[0]][self.args_max[1]]
                npred_others = npp['others'][self.args_max[0]][self.args_max[1]] if ip==0 else npred_others+npp['others'][self.args_max[0]][self.args_max[1]]
            logger.debug('Nobs: {0}'.format(nobs))
            logger.debug('Npred (GRB): {0}'.format(npred_target))
            logger.debug('Npred (others): {0}'.format(npred_others))
            npred_total =  npred_target+npred_others
            resid = (nobs-(npred_total))# / np.sqrt(npred_total) #_target
            if npred_total>100:
                ts = pow(nobs-npred_total,2) / npred_total
            elif nobs>0:
                ts = -2 * ( np.log(pow(npred_total, float(nobs))*np.exp(-npred_total)/float(np.math.factorial(nobs))) - np.log(pow(float(nobs), float(nobs))*np.exp(-float(nobs))/float(np.math.factorial(nobs))) )
            else:
                ts = -2 * ( np.log(np.exp(-npred_total)) )
            resid_best[grb['GRBNAME']] = resid
            ts_best[grb['GRBNAME']] = ts
            ts_sqrt_best[grb['GRBNAME']] = min(np.sqrt(ts_best[grb['GRBNAME']]), 5)
            str_info = """* GRB {name}
  Fractional residual: {resid}
  TS: {ts}
  Nobs: {obs}
  Npred (GRB): {pred1}
  Npred (others): {pred0}
""".format(name=grb['GRBNAME'], resid=resid_best[grb['GRBNAME']], ts=ts_best[grb['GRBNAME']], obs=nobs, pred1=npred_target, pred0=npred_others)
            if ts_best[grb['GRBNAME']] > 3:
                logger.warning(str_info)
            else:
                logger.info(str_info)

            # XRT
            grb_xrt = ReadSwiftCatalogueInfo.read_one_row(self.table_xrt, grb['GCNNAME'])
            #logger.info(grb_xrt)
            if len(grb_xrt.index)>0:
                flux_xrt = grb_xrt['XRT11HourFlux']#.values[0]
                if flux_xrt!='n/a':
                    if float(flux_xrt)>0:
                        fluxes_xrt[grb['GRBNAME']] = float(flux_xrt)
                        lst_fluxesxrt.append(float(flux_xrt)/grb['GBM']['FLUENCE'])
                        lst_resid_xrt.append(resid_best[grb['GRBNAME']])
                        logger.debug('XRT flux: {0}x10^-11erg/cm^2/s'.format(fluxes_xrt[grb['GRBNAME']]))
                        trefs[grb['GRBNAME']] = np.array([t[self.args_max[0]][self.args_max[1]] for t in self.load_values(grb['GRBNAME'], 'tref')])
                        nresids[grb['GRBNAME']] = np.zeros_like(nobss)
                        for kperiod, (obs, pred) in enumerate(zip(nobss, npred)):
                            nresids[grb['GRBNAME']][kperiod] = obs - (pred['target'][self.args_max[0]][self.args_max[1]]+pred['others'][self.args_max[0]][self.args_max[1]])

        ax[0].hist(resid_best.values(), bins=600, range=(-100, 500))
        ax[0].set_title('Fractional residual')
        ax[0].set_xlabel(r'$(N_{obs}-N_{pred})/ \sqrt{N_{pred}}$')
        ax[0].set_ylabel('[bursts]')
        ax[0].grid()

        ax[1].hist(ts_sqrt_best.values(), bins=51, range=(0, 5.1))
        ax[1].set_title(r'Deviation of $N_{pred}$ from $N_{obs}$')
        ax[1].set_xlabel(r'Min($\sqrt{TS}$, 5)')
        ax[1].set_ylabel('[bursts]')
        ax[1].grid()

        ax[2].scatter(lst_fluxesxrt, lst_resid_xrt)
        ax[2].set_xlim((2e2, 2e6))
        ax[2].set_xscale('log')
        ax[2].set_xlabel(r'XRT 11Hour Flux / GBM Fluence $\rm{[10^{-11}/s]}$')
        ax[2].set_ylabel(r'$(N_{obs}-N_{pred})/ \sqrt{N_{pred}}$')
        ax[2].grid()        

        for grbname, fluxxrt in fluxes_xrt.items():
            z = np.log10(fluxxrt/ReadLATCatalogueInfo.select_one_by_name(TABLE_LAT, grbname, TABLE_GBM)['GBM']['FLUENCE'])
            logger.info('{0}: {1}'.format(grbname, z))
            ax[3].scatter(trefs[grbname], nresids[grbname], s=z*z*2, alpha=0.5)
        #ax[3].set_xlim((2e2, 2e6))
        ax[3].set_xscale('log')
        ax[3].set_xlabel('Time [s]')
        ax[3].set_ylabel(r'$N_{obs}-N_{pred}$')
        ax[3].grid()        
        for ff in figforms:
            fig.savefig('{dire}/{name}.{form}'.format(dire=path_save, name=name_save, form=ff))
                


    def sum_conservative_efluence(self):
        fluence_sum_all = np.zeros(shape=(len(self.norms), len(self.lcindices)))
        fluence_sum_all_nominal = np.zeros(shape=(len(self.norms), len(self.lcindices)))
        efluences_nominal = {}
        conversion = -1
        logger.debug('Best indexes for calc fluence: ({0}, {1}'.format(self.args_max[0], self.args_max[1]))
        for grb in self.table_grb:
            logger.debug(grb['GRBNAME'])
            periods = self.load_values(grb['GRBNAME'], 'period')
            logger.debug('Loading time periods finished.')
            logger.debug('In the first period: {0} - {1} s'.format(periods[0][0], periods[0][1]))

            fluences = self.load_values(grb['GRBNAME'], 'efluence')
            logger.debug('Loading energy fluence values finished.')
            logger.debug('In the first period: {0}'.format(fluences[0][self.args_max[0]][self.args_max[1]]))

            eflux = self.load_values(grb['GRBNAME'], 'eflux')
            logger.debug('Loading energy flux values finished.')
            logger.debug('In the first period: {0}'.format(eflux[0][self.args_max[0]][self.args_max[1]]))

            prefactors = self.load_values(grb['GRBNAME'], 'normalization')
            logger.debug('Loading normalization values finished.')
            logger.debug('In the first period: {0}'.format(prefactors[0][self.args_max[0]][self.args_max[1]]))

            for n,e in zip(prefactors, eflux):
                if n[0][0]>0 and e[0][0]>0:
                    conversion = e[0][0]/n[0][0]
                    break
            if conversion==-1:
                logger.critical('Conversion factor has NOT been found!!!')
                sys.exit(1)
            logger.debug('Conversion factor from prefactor to efluence: {0}'.format(conversion))

            tnorm = 10.
            tstart = grb['GBM']['T90_START'] + grb['GBM']['T90']
            tstop = 100000.
            fluence_gbm = grb['GBM']['FLUENCE']
            prefactors_timeintegral = self.norms_mesh * tnorm / (self.lcindices_mesh+1) * (pow(tstop/tnorm, self.lcindices_mesh+1) - pow(tstart/tnorm, self.lcindices_mesh+1))

            #t_conversion = (self.lcindices_mesh!=-1) * pow((pow(periods[0][0], self.lcindices_mesh+1) + pow(periods[0][1], self.lcindices_mesh+1))/2., 1./((self.lcindices_mesh!=-1)*self.lcindices_mesh+1)) + (self.lcindices_mesh==-1) * np.exp((np.log(periods[0][1]) + np.log(periods[0][0]))/2. )
            logger.debug('Period: {0} - {1}'.format(periods[0][0], periods[-1][1]))

            #integral_nominal = prefactors[0] * t_conversion * ( ( (self.lcindices_mesh!=-1) / ((self.lcindices_mesh!=-1)*self.lcindices_mesh+1) * (pow(periods[-1][1]/t_conversion, self.lcindices_mesh+1) - pow(periods[0][0]/t_conversion, self.lcindices_mesh+1))) + ((self.lcindices_mesh==-1) *  (np.log(periods[-1][1]/t_conversion) - np.log(periods[0][0]/t_conversion)) ) )
            efluences_nominal[grb['GRBNAME']] = conversion *fluence_gbm * prefactors_timeintegral
            fluence_sum_all += np.sum(fluences, axis=0)
        self.efluences_sum = fluence_sum_all
        self.fluence_sum_best = fluence_sum_all[self.args_max[0]][self.args_max[1]]
        self.efluences_nominal = efluences_nominal
        self.efluence_sum_nominal = np.sum(efluences_nominal.values(), axis=0)
        self.efluence_sum_best_nominal = self.efluence_sum_nominal[self.args_max[0]][self.args_max[1]]

        logger.info('Energy fluence (conservative): {v:1.3E} erg/cm^2'.format(v=self.fluence_sum_best*MEVtoERG))
        logger.info('Energy fluence (nominal): {v:1.3E} erg/cm^2'.format(v=self.efluence_sum_best_nominal*MEVtoERG))

        return self.fluence_sum_best


    def sum_gbm_fluences(self):
        fluesum = 0
        for grb in self.table_grb:
            logger.debug(grb['GRBNAME'])
            fluesum += grb['GBM']['FLUENCE']
        return fluesum


    def plot(self, path_save, name_save, figforms=('png', 'pdf')):
        fig, ax = plt.subplots(3, 1, sharex=False, sharey=False, figsize=(5, 12))
        x, y = np.meshgrid(self.lcindices, self.norms)
        ax[0].grid()
        ax[0].set_xlabel('Temporal decay index')
        ax[0].set_ylabel(r'Normalization factor [a.u.]')

        #range_1sigma = find_range_shown(x=x, y=y, f=lambda u,v:self.dloglike_doubled[u][v]<=2.30)
        #logger.info('1 sigma range: {0}'.format(range_1sigma))
        range_shown = find_range_shown(x=self.lcindices_mesh, y=self.norms_mesh, f=lambda u,v:self.dloglike_doubled[u][v]<=max(self.cont_levels)*5)
                    
        ax[0].set_ylim((range_shown[0][0], range_shown[0][1]))
        ax[0].set_xlim((range_shown[1][0], range_shown[1][1]))
        cont = ax[0].contour(self.lcindices_mesh, self.norms_mesh, self.dloglike_doubled, levels=self.cont_levels, colors='black')
        cont.clabel(fmt='%1.2f', fontsize=12)
        logger.debug(self.dloglike_doubled[:, 1])
        ax[0].set_yscale('log')

        efluence_yrange_1sigma, efluence_xrange_1sigma = find_range_shown(x=self.lcindices_mesh, y=self.efluences_sum, f=lambda u,v:self.dloglike_doubled[u][v]<=2.30)
        logger.info('1 sigma range of temporal index: {0} - {1}'.format(efluence_xrange_1sigma[0], efluence_xrange_1sigma[1]))
        logger.info('1 sigma range of energy fluence: {0:1.2E} - {1:1.2E} erg/cm^2'.format(efluence_yrange_1sigma[0]*MEVtoERG, efluence_yrange_1sigma[1]*MEVtoERG))
        efluence_yrange_shown, efluence_xrange_shown = find_range_shown(x=self.lcindices_mesh, y=self.efluences_sum, f=lambda u,v:self.dloglike_doubled[u][v]<=max(self.cont_levels)*5)
        #logger.debug('{0}'.format(efluence_range_shown))
        cont_efluence = ax[1].contour(self.lcindices_mesh, self.efluences_sum*MEVtoERG, self.dloglike_doubled, levels=self.cont_levels, colors='black')
        cont_efluence_nominal = ax[1].contour(self.lcindices_mesh, self.efluence_sum_nominal*MEVtoERG, self.dloglike_doubled, levels=self.cont_levels, colors='gray')
        #cont_efluence = ax[1].contour(self.lcindices_mesh, self.efluence_sum_nominal*MEVtoERG, self.dloglike_doubled, levels=self.cont_levels, colors='black')
        cont_efluence.clabel(fmt='%1.2f', fontsize=10)
        ax[1].axhline(self.sum_gbm_fluences(), c='g', label='GBM (prompt)')
        ax[1].set_yscale('log')
        ax[1].set_xlabel('Temporal decay index')
        ax[1].set_ylabel(r'Stacked energy fluence $\mathrm{[erg/cm^{2}]}$')
        ax[1].grid()
        ax[1].set_xlim(efluence_xrange_shown[0], efluence_xrange_shown[1]) 
        #ax[1].set_ylim(3e-6, 10e-5) 
        ax[1].legend(loc=0)
        ax[1].set_ylim(efluence_yrange_shown[0]*MEVtoERG, efluence_yrange_shown[1]*MEVtoERG) 

        ax[2].grid()
        ax[2].grid(axis='both', which='major',color='black',linestyle='--')
        ax[2].grid(axis='y', which='minor',color='black',linestyle='-.', alpha=0.25)

        ax[2].set_xlabel('Temporal decay index')
        ax[2].set_ylabel(r'$\nu F_{\nu} / F_{GBM}\ \rm{[/s] \ at \ 10 \, s}$')

        vFv_at10s = self.get_e2dnde(10, ENERGY_RANGES[self.erange][0])
        yrange_shown_vFv, xrange_shown_vFv = find_range_shown(x=self.lcindices_mesh, y=vFv_at10s, f=lambda u,v:self.dloglike_doubled[u][v]<=max(self.cont_levels)*2)
                    
        ax[2].set_xlim((xrange_shown_vFv[0], xrange_shown_vFv[1]))
        ax[2].set_ylim((yrange_shown_vFv[0], yrange_shown_vFv[1]))
        cont_vFv = ax[2].contour(self.lcindices_mesh, vFv_at10s, self.dloglike_doubled, levels=self.cont_levels, colors='black')
        cont_vFv.clabel(fmt='%1.1f', fontsize=10)
        #if self.erange in ('high'):
        ax[2].set_yscale('log')

        vFv_yindexrange_1sigma, vFv_xindexrange_1sigma = find_indexrange(x=self.lcindices_mesh, y=vFv_at10s, f=lambda u,v:self.dloglike_doubled[u][v]<=2.30)
        #logger.debug(vFv_xindexrange_1sigma)

        # vFv_yrange_minAlpha, vFv_xrange_minAlpha = find_range_shown(x=self.lcindices_mesh, y=vFv_at10s, f=lambda u,v:self.dloglike_doubled[u][v]<=self.cont_levels[0], nx_restricted=[vFv_xindexrange_1sigma[0]])#[0]
        # logger.debug('E^2dN/dE with minimum alpha: {0}'.format(vFv_yrange_minAlpha))
        # logger.info('E^2dN/dE: {mi:1.2E} - {ma:1.2E} with minimum alpha={al:1.2f}'.format(mi=vFv_yrange_minAlpha[0], ma=vFv_yrange_minAlpha[1], al=self.lcindices[vFv_xindexrange_1sigma[0]]))
        # vFv_yrange_maxAlpha, vFv_xrange_maxAlpha = find_range_shown(x=self.lcindices_mesh, y=vFv_at10s, f=lambda u,v:self.dloglike_doubled[u][v]<=self.cont_levels[0], nx_restricted=[vFv_xindexrange_1sigma[1]])#[0]
        # logger.info('E^2dN/dE: {mi:1.2E} - {ma:1.2E} with maximum alpha={al:1.2f}'.format(mi=vFv_yrange_maxAlpha[0], ma=vFv_yrange_maxAlpha[1], al=self.lcindices[vFv_xindexrange_1sigma[1]]))

        fig.tight_layout() #subplots_adjust(hspace=0)
        for ff in figforms:
            fig.savefig('{dire}/{name}.{form}'.format(dire=path_save, name=name_save, form=ff))


    def eval_individual_norms(self):
        ngrb = len(self.table_grb)
        norms = {'best':np.ndarray((ngrb,)),
                 'poserr':np.ndarray((ngrb,)),
                 'negerr':np.ndarray((ngrb,)),
                 'ul':np.ndarray((ngrb,)),
                 'll':np.ndarray((ngrb,)),}
        for igrb, grb in enumerate(self.table_grb):
            logger.debug('Loading likelihood valus of {0}'.format(grb['GRBNAME']))
            loglikes = self.load_values(grb['GRBNAME'], 'loglike')
            loglike_tsum = sum(loglikes)
            loglike_bestlcindex = loglike_tsum[:, self.args_max[1]]
            loglike_bestlcindex_max = loglike_bestlcindex.max()
            dloglike_bestlcindex = - loglike_bestlcindex + loglike_bestlcindex_max
            dloglike_doubled = dloglike_bestlcindex * 2.
            #dloglike_doubled_smoothed = interp1d(self.norms, dloglike_doubled)
            logger.debug("""Doubled dloglikelihood: 
{0}""".format(dloglike_doubled))
            arg_max = np.where(loglike_bestlcindex==loglike_bestlcindex_max)[0]
            norms['best'][igrb] = self.norms[arg_max] #, self.args_max[1]]
            interpolated = interp1d(self.norms, dloglike_doubled-1.0)
            logger.debug('{0}: {1}'.format(self.norms[0], interpolated(self.norms[0])))
            logger.debug('{0}: {1}'.format(norms['best'][igrb], interpolated(norms['best'][igrb])))
            logger.debug('{0}: {1}'.format(self.norms[-1], interpolated(self.norms[-1])))
            roots_1sigma_hi = optimize.brentq(interpolated, norms['best'][igrb], self.norms[-1])
            if dloglike_doubled[0]>=1.0:
                roots_1sigma_lo = optimize.brentq(interpolated, self.norms[0], norms['best'][igrb])
                norms['poserr'][igrb] = roots_1sigma_hi - norms['best'][igrb] #np.max(roots_1sigma_hi) - norms['best'][igrb]
                norms['negerr'][igrb] = norms['best'][igrb] - roots_1sigma_lo #norms['best'][igrb] - np.min(roots_1sigma_lo)
            else:
                norms['best'][igrb] = norms['best'][igrb]+roots_1sigma_hi
                norms['poserr'][igrb] = 0
                norms['negerr'][igrb] = 0
            logger.debug('{0} + {1} - {2}'.format(norms['best'][igrb], norms['poserr'][igrb], norms['negerr'][igrb]))
            #roots_2sigma = optimize.fsolve(interp1d(self.norms, dloglike_doubled-4.0), norms['best'][igrb])
            #norms['ul'][igrb] = np.max(roots_2sigma)
            #norms['ll'][igrb] = np.min(roots_2sigma)
        self.norms_bestlcindex = norms


    def eval_individual_efluence(self):
        ngrb = len(self.table_grb)
        efluences = {'best':np.ndarray((ngrb,)),
                 'poserr':np.ndarray((ngrb,)),
                 'negerr':np.ndarray((ngrb,)),
                 'ul':np.ndarray((ngrb,)),
                 'll':np.ndarray((ngrb,)),}
        for igrb, grb in enumerate(self.table_grb):
            logger.debug('Loading likelihood valus of {0}'.format(grb['GRBNAME']))
            loglikes = self.load_values(grb['GRBNAME'], 'loglike')
            loglike_tsum = sum(loglikes)
            loglike_tsum_max = loglike_tsum.max()
            dloglike_tsum = - loglike_tsum + loglike_tsum_max
            dloglike_doubled = dloglike_tsum * 2.
            logger.debug("""Doubled dloglikelihood: 
{0}""".format(dloglike_doubled))
            args_max = zip(np.where(loglike_tsum==loglike_tsum_max)[0], np.where(loglike_tsum==loglike_tsum_max)[1])[0]
            logger.debug('Y best index: {0}'.format(args_max[0]))
            logger.debug('X best index: {0}'.format(args_max[1]))
            logger.debug('Best dloglikelihood: {0}'.format(dloglike_doubled[args_max[0]][args_max[1]]))
            efluences['best'][igrb] = self.efluences_nominal[grb['GRBNAME']][args_max[0]][args_max[1]]
            efluence_1sigma, lcindex_1sigma = find_range_shown(self.efluences_nominal[grb['GRBNAME']], self.lcindices_mesh, f=lambda u,v:dloglike_doubled[u][v]<=self.cont_levels[0])
            logger.info('Energy fluence best value: {0}'.format(efluences['best'][igrb]))
            logger.info('Energy fluence 1 sigma range: {0}'.format(efluence_1sigma))
            logger.info('Temporal index 1 sigma range: {0}'.format(lcindex_1sigma))
            efluence_2sigma, lcindex_2sigma = find_range_shown(self.efluences_nominal[grb['GRBNAME']], self.lcindices_mesh, f=lambda u,v:dloglike_doubled[u][v]<=self.cont_levels[1])
            if efluence_2sigma[0]>self.efluences_nominal[grb['GRBNAME']].min():
                efluences['poserr'][igrb] = efluence_1sigma[1] - efluences['best'][igrb]
                efluences['negerr'][igrb] = efluences['best'][igrb] - efluence_1sigma[0]
            else:
                efluences['best'][igrb] = efluence_2sigma[1]
                efluences['poserr'][igrb] = 0
                efluences['negerr'][igrb] = 0
            logger.info('{0} + {1} - {2}'.format(efluences['best'][igrb], efluences['poserr'][igrb], efluences['negerr'][igrb]))
        self.efluences_freelcindex = efluences


    def plot_individuals(self, path_save, name_save, figforms=('png', 'pdf')):
        flux_xrt11hr = []
        fluence_gbm = []
        igrbs_val = []
        igrbs_ul = []
        for igrb, grb in enumerate(self.table_grb):
            grb_xrt = ReadSwiftCatalogueInfo.read_one_row(self.table_xrt, grb['GCNNAME'])
            flux_xrt11hr.append(float(grb_xrt['XRT11HourFlux']))
            fluence_gbm.append(grb['GBM']['FLUENCE'])
            if self.efluences_freelcindex['negerr'][igrb]>0 and self.efluences_freelcindex['poserr'][igrb]>0:
                igrbs_val.append(igrb)
            else:
                igrbs_ul.append(igrb)
            if flux_xrt11hr[-1]/fluence_gbm[-1]>1E5:
                logger.info('GRB {0} XRT/GBM = {1}'.format(grb, flux_xrt11hr[-1]/fluence_gbm[-1]))

        logger.info('Value points: {0}'.format(igrbs_val))
        logger.info('UL points: {0}'.format(igrbs_ul))

        npa_flux_xrt11hr = np.array(flux_xrt11hr)
        npa_fluence_gbm = np.array(fluence_gbm)
        fig, ax = plt.subplots(1,2,sharex=False, sharey=False, figsize=(10, 5))
        
        ax[0].grid()
        ax[0].errorbar(flux_xrt11hr, npa_fluence_gbm*self.efluences_freelcindex['best'], xerr=0, yerr=[npa_fluence_gbm*self.efluences_freelcindex['negerr'], npa_fluence_gbm*self.efluences_freelcindex['poserr']], fmt='o', lw=1, markevery=igrbs_val)
        ax[0].plot(flux_xrt11hr, npa_fluence_gbm*(self.efluences_freelcindex['best']+self.efluences_freelcindex['poserr']), 'v', markevery=igrbs_ul)
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        ax[0].set_ylabel('LAT afterglow fluence')
        ax[0].set_xlabel('XRT flux at 11 hrs')

        ax[1].grid()
        ax[1].errorbar(npa_flux_xrt11hr/npa_fluence_gbm, self.efluences_freelcindex['best'], xerr=0, yerr=[self.efluences_freelcindex['negerr'], self.efluences_freelcindex['poserr']], fmt='o', lw=1, markevery=igrbs_val)
        ax[1].plot(flux_xrt11hr/npa_fluence_gbm, (self.efluences_freelcindex['best']+self.efluences_freelcindex['poserr']), 'v', markevery=igrbs_ul)
        ax[1].set_xscale('log')
        ax[1].set_yscale('log')
        ax[1].set_ylabel('LAT afterglow fluence / GBM prompt fluence')
        ax[1].set_xlabel('XRT flux at 11 hrs / GBM prompt fluence')

        for ff in figforms:
            fig.savefig('{dire}/{name}.{form}'.format(dire=path_save, name=name_save, form=ff))


    def plot_combined(self, ax, iplot=0, label=None):
        ax.contour(self.lcindices_mesh, self.efluences_sum*MEVtoERG, self.dloglike_doubled, levels=self.cont_levels[:2], colors=TPL_COLOR[iplot], ls=TPL_LINE[iplot], label=label)


@click.command()
@click.option('--low', '-l', type=str, default='E0000100-0001000MeV/r12deg')
@click.option('--mid', '-m', type=str, default='E0001000-0010000MeV/r03deg')
@click.option('--high', '-h', type=str, default='E0010000-0100000MeV/r01deg')
@click.option('--lowmid', type=str, default='E0000100-0010000MeV/r12deg')
@click.option('--midhigh', type=str, default='E0001000-0100000MeV/r03deg')
@click.option('--twindow', '-t', type=click.Choice(['afterglow', 'T95to01ks', '01ksto100ks', 'T95to03ks', '03ksto100ks']), default='afterglow')
@click.option('--namemin', type=str, default='000000000')
@click.option('--namemax', type=str, default='999999999')
@click.option('--normanchor', '-a', type=float, default=1.0)
@click.option('--combine', '-c', type=str, multiple=True, help='low, mid, high lowmid midhigh')
@click.option('--suffix', '-s', type=str, default='')
@click.option('--outdir', '-o', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/LongGRBs/Stacking/briefslots')
@click.option('--figform', type=str, default=('png','pdf'), multiple=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(low, mid, high, lowmid, midhigh, twindow, namemin, namemax, outdir, normanchor, combine, suffix, figform, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    ##### Combind plot ######
    fig_com = plt.figure()
    ax_com = fig_com.add_axes((0.2, 0.1, 0.75, 0.75))

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    dct_indir = {'low': low,
                 'mid': mid,
                 'high': high,
                 'lowmid':lowmid,
                 'midhigh':midhigh}

    ##### Suffix #####
    suffix_grb = ''
    if namemin==namemax:
        suffix_grb = 'GRB{0}'.format(namemin)
    elif namemin!='000000000' or namemax!='999999999':
        suffix_grb = 'GRB{0}-{1}'.format(namemin, namemax)
    suffix = suffix_grb + suffix

    ##### Catalogue information #####
    tb_lat = ReadLATCatalogueInfo.select_by_name(TABLE_LAT, namemin, namemax, TABLE_GBM)
    tb_lat = ReadLATCatalogueInfo.remove_by_name(tb_lat, excludes=('130427324', '150314205', '130702004'))
    tb_lat = ReadLATCatalogueInfo.select_gbm_exist(tb_lat)
    tb_lat = ReadLATCatalogueInfo.select_long(tb_lat)
    tb_lat = ReadLATCatalogueInfo.select_small_error(tb_lat, 0.3)
    tb_lat = ReadLATCatalogueInfo.select_zenith_angle(tb_lat, 100.)
    dframe_xrt = ReadSwiftCatalogueInfo.open_table()
    tb_lat = ReadLATCatalogueInfo.select_by_swift(tb_lat, dframe_xrt, ['XRT11HourFlux', 'XRT24HourFlux', 'XRTSpectralIndex'])

    ngrb = len(tb_lat)
    gbm_fluence_sum = 0
    for g in tb_lat:
        gbm_fluence_sum += g['GBM']['FLUENCE']
    gbm_fluence_ave = gbm_fluence_sum / ngrb

    for inum, enr in enumerate(combine):#('low', 'mid', 'high')):
        logger.info('===== {0} energy ====='.format(enr))
        indir = dct_indir[enr]
        pathout = '/'.join((outdir, indir))
        if not os.path.exists(pathout):
            os.makedirs(pathout)
        lognorm_anc = 2.5-2.*np.log10(ENERGY_RANGES[enr][0])
        norms = 10**np.linspace(lognorm_anc-3., lognorm_anc+3., 181)*normanchor
        lcindices = np.linspace(-0.655, -1.555, 91)
        model_ensemble = ModelEnsemble('{0} in {1} - {2} GeV'.format(twindow, ENERGY_RANGES[enr][0], ENERGY_RANGES[enr][1]), enr, twindow, indir, norms=norms, lcindices=lcindices, specindices=np.array([-2]), suffix=suffix, table_catalogue=tb_lat)
        model_ensemble.joint_likelihood()
        model_ensemble.delta_loglikelihood()
        model_ensemble.get_e2dnde(10, ENERGY_RANGES[enr][0])
        model_ensemble.sum_conservative_efluence()
        model_ensemble.plot(path_save=outdir, name_save='JointBriefSlotLikelihood{0}_{1}E_{2}'.format(model_ensemble.suffix, enr, twindow), figforms=figform)
        model_ensemble.eval_individual_efluence()
        model_ensemble.plot_individuals(path_save=outdir, name_save='IndividualLikelihood_XRT{0}_{1}E_{2}'.format(model_ensemble.suffix, enr, twindow), figforms=('png', 'pdf'))
        model_ensemble.plot_count(path_save=outdir, name_save='JointBriefSlotLikelihood_counts{0}_{1}E_{2}'.format(model_ensemble.suffix, enr, twindow), figforms=('png', 'pdf'))#, skip=('130427324', '150314205'))
        
        ax_com.contour(model_ensemble.lcindices_mesh, model_ensemble.efluence_sum_nominal*MEVtoERG, model_ensemble.dloglike_doubled, levels=model_ensemble.cont_levels[:2], colors=TPL_COLOR[inum], ls=TPL_LINE[inum], label='{0:.0f} - {1:.0f} GeV'.format(ENERGY_RANGES[enr][0], ENERGY_RANGES[enr][1]))
        gc.collect()
    ax_com.legend(loc=0, fancybox=True, framealpha=0.5)
    ax_com.set_xlim(-1.5, -0.7) 
    ax_com.set_ylim(1e-4, 0.5e-3) 
    ax_com.set_yscale('log')
    ax_com.set_xlabel('Temporal decay index')
    ax_com.set_ylabel(r'Stacked energy fluence $\mathrm{[erg/cm^{2}]}$')
    ax_com.grid()
    name_save = 'JointBriefSlotLikelihood_combined{0}_{1}'.format(suffix if suffix=='' else '_'+suffix, twindow)
    for ff in figform:
        path_plot = '{dire}/{name}.{form}'.format(dire=outdir, name=name_save, form=ff)
        fig_com.savefig(path_plot)
        logger.debug(path_plot)


if __name__ == '__main__':
    main()
