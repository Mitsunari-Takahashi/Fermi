#!/usr/bin/env python

import sys
import os
import os.path
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
import click
import itertools
import subprocess
import matplotlib as mpl
import matplotlib.pyplot as plt
import pLATLikelihoodConfig
import pickle_utilities
from STLikelihoodAnalysis import get_module_logger
import ReadLATCatalogueInfo
import ReadGBMCatalogueInfo
from STLikelihoodGRBAnalysis import likelihood_grb_analysis
from pCommon import MEVtoERG


##### Logger #####
logger = get_module_logger(__name__)
handler = StreamHandler()


##### Matplotlib #####
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
plt.rcParams["font.size"] = 15


def likelihood_grb_analysis_three_energies(name, mode, roi=None, spectraltype='PowerLaw', emax=100000., refit=True, force=False, suffix='', grbcatalogue=pLATLikelihoodConfig.GRB_CATALOGUE_LAT, modelonly=False, sptmin=None, sptmax=None, redshift=0., calonly=[None]*5, outdir='.', binned=False, masifps=False, scalefactor=1., gti_external=None):
    dct_emin = {'whole_energies':100, 'lower_energies':100, 'highest_energies':10000, 'middle_energies': 1000}
    dct_emax = {'whole_energies':emax, 'lower_energies':1000, 'highest_energies':emax, 'middle_energies': 10000}
    dct_eref = {'whole_energies':1000, 'lower_energies':300, 'highest_energies':30000, 'middle_energies': 3000}
    if roi is not None and roi>0:
        dct_roi = {'whole_energies':roi, 'lower_energies':roi, 'highest_energies':roi, 'middle_energies':roi}
    else:
        dct_roi = {'whole_energies':12., 'lower_energies':12., 'highest_energies':1., 'middle_energies':3.}
    dct_summary = {}
    for erange in ('highest_energies', 'whole_energies', 'lower_energies', 'middle_energies'):
        if dct_emin[erange]>=dct_emax[erange]:
            logger.critical('Emax={0} is NOT larger than Emin={1}!!!'.format(emax, emin))
            sys.exit(1)
        dct_summary[erange] = likelihood_grb_analysis(name=name, mode=mode, emin=dct_emin[erange], emax=dct_emax[erange], eref=dct_eref[erange], roi=dct_roi[erange], spectraltype=spectraltype, refit=refit, force=force, suffix=suffix, grbcatalogue=grbcatalogue, modelonly=modelonly, sptmin=sptmin, sptmax=sptmax, redshift=redshift, calonly=calonly, outdir=outdir, binned=binned, masifps=masifps, scalefactor=scalefactor, gti_external=gti_external)
        dct_summary[erange]['emin'] = dct_emin[erange]
        dct_summary[erange]['emax'] = dct_emax[erange]
        dct_summary[erange]['eref'] = dct_eref[erange]
    pickle_utilities.dump('{0}/{1}/three_eranges_{2}_{3}{4}.pickle'.format(outdir, name, mode, spectraltype, '' if suffix=='' else '_'+suffix), dct_summary)
    return dct_summary


def plot_spectrum(dict_info, name, mode, spectraltype='PowerLaw', suffix='', outdir='.'):
    """Make a plot of spectrum over the three energy bins.
Input a dictionary returned by likelihood_grb_analysis_three_energies.
"""

    fig = plt.figure(figsize=(15,7.5))
    ax = fig.add_axes((0.07, 0.15, 0.9, 0.75))
    flags_label = []
    for erange,dict_summary in dict_info.items():
        if ('likeratioordering' in dict_summary) and ('sed_bowtie' in dict_summary['likeratioordering']) and ('energies' in dict_summary['likeratioordering']['sed_bowtie']):
            seds = dict_summary['likeratioordering']['sed_bowtie']
        else:
            seds = dict_summary['dloglike']['sed_bowtie']
        for item_shown in seds['energies'].keys():
            if seds['energies'][item_shown] is not None:
                lim_markers = [ True if ie%3==0 and item_shown=='2sigma' else False for ie in range(len(seds['energies'][item_shown])) ]
                if erange in ('whole_energies'):
                    ax.fill_between(seds['energies'][item_shown], seds['curve_lo'][item_shown]*MEVtoERG, seds['curve_hi'][item_shown]*MEVtoERG, alpha=0.8 if item_shown is '1sigma' else 0.3, label=item_shown, facecolor='lime' if item_shown=='1sigma' else 'lightsage')
                else:
                    #ax.errorbar(seds['energies'][item_shown], seds['curve_lo'][item_shown]*MEVtoERG, xerr=0, yerr=[[0]*len(lim_markers), seds['curve_lo'][item_shown]*MEVtoERG*lim_markers], lolims=lim_markers, uplims=False, ls='-', fmt='', c='k', lw=1, marker='None')
                    ax.plot(seds['energies'][item_shown], seds['curve_lo'][item_shown]*MEVtoERG, '-^', c='k' if item_shown=='1sigma' else 'darkred', lw=1, markevery=3, alpha=1.0)
                    if not item_shown in flags_label:
                        #ax.errorbar(seds['energies'][item_shown], seds['curve_hi'][item_shown]*MEVtoERG, xerr=0, yerr=[seds['curve_hi'][item_shown]*MEVtoERG*lim_markers*0.5,[0]*len(lim_markers)],uplims=lim_markers, lolims=False, ls='-', fmt='', label=item_shown, c='k', lw=1, marker='None')
                        ax.plot(seds['energies'][item_shown], seds['curve_hi'][item_shown]*MEVtoERG, '-v', c='k' if item_shown=='1sigma' else 'darkred', lw=1, label=item_shown, markevery=3, alpha=1.0)
                        flags_label.append(item_shown)
                    else:
                        #ax.errorbar(seds['energies'][item_shown], seds['curve_hi'][item_shown]*MEVtoERG, xerr=0, yerr=[seds['curve_hi'][item_shown]*MEVtoERG*lim_markers*0.5,[0]*len(lim_markers)],uplims=lim_markers, lolims=False, ls='-', fmt=',', c='k', lw=1, marker='None')
                        ax.plot(seds['energies'][item_shown], seds['curve_hi'][item_shown]*MEVtoERG, '-v', c='k' if item_shown=='1sigma' else 'darkred', lw=1, markevery=3, alpha=1.0)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim((1E-15, 1E-7))
        ax.set_xlabel('Energy [MeV]')
        ax.set_ylabel(r'$\nu F_{\nu} \ \rm{[erg/cm^2 \cdot s]}$')
        ax.grid(color='k', linestyle='--', linewidth=1, alpha=0.5)
        ax.legend(loc=0, fontsize=12, fancybox=True, framealpha=0.5)

    for ff in ('pdf', 'png'):
        path_save = '{0}/{1}/sed_bowtie_combined_{2}_{3}{4}.{5}'.format(outdir, name, mode, spectraltype, '' if suffix=='' or suffix[0]=='_' else '_'+suffix, ff)#"{0}/{1}_sed_bowtie_{2}{3}.{4}".format(self.dir_work, name, self.target.name, self.suffix, ff)
        fig.savefig(path_save)
        logger.info('{0} has been saved.'.format(path_save))


@click.command()
@click.option('--namemin', type=str, default='0')
@click.option('--namemax', type=str, default='200000000')
@click.option('--grbcatalogue', '-c', type=str, default=pLATLikelihoodConfig.GRB_CATALOGUE_LAT)
@click.option('--mode', '-m', type=click.Choice(['unified', 'prompt', 'primary', 'intermittent', 'afterglow', 'earlyAG', 'lateAG', 'farAG', 'T95to01ks', 'T95to03ks', '01ksto10ks', '01ksto100ks', 'T95to03ks', '03ksto10ks', '03ksto100ks', 'lightcurve', 'special']))
@click.option('--roi', type=float, default=0)
@click.option('--spectraltype', type=click.Choice(['PowerLaw', 'PowerLaw2', 'EblAtten::PowerLaw2', 'ScaleFactor::PowerLaw2', 'ExpCutoff', 'BrokenPowerLaw']), default='PowerLaw')
@click.option('--emax', type=float, default=100000.)
@click.option('--suffix', '-s', type=str, default='')
@click.option('--force', '-f', is_flag=True)
@click.option('--modelonly', is_flag=True)
@click.option('--sptmin', type=float, default=None)
@click.option('--sptmax', type=float, default=None)
@click.option('--redshift', '-z', type=float, default=0.)
@click.option('--calonly', type=(str, str, str, int, str), default=[None]*5, help='path_onevt, path_onexp, path_offevt, path_offexp, nhtg, rclass')
@click.option('--masifps', is_flag=True)
@click.option('--refit', '-r', is_flag=True)
@click.option('--binned', '-b', is_flag=True)
@click.option('--plotonly', '-p', type=str, default=None)
@click.option('--outdir', '-o', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone')
@click.option('--quanta', '-q', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/QuantiledGRBs_longonly3_GTI.pickle')
@click.option('--bsub', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(namemin, namemax, mode, roi, spectraltype, emax, refit, force, suffix, grbcatalogue, modelonly, sptmin, sptmax, redshift, calonly, outdir, binned, plotonly, masifps, quanta, bsub, loglevel):
    if mode=='special':
        suffix  = 'T{0:0>6.0f}-{1:0>6.0f}s{2}'.format(sptmin, sptmax, suffix if suffix is '' or suffix[0]=='_' else '_'+suffix)
    ##### Logger ######
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    tb_lat = ReadLATCatalogueInfo.open_table()
    tb_gbm = ReadGBMCatalogueInfo.open_table()
    print 'Loading {0}...'.format(quanta)
    dct_quantile = pickle_utilities.load(quanta)
    lc_indices = dct_quantile['fluence_scaled_gbm']['indices']
    phases = dct_quantile['fluence_scaled_gbm']['phases']

    if bsub==False:
        tb = ReadLATCatalogueInfo.select_one_by_name(tb_lat, namemin, tb_gbm)
        if len(tb)<1:
            logger.error('GRB {0} has NOT been found.'.format(namemin))
            sys.exit(1)
        scalefactor = 1
        if spectraltype[:13]=='ScaleFactor::':
            scalefactor = tb['GBM']['FLUENCE'] #scalefactors = dct_quantile['fluence_scaled_gbm']['scaled'][namemin]
        if plotonly is None:
            results = likelihood_grb_analysis_three_energies(namemin, mode, roi, spectraltype, emax, refit, force, suffix, grbcatalogue, modelonly, sptmin, sptmax, redshift, calonly, outdir, binned, masifps, scalefactor=scalefactor)
        else:
            results = pickle_utilities.load(plotonly)
        plot_spectrum(results, namemin, mode, spectraltype, suffix, outdir)
    else:
        tb = ReadLATCatalogueInfo.select_by_name(tb_lat, namemin, namemax, tb_gbm)
        tb = ReadLATCatalogueInfo.select_gbm_exist(tb)
        if len(tb)<1:
            print 'No GRBs.'
            return 1
        for (irow , row) in enumerate(tb):
            name = row['GRBNAME']
            print '##### No.{0} GRB{1} #####'.format(irow, name)
            if not os.path.exists(name):
                os.mkdir(name)
            acmd = ['bsub', '-o','{0}/{1}/GRB{1}_{2}{3}.log'.format(outdir, name, mode, suffix if suffix=='' or suffix[0]=='_' else '_'+suffix), '-J','{0}{1}'.format(name[:6], mode[:2]), '-W','600', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/STLikelihoodAnalysis/STLikelihoodGRBAnalysis_ThreeEnergies.py', '-m', mode, '-s', suffix, '--roi', str(roi), '--spectraltype', spectraltype, '--emax', str(emax), '--sptmin', str(sptmin), '--sptmax', str(sptmax), '--redshift', str(redshift), '--calonly', str(calonly), '--outdir', '{0}'.format(outdir), '--namemin', name] #, '--quanta', quanta]
            if force==True:
                acmd.append('--force')
            if refit==True:
                acmd.append('--refit')
            if modelonly==True:
                acmd.append('--modelonly')
            if masifps==True:
                acmd.append('--masifps')
            if binned==True:
                acmd.append('--binned')
            if calonly != tuple([None]*5):
                acmd.append('--calonly')
                for i in range(5):
                    acmd.append(str(calonly[i]))
            if plotonly is not None:
                acmd.append('--plotonly {0}'.format(plotonly))
            print acmd
            subprocess.call(acmd)
        return 0



if __name__ == '__main__':
    main()
