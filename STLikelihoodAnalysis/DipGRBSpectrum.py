#!/usr/bin/env python

import sys
import os
import os.path
path_upstairs = os.path.join(os.path.dirname(__file__), '../')
import click
import subprocess
import itertools
import pLATLikelihoodConfig
from STLikelihoodAnalysis import get_module_logger
import ReadGBMCatalogueInfo
import ReadLATCatalogueInfo
import pickle_utilities

##### Logger #####
logger = get_module_logger(__name__)


def likelihood_grb_dip(name, mode, emin, emax, roi, refit, force, suffix, grbcatalogue, sptmin, sptmax, outdir, binned, masifps, gti_external=None):

    dct_summary_functions = {}

    logger.info('Fitting by PowerLaw2...')
    grb = pLATLikelihoodConfig.GRBTarget(name, grbcatalogue, spectraltype='PowerLaw2', spectralpars={'Integral':1e-5, 'Index':-1.1, 'LowerLimit':emin, 'UpperLimit':emax})
    ana = pLATLikelihoodConfig.GRBConfig(grb, mode, emin=emin, emax=emax, deg_roi=roi, tmin_special=sptmin, tmax_special=sptmax, binned=binned, psForce=masifps, ft2interval='1s', gti_external=gti_external)
    ana.setup(force={'download':False, 'filter':False, 'maketime':True, 'evtbin':False, 'livetime':False, 'exposure':False, 'model_3FGL_sources':True, 'diffuse_responses':True, 'srcmaps':True})
    nfreepars_pl = len(grb.spectralpars) - len(grb.spectralpars_fixed)
    fr = ana.fit(bredo=True)
    loglike_inv_pl = ana.loglike_inversed
    if fr!=1:
        ana.summarize_fit_results()
        ana.plot_countspectra_fitted()
    logger.info(ana.dct_summary_results)
    dct_summary_functions['PowerLaw2'] = ana.dct_summary_results
    logger.info('TS = {0}'.format(dct_summary_functions['PowerLaw2']['TS']))

    logger.info('Fitting by BrokenPowerLaw2...')
    #betas = np.linspace(-4, 0, 41)
    #deltas = np.linspace(0, 2, 20)
    #gamma1_idx = self.like.par_index(self.target.name, 'Index1')
    #gamma2_idx = self.like.par_index(self.target.name, 'Index2')

    grb = pLATLikelihoodConfig.GRBTarget(name, grbcatalogue, spectraltype='BrokenPowerLaw2', spectralpars={'Integral':max(1e-7,dct_summary_functions['PowerLaw2']['Integral']['value']), 'Index1':-2.2, 'Index2':-2./3., 'BreakValue':5000, 'LowerLimit':emin, 'UpperLimit':emax}, spectralpars_fixed=['LowerLimit', 'UpperLimit'])
    ana = pLATLikelihoodConfig.GRBConfig(grb, mode, emin=emin, emax=emax, deg_roi=roi, tmin_special=sptmin, tmax_special=sptmax, binned=binned, psForce=masifps, ft2interval='1s', gti_external=gti_external)
    ana.setup(force={'download':False, 'filter':False, 'maketime':True, 'evtbin':False, 'livetime':False, 'exposure':False, 'model_3FGL_sources':False, 'diffuse_responses':False, 'srcmaps':False})
    nfreepars_bpl = len(grb.spectralpars) - len(grb.spectralpars_fixed)

    #for b, d in zip(betas, deltas):
    fr = ana.fit(bredo=True)
    loglike_inv_bpl = ana.loglike_inversed
    if fr!=1:
        ana.summarize_fit_results()
        ana.plot_countspectra_fitted()
    logger.info(ana.dct_summary_results)
    dct_summary_functions['BrokenPowerLaw2'] = ana.dct_summary_results
    logger.info('TS = {0}'.format(dct_summary_functions['BrokenPowerLaw2']['TS']))
    
    ts_break = 2.*(loglike_inv_bpl - loglike_inv_pl)
    logger.info('TS of spectral break: {0}'.format(ts_break))
    logger.info('Delta-digrees of freedom: {0}'.format(nfreepars_bpl-nfreepars_pl))

    pickle_utilities.dump('{0}/SummaryDip{1}.pickle'.format(ana.dir_work, '' if suffix=='' else '_'+suffix), dct_summary_functions)
    return dct_summary_functions


@click.command()
@click.option('--namemin', type=str, default='0')
@click.option('--namemax', type=str, default='200000000')
@click.option('--grbcatalogue', '-c', type=str, default=pLATLikelihoodConfig.GRB_CATALOGUE_LAT)
@click.option('--mode', '-m', type=click.Choice(['unified', 'prompt', 'primary', 'intermittent', 'afterglow', 'earlyAG', 'lateAG', 'farAG', 'T95to01ks', 'T95to03ks', '01ksto10ks', '01ksto100ks', 'T95to03ks', '03ksto10ks', '03ksto100ks', 'lightcurve', 'special']))
@click.option('--emin', type=float, default=100.)
@click.option('--emax', type=float, default=100000.)
@click.option('--roi', type=float, default=12.)
@click.option('--suffix', '-s', type=str, default='')
@click.option('--force', '-f', is_flag=True)
@click.option('--masifps', is_flag=True)
@click.option('--refit', '-r', is_flag=True)
@click.option('--sptmin', type=float, default=-1, help='Only for mode special')
@click.option('--sptmax', type=float, default=-1, help='Only for mode special')
@click.option('--binned', '-b', is_flag=True)
@click.option('--outdir', '-o', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone')
@click.option('--bsub', is_flag=True)
def main(namemin, namemax, mode, emin, emax, roi, refit, force, suffix, grbcatalogue, sptmin, sptmax, outdir, binned, masifps, bsub):

    tb_lat = ReadLATCatalogueInfo.open_table()
    tb_gbm = ReadGBMCatalogueInfo.open_table()

    if bsub==False:
        tb = ReadLATCatalogueInfo.select_one_by_name(tb_lat, namemin, tb_gbm)
        if len(tb)<1:
            logger.error('GRB {0} has NOT been found.'.format(namemin))
            sys.exit(1)
        likelihood_grb_dip(namemin, mode, emin, emax, roi, refit, force, suffix, grbcatalogue, sptmin, sptmax, outdir, binned, masifps)
    else:
        tb = ReadLATCatalogueInfo.select_by_name(tb_lat, namemin, namemax, tb_gbm)
        tb = ReadLATCatalogueInfo.select_gbm_exist(tb)
        if len(tb)<1:
            logger.error('No GRBs.')
            sys.exit(1)
        for (irow , row) in enumerate(tb):
            name = row['GRBNAME']
            print '##### No.{0} GRB{1} #####'.format(irow, name)
            if not os.path.exists(name):
                os.mkdir(name)
            acmd = ['bsub', '-o','{0}/{1}/GRB{1}_Dip_{2}{3}.log'.format(outdir, name, mode, suffix if suffix=='' else '_'+suffix), '-J','{0}{1}'.format(name[:6], mode[:2]), '-W','300', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/STLikelihoodAnalysis/STLikelihoodGRBAnalysis.py', '-m', mode, '--emin', str(emin), '--emax', str(emax), '-s', suffix, '--roi', str(roi), '--sptmin', str(sptmin), '--sptmax', str(sptmax), '--outdir', '{0}/{1}'.format(outdir, name), '--namemin', name]
            if force==True:
                acmd.append('--force')
            if refit==True:
                acmd.append('--refit')
            if masifps==True:
                acmd.append('--masifps')
            if binned==True:
                acmd.append('--binned')
            print acmd
            subprocess.call(acmd)
        return 0


if __name__ == '__main__':
    main()
