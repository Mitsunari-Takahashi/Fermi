#!/usr/bin/env python

import sys
import os
import os.path
import click
import subprocess
import pLATLikelihoodConfig
from STLikelihoodAnalysis import get_module_logger
import ReadLTFCatalogueInfo

##### Logger #####
logger = get_module_logger(__name__)


def likelihood_grb_analysis(name, mode, emin, emax, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, binned, masifps):
    if spectraltype=='PowerLaw':
        spectralpars = {'Prefactor':1e-10, 'Index':-2.0, 'Scale':1000.}
    elif spectraltype=='PowerLaw2':
        spectralpars = {'Integral':1e-5, 'Index':-2.0, 'LowerLimit':emin, 'UpperLimit':emax}
    elif spectraltype=='ExpCutoff':
        spectralpars = {'Prefactor':1e-10, 'Index':-2.0, 'Scale':1000., 'Ebreak':10.0, 'P1':10000., 'P2':0, 'P3':0}
    elif spectraltype=='BrokenPowerLaw':
        spectralpars = {'Prefactor':1e-12, 'Index1':-2.0, 'Index2':-2.0, 'BreakValue':10000.}
    else:
        logger.critical("""{0} is NOT available!!! Use PowerLaw or ExpCutoff.""".format(spectraltype))
        sys.exit(1)

    grb = pLATLikelihoodConfig.GRBTarget(name, grbcatalogue, spectraltype=spectraltype, spectralpars=spectralpars)
    ana = pLATLikelihoodConfig.GRBConfig(grb, mode, emin=emin, emax=emax, deg_roi=roi, binned=binned, psForce=masifps)
    if ana.tmin == ana.tmax:
        logger.warning('Time range of GRB config is NOT valid!')
        sys.exit(1)
    nrough = ana.setup(force={'download':False, 'filter':force, 'maketime':force, 'evtbin':force, 'livetime':force, 'exposure':force, 'model_3FGL_sources':True, 'diffuse_responses':force, 'srcmaps':force}, skip_zero_data=False)
    #if nrough<1:
    #    ana.set_likelihood()
    #    logger.info("""Finishing the analysis because the number of events turned out to be zero.""")
        #return 0
    if modelonly==True:
        esl = ana.set_likelihood()
        sys.exit(esl)
    ana.fit(bredo=True)
    ana.plot_countspectra_fitted()
    ana.eval_flux_and_error()
    if spectraltype in ('PowerLaw', 'PowerLaw2'):
        ana.eval_limits_powerlaw(str_index_fixed=['best'])
    logger.info(ana.dct_summary_results)
    return ana.dct_summary_results


@click.command()
#@click.argument('name', type=str)
@click.option('--namemin', type=str, default='0')
@click.option('--namemax', type=str, default='200000000')
@click.option('--grbcatalogue', '-c', type=str, default=pLATLikelihoodConfig.GRB_CATALOGUE_LTF)
@click.option('--mode', '-m', type=click.Choice(['unified', 'prompt', 'primary', 'intermittent', 'afterglow', 'earlyAG', 'lateAG', 'farAG', 'T95to01ks', '01ksto10ks', 'T95to03ks', '03ksto10ks', 'lightcurve', 'special']))
@click.option('--emin', type=float, default=100.)
@click.option('--emax', type=float, default=100000.)
@click.option('--roi', type=float, default=12.)
@click.option('--spectraltype', type=click.Choice(['PowerLaw', 'PowerLaw2', 'ExpCutoff', 'BrokenPowerLaw']), default='PowerLaw')
@click.option('--suffix', '-s', type=str, default='')
@click.option('--force', '-f', is_flag=True)
@click.option('--modelonly', is_flag=True)
@click.option('--masifps', is_flag=True)
@click.option('--refit', '-r', is_flag=True)
@click.option('--binned', '-b', is_flag=True)
@click.option('--outdir', '-o', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone')
@click.option('--bsub', is_flag=True)
def main(namemin, namemax, mode, emin, emax, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, binned, masifps, bsub):

    if bsub==False:
        likelihood_grb_analysis(namemin, mode, emin, emax, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, binned, masifps)
    else:
        tb_ltf = ReadLTFCatalogueInfo.open_table()
        tb_gbm = ReadLTFCatalogueInfo.select_gbm_exist(tb_ltf)
        tb = ReadLTFCatalogueInfo.select_by_name(tb_gbm, namemin, namemax)
        if len(tb)<1:
            print 'No GRBs.'
            return 1
        for (irow , row) in enumerate(tb):
            name = row['GRBNAME']
            print '##### No.{0} GRB{1} #####'.format(irow, name)
            if not os.path.exists(name):
                os.mkdir(name)
            acmd = ['bsub', '-o','{0}/{1}/GRB{1}_{2}{3}.log'.format(outdir, name, mode, suffix if suffix=='' else '_'+suffix), '-J','{0}{1}'.format(name[:6], mode[:2]), '-W','300', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/STLikelihoodAnalysis/STLikelihoodGRBAnalysis.py', '-m', mode, '--emin', str(emin), '--emax', str(emax), '-s', suffix, '--roi', str(roi), '--spectraltype', spectraltype, '--outdir', '{0}/{1}'.format(outdir, name), '--namemin', name]
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
            print acmd
            subprocess.call(acmd)
        return 0



if __name__ == '__main__':
    main()
