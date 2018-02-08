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


def likelihood_grb_analysis(name, mode, emin, emax, eref, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, sptmin, sptmax, redshift, calonly, outdir, binned, masifps, scalefactor=1., gti_external=None):
    if spectraltype=='PowerLaw':
        spectralpars = {'Prefactor':1e-10, 'Index':-2.0, 'Scale':eref}
    elif spectraltype=='PowerLaw2':
        spectralpars = {'Integral':1e-5, 'Index':-2.0, 'LowerLimit':emin, 'UpperLimit':emax}
    elif spectraltype=='ScaleFactor::PowerLaw2':
        spectralpars = {'Integral':1e-5, 'Index':-2.0, 'LowerLimit':emin, 'UpperLimit':emax, 'ScaleFactor':scalefactor}
    elif spectraltype=='EblAtten::PowerLaw2':
        spectralpars = {'Integral':1e-5, 'Index':-2.0, 'LowerLimit':emin, 'UpperLimit':emax, 'tau_norm':1., 'redshift':redshift, 'ebl_model':4}
    elif spectraltype=='ExpCutoff':
        spectralpars = {'Prefactor':1e-10, 'Index':-2.0, 'Scale':eref, 'Ebreak':10.0, 'P1':10000., 'P2':0, 'P3':0}
    elif spectraltype=='BrokenPowerLaw':
        spectralpars = {'Prefactor':1e-12, 'Index1':-2.0, 'Index2':-2.0, 'BreakValue':eref}
    elif spectraltype=='BrokenPowerLaw2':
        spectralpars = {'Integral':1e-5, 'Index1':-2.0, 'Index2':-2.0, 'BreakValue':eref, 'LowerLimit':emin, 'UpperLimit':emax}
    else:
        logger.critical("""{0} is NOT available!!! Use PowerLaw or ExpCutoff.""".format(spectraltype))
        sys.exit(1)

    grb = pLATLikelihoodConfig.GRBTarget(name, grbcatalogue, spectraltype=spectraltype, spectralpars=spectralpars)
    ana = pLATLikelihoodConfig.GRBConfig(grb, mode, emin=emin, emax=emax, deg_roi=roi, suffix=suffix, tmin_special=sptmin, tmax_special=sptmax, binned=binned, psForce=masifps, ft2interval='1s', gti_external=gti_external)

    if ana.tmin == ana.tmax:
        logger.warning('Time range of GRB config is NOT valid!')
        sys.exit(1)
    nrough = ana.setup(force={'download':False, 'filter':force, 'maketime':True, 'evtbin':force, 'livetime':force, 'exposure':force, 'model_3FGL_sources':True, 'diffuse_responses':True, 'srcmaps':force}, skip_zero_data=True)
    if ana.duration<=0:
        ana.dct_summary_results['TS'] = 0
        pickle_utilities.dump('{0}/Summary{1}.pickle'.format(ana.dir_work, '' if suffix=='' else '_'+suffix), ana.dct_summary_results)
        return ana.dct_summary_results
    #if nrough<1:
    #    ana.set_likelihood()
    #    logger.info("""Finishing the analysis because the number of events turned out to be zero.""")
        #return 0

    if calonly != tuple([None]*5):
        calonly_provided = True
    else:
        calonly_provided = False
    if calonly_provided==True:
        logger.info('CalOnly data are laoded...')
        logger.info(calonly)
        ana.add_calonly(path_onevt=calonly[0], path_onexp=calonly[1], path_bkgevt=calonly[2], nhistogram=calonly[3], rclass=calonly[4])
        if ana.calonly.ne_bins<1:
            calonly_provided = False            
    else:
        logger.info('CalOnly data are skipped.')

    if modelonly==True:
        esl = ana.set_likelihood()
        sys.exit(esl)
    fr = ana.fit(bredo=True)
    if fr!=1:
        ana.summarize_fit_results(use_calonly=calonly_provided)
        ana.plot_countspectra_fitted()
        ana.eval_flux_and_error()
        if spectraltype in ('PowerLaw', 'PowerLaw2', 'ScaleFactor::PowerLaw2', 'EblAtten::PowerLaw2'):
            ana.eval_limits_powerlaw(eref=eref, str_index_fixed=['best'])
            ana.set_likelihood_external_model(ana.path_model_xml_new)
            ana.scan_norm_and_index(eref=eref, use_calonly=calonly_provided)
            #ana.eval_limits_powerlaw_index()
    logger.info(ana.dct_summary_results)

    ##### Likelihood ratio ordering #####
    if nrough<10:
        ana.order_likeratio(eref=eref, use_calonly=calonly_provided)
    
    pickle_utilities.dump('{0}/Summary{1}.pickle'.format(ana.dir_work, '' if suffix=='' else '_'+suffix), ana.dct_summary_results)
    return ana.dct_summary_results


@click.command()
#@click.argument('name', type=str)
@click.option('--namemin', type=str, default='0')
@click.option('--namemax', type=str, default='200000000')
@click.option('--grbcatalogue', '-c', type=str, default=pLATLikelihoodConfig.GRB_CATALOGUE_LAT)
@click.option('--mode', '-m', type=click.Choice(['unified', 'prompt', 'primary', 'intermittent', 'afterglow', 'earlyAG', 'lateAG', 'farAG', 'T95to01ks', 'T95to03ks', '01ksto10ks', '01ksto100ks', 'T95to03ks', '03ksto10ks', '03ksto100ks', 'lightcurve', 'special']))
@click.option('--emin', type=float, default=100.)
@click.option('--emax', type=float, default=100000.)
@click.option('--eref', type=float, default=1000.)
@click.option('--roi', type=float, default=12.)
@click.option('--spectraltype', type=click.Choice(['PowerLaw', 'PowerLaw2', 'EblAtten::PowerLaw2', 'ScaleFactor::PowerLaw2', 'ExpCutoff', 'BrokenPowerLaw']), default='PowerLaw')
@click.option('--suffix', '-s', type=str, default='')
@click.option('--force', '-f', is_flag=True)
@click.option('--modelonly', is_flag=True)
@click.option('--masifps', is_flag=True)
@click.option('--refit', '-r', is_flag=True)
@click.option('--sptmin', type=float, default=-1, help='Only for mode special')
@click.option('--sptmax', type=float, default=-1, help='Only for mode special')
@click.option('--redshift', '-z', type=float, default=0.)
@click.option('--calonly', type=(str, str, str, int, str), default=None, help='path_onevt, path_onexp, path_bkgevt, nhtg, rclass')
@click.option('--binned', '-b', is_flag=True)
@click.option('--outdir', '-o', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone')
@click.option('--quanta', '-q', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/QuantiledGRBs_longonly3_GTI.pickle')
@click.option('--bsub', is_flag=True)
def main(namemin, namemax, mode, emin, emax, eref, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, sptmin, sptmax, redshift, calonly, outdir, binned, masifps, quanta, bsub):

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
            scalefactors = dct_quantile['fluence_scaled_gbm']['scaled'][namemin]
        else:
            scalefactors = 1
        # for idx, pha in itertools.product(range(len(lc_indices)), range(len(phases))):
        #     if lc_indices[idx]==-1 and phases[pha]==mode:
        #         logger.info('ScaleFactor = {0}'.format(scalefactors[idx][pha]))
        likelihood_grb_analysis(namemin, mode, emin, emax, eref, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, sptmin, sptmax, redshift, calonly, outdir, binned, masifps) #, scalefactor=scalefactors[idx][pha])
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
            acmd = ['bsub', '-o','{0}/{1}/GRB{1}_{2}{3}.log'.format(outdir, name, mode, suffix if suffix=='' else '_'+suffix), '-J','{0}{1}'.format(name[:6], mode[:2]), '-W','300', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/STLikelihoodAnalysis/STLikelihoodGRBAnalysis.py', '-m', mode, '--emin', str(emin), '--emax', str(emax), '--eref', str(eref), '-s', suffix, '--roi', str(roi), '--spectraltype', spectraltype, '--sptmin', str(sptmin), '--sptmax', str(sptmax), '--redshift', str(redshift), '--calonly', str(calonly), '--outdir', '{0}/{1}'.format(outdir, name), '--namemin', name, '--quanta', quanta]
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
