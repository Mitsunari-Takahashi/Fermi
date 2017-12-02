#!/usr/bin/env python

import sys
import os
import os.path
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
import click
import itertools
import subprocess
import pLATLikelihoodConfig
import pickle_utilities
from STLikelihoodAnalysis import get_module_logger
import ReadLATCatalogueInfo
import ReadGBMCatalogueInfo
from STLikelihoodGRBAnalysis import likelihood_grb_analysis

##### Logger #####
logger = get_module_logger(__name__)
handler = StreamHandler()


def likelihood_grb_analysis_three_energies(name, mode, roi=None, spectraltype='PowerLaw', refit=True, force=False, suffix='', grbcatalogue=pLATLikelihoodConfig.GRB_CATALOGUE_LAT, modelonly=False, outdir='.', binned=False, masifps=False, scalefactor=1., gti_external=None):
    dct_emin = {'whole_energies':100, 'lower_energies':100, 'highest_energies':10000, 'middle_energies': 1000}
    dct_emax = {'whole_energies':100000, 'lower_energies':1000, 'highest_energies':100000, 'middle_energies': 10000}
    #dct_emin = {'whole_energies':100, 'lower_energies':100, 'highest_energies':10000}
    #dct_emax = {'whole_energies':100000, 'lower_energies':10**3.75, 'highest_energies':100000}
    if roi is not None and roi>0:
        dct_roi = {'whole_energies':roi, 'lower_energies':roi, 'highest_energies':roi, 'middle_energies':roi}
    else:
        dct_roi = {'whole_energies':12., 'lower_energies':12., 'highest_energies':1., 'middle_energies':3.}
    dct_summary = {}
    for erange in ('whole_energies', 'lower_energies', 'highest_energies', 'middle_energies'):
        # if erange in ('lower_energies', 'highest_energies', 'middle_energies'):
        #     gti_external = 'E{emin:0>7.0f}-{emax:0>7.0f}MeV/r{roi:0>2.0f}deg'.format(emin=dct_emin['whole_energies'], emax=dct_emax['whole_energies'], roi=dct_roi['whole_energies'])
        # else:
        #gti_external = None
        dct_summary[erange] = likelihood_grb_analysis(name=name, mode=mode, emin=dct_emin[erange], emax=dct_emax[erange], roi=dct_roi[erange], spectraltype=spectraltype, refit=refit, force=force, suffix=suffix, grbcatalogue=grbcatalogue, modelonly=modelonly, outdir=outdir, binned=binned, masifps=masifps, scalefactor=scalefactor, gti_external=gti_external)
        dct_summary[erange]['emin'] = dct_emin[erange]
        dct_summary[erange]['emax'] = dct_emax[erange]
    pickle_utilities.dump('{0}/{1}/three_eranges_{2}_{3}{4}.pickle'.format(outdir, name, mode, spectraltype, '' if suffix=='' else '_'+suffix), dct_summary)


@click.command()
@click.option('--namemin', type=str, default='0')
@click.option('--namemax', type=str, default='200000000')
@click.option('--grbcatalogue', '-c', type=str, default=pLATLikelihoodConfig.GRB_CATALOGUE_LAT)
@click.option('--mode', '-m', type=click.Choice(['unified', 'prompt', 'primary', 'intermittent', 'afterglow', 'earlyAG', 'lateAG', 'farAG', 'T95to01ks', 'T95to03ks', '01ksto10ks', '01ksto100ks', 'T95to03ks', '03ksto10ks', '03ksto100ks', 'lightcurve', 'special']))
@click.option('--roi', type=float, default=0)
@click.option('--spectraltype', type=click.Choice(['PowerLaw', 'PowerLaw2', 'ScaleFactor::PowerLaw2', 'ExpCutoff', 'BrokenPowerLaw']), default='PowerLaw')
@click.option('--suffix', '-s', type=str, default='')
@click.option('--force', '-f', is_flag=True)
@click.option('--modelonly', is_flag=True)
@click.option('--masifps', is_flag=True)
@click.option('--refit', '-r', is_flag=True)
@click.option('--binned', '-b', is_flag=True)
@click.option('--outdir', '-o', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone')
@click.option('--quanta', '-q', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone/QuantiledGRBs_longonly3_GTI.pickle')
@click.option('--bsub', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(namemin, namemax, mode, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, binned, masifps, quanta, bsub, loglevel):
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
        likelihood_grb_analysis_three_energies(namemin, mode, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, binned, masifps, scalefactor=scalefactor)
        #for idx, pha in itertools.product(range(len(lc_indices)), range(len(phases))):
            #print lc_indices[idx], phases[pha]
         #   if lc_indices[idx]==-1 and phases[pha]==mode:
          #      logger.info('ScaleFactor = {0}'.format(scalefactors[idx][pha]))
           #     likelihood_grb_analysis_three_energies(namemin, mode, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, binned, masifps, scalefactor=scalefactors[idx][pha])
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
            acmd = ['bsub', '-o','{0}/{1}/GRB{1}_{2}{3}.log'.format(outdir, name, mode, suffix if suffix=='' else '_'+suffix), '-J','{0}{1}'.format(name[:6], mode[:2]), '-W','600', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/STLikelihoodAnalysis/STLikelihoodGRBAnalysis_ThreeEnergies.py', '-m', mode, '-s', suffix, '--roi', str(roi), '--spectraltype', spectraltype, '--outdir', '{0}'.format(outdir), '--namemin', name, '--quanta', quanta]
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
