#!/usr/bin/env python

import sys
import os
import os.path
import click
import subprocess
import pLATLikelihoodConfig
import pickle_utilities
from STLikelihoodAnalysis import get_module_logger
import ReadLTFCatalogueInfo
from STLikelihoodGRBAnalysis import likelihood_grb_analysis

##### Logger #####
logger = get_module_logger(__name__)


def likelihood_grb_analysis_three_energies(name, mode, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, binned, masifps):
    dct_emin = {'whole_energies':100, 'lower_energies':100, 'highest_energies':10000}
    dct_emax = {'whole_energies':100000, 'lower_energies':10**3.75, 'highest_energies':100000}
    dct_summary = {}
    for erange in ('whole_energies', 'lower_energies', 'highest_energies'):
        dct_summary[erange] = likelihood_grb_analysis(name, mode, dct_emin[erange], dct_emax[erange], roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, binned, masifps)
        dct_summary[erange]['emin'] = dct_emin[erange]
        dct_summary[erange]['emax'] = dct_emax[erange]
    pickle_utilities.dump('{0}/{1}/three_eranges_{2}{3}.pickle'.format(outdir, name, mode, suffix), dct_summary)


@click.command()
#@click.argument('name', type=str)
@click.option('--namemin', type=str, default='0')
@click.option('--namemax', type=str, default='200000000')
@click.option('--grbcatalogue', '-c', type=str, default=pLATLikelihoodConfig.GRB_CATALOGUE_LTF)
@click.option('--mode', '-m', type=click.Choice(['unified', 'prompt', 'primary', 'intermittent', 'afterglow', 'earlyAG', 'lateAG', 'farAG', 'T95to01ks', '01ksto10ks', 'T95to03ks', '03ksto10ks', 'lightcurve', 'special']))
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
def main(namemin, namemax, mode, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, binned, masifps, bsub):

    if bsub==False:
        likelihood_grb_analysis_three_energies(namemin, mode, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, binned, masifps)
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
            acmd = ['bsub', '-o','{0}/{1}/GRB{1}_{2}{3}.log'.format(outdir, name, mode, suffix if suffix=='' else '_'+suffix), '-J','{0}{1}'.format(name[:6], mode[:2]), '-W','300', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/STLikelihoodAnalysis/STLikelihoodGRBAnalysis_ThreeEnergies.py', '-m', mode, '-s', suffix, '--roi', str(roi), '--spectraltype', spectraltype, '--outdir', '{0}'.format(outdir), '--namemin', name]
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
