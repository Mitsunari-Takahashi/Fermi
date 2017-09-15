#!/usr/bin/env python
"""Submit jobs of AnalyzeGRB_fermipy.py to SLAC server.
"""
import sys
import os
path_upstairs = os.path.join(os.path.dirname(__file__), '../')
sys.path.append(path_upstairs)
import os.path
import subprocess
import click
import csv
import ReadLTFCatalogueInfo


@click.command()
@click.option('--namemin', type=str, default='0')
@click.option('--namemax', type=str, default='200000000')
@click.option('--emin', type=float, default=100.0)
@click.option('--emax', type=float, default=100000.)
@click.option('--mode', '-m', type=click.Choice(['unified', 'prompt', 'afterglow', 'earlyAG', 'lateAG', 'lightcurve', 'special']))
@click.option('--roi', type=float, default=12.)
@click.option('--spectraltype', type=click.Choice(['PowerLaw', 'ExpCutoff']))
@click.option('--refit', '-r', is_flag=True)
@click.option('--force', '-f', is_flag=True)
@click.option('--directory', '-d', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone')
@click.option('--suffix', '-s', type=str, default='')
@click.option('--modelonly', is_flag=True)
@click.option('--masifps', is_flag=True)
def main(namemin, namemax, emin, emax, mode, roi, spectraltype, directory, suffix, refit, force, modelonly, masifps):
    tb_ltf = ReadLTFCatalogueInfo.open_table()
    tb_gbm = ReadLTFCatalogueInfo.select_gbm_exist(tb_ltf)
    tb_selected = ReadLTFCatalogueInfo.select_by_name(tb_gbm, namemin, namemax)
    tb = tb_selected
    if len(tb)<1:
        print 'No GRBs.'
        return 1
    for (irow , row) in enumerate(tb):
        name = row['GRBNAME']
        print '##### No.{0} GRB{1} #####'.format(irow, name)
        if not os.path.exists(name):
            os.mkdir(name)
        acmd = ['bsub', '-o','{0}/{1}/GRB{1}_{2}{3}.log'.format(directory, name, mode, suffix if suffix=='' else '_'+suffix), '-J','{0}{1}'.format(name[:6], mode[:2]), '-W','300', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/STLikelihoodAnalysis/UnbinnedGRBAnalysis.py', '-m', mode, '--emin', str(emin), '--emax', str(emax), '-s', suffix, '--roi', str(roi), '--spectraltype', spectraltype, '--outdir', '{0}/{1}'.format(directory, name), name]
        if force==True:
            acmd.append('--force')
        if refit==True:
            acmd.append('--refit')
        if modelonly==True:
            acmd.append('--modelonly')
        if masifps==True:
            acmd.append('--masifps')
        print acmd
        subprocess.call(acmd)


if __name__ == '__main__':
    main()
