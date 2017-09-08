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
@click.option('--eminfit', type=float, default=177.828)
@click.option('--emaxfit', type=float, default=5623.41)
@click.option('--eminextrapolate', type=float, default=10000.)
@click.option('--emaxextrapolate', type=float, default=100000.)
@click.option('--mode', type=str, default='unified')
@click.option('--directory', '-d', type=str, default='/u/gl/mtakahas/work/FermiAnalysis/GRB/Regualr/HighestFluenceGRBs/LatAlone')
@click.option('--suffix', type=str, default='')
@click.option('--force', '-f', is_flag=True)
@click.option('--refit', '-r', is_flag=True)
@click.option('--roi', '-r', type=float, default=7.)
def main(namemin, namemax, eminfit, emaxfit, eminextrapolate, emaxextrapolate, mode, roi, directory, suffix, refit, force):
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
        acmd = ['bsub', '-o','{0}/{1}/GRB{1}_{2}{3}.log'.format(directory, name, mode, suffix if suffix=='' else '_'+suffix), '-J','{0}{1}'.format(name[:6], mode[:2]), '-W','300', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/STLikelihoodAnalysis/ExtrapolateGRBSpectrum.py', '-m', mode, '--eminfit', str(eminfit), '--emaxfit', str(emaxfit), '--eminextrapolate', str(eminextrapolate), '--emaxextrapolate', str(emaxextrapolate), '-s', suffix, '--roi', str(roi), '--outdir', '{0}/{1}'.format(directory, name), name]
        if force==True:
            acmd.append('--force')
        if refit==True:
            acmd.append('--refit')
        print acmd
        subprocess.call(acmd)


if __name__ == '__main__':
    main()
