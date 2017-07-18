#!/usr/bin/env python
"""Submit jobs of AnalyzeGRB_fermipy.py to SLAC server.
"""
import sys
import os
import subprocess
import click
import csv


@click.command()
@click.argument('name', type=str)
@click.argument('timefile', type=str)
@click.argument('tmin', type=float)
@click.argument('tmax', type=float)
@click.option('--outpath', '-o', default='.')
@click.option('--mode', '-m', type=click.Choice(['prompt', 'afterglow']))
@click.option('--catalogues', '-c', multiple=True, default=None, type=str)
@click.option('--goodstat', '-g', type=int, default=0)
#@click.option('--thirtygev', is_flag=True)
@click.option('--suffix', type=str, default='')
@click.option('--force', '-f', is_flag=True)
def main(name, timefile, tmin, tmax, outpath, mode, catalogues, goodstat, suffix, force):
    with open(timefile, 'r') as csvfile:
        reader = csv.reader(csvfile) #, delimiter=' ')
        header = next(reader)
        for (irow,row) in enumerate(reader):
            start = float(row[0])
            stop = float(row[1])
            print 'Time:', start, stop
            acmd = ['bsub', '-o','GRB{0}_{1}_T{2:0>6}-{3:0>6}.log'.format(name, mode, int(start), int(stop)), '-J','{0}{1}{2}'.format(name[:4], mode[:2], irow), '-W','300', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/AnalyzeGRB_fermipy.py', '--emin', '100.0', '--emax', '316228.0',  '-o', outpath, '--skipsed', '--skipts', '--skipresid', '-m', mode, '-s', suffix, name, str(start), str(stop)]
            if catalogues is not None and len(catalogues)>0:
                for cat in catalogues:
                    acmd.append('-c')
                    acmd.append(cat)
            if goodstat>0:
                acmd.append('--goodstat')
                acmd.append(goodstat)
            #if thirtygev is True:
             #   acmd.append('--thirtygev')
            if force is True:
                acmd.append('--force')
            print acmd
            subprocess.call(acmd)


if __name__ == '__main__':
    main()
