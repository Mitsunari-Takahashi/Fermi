#!/usr/bin/env python
"""Submit jobs of AnalyzeGRB_fermipy.py to SLAC server.
"""
import sys
import os
import os.path
import subprocess
import click
import csv
import ReadLTFCatalogueInfo


@click.command()
@click.option('--namemin', type=str, default='0')
@click.option('--namemax', type=str, default='200000000')
#@click.argument('timefile', type=str)
@click.option('--tmin', type=float, default=None)
@click.option('--tmax', type=float, default=None)
@click.option('--emin', type=float, default=0)
@click.option('--emax', type=float, default=0)
@click.option('--outpath', '-o', default=None)
@click.option('--mode', '-m', type=click.Choice(['prompt', 'afterglow', 'unified', 'earlyAG', 'lateAG', 'special']))
@click.option('--catalogues', '-c', multiple=True, default=None, type=str)
@click.option('--func', multiple=True, default=None, type=str)
@click.option('--goodstat', '-g', type=int, default=0)

@click.option('--suffix', type=str, default='')
@click.option('--force', '-f', is_flag=True)
@click.option('--edisp', is_flag=True)
@click.option('--download', '-d', is_flag=True)
@click.option('--shiftenergies', is_flag=True)
@click.option('--roi', '-r', type=float, default=12)
@click.option('--sedadjusted', is_flag=True)
def main(namemin, namemax, tmin, tmax, emin, emax, outpath, mode, catalogues, goodstat, edisp, func, suffix, force, download, shiftenergies, roi, sedadjusted):
    #with open(timefile, 'r') as csvfile:
     #   reader = csv.reader(csvfile) #, delimiter=' ')
     #   header = next(reader)
        #for (irow,row) in enumerate(reader):
    tb_ltf = ReadLTFCatalogueInfo.open_table()
    tb_gbm = ReadLTFCatalogueInfo.select_gbm_exist(tb_ltf)
    tb_selected = ReadLTFCatalogueInfo.select_by_name(tb_gbm, namemin, namemax)
    #tb_zknown = ReadLTFCatalogueInfo.select_redshift_known(tb_selected)
    tb = tb_selected
    if len(tb)<1:
        print 'No GRBs.'
        return 1
    for (irow , row) in enumerate(tb):
        name = row['GRBNAME']
        print '##### No.{0} GRB{1} #####'.format(irow, name)
        # if mode=='prompt':
        #     start = row['T90_START']
        #     stop = row['T90_START'] + row['T90']
        # elif mode=='afterglow':
        #     start = row['T90_START'] + row['T90']
        #     stop = 10000.
        # elif mode=='unified':
        #     start = row['T90_START']
        #     stop = 10000.
        # if tmin is not None:
        #     start = tmin
        # if tmax is not None:
        #     stop = tmax
        # print 'Time:', start, '-', stop
        if not os.path.exists('logs'):
            os.mkdir('logs')
#        acmd = ['bsub', '-o','logs/GRB{0}_{1}_T{2:0>6}-{3:0>6}.log'.format(name, mode, int(start), int(stop)), '-J','{0}{1}'.format(name[:6], mode[:2]), '-W','300', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/AnalyzeGRB_fermipy.py', '--emin', str(emin), '--emax', str(emax), '--skipts', '--skipresid', '-m', mode, '-s', suffix, '--tmin', str(start), '--tmax', str(stop), '--roi', str(roi), name]
        acmd = ['bsub', '-o','logs/GRB{0}_{1}{2}.log'.format(name, mode, suffix), '-J','{0}{1}'.format(name[:6], mode[:2]), '-W','300', 'python', '/u/gl/mtakahas/work/PythonModuleMine/Fermi/AnalyzeGRB_fermipy.py', '--emin', str(emin), '--emax', str(emax), '--skipts', '--skipresid', '-m', mode, '-s', suffix, '--roi', str(roi), name]
        if catalogues is not None and len(catalogues)>0:
            for cat in catalogues:
                acmd.append('-c')
                acmd.append(cat)
        if func is not None and len(func)>0:
            for fc in func:
                acmd.append('--func')
                acmd.append(fc)
        if goodstat>0:
            acmd.append('--goodstat')
            acmd.append(goodstat)
        if outpath is not None:
            acmd.append('-o')
            acmd.append(outpath)
        if download==True:
            acmd.append('--download')
        if force is True:
            acmd.append('--force')
        if edisp is True:
            acmd.append('--edisp')
        if shiftenergies is True:
            acmd.append('--shiftenergies')
        if sedadjusted is True:
            acmd.append('--sedadjusted')
        print acmd
        subprocess.call(acmd)


if __name__ == '__main__':
    main()
