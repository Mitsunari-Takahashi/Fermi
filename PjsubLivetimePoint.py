#!/usr/bin/env python

import sys
import subprocess
import click
from pLsList import ls_list
from pReadGBMCatalogueInfo import ReadGBMCatalogueOneLine


@click.command()
@click.argument('srcname', type=str)
@click.option('--ra', type=float, default=None)
@click.option('--dec', type=float, default=None)
@click.option('--start', type=float, default=239557417.)
@click.option('--stop', type=float, default=501033396.)
@click.option('--exstart', type=float, default=0)
@click.option('--exstop', type=float, default=0)
@click.option('--fixenergy', type=float, default=0.)
@click.option('--fixinclin', type=float, default=0.)
@click.option('--truepoint', is_flag=True)
@click.option('--rbkg', '-b', type=click.Choice(['100', '030', '010', '003']))
@click.option('--rscgrp', '-g', default="B", type=str)
@click.option('--gbmtimeoff', '-m', is_flag=True)
def main(srcname, ra, dec, start, stop, fixenergy, fixinclin, truepoint, rscgrp, rbkg, exstart, exstop, gbmtimeoff):
    if gbmtimeoff is True:
        dct_grb = ReadGBMCatalogueOneLine(srcname)
        ra = dct_grb['ra']
        dec = dct_grb['dec']
        exstart = dct_grb['trigger_time'] - 2.*86400.
        exstop = dct_grb['trigger_time'] + 20.*86400.
    elif ra is None or dec is None:
        print 'Please use --gbmcatalogue or both of --ra and --dec.'
    perf = '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28bin_Cth40bins_axisObs_CalOnly_R{0}_perf.root'.format(rbkg)
    print 'IRF:', perf
    if start<234316801:
        raise click.BadParameter('start MET is wrong!')
    if start>=stop:
        raise click.BadParameter('stop MET is wrong!')
    tduration = stop - start
    TWEEK = 7.*86400.
    nweek = int(tduration/TWEEK+1)
    li_time_set = [start]
    for iweek in range(nweek-1):
        li_time_set.append(li_time_set[-1]+TWEEK)
    li_time_set.append(stop)    
    print "Time division:"
    print li_time_set

    str_truepoint = ''
    if truepoint==True:
        str_truepoint = '--truepoint'
#    if suffix!="":
#        suffix = "_" + suffix

    for jweek in range(nweek):
        str_script = """#!/bin/sh
#------ pjsub option --------#
#PJM -L "rscunit=common"
#PJM -L "vnode=1"
#PJM -L "rscgrp={0}"
#------- Program execution -------#
alias date='\date +%Y%m%d%H%M'
export PATH=/home/takhsm/app/anaconda2/bin:/home/takhsm/app/Python2/bin:/usr/local/gcc473/bin:${{ROOTSYS}}/bin:$PATH
export LD_LIBRARY_PATH=/home/takhsm/lib/lib:/usr/local/gcc473/lib64:/usr/local/gcc473/lib:${{ROOTSYS}}/lib:${{ROOTSYS}}/lib/root:/home/takhsm/app/lib:${{ROOTSYS}}:${{MARSSYS}}
export PYTHONPATH=/usr/local/gcc473/ROOT/lib/root:/home/takhsm/FermiMVA/eventSelect/python:/home/takhsm/app/yaml/lib64/python:/home/takhsm/app/lib/python2.7/site-packages:/home/takhsm/app/Python2/lib/python2.7/site-packages:/home/takhsm/PythonModuleMine
export PATH_PYTHON_MODULE_MINE=/home/takhsm/PythonModuleMine
export ROOTSYS=/usr/local/gcc473/root_v5.34.14/
export LD_LIBRARY_PATH=/usr/local/gcc473/lib64:/usr/local/gcc473/lib:${{ROOTSYS}}/lib:${{ROOTSYS}}/lib/root:${{ROOTSYS}}
export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH

#python ~/PythonModuleMine/Fermi/pLivetimePoint.py {1}_{2:0>3} {3} {4} {5} {6}
python ~/PythonModuleMine/Fermi/LivetimePointModel.py --perf {10} --suffix {2:0>3} {7} --start {3} --stop {4} --energy {5} --inclin {6} --exstart {11} --exstop {12} -- GRB{1} {8} {9}
""".format(rscgrp, srcname, jweek, li_time_set[jweek], li_time_set[jweek+1], fixenergy, fixinclin, str_truepoint, ra, dec, perf, exstart, exstop)
        name_pjscript = "LtP_{0}{1}.sh".format(srcname, jweek)
        f = open(name_pjscript,"w")
        f.write(str_script)
        f.close()
        #subprocess.call(['cat', name_pjscript])
        subprocess.call(['pjsub', name_pjscript])


if __name__ == '__main__':
    main()
