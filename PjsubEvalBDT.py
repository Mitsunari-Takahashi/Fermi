#!/usr/bin/env python

import os
import sys
import subprocess
import click
from pLsList import ls_list


@click.command()
@click.argument('pathfiles', type=str)
@click.argument('pathweight', default='/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/weights/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060.weights.xml')
@click.option('--rscgrp', '-g', default="B", help="Remove the current friends")
@click.option('--priority', '-p', type=int, default="127", help="Priority of batch request. Set larger number for higher priority from 0 to 255.")
def main(pathfiles, rscgrp, priority):
    if priority<0 or priority>255:
        raise click.BadParameter('Priority is wrong!')
    li_path_files = ls_list(pathfiles)
    for path_file in li_path_files:
        path_file_base = os.path.basename(path_file)
        file_name, file_ext = os.path.splitext(path_file_base)
        index_week = file_name.find('_week')

        # BDT
        str_script = """#!/bin/sh
#------ pjsub option --------#
#PJM -L "rscunit=common"
#PJM -L "vnode=1"
#PJM -L "rscgrp={0}"
#------- Program execution -------#
alias date='\date +%Y%m%d%H%M'
export ROOTSYS=/home/takhsm/app/ROOT_TMVA1
export LD_LIBRARY_PATH=/usr/local/gcc473/lib64:/usr/local/gcc473/lib:${{ROOTSYS}}/lib:${{ROOTSYS}}/lib/root:${{ROOTSYS}}
export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH
export PATH=/home/takhsm/app/anaconda2/bin:/home/takhsm/app/Python2/bin:${{ROOTSYS}}/bin:/usr/local/gcc473/bin:$PATH
export LD_LIBRARY_PATH=${{ROOTSYS}}/lib:${{ROOTSYS}}/lib/root:/home/takhsm/app/lib:${{ROOTSYS}}:/home/takhsm/lib/lib:/usr/local/gcc473/lib64:/usr/local/gcc473/lib
export PYTHONPATH=/home/takhsm/app/ROOT_TMVA1/lib/root:/home/takhsm/FermiMVA/eventSelect/python:/home/takhsm/PythonModuleMine
export PATH_PYTHON_MODULE_MINE=/home/takhsm/PythonModuleMine

cd /disk/gamma/cta/store/takhsm/FermiMVA/AllSky/merit
python ~/FermiMVA/eventSelect/scripts/EvalClassifier.py {1} --weights={2}
""".format(rscgrp, path_file, path_weight)
        name_pjscript = "pjsub/EvalBDT_{0}.sh".format(file_name)
        f = open(name_pjscript,"w")
        f.write(str_script)
        f.close()
        if index_week>=2:
            job_name = 'D' + file_name[index_week-2:index_week] + 'w' + file_name[index_week+5:index_week+8]
        else:
            job_name = file_name
        aCmd = ['pjsub', name_pjscript, '-N', job_name, '-p', '{0}'.format(priority), '-o', 'pjsub/EvalBDT_{0}.log'.format(file_name), '-j']
        print aCmd
        subprocess.call(aCmd)


if __name__ == '__main__':
    main()
