#!/usr/bin/env python

import sys
import os
import os.path
import subprocess


files_input = sys.argv[1:]

for filein in files_input:
    print filein
    filebase = os.path.splitext(os.path.basename(filein))[0]
    subprocess.call(['bsub', '-o/nfs/farm/g/glast/u/mtakahas/data/LPA/AllSky/logs/AddNewVariables_{0}.log'.format(filebase), '-J NV{0}'.format(filebase[-10:]), '-W 20', '-We 2', '~/work/PythonModuleMine/Fermi/AddNewVariables', filein])
