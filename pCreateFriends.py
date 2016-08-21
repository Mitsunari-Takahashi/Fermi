#!/usr/bin/env python

import sys
import yaml
from array import array
par = sys.argv
import subprocess
from datetime import datetime


pathWeight = "/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1_BEP/Carmelo_caseE_new/WP8CalOnlyBEPCaseE_myBDT.weights.xml"
print pathWeight
timeTmp = datetime.now()
strTemp = "{0}{1}{2}{3}{4}".format(timeTmp.year, timeTmp.month, timeTmp.day, timeTmp.hour, timeTmp.minute)
strSuffix = par[1]
aPathFile = par[2:]

for (iFile, pathFile) in enumerate(aPathFile):
    aCmd = ['bsub', '-o' 'logs/{0}evalBEP_{1}_{2}.log'.format(strTemp, strSuffix, iFile), "-J", "evalBEP_{0}_{1}".format(strSuffix, iFile), "-W 1000", "python", "~/eventSelect/scripts/EvalClassifier.py", pathFile, "--weights={0}".format(pathWeight)]
    print aCmd
    subprocess.call(aCmd)
    bCmd = ['bsub', '-o' 'logs/{0}addNewVar_{1}_{2}.log'.format(strTemp, strSuffix, iFile), "-J", "addNewVar_{0}_{1}".format(strSuffix, iFile), "-W 1000", "python", "/nfs/farm/g/glast/u/mtakahas/PythonModuleMine/Fermi/pAddNewVariable.py", pathFile]
    print bCmd
    subprocess.call(bCmd)
