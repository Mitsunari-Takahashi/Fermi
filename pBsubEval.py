#!/usr/bin/env python

import sys
import yaml
from array import array
par = sys.argv
import subprocess
from datetime import datetime

#if len(sys.argv)<2:
pathWeight = "/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/weights/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060.weights.xml" #"/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatThree_15/weights/S18V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatThree_15_2_BDTG500D06.weights.xml"
#else:
#pathWeight = par[1]
print pathWeight
timeTmp = datetime.now()
strTemp = "{0}{1}{2}{3}{4}".format(timeTmp.year, timeTmp.month, timeTmp.day, timeTmp.hour, timeTmp.minute)

strSuffix = par[1]
nStart = int(par[2])
aPathFile = par[3:]

for (iFile, pathFile) in enumerate(aPathFile):
    aCmd = ['bsub', '-o' 'logs/{0}evalBDT_{1}_{2}.log'.format(strTemp, strSuffix, nStart+iFile), "-J", "evalBDT_{0}_{1}".format(strSuffix, nStart+iFile), "-W 1000", "python", "~/eventSelect/scripts/EvalClassifier.py", pathFile, "--weights={0}".format(pathWeight)]
    print aCmd
    subprocess.call(aCmd)
