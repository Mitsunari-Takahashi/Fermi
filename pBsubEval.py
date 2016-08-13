#!/usr/bin/env python

import sys
import yaml
from array import array
par = sys.argv
import subprocess
from datetime import datetime

if len(sys.argv)<2:
    pathWeight = "/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatThree_15/weights/S18V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatThree_15_2_BDTG500D06.weights.xml"
else:
    pathWeight = par[1]
print pathWeight
timeTmp = datetime.now()
strTemp = "{0}{1}{2}{3}{4}".format(timeTmp.year, timeTmp.month, timeTmp.day, timeTmp.hour, timeTmp.minute)
dictFile = {"AG":"/nfs/farm/g/glast/u/mtakahas/data/MC/AG200909_62_2016Jun.root",
"BKG_0":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_0.root",
"BKG_1":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_1.root",
"BKG_2":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_2.root",
"BKG_3":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_3.root",
"BKG_4":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_4.root",
"BKGp2":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62p2MCE2e4.root",
"BKGp3b":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62p3bMCE2e4.root",
"BKGb8":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62b8MCE2e4.root",
"BKGb9":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62b9MCE2e4.root",
"BKGb10":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62b10MCE2e4.root",
"BKGb11":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62b11MCE2e4.root",
"BKGb12":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62b12MCE2e4.root",
"BKGb13":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62b13MCE2e4.root"
}
liKyes = dictFile.keys()
liValues = dictFile.values()
for iEl in range(len(liKyes)):
    strNN = liKyes[iEl]
    strFile = liValues[iEl]
    aCmd = ['bsub', '-o' '{0}eval{1}.log'.format(strTemp, strNN), "-J", "eval{0}".format(strNN), "-W 1000", "python", "~/eventSelect/scripts/EvalClassifier.py", strFile, "--weights={0}".format(pathWeight)]
    print aCmd
    subprocess.call(aCmd)
