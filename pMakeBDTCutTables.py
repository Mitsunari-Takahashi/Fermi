#!/usr/bin/env python

import sys
import ROOT
from ROOT import TTree
from ROOT import TChain
from ROOT import TH2D
import numpy as np
import yaml
import xml.etree.ElementTree as ET
import datetime
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
from pColor import *

ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)

#from pCutBDT import cutBDT
from pAnalysisConfig import *

# ----- Event class setup -----
par = sys.argv
cfg = ClassConfig('Both', [10, 3, 1], 1)
aCutEGB = cfg.aCutEGB

nameFileRoc = ["/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_ZDIR020to060_S18ZDIR020catTwoZDIR060Log_ZDIR020to060_roc.root", "/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_ZDIR060to100_S18ZDIR020catTwoZDIR060Log_ZDIR060to100_roc.root"]
nameVarBDT = "S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060"
aCutMVA = []
for iRocFile in range(len(nameFileRoc)):
    aCutMVA.append(CutBDT(nameFileRoc[iRocFile], aCutEGB))
aCosZ=[0.2, 0.6, 1.0]
cutMVA_comb = CutBDT_2D(aCutMVA, aCosZ, nameVarBDT)
print cutMVA_comb
aaaValCutBDT = cutMVA_comb.aaValCutBDT
print aaaValCutBDT
