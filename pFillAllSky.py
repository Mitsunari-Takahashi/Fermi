#!/usr/bin/env python

import sys
from array import array
import math
import yaml
import datetime
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle
# Original modules
from pAnalysisConfig import *
import pObjectCatalogue
from pTarget import *
from pColor import *

#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()

par = sys.argv
listFileTr = par[1:]
ch = ROOT.TChain("trGammas")
for nameFileTr in listFileTr:
    ch.Add(nameFileTr)
print "Total number of events:", ch.GetEntries()

nameFileOut = "PlotsAllSky.root"
fileOut = ROOT.TFile(nameFileOut, "RECREATE")

cfg = ClassConfig('Both', [10, 2, 1])
er = EnergyLogRegion(3, 4.75, 0.25)
aaStrSelect = cfg.aaStrSelect

fileOut.cd()
aTgt = []
print "Target list: "
print "*", "Galactic ridge"
ridge = GalacticRidge("GalacticRidge")
for obj in pObjectCatalogue.aObjectDict:
    print "*", obj["Name"]
    aTgt.append(PointSource(obj["Name"], obj["L"], obj["B"], 12., [12., 12.], [1., 8.], [0.45, 4.], cfg, er))
    #tgt = PointSource(obj["Name"], obj["L"], obj["B"], 12., [12., 12.], [1., 8.], [0.45, 4.], cfg, er)
print ""
timeStart = datetime.datetime.now()
print "Event filling started at", timeStart
step1 = ch.GetEntries()/200/100
step2 = step1*100
iStep=0
for iEve in range(ch.GetEntries()):
    ch.GetEntry(iEve)
    ridge.fill(ch.l, ch.b, ch.e, ch.c, ch.s)
    for tgt in aTgt:
        tgt.fill(ch.l, ch.b, ch.e, ch.c, ch.s)
    meter1 = ""
    if iEve%(step1)==0:
        iStep=iStep+1
        if iEve%(step2)==0:
            rate = int((iEve*100.)/ch.GetEntries()+0.5)
            if rate>0:
                nt = (datetime.datetime.now() - timeStart).seconds * (100.-rate)/rate
                meter0 = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
                meter1 = " Wait {0} hr {1} min".format(int(nt/3600), int((nt%3600)/60)+1)
            else:
                meter0 = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
            iStep=0
        meter2 = meter0 + "["
        jStep = iStep
        for iCount in [6,5,4,3,2,1,0]:
            if jStep/(2**iCount):
                meter2 = meter2 + "|"
                jStep = jStep - 2**iCount
            else:
                meter2 = meter2 + " "
        meter2 = meter2 + "]" + meter1
        # if iStep%2==0:
        #     meter2 = "\r" + "[+]" + meter
        # else:
        #     meter2 = "\r" + "[x]" + meter
        sys.stdout.write(meter2)
        sys.stdout.flush()
print ""
