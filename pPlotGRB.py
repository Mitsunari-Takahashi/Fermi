#!/usr/bin/env python

import sys
from array import array
import math
#import yaml
import datetime
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
# Original modules
from pAnalysisConfig import *
import pGrbCatalogue
from pTarget import *
#from pTarget_GRB import *
from pColor import *
from pMETandMJD  import *
par = sys.argv
print par
if not len(par)==4:
    print "pPlotGRB.py [name of GRB] [suffix] [tree file]"
    sys.exit(0)
suffix = par[2]
nameGRB = par[1]
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()

ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)

psfCalOnly = 3.0
listFileTr = par[3]
print listFileTr
ch = ROOT.TChain("trGammas")
#for nameFileTr in listFileTr:
#ch.Add(nameFileTr)
ch.Add(listFileTr)
print "Total number of events:", ch.GetEntries()

nameFileOut = "PlotsGRB_" + nameGRB + "_" + suffix + ".root"
#nameFileOut = "PlotsAllSky_" + suffix + "_" + str(int(psfCalOnly*10)) + ".root"
fileOut = ROOT.TFile(nameFileOut, "RECREATE")
print "Output file:", fileOut.GetName()

#listPathFilePerf = ["/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R100_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R30_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R10_perf.root"]
listPathFilePerf = [["/home/takhsm/FermiMVA/CalTkr/v20r9p9_TRANSIENT_P8R1_TRANSIENT_R100_perf.root", "/home/takhsm/FermiMVA/CalTkr/v20r9p9_SOURCE_P8R1_SOURCE_perf.root"], ["/home/takhsm/FermiMVA/S11/S11V200909_020RAWE20ZCS000wwoTRKwBKwoMCZDIR00woRWdivCALE_15/v20r9p9_S11V200909_020RAWE20ZCS000wwoTRKwBKwoMCZDIR00woRWdivCALE_15_CalOnly_R100_perf.root", "/home/takhsm/FermiMVA/S11/S11V200909_020RAWE20ZCS000wwoTRKwBKwoMCZDIR00woRWdivCALE_15/v20r9p9_S11V200909_020RAWE20ZCS000wwoTRKwBKwoMCZDIR00woRWdivCALE_15_CalOnly_R30_perf.root", "/home/takhsm/FermiMVA/S11/S11V200909_020RAWE20ZCS000wwoTRKwBKwoMCZDIR00woRWdivCALE_15/v20r9p9_S11V200909_020RAWE20ZCS000wwoTRKwBKwoMCZDIR00woRWdivCALE_15_CalOnly_R10_perf.root"]]
cfg = ClassConfig('Both', [10, 3, 1])
er = EnergyLogRegion(5, 4.5, 0.25)
aaStrSelect = cfg.aaStrSelect
aStrSelect = cfg.aStrSelect

htgPerf = CutPerformanceHtg(listPathFilePerf)
fileOut.cd()

rOffMax=[15., 15.]#[1.7, 12.]
rAppa=rOffMax[1]
rOnMax=[0.45, psfCalOnly]
rOffMin=[2.*rOnMax[0], 2.*rOnMax[1]]

aGrb = []
print "Target list: "
#for obj in pGrbCatalogue.aObjectDict:
#    if obj["Name"]==nameGRB:
#        aGrb.append(PointSource(obj["Name"], obj["RA"], obj["DEC"], 0, 0, obj["Z"], rAppa, rOffMax, rOffMin, rOnMax, cfg, er, er, htgPerf, obj["Start"], obj["Start"]+10000))
        #aGrb.append(PointSource(obj["Name"], obj["RA"], obj["DEC"], obj["Z"], rAppa, rOffMax, rOffMin, rOnMax, cfg, er, htgPerf, ch.GetMinimum("Time"), ch.GetMaximum("Time")))

nameFileGrb = "/home/takhsm/FermiMVA/python/GBM_fluenceLT100micro_summary.root"
fileGrb = ROOT.TFile(nameFileGrb, "READ")
trGrb = fileGrb.Get("grb")
for iGrb in range(trGrb.GetEntries()):
    trGrb.GetEntry(iGrb)
    if trGrb.name==int(nameGRB):
        aGrb.append(PointSource(trGrb.name, trGrb.ra, trGrb.dec, trGrb.lii, trGrb.bii, "N/A", rAppa, rOffMax, rOffMin, rOnMax, cfg, er, er, htgPerf, ConvertMjdToMet(trGrb.trigger_time), ConvertMjdToMet(trGrb.trigger_time)+10000))

print ""
print "================"
print "Filling events."
print "================"
timeStart = datetime.datetime.now()
print "Started at", timeStart
step1 = ch.GetEntries()/100
if step1>=1:
    iStep=0
    for iEve in range(ch.GetEntries()):
        ch.GetEntry(iEve)
        for tgt in aGrb:
            tgt.fill(ch.ra, ch.dec, ch.l, ch.b, ch.e, ch.c, ch.s, ch.z, ch.t, ch.cth)
        meter1 = ""
        if iEve%(step1)==0:
            iStep=iStep+1
            rate = int((iEve*100.)/ch.GetEntries()+0.5)
            if rate>0:
                nt = (datetime.datetime.now() - timeStart).seconds * (100.-rate)/rate
                meter0 = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
                meter1 = " Wait {0} hr {1} min".format(int(nt/3600), int((nt%3600)/60)+1)
            else:
                meter0 = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
            iStep=0
        sys.stdout.flush()
print ""

print "================"
print "Making plots."
print "================"

for tgt in aGrb:
    print "*", tgt.name
    fileOut.cd()
#    fileOut.mkdir(tgt.name)
#    fileOut.cd(tgt.name)
    tgt.calc()
#    fileOut.cd(tgt.name)
    tgt.draw()
    tgt.writeObjects()
print "Plots were saved."
print "All analysis finished."
