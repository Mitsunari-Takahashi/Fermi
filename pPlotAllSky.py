#!/usr/bin/env python

import sys
from astropy.io import fits
from array import array
import math
import numpy as np
import datetime
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
import healpy as hp
from healpy import pixelfunc as hppf
# Original modules
from pAnalysisConfig import *
import pCandleCatalogue
import pTransientCatalogue
from pTarget import *
from pHealplot import Healcube, Setdistance
from pColor import *


par = sys.argv
print par
if len(par)<2:
    print "pPlotAllSky.py [suffix] [list of tree file]"
    sys.exit(0)
suffix = par[1]

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)
ROOT.gStyle.SetPalette(53)

psfCalOnly = 3.0
listFileTr = par[2:]
print listFileTr
ch = ROOT.TChain("EVENTS")
for nameFileTr in listFileTr:
    ch.Add(nameFileTr)
print "Total number of events:", ch.GetEntries()

nameFileOut = "PlotsAllSky_" + suffix + ".root"

fileOut = ROOT.TFile(nameFileOut, "RECREATE")
print "Output file:", fileOut.GetName()

listPathFilePerf = [
["/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_TRANSIENT100_P8R2_TRANSIENT100_perf.root", 
"/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_SOURCE_P8R2_SOURCE_perf.root"],
["/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R100_perf.root", 
"/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R30_perf.root", 
"/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R10_perf.root"]]
# listPathFilePerf = [
# ["/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/CalTkr/v20r09p09_P8R2_SOURCE_P8R2_TRANSIENT100_perf.root", 
# "/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/CalTkr/v20r09p09_P8R2_SOURCE_P8R2_SOURCE_perf.root"],
# ["/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R100_perf.root", 
# "/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R30_perf.root", 
# "/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R10_perf.root"]]

cfg = ClassConfig('Both', [10, 3, 1])
er = EnergyLogRegion(7, 4.35, 0.2)
erplot = er #EnergyLogRegion(4, 4.75, 0.25)
cthr = EnergyLogRegion(10, 0.0, 0.1)
aaStrSelect = cfg.aaStrSelect
aStrSelect = cfg.aStrSelect

htgPerf = CutPerformanceHtg(listPathFilePerf)
fileOut.cd()

rOffMax=[15., 15.]#[1.7, 12.]
rAppa=20#rOffMax[1]
rOnMax=[0.45, psfCalOnly]
rOffMin=[2.*rOnMax[0], 2.*rOnMax[1]]

aTgt = []
aTrs = []
print "Target list: "
print "*", "Galactic ridge" 
ridge = GalacticRidge("GalacticRidge", config=cfg, perf=htgPerf, eRegion=er, ePlotRegion=erplot)
limb = EarthLimb("EarthLimb", config=cfg, perf=htgPerf, eRegion=er, ePlotRegion=erplot)
#aInnerGal = [InnerGalaxy("InnerGalaxyR16", 16, config=cfg, perf=htgPerf, eRegion=er, ePlotRegion=erplot), InnerGalaxy("InnerGalaxyR41", 41, config=cfg, perf=htgPerf, eRegion=er, ePlotRegion=erplot)]
aInnerGal = [ InnerGalaxy("InnerGalaxyR41", 41, config=cfg, perf=htgPerf, eRegion=er, ePlotRegion=erplot)]
for obj in pCandleCatalogue.aObjectDict:
    aTgt.append(PointSource(obj["Name"], obj["RA"], obj["DEC"], obj["L"], obj["B"], obj["Z"], rAppa, rOffMax, rOffMin, rOnMax, cfg, er, erplot, htgPerf))

#HEALPix plot
NHPSIDE = 64
NHPPIX = hppf.nside2npix(NHPSIDE)
print 'HEALPix resolution:', degrees(hppf.nside2resol(NHPSIDE)), 'deg'
lstt_hp_htg = []
#li_hp_count_tel = []
for (icat, cat) in enumerate(aStrSelect):
    lstt_hp_htg.append([])
    #li_hp_count_tel.append([])
    for (icla, cla) in enumerate(aaStrSelect[icat]):
        lstt_hp_htg[-1].append(ROOT.TH3D('htg3D_{0}'.format(cla), cla, erplot.nBin, erplot.edgeLow, erplot.edgeUp, cthr.nBin, cthr.edgeLow, cthr.edgeUp, NHPPIX, 0, NHPPIX))
       #li_hp_count_tel[-1].append(np.zeros(NHPPIX))

# Smearing
Setdistance(NHPSIDE)

print ""
print "================"
print "Filling events."
print "================"
timeStart = datetime.datetime.now()
print "Started at", timeStart
step1 = ch.GetEntries()/100
step2 = step1*100
iStep=0
for iEve in range(ch.GetEntries()):
    ch.GetEntry(iEve)
    cls = 0
    if ch.s==128 or ch.s==16384:
        cls = 3
    elif ch.s==8192:
        cls = 2
    elif ch.s==4 or ch.s==4096:
        cls = 1
    energybin = erplot.findBin(ch.e)
    if ch.z<90 and energybin>=0 and energybin<erplot.nBin:
        for clsPlus in range(cls-int(ch.c==1 and cls==3)):
            lstt_hp_htg[ch.c-1][clsPlus].Fill(ch.e, ch.cth, hppf.ang2pix(NHPSIDE, math.pi/2.-math.radians(ch.dec), math.radians(ch.ra))+0.5)
           #li_hp_count_tel[ch.c-1][clsPlus][hppf.ang2pix(NHPSIDE, pi/2.-math.radians(ch.dec), math.radians(ch.ra))]=+1
    ridge.fill(ch.ra, ch.dec, ch.l, ch.b, ch.e, ch.c, cls, ch.z, ch.t, ch.cth)
    limb.fill(ch.ra, ch.dec, ch.l, ch.b, ch.e, ch.c, cls, ch.z, ch.t, ch.cth)
    for ingal in aInnerGal:
        ingal.fill(ch.ra, ch.dec, ch.l, ch.b, ch.e, ch.c, cls, ch.z, ch.t)
    for tgt in aTgt:
        tgt.fill(ch.ra, ch.dec, ch.l, ch.b, ch.e, ch.c, cls, ch.z, ch.t, ch.cth)
    for trs in aTrs:
        trs.fill(ch.ra, ch.dec, ch.l, ch.b, ch.e, ch.c, cls, ch.z, ch.t, ch.cth)
    meter1 = ""
    if iEve%(step1)==0:
#        print ch.ra, ch.dec, ch.l, ch.b, ch.e, ch.c, cls, ch.z, ch.t, ch.cth
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

print "================"
print "Making plots."
print "================"

print "*", ridge.name
fileOut.cd()
fileOut.mkdir(ridge.name)
fileOut.cd(ridge.name)
ridge.calc()
fileOut.cd(ridge.name)
ridge.draw(lstt_hp_htg)
ridge.writeObjects()
for lstcat in lstt_hp_htg:
    for htgcla in lstcat:
        htgcla.Write()
print "*", limb.name
fileOut.cd()
fileOut.mkdir(limb.name)
fileOut.cd(limb.name)
limb.calc()
fileOut.cd(limb.name)
limb.draw()
limb.writeObjects()
for ingal in aInnerGal:
    print "*", ingal.name
    fileOut.cd()
    fileOut.mkdir(ingal.name)
    fileOut.cd(ingal.name)
    ingal.calc()
    fileOut.cd(ingal.name)
    ingal.draw()
    ingal.writeObjects()
for tgt in aTgt:
    print "*", tgt.name
    fileOut.cd()
    fileOut.mkdir(tgt.name)
    fileOut.cd(tgt.name)
    tgt.calc()
    fileOut.cd(tgt.name)
    tgt.draw(lstt_hp_htg)
    tgt.writeObjects()
for trs in aTrs:
    print "*", trs.name
    fileOut.cd()
    fileOut.mkdir(trs.name)
    fileOut.cd(trs.name)
    trs.calc()
    fileOut.cd(trs.name)
    trs.draw()
    trs.writeObjects()
print "Plots were saved."

print "================="
print "Stacking analysis"
print "================="
fileOut.mkdir("Stacking")
print "Drawing maps."
aaaGrpMap = []
aaMgrMap = []
cMap = ROOT.TCanvas("cMapStack", "%s sources stacked map" % len(pCandleCatalogue.aObjectDict), 900, 600)
nDown = len(aaStrSelect)
nAcross = 0
for aS in aaStrSelect:    
    nAcross = max(nAcross, len(aS))
cMap.Divide(nAcross, nDown)
legGrpMap = ROOT.TLegend(0.1,0.1,0.9,0.9,"Event energy", 'NDC')

aaTitleMap = []
for ctEvent in range(len(aStrSelect)):
    aaTitleMap.append([])
    for clEvent in range(len(aaStrSelect[ctEvent])):
        aaTitleMap[-1].append(ROOT.TPaveText(0.1, 1.05, 0.9, 0.9, "NDC"))
        aaTitleMap[-1][-1].AddText("%s map of %s stacked sources" % (aaStrSelect[ctEvent][clEvent], len(pCandleCatalogue.aObjectDict)))
        aaTitleMap[-1][-1].SetFillStyle(1001)
        aaTitleMap[-1][-1].SetFillColor(kWhite)

for pS in range(len(aaStrSelect)):
    aaaGrpMap.append([])
    aaMgrMap.append([])
    for qS in range(len(aaStrSelect[pS])):
        aaaGrpMap[-1].append([])
        aaMgrMap[-1].append(ROOT.TMultiGraph("mgrMap%s_%s" % (pS,qS),"%s map (%s sources stacked)" % (aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict))))
        for rE in range(erplot.nBin):
            aaaGrpMap[-1][-1].append(ROOT.TGraphPolar())            
            aaaGrpMap[-1][-1][-1].SetName("grpMapStack%s_%s_%s" % (pS,qS,rE))
            aaaGrpMap[-1][-1][-1].SetMarkerStyle(6)
            aaaGrpMap[-1][-1][-1].SetTitle("{0} map ({1} sources stacked) in {2:.1f} - {3:.1f} GeV".format(aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict), 10**(erplot.getBin(rE)[0]-3), 10**(erplot.getBin(rE)[1]-3)))
            aaaGrpMap[-1][-1][-1].SetMinRadial(0.)
            aaaGrpMap[-1][-1][-1].SetMaxRadial(15.)
            aaaGrpMap[-1][-1][-1].SetMarkerColor(pColor.akColor(rE))
            aaaGrpMap[-1][-1][-1].SetMarkerStyle(6)
            for pbj in pCandleCatalogue.aObjectDict:
                fileOut.cd(pbj["Name"])
                grp = gDirectory.Get("grpMap%s_%s_%s_%s" % (pbj["Name"], pS, qS, rE))
                ps = [[ROOT.Double(), ROOT.Double()] for idou in range(grp.GetN())]
                for iEvent in range(grp.GetN()):
                    grp.GetPoint(iEvent, ps[iEvent][0], ps[iEvent][1])
                    aaaGrpMap[-1][-1][-1].SetPoint(aaaGrpMap[-1][-1][-1].GetN(), ps[iEvent][0], ps[iEvent][1])
            aaMgrMap[-1][-1].Add(aaaGrpMap[-1][-1][-1])
            if pS==0 and qS==0:
                legGrpMap.AddEntry(aaaGrpMap[pS][qS][rE], "{0:.1f} - {1:.1f} GeV".format(10**(erplot.getBin(rE)[0]-3), 10**(erplot.getBin(rE)[1]-3)), 'p')
            fileOut.cd("Stacking")
            aaaGrpMap[-1][-1][-1].Write("", ROOT.TObject.kOverwrite)
            aaaGrpMap[-1][-1][-1].SetTitle("")
        cMap.cd(1+pS*nAcross+qS)
        aaMgrMap[-1][-1].Draw("NP")
        aaTitleMap[pS][qS].Draw('same')
        gPad.Update()
        for rE in range(len(aaaGrpMap[pS][qS])):
            try:
                aaaGrpMap[pS][qS][rE].GetPolargram().SetRangeRadial(0.0, 12.0)
                aaaGrpMap[pS][qS][rE].GetPolargram().SetNdivPolar(104)
                aaaGrpMap[pS][qS][rE].GetPolargram().SetNdivRadial(101)
            except Exception:
                print aaaGrpMap[pS][qS][rE].GetName(), "Polargram failed."
        fileOut.cd("Stacking")
        aaMgrMap[-1][-1].Write("", ROOT.TObject.kOverwrite)
cMap.cd(3)
legGrpMap.Draw()
fileOut.cd("Stacking")
cMap.Write("", ROOT.TObject.kOverwrite)

mDown = int(math.sqrt(erplot.nBin+1))
mAcross = int(math.ceil((erplot.nBin+1) / mDown))
lDown = len(aaStrSelect)
lAcross = 0
aaaHtgMapAllSky = []
for aS in aaStrSelect:    
    lAcross = max(lAcross, len(aS))

print "Plotting theta square and calculating some values."
rOnMax = []
rOffMin = []
rOffMax = []
fRadiusOnMax = []
fRadiusOffMin = []
fRadiusOffMax = []
for ctEvent in range(len(aaStrSelect)):
    rOnMax.append([])
    rOffMin.append([])
    rOffMax.append([])
    fRadiusOnMax.append([])
    fRadiusOffMin.append([])
    fRadiusOffMax.append([])
    for clEvent in range(len(aaStrSelect[ctEvent])):
        rOnMax[ctEvent].append([])
        rOffMin[ctEvent].append([])
        rOffMax[ctEvent].append([])
        fRadiusOnMax[ctEvent].append([])
        fRadiusOffMin[ctEvent].append([])
        fRadiusOffMax[ctEvent].append([])
        for binE in range(er.nBin):
            rOnMax[ctEvent][clEvent].append(htgPerf.getPSF68(ctEvent, clEvent, er.aBin[binE]+er.wBin/2.0))
            rOffMin[ctEvent][clEvent].append(htgPerf.getPSF68(ctEvent, clEvent, er.aBin[binE]+er.wBin/2.0))
            rOffMax[ctEvent][clEvent].append(rAppa)

            fRadiusOnMax[ctEvent][clEvent].append(ROOT.TF2("fRadiusOnMax%s_%s_%s" % (ctEvent, clEvent, binE), "x**2+y**2 - %s**2" % rOnMax[ctEvent][clEvent][binE], -rAppa, rAppa, -rAppa, rAppa))
            fRadiusOnMax[-1][-1][-1].SetMinimum(0)
            fRadiusOnMax[-1][-1][-1].SetMaximum(0)
            fRadiusOnMax[-1][-1][-1].SetLineWidth(1)
            fRadiusOnMax[-1][-1][-1].SetLineColor(kGray)
            fRadiusOffMin[ctEvent][clEvent].append(ROOT.TF2("fRadiusOffMin%s_%s_%s" % (ctEvent, clEvent, binE), "x**2+y**2 - %s**2" % rOffMin[ctEvent][clEvent][binE], -rAppa, rAppa, -rAppa, rAppa))
            fRadiusOffMin[-1][-1][-1].SetMinimum(0)
            fRadiusOffMin[-1][-1][-1].SetMaximum(0)
            fRadiusOffMin[-1][-1][-1].SetLineWidth(1)
            fRadiusOffMin[-1][-1][-1].SetLineColor(kGray)
            fRadiusOffMax[ctEvent][clEvent].append(ROOT.TF2("fRadiusOffMax%s_%s_%s" % (ctEvent, clEvent, binE), "x**2+y**2 - %s**2" % rOffMax[ctEvent][clEvent][binE], -rAppa, rAppa, -rAppa, rAppa))
            fRadiusOffMax[-1][-1][-1].SetMinimum(0)
            fRadiusOffMax[-1][-1][-1].SetMaximum(0)
            fRadiusOffMax[-1][-1][-1].SetLineWidth(1)
            fRadiusOffMax[-1][-1][-1].SetLineColor(kGray)

aCanSpacial = []
aaaHtgMap = []
aaaHtgDummySpacial = []
for rE in range(erplot.nBin):
    aCanSpacial.append(ROOT.TCanvas("cSpacialStack_%s" % rE, "Spacial distribution of stacked sources in {0:.1f} - {1:.1f} GeV".format(10**(erplot.getBin(rE)[0]-3), 10**(erplot.getBin(rE)[1]-3)), 1200, 800))
    aCanSpacial[-1].Divide(lAcross, lDown)
    aaaHtgMap.append([])
    aaaHtgDummySpacial.append([])
    for pS in range(len(aaStrSelect)):
        aaaHtgMap[-1].append([])
        aaaHtgDummySpacial[-1].append([])
        for qS in range(len(aaStrSelect[pS])):
            aaaHtgMap[-1][-1].append(ROOT.TH2D("hMapStack_%s_%s_%s" % (pS, qS, rE), "{0} spacial distribution of {1} in {2:.1f} - {3:.1f} GeV;;;[counts/sr]".format(aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict), 10**(erplot.getBin(rE)[0]-3), 10**(erplot.getBin(rE)[1]-3)), 12, 0., 2.*math.pi, 12, 0, 12))
            aaaHtgDummySpacial[-1][-1].append(ROOT.TH2D("hDummySpacial_%s_%s_%s" % (rE, pS, qS), "{0} pacial distribution of {1} stacked sources in {2:.1f} - {3:.1f} GeV;[#circ];[#circ]".format(aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict), 10**(erplot.getBin(rE)[0]-3), 10**(erplot.getBin(rE)[1]-3)), 2*12, -12, 12, 2*12, -12, 12))
            listHtgMap = ROOT.TList()
            for rbj in pCandleCatalogue.aObjectDict:
                fileOut.cd(rbj["Name"])
                listHtgMap.Add(gDirectory.Get("hMap%s_%s_%s_%s" % (rbj["Name"], pS, qS, rE)))
            aaaHtgMap[rE][pS][qS].Merge(listHtgMap)
            aCanSpacial[rE].cd(1+pS*lAcross+qS)
            aCanSpacial[rE].cd(1+pS*lAcross+qS).SetGridx(0)
            aCanSpacial[rE].cd(1+pS*lAcross+qS).SetGridy(0)
            #aaaHtgDummySpacial[rE][pS][qS].SetStats(kFALSE)
            aaaHtgDummySpacial[rE][pS][qS].Draw()
            aaaHtgMap[rE][pS][qS].SetStats(kFALSE)
            aaaHtgMap[rE][pS][qS].Draw("POL COLZ SAMES")
            fRadiusOnMax[pS][qS][rE].Draw('same')
            fRadiusOffMin[pS][qS][rE].Draw('same')
            fRadiusOffMax[pS][qS][rE].Draw('same')
            gPad.Update()
            palette = aaaHtgMap[rE][pS][qS].GetListOfFunctions().FindObject("palette")
            try:
                palette.SetX1NDC(0.88)
                palette.SetX2NDC(0.93)
            except Exception:
                print "Setting of", palette, "failed."

            fileOut.cd("Stacking")
            aaaHtgMap[rE][pS][qS].Write("", ROOT.TObject.kOverwrite)
    aCanSpacial[rE].Write("", ROOT.TObject.kOverwrite)

cTheta = ROOT.TCanvas("cThetaStack", "%s sources stacked #theta^{2} plot" % len(pCandleCatalogue.aObjectDict), 900, 600)
cTheta.Divide(mAcross, mDown)
legTheta = ROOT.TLegend(0.1,0.1,0.9,0.9,"Event class",'NDC')
aaaHtgTheta = []
aHStackTheta = []

legHtg = ROOT.TLegend(0.7, 0.85, 0.95, 0.55, "Event class", 'NDC')

cNumOn = ROOT.TCanvas("cNumOnStack", "%s sources stacked number of ON events" % len(pCandleCatalogue.aObjectDict), 900, 600)
aaHtgNumOn = []
hsNumOn = ROOT.THStack("hsNumOn", "Number of ON events of %s stacked sources;log_{10}Energy[MeV];[events]" % len(pCandleCatalogue.aObjectDict))

cNumOff = ROOT.TCanvas("cNumOffStack", "%s sources stacked number of OFF events;log_{10}Energy[MeV];[events]" % len(pCandleCatalogue.aObjectDict), 900, 600)
aaHtgNumOff = []
hsNumOff = ROOT.THStack("hsNumOff", "Number of OFF events of %s stacked sources;log_{10}Energy[MeV];[events]" % len(pCandleCatalogue.aObjectDict))

cNumSig = ROOT.TCanvas("cNumSigStack", "%s sources stacked number of Signal events;log_{10}Energy[MeV];[events]" % len(pCandleCatalogue.aObjectDict), 900, 600)
aaHtgNumSig = []
hsNumSig = ROOT.THStack("hsNumSig", "Number of Signal events of %s stacked sources;log_{10}Energy[MeV];[events]" % len(pCandleCatalogue.aObjectDict))

cNumBkg = ROOT.TCanvas("cNumBkgStack", "%s sources stacked number of Background events;log_{10}Energy[MeV];[events]" % len(pCandleCatalogue.aObjectDict), 900, 600)
aaHtgNumBkg = []
hsNumBkg = ROOT.THStack("hsNumBkg", "Number of Background events of %s stacked sources;log_{10}Energy[MeV];[events]" % len(pCandleCatalogue.aObjectDict))

cSBratio = ROOT.TCanvas("cSBratioStack", "%s sources stacked Signal/Background ratio;log_{10}Energy[MeV];[events]" % len(pCandleCatalogue.aObjectDict), 900, 600)
cNumSig.cd()
aaHtgSBratio = []
hsSBratio = ROOT.THStack("hsSBratio", "%s sources stacked Signal/Background ratio;log_{10}Energy[MeV];[events]" % len(pCandleCatalogue.aObjectDict))

#cSgnf = ROOT.TCanvas("cSignificanceStack", "%s sources stacked Significance;log_{10}Energy[MeV];[events]" % len(pCandleCatalogue.aObjectDict), 900, 600)
# hSgnfDummy = ROOT.TH1D("hSignificanceDummy", "Significance of %s stacked sources;log_{10}Energy[MeV];[#sigma]" % len(pCandleCatalogue.aObjectDict), er.nBin, er.edgeLow, er.edgeUp)
# cNumSig.cd()
# hNumSigDummy.Draw()
#aaHtgSgnf = []

cEnergy = ROOT.TCanvas("cEnergyStack", "%s sources stacked Energy plot" % len(pCandleCatalogue.aObjectDict), 800, 500)
aaHtgEnergy = []
hsEnergy = ROOT.THStack("hsEnergy", "%s sources stacked energy plot" % len(pCandleCatalogue.aObjectDict))

print "Drawing histograms."
for pS in range(len(aaStrSelect)):
    aaaHtgTheta.append([])
    aaHtgNumOn.append([])
    aaHtgNumOff.append([])
    aaHtgNumSig.append([])
    aaHtgNumBkg.append([])
    aaHtgSBratio.append([])
#    aaHtgSgnf.append([])
    aaHtgEnergy.append([])
    for qS in range(len(aaStrSelect[pS])):
        aaaHtgTheta[-1].append([])

        aaHtgNumOn[-1].append(ROOT.TH1D("hNumOnStack%s_%s" % (pS,qS), "Number of %s ON events (%s sources stacked);log_{10}Energy[MeV];[events]" % (aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict)), erplot.nBin, erplot.edgeLow, erplot.edgeUp))
        aaHtgNumOn[-1][-1].SetFillStyle(0)
        aaHtgNumOn[-1][-1].SetLineWidth(3-pS)
        aaHtgNumOn[-1][-1].SetLineStyle(2-pS)
        aaHtgNumOn[-1][-1].SetLineColor(pColor.akColor(qS+(1-pS)*qS))
        aaHtgNumOn[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(pS,qS))
        aaHtgNumOn[-1][-1].SetMarkerSize(pColor.aakMarkerSize(pS,qS))
        aaHtgNumOn[-1][-1].SetMarkerColor(pColor.akColor(qS+(1-pS)*qS))

        aaHtgNumOff[-1].append(ROOT.TH1D("hNumOffStack%s_%s" % (pS,qS), "Number of %s OFF events (%s sources stacked);log_{10}Energy[MeV];[events]" % (aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict)), erplot.nBin, erplot.edgeLow, erplot.edgeUp))
        aaHtgNumOff[-1][-1].SetFillStyle(0)
        aaHtgNumOff[-1][-1].SetLineWidth(3-pS)
        aaHtgNumOff[-1][-1].SetLineStyle(2-pS)
        aaHtgNumOff[-1][-1].SetLineColor(pColor.akColor(qS+(1-pS)*qS))
        aaHtgNumOff[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(pS,qS))
        aaHtgNumOff[-1][-1].SetMarkerSize(pColor.aakMarkerSize(pS,qS))
        aaHtgNumOff[-1][-1].SetMarkerColor(pColor.akColor(qS+(1-pS)*qS))

        aaHtgNumSig[-1].append(ROOT.TH1D("hNumSigStack%s_%s" % (pS,qS), "Number of %s Signal events (%s sources stacked);log_{10}Energy[MeV];[events]" % (aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict)), erplot.nBin, erplot.edgeLow, erplot.edgeUp))
        aaHtgNumSig[-1][-1].SetFillStyle(0)
        aaHtgNumSig[-1][-1].SetLineWidth(3-pS)
        aaHtgNumSig[-1][-1].SetLineStyle(2-pS)
        aaHtgNumSig[-1][-1].SetLineColor(pColor.akColor(qS+(1-pS)*qS))
        aaHtgNumSig[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(pS,qS))
        aaHtgNumSig[-1][-1].SetMarkerSize(pColor.aakMarkerSize(pS,qS))
        aaHtgNumSig[-1][-1].SetMarkerColor(pColor.akColor(qS+(1-pS)*qS))

        aaHtgNumBkg[-1].append(ROOT.TH1D("hNumBkgStack%s_%s" % (pS,qS), "Number of %s Background events (%s sources stacked);log_{10}Energy[MeV];[events]" % (aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict)), erplot.nBin, erplot.edgeLow, erplot.edgeUp))
        aaHtgNumBkg[-1][-1].SetFillStyle(0)
        aaHtgNumBkg[-1][-1].SetLineWidth(3-pS)
        aaHtgNumBkg[-1][-1].SetLineStyle(2-pS)
        aaHtgNumBkg[-1][-1].SetLineColor(pColor.akColor(qS+(1-pS)*qS))
        aaHtgNumBkg[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(pS,qS))
        aaHtgNumBkg[-1][-1].SetMarkerSize(pColor.aakMarkerSize(pS,qS))
        aaHtgNumBkg[-1][-1].SetMarkerColor(pColor.akColor(qS+(1-pS)*qS))

        aaHtgSBratio[-1].append(ROOT.TH1D("hSBratioStack%s_%s" % (pS,qS), "%s Signal/Background ratio (%s sources stacked);log_{10}Energy[MeV];[events]" % (aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict)), erplot.nBin, erplot.edgeLow, erplot.edgeUp))
        aaHtgSBratio[-1][-1].SetFillStyle(0)
        aaHtgSBratio[-1][-1].SetLineWidth(3-pS)
        aaHtgSBratio[-1][-1].SetLineStyle(2-pS)
        aaHtgSBratio[-1][-1].SetLineColor(pColor.akColor(qS+(1-pS)*qS))
        aaHtgSBratio[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(pS,qS))
        aaHtgSBratio[-1][-1].SetMarkerSize(pColor.aakMarkerSize(pS,qS))
        aaHtgSBratio[-1][-1].SetMarkerColor(pColor.akColor(qS+(1-pS)*qS))

        # aaHtgSgnf[-1].append(ROOT.TH1D("hSignificanceStack%s_%s" % (pS,qS), "%s Significance (%s sources stacked);log_{10}Energy[MeV];[events]" % (aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict)), er.nBin, er.edgeLow, er.edgeUp))
        # aaHtgSgnf[-1][-1].SetFillStyle(0)
        # aaHtgSgnf[-1][-1].SetLineWidth(3-pS)
        # aaHtgSgnf[-1][-1].SetLineStyle(2-pS)
        # aaHtgSgnf[-1][-1].SetLineColor(pColor.akColor(qS+(1-pS)*qS))

        aaHtgEnergy[-1].append(ROOT.TH1D("hEnergyStack%s_%s" % (pS,qS), "%s Energy plot (%s sources stacked);log_{10}Energy[MeV];[events]" % (aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict)), erplot.nBin*5, erplot.edgeLow, erplot.edgeUp))
        aaHtgEnergy[-1][-1].SetFillStyle(0)
        aaHtgEnergy[-1][-1].SetLineWidth(2-pS)
        aaHtgEnergy[-1][-1].SetLineStyle(2-pS)
        aaHtgEnergy[-1][-1].SetLineColor(pColor.akColor(qS+(1-pS)*qS))
        aaHtgEnergy[-1][-1].SetLineColor(pColor.akColor(qS+(1-pS)*qS))
        aaHtgEnergy[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(pS,qS))
        aaHtgEnergy[-1][-1].SetMarkerSize(pColor.aakMarkerSize(pS,qS))
        aaHtgEnergy[-1][-1].SetMarkerColor(pColor.akColor(qS+(1-pS)*qS))

        print "Merging histograms."
        
        print "Saving histograms."
        fileOut.cd("Stacking")
        # aaHtgNumOn[pS][qS].Sumw2()
        # aaHtgNumOff[pS][qS].Sumw2()
        # aaHtgNumSig[pS][qS].Sumw2()
        # aaHtgNumBkg[pS][qS].Sumw2()
        for rE in range(erplot.nBin):
            #print pS,qS, er.getBin(rE)
            neOn = 0
            neOff = 0
            neSig = 0
            neBkg = 0
            neOnErrSq = 0
            neOffErrSq = 0
            neSigErrSq = 0
            neBkgErrSq = 0
            for rbj in pCandleCatalogue.aObjectDict:
                fileOut.cd(rbj["Name"])
                hNumOn = gDirectory.Get("hNumOn%s_%s_%s" % (rbj["Name"], pS, qS))
                neOn = neOn + hNumOn.GetBinContent(rE+1)
                neOnErrSq = neOnErrSq + hNumOn.GetBinError(rE+1)**2
                hNumOff = gDirectory.Get("hNumOff%s_%s_%s" % (rbj["Name"], pS, qS))
                neOff = neOff + hNumOff.GetBinContent(rE+1)
                neOffErrSq = neOffErrSq + hNumOff.GetBinError(rE+1)**2
                hNumSig = gDirectory.Get("hNumSig%s_%s_%s" % (rbj["Name"], pS, qS))
                neSig = neSig + hNumSig.GetBinContent(rE+1)
                #print "Number of signal events until", rbj["Name"], " :", neSig
                neSigErrSq = neSigErrSq + hNumSig.GetBinError(rE+1)**2
                hNumBkg = gDirectory.Get("hNumBkg%s_%s_%s" % (rbj["Name"], pS, qS))
                neBkg = neBkg + hNumBkg.GetBinContent(rE+1)
                neBkgErrSq = neBkgErrSq + hNumBkg.GetBinError(rE+1)**2
            aaHtgNumOn[pS][qS].SetBinContent(rE+1,neOn)
            aaHtgNumOn[pS][qS].SetBinError(rE+1,math.sqrt(neOnErrSq))
            aaHtgNumOff[pS][qS].SetBinContent(rE+1,neOff)
            aaHtgNumOff[pS][qS].SetBinError(rE+1,math.sqrt(neOffErrSq))
            aaHtgNumSig[pS][qS].SetBinContent(rE+1,neSig)
            #print "Total number of signal events :", neSig
            aaHtgNumSig[pS][qS].SetBinError(rE+1,math.sqrt(neSigErrSq))
            aaHtgNumBkg[pS][qS].SetBinContent(rE+1,neBkg)
            aaHtgNumBkg[pS][qS].SetBinError(rE+1,math.sqrt(neBkgErrSq))

            aaaHtgTheta[pS][qS].append(ROOT.TH1D("hThetaStack%s_%s_%s" % (pS,qS,rE), "{0} #theta^2 plot ({1} sources stacked) in {2:.1f} - {3:.1f} GeV;#theta^2 [#circ];[events]".format( aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict), 10**(erplot.getBin(rE)[0]-3), 10**(erplot.getBin(rE)[1]-3)), int((12**2)*10), 0, (12.**2)))
            aaaHtgTheta[pS][qS][-1].SetLineColor(pColor.akColor(qS+(1-pS)*qS))
            aaaHtgTheta[pS][qS][-1].SetLineWidth(2-pS)
            aaaHtgTheta[pS][qS][-1].SetLineStyle(2-pS)
            aaaHtgTheta[pS][qS][-1].Rebin(40)
            #listHtgTheta = ROOT.TList()
            for iThe in range(aaaHtgTheta[pS][qS][-1].GetNbinsX()):
                neThe = 0
                for qbj in pCandleCatalogue.aObjectDict:
                    #  print qbj["Name"]
                    fileOut.cd(qbj["Name"])
                    hThe = gDirectory.Get("hTheta%s_%s_%s_%s" % (qbj["Name"], pS, qS, rE))
                    neThe = neThe + hThe.GetBinContent(iThe+1)
                #aaaHtgTheta[pS][qS][rE].Add(gDirectory.Get("hTheta%s_%s_%s_%s" % (qbj["Name"], pS, qS, rE)))
                #listHtgTheta.Add(gDirectory.Get("hTheta%s_%s_%s_%s" % (qbj["Name"], pS, qS, rE)))
                aaaHtgTheta[pS][qS][-1].SetBinContent(iThe+1, neThe)
            #print "Merging",
            #listHtgTheta.Print()
            #aaaHtgTheta[pS][qS][rE].Merge(listHtgTheta)
            #del listHtgTheta
            fileOut.cd("Stacking")
            print "Saving."
            aaaHtgTheta[pS][qS][-1].Write("", ROOT.TObject.kOverwrite)
            print "Make THStack."
            if pS==0 and qS==0:
                aHStackTheta.append(ROOT.THStack("hsTheta%s" % rE, "{0} sources {1} stacked #theta^2 plot in {2:.1f} - {3:.1f} GeV;[#circ]^2;[counts]".format(aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict), 10**(erplot.getBin(rE)[0]-3), 10**(erplot.getBin(rE)[1]-3)) ))
            aHStackTheta[rE].Add(aaaHtgTheta[pS][qS][rE])
            if rE==0:
                legTheta.AddEntry(aaaHtgTheta[pS][qS][rE], aaStrSelect[pS][qS], 'l')
            if pS==len(aaStrSelect)-1 and qS==len(aaStrSelect[pS])-1:
                cTheta.cd(1+rE)
                aHStackTheta[rE].Draw("nostack")
                cTheta.cd(4)
                legTheta.Draw()
            # if pS==0 and qS==0:
            #     aaaHtgTheta[-1][-1][-1].Draw()
            # else:
            #     aaaHtgTheta[-1][-1][-1].Draw("same")
            fileOut.cd("Stacking")
            aHStackTheta[rE].Write("",ROOT.TObject.kOverwrite)
        aaHtgSBratio[pS][qS] = aaHtgNumSig[pS][qS].Clone("hSBratioStack%s_%s" % (pS,qS))
        aaHtgSBratio[pS][qS].SetTitle("%s Signal/Background ratio (%s sources stacked);log_{10}Energy[MeV];[events]" % (aaStrSelect[pS][qS], len(pCandleCatalogue.aObjectDict)))
        aaHtgSBratio[pS][qS].Divide(aaHtgNumBkg[pS][qS])
        aaHtgNumOn[pS][qS].Write("", ROOT.TObject.kOverwrite)
        aaHtgNumOff[pS][qS].Write("", ROOT.TObject.kOverwrite)
        aaHtgNumSig[pS][qS].Write("", ROOT.TObject.kOverwrite)
        aaHtgNumBkg[pS][qS].Write("", ROOT.TObject.kOverwrite)
        aaHtgSBratio[pS][qS].Write("", ROOT.TObject.kOverwrite)

        for sE in range(erplot.nBin*10):
            neEne = 0
            neEneErrSq = 0
            for rbj in pCandleCatalogue.aObjectDict:
                fileOut.cd(rbj["Name"])
                hEnergy = gDirectory.Get("hEnergy%s_%s_%s" % (rbj["Name"], pS, qS))
                neEne = neEne + hEnergy.GetBinContent(sE+1)
                neEneErrSq = neEneErrSq + hEnergy.GetBinError(sE+1)**2
            aaHtgEnergy[pS][qS].SetBinContent(sE+1,neEne)
            aaHtgEnergy[pS][qS].SetBinError(sE+1,math.sqrt(neEneErrSq))
        aaHtgEnergy[pS][qS].Write("", ROOT.TObject.kOverwrite)

        hsNumOn.Add(aaHtgNumOn[pS][qS])
        hsNumOff.Add(aaHtgNumOff[pS][qS])
        hsNumSig.Add(aaHtgNumSig[pS][qS])
        hsNumBkg.Add(aaHtgNumBkg[pS][qS])
        hsSBratio.Add(aaHtgSBratio[pS][qS])
        hsEnergy.Add(aaHtgEnergy[pS][qS])
        legHtg.AddEntry(aaHtgNumOn[pS][qS], aaStrSelect[pS][qS], 'lp')

cNumOn.cd()
hsNumOn.Draw("nostack E1")
legHtg.Draw('same')
cNumOff.cd()
hsNumOff.Draw("nostack E1")
legHtg.Draw('same')
cNumSig.cd()
hsNumSig.Draw("nostack E1")
legHtg.Draw('same')
cNumBkg.cd()
hsNumBkg.Draw("nostack E1")
legHtg.Draw('same')
cSBratio.cd()
hsSBratio.Draw("nostack E1")
legHtg.Draw('same')
cEnergy.cd()
hsEnergy.Draw("nostack E1")
legHtg.Draw('same')
print "Save plots."
fileOut.cd("Stacking")
cTheta.Write("", ROOT.TObject.kOverwrite)
cNumOn.Write("", ROOT.TObject.kOverwrite)
hsNumOn.Write("", ROOT.TObject.kOverwrite)
cNumOff.Write("", ROOT.TObject.kOverwrite)
hsNumOff.Write("", ROOT.TObject.kOverwrite)
cNumSig.Write("", ROOT.TObject.kOverwrite)
hsNumSig.Write("", ROOT.TObject.kOverwrite)
cNumBkg.Write("", ROOT.TObject.kOverwrite)
hsNumBkg.Write("", ROOT.TObject.kOverwrite)
cSBratio.Write("", ROOT.TObject.kOverwrite)
hsSBratio.Write("", ROOT.TObject.kOverwrite)
cEnergy.Write("", ROOT.TObject.kOverwrite)
hsEnergy.Write("", ROOT.TObject.kOverwrite)

aaGrPerformance = []
aMgrPerformance = []
cPerformance = ROOT.TCanvas("cPerformanceStack", "Performance plots with %s stacked sources" % len(pCandleCatalogue.aObjectDict), 1200, 800)
oDown = int(math.sqrt(erplot.nBin+1))
oAcross = int(math.ceil((erplot.nBin+1) / oDown))
cPerformance.Divide(oAcross, oDown)
legPerformance = ROOT.TLegend(0.1,0.1,0.9,0.9,"Event category",'NDC')
for tE in range(erplot.nBin):
    aaGrPerformance.append([])
    aMgrPerformance.append(ROOT.TMultiGraph("mgrPerformanceStack_%s" % tE, "Performance plots with {0} stacked sources in {1:.1f} - {2:.1f} GeV".format(len(pCandleCatalogue.aObjectDict), 10**(erplot.getBin(tE)[0]-3), 10**(erplot.getBin(tE)[1]-3))))
    for kAS in range(len(aStrSelect)):
        aaGrPerformance[-1].append(ROOT.TGraphErrors(len(aaStrSelect[kAS])))
        aaGrPerformance[-1][-1].SetName("grPerformanceStack_%s_%s_%s" % (tE, kAS, suffix))
        aaGrPerformance[-1][-1].SetTitle("{0} performance plot with {1} stacked sources  in {2:.1f} - {3:.1f} GeV".format(aStrSelect[kAS], len(pCandleCatalogue.aObjectDict), 10**(erplot.getBin(tE)[0]-3), 10**(erplot.getBin(tE)[1]-3)))
        aaGrPerformance[-1][-1].GetXaxis().SetTitle("Number of signal events[counts]")
        aaGrPerformance[-1][-1].GetYaxis().SetTitle("S/B ratio")
        aaGrPerformance[-1][-1].SetMarkerStyle(25-kAS*5)
        aaGrPerformance[-1][-1].SetLineStyle(2-kAS)
        aaGrPerformance[-1][-1].SetLineWidth(2)
        for kSS in range(len(aaStrSelect[kAS])):
            aaGrPerformance[tE][kAS].SetPoint(kSS, aaHtgNumSig[kAS][kSS].GetBinContent(tE+1), aaHtgSBratio[kAS][kSS].GetBinContent(tE+1))
            aaGrPerformance[tE][kAS].SetPointError(kSS, aaHtgNumSig[kAS][kSS].GetBinError(tE+1), aaHtgSBratio[kAS][kSS].GetBinError(tE+1))
        aaGrPerformance[tE][kAS].Write("", ROOT.TObject.kOverwrite)
        aMgrPerformance[tE].Add(aaGrPerformance[tE][kAS])
        if tE==0:
            legPerformance.AddEntry(aaGrPerformance[tE][kAS], aStrSelect[kAS], 'lp')
    cPerformance.cd(tE+1)
    aMgrPerformance[tE].Draw("APL")
    aMgrPerformance[tE].GetXaxis().SetTitle("Number of signal events [counts]")
    aMgrPerformance[tE].GetYaxis().SetTitle("S/B ratio")
cPerformance.cd(4)
legPerformance.Draw()
cPerformance.Write("", ROOT.TObject.kOverwrite)
print "All analysis finished."
