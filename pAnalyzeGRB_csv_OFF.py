#!/usr/bin/env python

import sys
import ROOT
ROOT.gROOT.SetBatch()
from ROOT import TTree, TChain
import numpy as np
import yaml
import datetime
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
import pandas
from pMETandMJD import *
from pColor import *

ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)

from pAnalysisConfig import *

# ----- Event class setup -----
par = sys.argv
cfg = ClassConfig('Both', [10, 3, 1], 1)
aCutEGB = cfg.aCutEGB
aaStrSelect = cfg.aaStrSelect
nStartBin = cfg.nStartBin

nameFileRoc = "/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZCS000wwoTRKwoMCZDIR00woRW_15_S11D_catTwoZDIR050Log_roc.root" #par[2]
nameVarBDT = "S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_BDTG1000D06"
nameFileSuffix = par[1]
cutMVA = CutBDT(nameFileRoc, aCutEGB)
aaValCutBDT = cutMVA.aaValCutBDT[0:]
print aaValCutBDT
nEnergyBin = cutMVA.aaValCutBDT[0]['numBin'] - nStartBin
vEnergyBinWidth = cutMVA.aaValCutBDT[0]['widthBin']
vEnergyLow = cutMVA.aaValCutBDT[0]['edgeLow'] + nStartBin*vEnergyBinWidth
vEnergyUp = vEnergyLow + nEnergyBin*vEnergyBinWidth
aaNumEventClass=[]
for hS in range(len(aaStrSelect)):
    aaNumEventClass.append([])
    for iS in range(len(aaStrSelect[hS])):
        aaNumEventClass[hS].append(0)

#IRF
listPathFilePerf = [['/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_TRANSIENT100_P8R2_TRANSIENT100_perf.root', '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_SOURCE_P8R2_SOURCE_perf.root'], 
                    ['/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_CalOnly_R100_perf.root', '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_CalOnly_R30_perf.root', '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_CalOnly_R10_perf.root']]
htgPerf = CutPerformanceHtg(listPathFilePerf)
aSepCut = [68, 95]
zenithCut = 90

# Data
pathList = "/nfs/farm/g/glast/u/mtakahas/data/lists/GBMcatalogue_FluenceLt1Micro_20160514.csv"
csv = pandas.read_csv(pathList)
num_lines = sum(1 for line in open(pathList))

# OFF regions
nOff = 0 #4;
degOffOffset = 14.0

print "===================="
# Making all sky map
listFileIn = par[2]
trGammas = ROOT.TChain("trGammas", "Gamma-like events")
trGammas.Add(listFileIn)
nEvent = trGammas.GetEntries()
if not nEvent>0:
    print "No events."
    exit
print "Total number of events:", nEvent
tFirst = trGammas.GetMinimum("t")
tLast = trGammas.GetMaximum("t")
tTotalTime = tLast - tFirst
listGrb = ['100325275', '160509374', '081009140']
listPhCosTh = [0.956564320626, 0.985249156718, 0.567050376683]
listPhEnergy = [5.106870174, 5.063819408, 4.526272297]
for iGrb in range(len(listGrb)):
    nameGrb = listGrb[iGrb]
    costhPh = listPhCosTh[iGrb]
    enePh = listPhEnergy[iGrb]
    cEvent = [ ROOT.TCanvas("cEvent68", "Gamma-like events within 68% from the GRB{0}".format(nameGrb)), ROOT.TCanvas("cEvent95", "Gamma-like events within 95% from the GRB{0}".format(nameGrb))]
    cZenith = [ ROOT.TCanvas("cZenith68", "Zenith angle of ON/OFF events within 68% from the GRB{0}".format(nameGrb)), ROOT.TCanvas("cZenith95", "Zenith angle of ON/OFF events within 95% from the GRB{0}".format(nameGrb))]
    print ""
    nameFileOut = "GRB" + nameGrb + "_" + nameFileSuffix + ".root"
    fileOut = ROOT.TFile(nameFileOut, 'UPDATE')

    #------ Source data -----
    for jGrb in range(num_lines-1):
        if int(nameGrb) == int(csv.ix[jGrb,'name']): #grb.findtext("./GRBNAME")==nameGrb: 
            raSrc = float(csv.ix[jGrb,'ra'])
            decSrc = float(csv.ix[jGrb,'dec']) 
            trigger_time = ConvertMjdToMet(float(csv.ix[jGrb,'trigger_time']))
            err_rad = float(csv.ix[jGrb,'error_radius'])
    print ""
    print "==============="
    print "GRB", nameGrb
    print "==============="
    print "(", raSrc, ",", decSrc, "), Error radius:", err_rad, "Trigger MET:", trigger_time 

    # Plot
    hEvent = []
    hZenith = []
    hZenithTime = []
    mgr = []
    greOn = []
    greOff=[]
    mgrZenith = []
    grZenith = []
    for vSepCut in aSepCut:
        mgr.append(ROOT.TMultiGraph("mgrPSF{0}".format(vSepCut), "Gamma-like events within PSF{0}% from the GRB{1}".format(vSepCut, nameGrb)))
        greOn.append([])
        greOff.append([])
        hEvent.append([])
        hZenith.append([])
        hZenithTime.append([])
        for pC in range(len(aaStrSelect)):
            greOn[-1].append([])
            greOff[-1].append([])
            hEvent[-1].append([])
            hZenith[-1].append([])
            hZenithTime[-1].append([])
            for qC in range(len(aaStrSelect[pC])):
                hEvent[-1][-1].append(ROOT.TH2D("hEvent{0}_{1}_{2}".format(vSepCut, pC, qC), "{0} ON PSF{1}%".format(aaStrSelect[pC][qC], vSepCut), (int)(tTotalTime/1000), tFirst, tLast, 7, 4.35, 5.75))
                hZenith[-1][-1].append(ROOT.TH2D("hZenith{0}_{1}_{2}".format(vSepCut, pC, qC), "{0} ON PSF{1}%".format(aaStrSelect[pC][qC], vSepCut), 7, 4.35, 5.75, 140, 0, 140))
                hZenithTime[-1][-1].append(ROOT.TH2D("hZenithTime{0}_{1}_{2}".format(vSepCut, pC, qC), "{0} ON PSF{1}%".format(aaStrSelect[pC][qC], vSepCut), (int)(tTotalTime/1000), tFirst, tLast, 140, 0, 140))
                greOn[-1][-1].append(ROOT.TGraphErrors())
                greOn[-1][-1][-1].SetName("greOn_{0}_{1}".format(pC, qC))
                greOn[-1][-1][-1].SetTitle("{0} ON".format(aaStrSelect[pC][qC]))
                greOn[-1][-1][-1].SetMarkerStyle(20)
                if pC==0:
                    greOn[-1][-1][-1].SetMarkerColor(13-12*qC)
                elif pC==1:
                    greOn[-1][-1][-1].SetMarkerColor(kRed+3*(qC-2))
                if pC==1:
                    mgr[-1].Add(greOn[-1][-1][-1])
                greOff[-1][-1].append([])
                for hRegio in range(nOff):
                    greOff[-1][-1][-1].append(ROOT.TGraphErrors())
                    greOff[-1][-1][-1][-1].SetName("greOff_{0}_{1}_{2}".format(pC, qC, hRegio+1))
                    greOff[-1][-1][-1][-1].SetTitle("{0} Off{1} events".format(aaStrSelect[pC][qC], hRegio+1))
                    if pC==0:
                        greOff[-1][-1][-1][-1].SetMarkerColor(13-12*qC)
                    elif pC==1:
                        greOff[-1][-1][-1][-1].SetMarkerColor(kRed+3*(qC-2))
                    greOff[-1][-1][-1][-1].SetMarkerStyle(25+hRegio)
                    if pC==1:
                        mgr[-1].Add(greOff[-1][-1][-1][-1])
        mgrZenith.append(ROOT.TMultiGraph("mgrZenith{0}".format(vSepCut), "Zenith angle within ON/OFF regions of PSF {0}%".format(vSepCut)))
        grZenith.append([])
        for gRegio in range(nOff+1):
            grZenith[-1].append(ROOT.TGraph())
            grZenith[-1][-1].SetName("grZenith{0}_{1}".format(vSepCut, gRegio))
            if gRegio==0:
                grZenith[-1][0].SetTitle("ON")
            else:
                grZenith[-1][gRegio].SetTitle("OFF{0}".format(gRegio))
                grZenith[-1][gRegio].SetMarkerStyle(7)
                grZenith[-1][gRegio].SetMarkerColor(akColor(gRegio))
            mgrZenith[-1].Add(grZenith[-1][-1])

    vecTgt = []
    vecTgt.append(np.array([cos(radians(decSrc))*cos(radians(raSrc)), cos(radians(decSrc))*sin(radians(raSrc)), sin(radians(decSrc))]))
    vecTgt.append(np.array([cos(radians(decSrc-degOffOffset))*cos(radians(raSrc)), cos(radians(decSrc-degOffOffset))*sin(radians(raSrc)), sin(radians(decSrc-degOffOffset))]))
    vecTgt.append(np.array([cos(radians(decSrc))*cos(radians(raSrc-degOffOffset/cos(radians(decSrc)))), cos(radians(decSrc))*sin(radians(raSrc-degOffOffset/cos(radians(decSrc)))), sin(radians(decSrc))]))
    vecTgt.append(np.array([cos(radians(decSrc+degOffOffset))*cos(radians(raSrc)), cos(radians(decSrc+degOffOffset))*sin(radians(raSrc)), sin(radians(decSrc+degOffOffset))]))
    vecTgt.append(np.array([cos(radians(decSrc))*cos(radians(raSrc+degOffOffset/cos(radians(decSrc)))), cos(radians(decSrc))*sin(radians(raSrc+degOffOffset/cos(radians(decSrc)))), sin(radians(decSrc))]))

    for iEvent in range(nEvent):
        trGammas.GetEntry(iEvent)
        vecEvt = np.array([cos(radians(trGammas.dec))*cos(radians(trGammas.ra)), cos(radians(trGammas.dec))*sin(radians(trGammas.ra)), sin(radians(trGammas.dec))])
        aDist = []
        distSepCut = []
        distSepCut.append(htgPerf.getPSF68_cth(trGammas.c-1, 0, trGammas.e, costhPh) + err_rad)
        distSepCut.append(htgPerf.getPSF95_cth(trGammas.c-1, 0, trGammas.e, costhPh) + err_rad)
        for iRegio in range(1+nOff):
            radTheta = acos(np.dot(vecTgt[iRegio], vecEvt))
            aDist.append(degrees(radTheta))
            for kCut in range(len(aSepCut)):
                if aDist[iRegio] < distSepCut[kCut]:
                    grZenith[kCut][iRegio].SetPoint(grZenith[kCut][iRegio].GetN(), trGammas.t, trGammas.z)
                    if iRegio==0:
                        if trGammas.c == 1:
                            if trGammas.s == 4:
                                hZenith[kCut][0][0].Fill(trGammas.e, trGammas.z)
                                hZenithTime[kCut][0][0].Fill(trGammas.t, trGammas.z)
                                if trGammas.z<zenithCut:
                                    hEvent[kCut][0][0].Fill(trGammas.t, trGammas.e)
                                    greOn[kCut][0][0].SetPoint(greOn[kCut][0][0].GetN(), trGammas.t, pow(10, trGammas.e-3))
                            elif trGammas.s == 128:
                                hZenith[kCut][0][0].Fill(trGammas.e, trGammas.z)
                                hZenithTime[kCut][0][0].Fill(trGammas.t, trGammas.z)
                                hZenith[kCut][0][1].Fill(trGammas.e, trGammas.z)
                                hZenithTime[kCut][0][1].Fill(trGammas.t, trGammas.z)
                                if trGammas.z<zenithCut:
                                    hEvent[kCut][0][0].Fill(trGammas.t, trGammas.e)
                                    hEvent[kCut][0][1].Fill(trGammas.t, trGammas.e)
                                    greOn[kCut][0][1].SetPoint(greOn[kCut][0][1].GetN(), trGammas.t, pow(10, trGammas.e-3))
                        elif trGammas.c == 2:
                            if trGammas.s == 4096:
                                hZenith[kCut][1][0].Fill(trGammas.e, trGammas.z)
                                hZenithTime[kCut][1][0].Fill(trGammas.t, trGammas.z)
                                if trGammas.z<zenithCut:
                                    hEvent[kCut][1][0].Fill(trGammas.t, trGammas.e)
                                    greOn[kCut][1][0].SetPoint(greOn[kCut][1][0].GetN(), trGammas.t, pow(10, trGammas.e-3))
                            elif trGammas.s == 8192:
                                hZenith[kCut][1][0].Fill(trGammas.e, trGammas.z)
                                hZenithTime[kCut][1][0].Fill(trGammas.t, trGammas.z)
                                hZenith[kCut][1][1].Fill(trGammas.e, trGammas.z)
                                hZenithTime[kCut][1][1].Fill(trGammas.t, trGammas.z)
                                if trGammas.z<zenithCut:
                                    hEvent[kCut][1][0].Fill(trGammas.t, trGammas.e)
                                    hEvent[kCut][1][1].Fill(trGammas.t, trGammas.e)
                                    greOn[kCut][1][1].SetPoint(greOn[kCut][1][1].GetN(), trGammas.t, pow(10, trGammas.e-3))
                            elif trGammas.s == 16384:
                                hZenith[kCut][1][0].Fill(trGammas.e, trGammas.z)
                                hZenithTime[kCut][1][0].Fill(trGammas.t, trGammas.z)
                                hZenith[kCut][1][1].Fill(trGammas.e, trGammas.z)
                                hZenithTime[kCut][1][1].Fill(trGammas.t, trGammas.z)
                                hZenith[kCut][1][2].Fill(trGammas.e, trGammas.z)
                                hZenithTime[kCut][1][2].Fill(trGammas.t, trGammas.z)
                                if trGammas.z<zenithCut:
                                    hEvent[kCut][1][0].Fill(trGammas.t, trGammas.e)
                                    hEvent[kCut][1][1].Fill(trGammas.t, trGammas.e)
                                    hEvent[kCut][1][2].Fill(trGammas.t, trGammas.e)
                                    greOn[kCut][1][2].SetPoint(greOn[kCut][1][2].GetN(), trGammas.t, pow(10, trGammas.e-3))
                    elif iRegio>0 and trGammas.z<zenithCut:
                        if trGammas.s>0:
                            if trGammas.c == 1:
                                if trGammas.s == 4:
                                    greOff[kCut][0][0][iRegio-1].SetPoint(greOff[kCut][0][0][iRegio-1].GetN(), trGammas.t, pow(10, trGammas.e-3))
                                elif trGammas.s == 128:
                                    greOff[kCut][0][1][iRegio-1].SetPoint(greOff[kCut][0][1][iRegio-1].GetN(), trGammas.t, pow(10, trGammas.e-3))
                            elif trGammas.c == 2:
                                if trGammas.s == 4096:
                                    greOff[kCut][1][0][iRegio-1].SetPoint(greOff[kCut][1][0][iRegio-1].GetN(), trGammas.t, pow(10, trGammas.e-3))
                                elif trGammas.s == 8192:
                                    greOff[kCut][1][1][iRegio-1].SetPoint(greOff[kCut][1][1][iRegio-1].GetN(), trGammas.t, pow(10, trGammas.e-3))
                                elif trGammas.s == 16384:
                                    greOff[kCut][1][2][iRegio-1].SetPoint(greOff[kCut][1][2][iRegio-1].GetN(), trGammas.t, pow(10, trGammas.e-3))
    fileOut.cd()
    leg = []
    legZenith = []
    for jSepCut in range(len(aSepCut)):
        cEvent[jSepCut].cd()
        mgr[jSepCut].Draw("AP")
        mgr[jSepCut].GetXaxis().SetTitle("Time [s]")
        mgr[jSepCut].GetYaxis().SetTitle("Energy [GeV]")
        mgr[jSepCut].GetYaxis().SetRangeUser(10, 10000)
        cEvent[jSepCut].SetLogy()
        leg.append(ROOT.TLegend(0.67, 0.5, 0.88, 0.88))
        for pD in range(len(aaStrSelect)):
            for qD in range(len(aaStrSelect[pD])):
                hEvent[jSepCut][pD][qD].Write()
                hZenith[jSepCut][pD][qD].Write()
                hZenithTime[jSepCut][pD][qD].Write()
                if pD==1:
                    leg[jSepCut].AddEntry(greOn[jSepCut][pD][qD], greOn[jSepCut][pD][qD].GetTitle(), "p")
        for rr in range(nOff):
            leg[jSepCut].AddEntry(greOff[jSepCut][1][0][rr], "OFF{0}".format(rr+1), "p")
        leg[jSepCut].Draw("same")
        cZenith[jSepCut].cd()
        mgrZenith[jSepCut].Draw("AP")
        legZenith.append(ROOT.TLegend(0.67, 0.5, 0.88, 0.88))
        for rz in range(nOff+1):
            legZenith[jSepCut].AddEntry(grZenith[jSepCut][rz], grZenith[jSepCut][rz].GetTitle(), "p")
        legZenith[jSepCut].Draw("same")
        print ""
        fileOut.cd()
        cEvent[jSepCut].Write()
        cZenith[jSepCut].Write()
    print "Analysis for GRB", nameGrb, "finished!"
