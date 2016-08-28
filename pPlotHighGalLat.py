#!/usr/bin/env python
"""For evaluation of CR background rejection with high-b flight data.
"""
import sys
from astropy.io import fits
from array import array
import math
import numpy as np
#import yaml
from ctypes import *
import datetime
import ROOT
ROOT.gROOT.SetBatch()
ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)
ROOT.gStyle.SetPalette(53)
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TMath, TTree, TChain, TH1, TH2, TH3
# Original modules
from pAnalysisConfig import *
import pTransientCatalogue
from pTarget import *
from pColor import *
from pFindHEALPix import *
import healpy as hp
from healpy import pixelfunc as hppf


par = sys.argv
print par
if len(par)<2:
    print "pPlotHighGalLat.py [nameVarBDT] [suffix] [fileCutBDT] [list of tree file]"
    sys.exit(0)
nameVarBDT = par[1]
strSuffix = par[2]
fileCutBDT = ROOT.TFile(par[3], 'READ')
aPathFileData = par[4:]

chData = ROOT.TChain('MeritTuple')
chBDT = ROOT.TChain('MeritTuple')

h3 = ROOT.TH3D('h3', 'High Galactic latitude events', 7, 4.35, 5.75, 8, 0.2, 1.0, 2000, 0, 5)

bdtR = c_float() 

# Region setup
NHPSIDE_OFF = 16
pathCatalogue = "/disk/gamma/cta/store/takhsm/FermiData/catalogue/gll_psch_v09.fit"
aHpxOFF = find_galoff_healpxs(NHPSIDE_OFF, 0, pathCatalogue)


for pathFileData in aPathFileData:
#    chData.Add(pathFileData)
    fileData = ROOT.TFile(pathFileData, 'READ')
    print fileData.GetName()
    trData = fileData.Get('MeritTuple')
    pathFileBDT = pathFileData[:-5] + '_' + nameVarBDT + '.root'
    #chBDT.Add(pathFileBDT)
    trData.AddFriend('{0} = MeritTuple'.format(nameVarBDT), pathFileBDT)
#    trData.AddFriend(chBDT)
    nEvt = trData.GetEntries()
    print "This file has", nEvt, "events"
    trData.SetBranchAddress(nameVarBDT,bdtR)

    for iEvt in range(nEvt):
        trData.GetEntry(iEvt)
        npix = hppf.ang2pix(NHPSIDE_OFF, math.pi/2.-math.radians(trData.FT1CalB), math.radians(trData.FT1CalL))
        if (npix in aHpxOFF) and trData.FT1CalZenithTheta<90 and trData.Cal1RawEnergySum>=20000 and (trData.TkrNumTracks==0 or (math.log10(max(trData.CalTrackAngle,1E-4)) > (0.529795)*(trData.EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(math.pow((trData.EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(math.pow((trData.EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(math.pow((trData.EvtJointLogEnergy-3.000000)/0.916667,3))))*(trData.EvtJointLogEnergy >= 3.000000 and trData.EvtJointLogEnergy <= 5.750000) + (-0.398962)*(trData.EvtJointLogEnergy >  5.750000))) and trData.FswGamState==0: #and (trData.GltGemSummary&0x20)==0):
            if trData.Cal1MomZCrossSide840>=0.0 and trData.Cal1MomNumIterations>0 and trData.Cal1TransRms>=10 and trData.Cal1TransRms<70 and trData.Acd2Cal1VetoSigmaHit>0 and trData.Cal1MomZDir>=0.2: #and (trData.WP8CalOnlyBEPCaseE_myBDT+1.)/2.>=0.06:
                h3.Fill(trData.EvtJointLogEnergy, trData.Cal1MomZDir, -math.log10(1.0-(1.0+bdtR.value)/2.0))
            else:
                h3.Fill(trData.EvtJointLogEnergy, trData.Cal1MomZDir, 0.0)
    print "Cumulative event number for plot:", h3.GetEntries(), "events"

fileOut = ROOT.TFile('Htg_rejection_GalOff_{0}_{1}.root'.format(nameVarBDT, strSuffix), 'UPDATE')
fileOut.cd()
h3.Write()

h2Raw = h3.Project3D("yx")
h2Raw.SetName("h2Raw")
h2Raw.SetTitle("Raw event distribution in high-b area")
h2Raw.Write()

# ----- Event class setup -----
cfg = ClassConfig('Both', [10, 3, 1], 1)
aCutEGB = cfg.aCutEGB
aaStrSelect = cfg.aaStrSelect

aHtgCut = []
aHtgRej = []
for ict in range(len(aaStrSelect[1])):
    print aaStrSelect[1][ict]
    htgCutBDT = fileCutBDT.Get('htgCutBDT{0}'.format(ict))
#    aHtgCut.append(h2Raw.Clone('h2Residusl{0}'.format(ict)))
    aHtgCut.append(ROOT.TH2D('h2Residusl{0}'.format(ict), 'Residual event distribution in {0}'.format(aaStrSelect[1][ict]), h2Raw.GetNbinsX(), h2Raw.GetXaxis().GetBinLowEdge(1), h2Raw.GetXaxis().GetBinUpEdge(h2Raw.GetNbinsX()), h2Raw.GetNbinsY(), h2Raw.GetYaxis().GetBinLowEdge(1), h2Raw.GetYaxis().GetBinUpEdge(h2Raw.GetNbinsY())))
#    aHtgCut[-1].SetTitle('Residual event distribution in {0}'.format(aaStrSelect[1][ict]))
    for ix in range(aHtgCut[-1].GetNbinsX()):
        for iy in range(aHtgCut[-1].GetNbinsY()):
            vCut = htgCutBDT.GetBinContent(htgCutBDT.GetXaxis().FindBin(aHtgCut[-1].GetXaxis().GetBinCenter(ix+1)), htgCutBDT.GetYaxis().FindBin(aHtgCut[-1].GetYaxis().GetBinCenter(iy+1)))
            print "BDT cut:", vCut
            aHtgCut[-1].SetBinContent(ix+1, iy+1, h3.Integral(ix+1, ix+1, iy+1, iy+1, h3.GetZaxis().FindBin(vCut), h3.GetNbinsZ()))
    aHtgCut[-1].Write()
    aHtgRej.append(aHtgCut[-1].Clone('h2Rejection{0}'.format(ict)))
    aHtgRej[-1].SetTitle('Rejection rate in {0}'.format(aaStrSelect[1][ict]))
    aHtgRej[-1].Divide(h2Raw)
    aHtgRej[-1].Write()
