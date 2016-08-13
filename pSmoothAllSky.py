#!/usr/bin/env python

import sys
from array import array
import math
import datetime
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
# Original modules
from pAnalysisConfig import *
import pCandleCatalogue
import pTransientCatalogue
from pTarget import *
from pColor import *
par = sys.argv
print par
if len(par)<2:
    print "pSmoothAllSky.py [suffix] [map file]"
    sys.exit(0)
suffix = par[1]
pathFile = par[2]
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()

ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)
ROOT.gStyle.SetPalette(53)

cfg = ClassConfig('Both', [10, 3, 1])
er = EnergyLogRegion(4, 4.75, 0.25)
aaStrSelect = cfg.aaStrSelect
aStrSelect = cfg.aStrSelect

htgPerf = CutPerformanceHtg(listPathFilePerf)

file = ROOT.TFile(pathFile, "UPDATE")
print "In/Output file:", file.GetName()

file.cd("GalacticRidge")
aaaHtgMapAllGal = []
aaaHtgMapAllCel = []
aaaHtgMapAllGalEx = []
aaaHtgMapAllCelEx = []

for pS in range(len(aaStrSelect)):
    aaaHtgMapGal.append([])
    aaaHtgMapCel.append([])
    aaaHtgMapGalEx.append([])
    aaaHtgMapCelEx.append([])
    for qS in range(len(aaStrSelect[pS])):
        aaaHtgMapGal[-1].append([])
        aaaHtgMapCel[-1].append([])
        aaaHtgMapGaExl[-1].append([])
        aaaHtgMapCelEx[-1].append([])
        for rE in range(er.nBin):
            aaaHtgMapGal[-1][-1].append(gDirectory.Get("hMapAllSky%s_%s_%s" %(pS, qS, rE)))
            aaaHtgMapCel[-1][-1].append(gDirectory.Get("hMapAllCel%s_%s_%s" %(pS, qS, rE)))
            aaaHtgMapAllGalEx[-1][-1].append(aaaHtgMapGal[-1][-1][-1].Clone("hMapAllSkyEx%s_%s_%s" %(pS, qS, rE)))
            aaaHtgMapAllCelEx[-1][-1].append(aaaHtgMapCel[-1][-1][-1].Clone("hMapAllCelEx%s_%s_%s" %(pS, qS, rE)))
            for igy in range(1, aaaHtgMapAllGal[-1][-1][-1].GetNbinsY()):
                for igx in range(1, aaaHtgMapAllGal[-1][-1][-1].GetXaxis().FindBin(-180)):
                    aaaHtgMapAllGalEx[-1][-1][-1].SetBinContent(igx, igy, aaaHtgMapGal[-1][-1][-1].GetBinContent(aaaHtgMapGal[-1][-1][-1].FindBin(aaaHtgMapGal[-1][-1][-1].GetBinCenter(igx)+360), igy))
                for igx in range(aaaHtgMapAllGal[-1][-1][-1].GetXaxis().FindBin(180)+1, aaaHtgMapAllGal[-1][-1][-1].GetNbinsX()+1):
                    aaaHtgMapAllGalEx[-1][-1][-1].SetBinContent(igx, igy, aaaHtgMapGal[-1][-1][-1].GetBinContent(aaaHtgMapGal[-1][-1][-1].FindBin(aaaHtgMapGal[-1][-1][-1].GetBinCenter(igx)-360), igy))
            for icy in range(1, aaaHtgMapAllCel[-1][-1][-1].GetNbinsY()):
                for icx in range(1, aaaHtgMapAllCel[-1][-1][-1].GetXaxis().FindBin(-180)):
                    aaaHtgMapAllCelEx[-1][-1][-1].SetBinContent(icx, icy, aaaHtgMapCel[-1][-1][-1].GetBinContent(aaaHtgMapCel[-1][-1][-1].FindBin(aaaHtgMapCel[-1][-1][-1].GetBinCenter(icx)+360), icy))
                for icx in range(aaaHtgMapAllCel[-1][-1][-1].GetXaxis().FindBin(180)+1, aaaHtgMapAllCel[-1][-1][-1].GetNbinsX()+1):
                    aaaHtgMapAllCelEx[-1][-1][-1].SetBinContent(icx, icy, aaaHtgMapCel[-1][-1][-1].GetBinContent(aaaHtgMapCel[-1][-1][-1].FindBin(aaaHtgMapCel[-1][-1][-1].GetBinCenter(icx)-360), icy))

                    #for lB in range(self.aaaHtgMapAllSky[kAS][kSS][tE].GetNbinsY()):
                    #    bPoint = self.aaaHtgMapAllSky[kAS][kSS][tE].GetYais().GetBinCenter(lB+1)
                    #    bRoiLow = bPoint-degPSF*4
                    #    if bRoiLow<self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetYaxis().GetBinLowEdge(1):
                    #        bRoiLow = self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetYaxis().GetBinUpEdge(self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetNbinsX()) + lRoiLow
                    #         lRoiUp = lPoint + degPSF*4
                    #         if lRoiUp>=self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetXaxis().GetBinUpEdge(self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetNbinsX()):
                    #             lRoiUp = self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetXaxis().GetBinLowEdge(1) + lRoiUp
                    #     for lL in range(self.aaaHtgMapAllSky[kAS][kSS][tE].GetNbinsX()):
                    #         lPoint = self.aaaHtgMapAllSky[kAS][kSS][tE].GetXais().GetBinCenter(lL+1)
                    #         lRoiLow = lPoint-degPSF*4
                    #         if lRoiLow<self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetXaxis().GetBinLowEdge(1):
                    #             lRoiLow = self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetXaxis().GetBinUpEdge(self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetNbinsX()) + lRoiLow
                    #         lRoiUp = lPoint + degPSF*4
                    #         if lRoiUp>=self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetXaxis().GetBinUpEdge(self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetNbinsX()):
                    #             lRoiUp = self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetXaxis().GetBinLowEdge(1) + lRoiUp

                    #         nCount = self.aaaHtgMapAllSky[kAS][kSS][tE].GetBinContent(lL+1, lB+1)
                    #         if nCount>0:
                    #             bPoint = self.aaaHtgMapAllSky[kAS][kSS][tE].GetXais().GetBinCenter(lB+1)
                    #             for mL in range(self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetNbinsX()):
                    #                 # for mL in range(self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetNbinsX()):
                    #                 for mB in range(self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetNbinsY()):
                    #                     degDist = pCommon.anglePointsDegToDeg(self.aaaHtgMapAllSky[kAS][kSS][tE].GetXaxis().GetBinCenter(lL+1), self.aaaHtgMapAllSky[kAS][kSS][tE].GetYaxis().GetBinCenter(lB+1), self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetXaxis().GetBinCenter(mL+1), self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetYaxis().GetBinCenter(mB+1))
                    #                     self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].Fill(mL+1, mB+1, fPSF.Eval(degDist) * cos(radians(90+self.aaaHtgMapAllSky[kAS][kSS][tE].GetYaxis().GetBinLowEdge(mB+1)))-cos(radians(90+self.aaaHtgMapAllSky[kAS][kSS][tE].GetYaxis().GetBinLowEdge(mB+2))) / self.aaaHtgMapAllSky[kAS][kSS][tE].GetNbinsX())

