import sys
import ROOT
import os
from astropy.io import fits
from array import array
import numpy as np
from pAnalysisConfig import *
import pColor
import pCommon
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import healpy as hp
from healpy import pixelfunc as hppf
from pMETandMJD import *
from math import cos, sin, tan, acos, asin, atan, radians, degrees
from ROOT import kWhite, kBlack, kGray, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TH1, TH2, TH3
from pHealplot import Healcube


class Target:
    def __init__(self, strName, config = ClassConfig(), eRegion=EnergyLogRegion(3,4.75,0.25), ePlotRegion=EnergyLogRegion(3, 4.75, 0.25), tStart=252460800.0, tEnd=504921604.0, nhpside=64):
        self.name = strName
        self.aaClass = config.aaStrSelect
        self.aClass = config.aStrSelect
        self.NSIDE = nhpside
        self.NPIX = hppf.nside2npix(self.NSIDE)
        self.npix_valid = self.NPIX
        self.tpl_pix_valid = tuple(range(self.npix_valid))
        self.HEALPLOT_MOLLWEIDE = False
        self.loncntr = 0.0
        self.latcntr = 0.0
        self.coord = 'E'
        self.energySet = eRegion
        self.energyPlot = ePlotRegion
        self.costhetaPlot = EnergyLogRegion(10, 0.0, 0.1)
        self.saOn = [1.,1.]
        self.saOff = [1.,1.]
        self.factorOn = 1.5
        self.ulOn = 4.5
        self.zCut = config.zCut
        self.tStart = tStart
        self.tEnd = tEnd
        self.tDuration = tEnd - tStart
        # numpy for FITS
        self.npaaENERGY = []
        self.npaaRA = []
        self.npaaDEC = []
        self.npaaL = []
        self.npaaB = []
        self.npaaTHETA = []
        self.npaaZENITH_ANGLE = []
        self.npaaTIME = []
        self.npaaEVENT_CLASS = []
        self.npaaEVENT_TYPE = []
        for aCat in self.aaClass:
            self.npaaENERGY.append([])
            self.npaaRA.append([])
            self.npaaDEC.append([])
            self.npaaL.append([])
            self.npaaB.append([])
            self.npaaTHETA.append([])
            self.npaaZENITH_ANGLE.append([])
            self.npaaTIME.append([])
            self.npaaEVENT_CLASS.append([])
            self.npaaEVENT_TYPE.append([])
            for aCla in aCat:
                self.npaaENERGY[-1].append(np.empty(0, dtype=np.float))
                self.npaaRA[-1].append(np.empty(0, dtype=np.float))
                self.npaaDEC[-1].append(np.empty(0, dtype=np.float))
                self.npaaL[-1].append(np.empty(0, dtype=np.float))
                self.npaaB[-1].append(np.empty(0, dtype=np.float))
                self.npaaTHETA[-1].append(np.empty(0, dtype=np.float))
                self.npaaZENITH_ANGLE[-1].append(np.empty(0, dtype=np.float))
                self.npaaTIME[-1].append(np.empty(0, dtype=np.double))
                self.npaaEVENT_CLASS[-1].append(np.empty(0, dtype=np.int32))
                self.npaaEVENT_TYPE[-1].append(np.empty(0, dtype=np.int32))

        self.trObj = ROOT.TTree('trFriend{0}'.format(self.name), 'Friend tree for {0}'.format(self.name))
        self.flagOn = np.zeros(1, dtype=bool)
        self.flagOff = np.zeros(1, dtype=bool)
        self.trObj.Branch('FLAG_ON', self.flagOn,'FLAG_ON/O')
        self.trObj.Branch('FLAG_OFF', self.flagOff,'FLAG_OFF/O')

        print "*", self.name
        print "Period (MET):", self.tStart, " - ", self.tEnd
        aaHtgNumOn = []
        aaHtgNumOff = []
        aaHtgNumSig = []
        aaHtgNumBkg = []
        aaHtgSBratio = []
        aaHtgSgnf = []
        aaHtgEnergy = []
        aaHtgLCCount = []
        aaHtgEnergyOn = []
        aaHtgLCCountOn = []
        aaHtgEnergyOff = []
        aaHtgLCCountOff = []
        for iS in range(len(self.aaClass)):
            aaHtgNumOn.append([])
            aaHtgNumOff.append([])
            aaHtgNumSig.append([])
            aaHtgNumBkg.append([])
            aaHtgSBratio.append([])
            aaHtgSgnf.append([])
            aaHtgEnergy.append([])
            aaHtgLCCount.append([])
            aaHtgEnergyOn.append([])
            aaHtgLCCountOn.append([])
            aaHtgEnergyOff.append([])
            aaHtgLCCountOff.append([])
            for jS in range(len(self.aaClass[iS])):
                aaHtgNumOn[-1].append(ROOT.TH1D("hNumOn%s_%s_%s" % (self.name,iS,jS), "Number of %s ON events of %s;log_{10}Energy[MeV];[counts]" % (self.aaClass[iS][jS], self.name), self.energyPlot.nBin, self.energyPlot.edgeLow, self.energyPlot.edgeUp))
                aaHtgNumOn[-1][-1].SetFillStyle(0)
                aaHtgNumOn[-1][-1].SetLineWidth(3-iS)
                aaHtgNumOn[-1][-1].SetLineStyle(2-iS)
                aaHtgNumOn[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumOn[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgNumOn[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgNumOn[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumOff[-1].append(ROOT.TH1D("hNumOff%s_%s_%s" % (self.name,iS,jS), "Number of %s OFF events of %s;log_{10}Energy[MeV];[counts]" % (self.aaClass[iS][jS], self.name), self.energyPlot.nBin, self.energyPlot.edgeLow, self.energyPlot.edgeUp))
                aaHtgNumOff[-1][-1].SetFillStyle(0)
                aaHtgNumOff[-1][-1].SetLineWidth(3-iS)
                aaHtgNumOff[-1][-1].SetLineStyle(2-iS)
                aaHtgNumOff[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumOff[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgNumOff[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgNumOff[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumSig[-1].append(ROOT.TH1D("hNumSig%s_%s_%s" % (self.name,iS,jS), "Number of %s Signal events of %s;log_{10}Energy[MeV];[counts]" % (self.aaClass[iS][jS], self.name), self.energyPlot.nBin, self.energyPlot.edgeLow, self.energyPlot.edgeUp))
                aaHtgNumSig[-1][-1].SetFillStyle(0)
                aaHtgNumSig[-1][-1].SetLineWidth(3-iS)
                aaHtgNumSig[-1][-1].SetLineStyle(2-iS)
                aaHtgNumSig[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumSig[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgNumSig[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgNumSig[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumBkg[-1].append(ROOT.TH1D("hNumBkg%s_%s_%s" % (self.name,iS,jS), "Number of %s Background events of %s;log_{10}Energy[MeV];[counts]" % (self.aaClass[iS][jS], self.name), self.energyPlot.nBin, self.energyPlot.edgeLow, self.energyPlot.edgeUp))
                aaHtgNumBkg[-1][-1].SetFillStyle(0)
                aaHtgNumBkg[-1][-1].SetLineWidth(3-iS)
                aaHtgNumBkg[-1][-1].SetLineStyle(2-iS)
                aaHtgNumBkg[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumBkg[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgNumBkg[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgNumBkg[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgSBratio[-1].append(ROOT.TH1D("hSBratio%s_%s_%s" % (self.name,iS,jS), "%s Signal/Background ratio of %s;log_{10}Energy[MeV];Ratio" % (self.aaClass[iS][jS], self.name), self.energyPlot.nBin, self.energyPlot.edgeLow, self.energyPlot.edgeUp))
                aaHtgSBratio[-1][-1].SetFillStyle(0)
                aaHtgSBratio[-1][-1].SetLineWidth(3-iS)
                aaHtgSBratio[-1][-1].SetLineStyle(2-iS)
                aaHtgSBratio[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgSBratio[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgSBratio[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgSBratio[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgSgnf[-1].append(ROOT.TH1D("hSignificance%s_%s_%s" % (self.name,iS,jS), "%s Significance of %s;log_{10}Energy[MeV];[#sigma]" % (self.aaClass[iS][jS], self.name), self.energyPlot.nBin, self.energyPlot.edgeLow, self.energyPlot.edgeUp))
                aaHtgSgnf[-1][-1].SetFillStyle(0)
                aaHtgSgnf[-1][-1].SetLineWidth(3-iS)
                aaHtgSgnf[-1][-1].SetLineStyle(2-iS)
                aaHtgSgnf[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgSgnf[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgSgnf[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgSgnf[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgEnergy[-1].append(ROOT.TH1D("hEnergy%s_%s_%s" % (self.name,iS,jS), "%s Energy plot of %s;log_{10}Energy[MeV];[counts]" % (self.aaClass[iS][jS], self.name), self.energyPlot.nBin*10, self.energyPlot.edgeLow, self.energyPlot.edgeUp))
                aaHtgEnergy[-1][-1].SetFillStyle(0)
                aaHtgEnergy[-1][-1].SetLineWidth(3-iS)
                aaHtgEnergy[-1][-1].SetLineStyle(2-iS)
                aaHtgEnergy[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgEnergy[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgEnergy[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgEnergy[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgEnergyOn[-1].append(aaHtgEnergy[-1][-1].Clone("hEnergyOn%s_%s_%s" % (self.name,iS,jS)))
                aaHtgEnergyOff[-1].append(aaHtgEnergy[-1][-1].Clone("hEnergyOff%s_%s_%s" % (self.name,iS,jS)))
                # Light curve in counts
                #if self.tDuration >= 12 * pCommon.JulYrInSec: # Year bin
                #hLCCountBinWidth = pCommon.JulYrInSec
                if self.tDuration >= 12 * pCommon.JulYrInSec/12.: # Month bin
                    hLCCountBinWidth = pCommon.JulYrInSec / 12.
                elif self.tDuration>=12 * 24*60*60: # Day bin
                    hLCCountBinWidth = 24.*60.*60.
                else: # Hour bin
                    hLCCountBinWidth = 60.*60.
                hLCCountNumBin = int(math.ceil(self.tDuration / hLCCountBinWidth))
                aaHtgLCCount[-1].append(ROOT.TH1D("hLCCount%s_%s_%s" % (self.name,iS,jS), "%s Light curve (in counts) of %s;MJD - %s [day;[counts]" % (self.aaClass[iS][jS], self.name, ConvertMetToMjd(self.tStart)), hLCCountNumBin, 0, self.tDuration/86400.))
                #aaHtgLCCount[-1].append(ROOT.TH1D("hLCCount%s_%s_%s" % (self.name,iS,jS), "%s Light curve (in counts) of %s;After MET%s [sec];[counts]" % (self.aaClass[iS][jS], self.name, self.tStart), hLCCountNumBin, 0, self.tDuration))
                aaHtgLCCount[-1][-1].SetFillStyle(3002-2001*iS)
                aaHtgLCCount[-1][-1].SetLineWidth(3-iS)
                aaHtgLCCount[-1][-1].SetLineStyle(2-iS)
                aaHtgLCCount[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgLCCount[-1][-1].SetFillColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgLCCount[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgLCCount[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgLCCount[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgLCCountOn[-1].append(aaHtgLCCount[-1][-1].Clone("hLCCountOn%s_%s_%s" % (self.name,iS,jS)))
                aaHtgLCCountOff[-1].append(aaHtgLCCount[-1][-1].Clone("hLCCountOff%s_%s_%s" % (self.name,iS,jS)))
        self.aaHtgNumOn = aaHtgNumOn
        self.aaHtgNumOff = aaHtgNumOff
        self.aaHtgNumSig = aaHtgNumSig
        self.aaHtgNumBkg = aaHtgNumBkg
        self.aaHtgSBratio = aaHtgSBratio
        self.aaHtgSgnf = aaHtgSgnf
        self.aaHtgEnergy = aaHtgEnergy
        self.aaHtgLCCount = aaHtgLCCount
        self.aaHtgEnergyOn = aaHtgEnergyOn
        self.aaHtgLCCountOn = aaHtgLCCountOn
        self.aaHtgEnergyOff = aaHtgEnergyOff
        self.aaHtgLCCountOff = aaHtgLCCountOff

        # Histograms for HEALPix plots
        self.lstt_hp_htg = []
#       self.lstt_npar_hp = []
        for (icat, lstcat) in enumerate(self.aaClass):
            self.lstt_hp_htg.append([])
           #self.lstt_npar_hp.append([])
            for (icla, cla) in enumerate(lstcat):
                npar_xaxis = array('d', self.energyPlot.aBin)
                npar_yaxis = array('d', self.costhetaPlot.aBin)
                npar_zaxis = array('d', range(self.NPIX+1))
                #self.lstt_hp_htg[-1].append(ROOT.TH3I('htg3D_{0}_{1}'.format(self.name, cla), '{0} {1}'.format(cla, self.name), int(self.energyPlot.nBin), float(self.energyPlot.edgeLow), float(self.energyPlot.edgeUp), int(self.costhetaPlot.nBin), float(self.costhetaPlot.edgeLow), float(self.costhetaPlot.edgeUp), int(self.NPIX), 0, float(self.NPIX)))
                self.lstt_hp_htg[-1].append(ROOT.TH3I('htg3D_{0}_{1}'.format(self.name, cla), '{0} {1}'.format(cla, self.name), int(self.energyPlot.nBin), npar_xaxis, int(self.costhetaPlot.nBin), npar_yaxis, int(self.NPIX), npar_zaxis))
                #self.lstt_npar_hp[-1].append()


    def calc(self):
        for iS in range(len(self.aaClass)):
            for jS in range(len(self.aaClass[iS])):
                self.aaHtgNumSig[iS][jS].Reset()
                self.aaHtgNumBkg[iS][jS].Reset()
                self.aaHtgSBratio[iS][jS].Reset()
                self.aaHtgSgnf[iS][jS].Reset()
                self.aaHtgEnergy[iS][jS].Reset()
                self.aaHtgLCCount[iS][jS].Reset()
                self.aaHtgEnergy[iS][jS].Sumw2()
                self.aaHtgLCCount[iS][jS].Sumw2()
                self.aaHtgEnergyOn[iS][jS].Sumw2()
                self.aaHtgLCCountOn[iS][jS].Sumw2()
                self.aaHtgEnergyOff[iS][jS].Sumw2()
                self.aaHtgLCCountOff[iS][jS].Sumw2()

        for iS in range(len(self.aaClass)):
            print "==========", self.aaClass[iS], "=========="
            for jS in range(len(self.aaClass[iS])):
                print "----------", self.aaClass[iS][jS], "----------"
                for iE in range(self.energyPlot.nBin):
                    print self.energyPlot.getBin(iE)
                    nOn=0
                    nOff=0
                    nSig=0
                    nBkg=0
                    rSB=0
                    sgnf=0
                    nOn = self.aaHtgNumOn[iS][jS].GetBinContent(iE+1)
                    nOnErr = self.aaHtgNumOn[iS][jS].GetBinError(iE+1)
                    print "Number of ON events:", nOn, "+/-", nOnErr
                    nOff = self.aaHtgNumOff[iS][jS].GetBinContent(iE+1)
                    nOffErr = self.aaHtgNumOff[iS][jS].GetBinError(iE+1)
                    #if nOff==0:
                    #    nOff = 1
                    #    nOffErr = max(nOffErr, 1)
                    print "Number of OFF events:", nOff, "+/-", nOffErr
                    saOn = self.saOn[iS][iS][iE]
                    saOff = self.saOff[iS][iS][iE]
                    nSig = nOn - nOff*saOn/saOff
                    nSigErr = math.sqrt( nOnErr**2 + (saOn/saOff*nOffErr)**2 )
                    print "Number of Signal events:", nSig, "+/-", nSigErr
                    self.aaHtgNumSig[iS][jS].SetBinContent(iE+1, nSig)
                    self.aaHtgNumSig[iS][jS].SetBinError(iE+1, nSigErr)
                    nBkg = nOff*saOn/saOff
                    nBkgErr = nOffErr*saOn/saOff
                    print "Number of Background events:", nBkg, "+/-", nBkgErr
                    self.aaHtgNumBkg[iS][jS].SetBinContent(iE+1, nBkg)
                    self.aaHtgNumBkg[iS][jS].SetBinError(iE+1, nBkgErr)
                    if not nBkg==0:
                        rSB = nSig/nBkg
                        rSBErr = rSB * math.sqrt((nSigErr/nSig)**2 + (nBkgErr/nBkg)**2)
                    elif nSig==0:
                        rSB = 0
                        rSBErr = rSB
                    else:
                        rSB = nSig
                        rSBErr = max(rSB * math.sqrt((nSigErr/nSig)**2), nSigErr)
                    print "S/B ratio:", rSB, "+/-", rSBErr
                    self.aaHtgSBratio[iS][jS].SetBinContent(iE+1, rSB)
                    self.aaHtgSBratio[iS][jS].SetBinError(iE+1, rSBErr)
                    if nOn==0:
                        sgnf = 0
                    elif nOff==0:
                        sgnf = math.sqrt( 2 * (nOn*math.log((saOff/saOn+1)*nOn/(nOn+1)) + 1*math.log((1+saOn/saOff)*1/(nOn+1)) ) )
                    else:
                        sgnf = math.sqrt( 2 * (nOn*math.log((saOff/saOn+1)*nOn/(nOn+nOff)) + nOff*math.log((1+saOn/saOff)*nOff/(nOn+nOff)) ) )
                    print "Significance:", sgnf
                    self.aaHtgSgnf[iS][jS].SetBinContent(iE+1, sgnf)
                self.aaHtgEnergy[iS][jS].Add(self.aaHtgEnergyOn[iS][jS], 1)
                self.aaHtgEnergy[iS][jS].Add(self.aaHtgEnergyOff[iS][jS], -saOn/saOff)
                self.aaHtgLCCount[iS][jS].Add(self.aaHtgLCCountOn[iS][jS], 1)
                self.aaHtgLCCount[iS][jS].Add(self.aaHtgLCCountOff[iS][jS], -saOn/saOff)
        print ""

    def draw(self):
        legHtg = ROOT.TLegend(0.7, 0.85, 0.95, 0.55, "Event class", 'NDC')

        cNumOn = ROOT.TCanvas("cNumOn_%s" % self.name, "Number of ON events of %s" % self.name, 800, 600)
        hsNumOn = ROOT.THStack("hsNumOn_%s" % self.name, "Number of ON events of %s;log_{10}Energy[MeV];[counts]" % self.name)
        for iaHtg in range(len(self.aaHtgNumOn)):
            for iHtg in range(len(self.aaHtgNumOn[iaHtg])):
                hsNumOn.Add(self.aaHtgNumOn[iaHtg][iHtg])
                legHtg.AddEntry(self.aaHtgNumOn[iaHtg][iHtg], self.aaClass[iaHtg][iHtg], 'lp')
        cNumOn.cd()
        hsNumOn.Draw("nostack E1")
        legHtg.Draw('same')
        self.hsNumOn = hsNumOn
        self.legHtg = legHtg
        self.cNumOn = cNumOn

        cNumOff = ROOT.TCanvas("cNumOff_%s" % self.name, "Number of OFF events of %s" % self.name, 800, 600)
        hsNumOff = ROOT.THStack("hsNumOff_%s" % self.name, "Number of OFF events of %s;log_{10}Energy[MeV];[counts]" % self.name)
        for aHtg in self.aaHtgNumOff:
            for htg in aHtg:
                hsNumOff.Add(htg)
        cNumOff.cd()
        hsNumOff.Draw("nostack E1")
        legHtg.Draw('same')
        self.hsNumOff = hsNumOff
        self.cNumOff = cNumOff

        cNumSig = ROOT.TCanvas("cNumSig_%s" % self.name, "Excess counts of %s" % self.name, 800, 600)
        hsNumSig = ROOT.THStack("hsNumSig_%s" % self.name, "Excess counts of %s;log_{10}Energy[MeV];[counts]" % self.name)
        for aHtg in self.aaHtgNumSig:
            for htg in aHtg:
                hsNumSig.Add(htg)
        cNumSig.cd()
        hsNumSig.Draw("nostack E1")
        legHtg.Draw('same')
        self.hsNumSig = hsNumSig
        self.cNumSig = cNumSig

        cNumBkg = ROOT.TCanvas("cNumBkg_%s" % self.name, "Expected background counts of %s" % self.name, 800, 600)
        hsNumBkg = ROOT.THStack("hsNumBkg_%s" % self.name, "Expected background counts of %s;log_{10}Energy[MeV];[counts]" % self.name)
        for aHtg in self.aaHtgNumBkg:
            for htg in aHtg:
                hsNumBkg.Add(htg)
        cNumBkg.cd()
        hsNumBkg.Draw("nostack E1")
        legHtg.Draw('same')
        self.hsNumBkg = hsNumBkg
        self.cNumBkg = cNumBkg

        cSBratio = ROOT.TCanvas("cSBratio_%s" % self.name, "Excess/Background ratio of %s" % self.name, 800, 600)
        hsSBratio = ROOT.THStack("hsSBratio_%s" % self.name, "Excess/Background ratio of %s;log_{10}Energy[MeV];Ratio" % self.name)
        for aHtg in self.aaHtgSBratio:
            for htg in aHtg:
                hsSBratio.Add(htg)
        cSBratio.cd()
        hsSBratio.Draw("nostack E1")
        legHtg.Draw('same')
        self.hsSBratio = hsSBratio
        self.cSBratio = cSBratio

        cSgnf = ROOT.TCanvas("cSgnf_%s" % self.name, "Significance of %s" % self.name, 800, 600)
        hsSgnf = ROOT.THStack("hsSgnf_%s" % self.name, "Significance(Li-Ma) of %s;log_{10}Energy[MeV];[#sigma]" % self.name)
        for aHtg in self.aaHtgSgnf:
            for htg in aHtg:
                hsSgnf.Add(htg)
        cSgnf.cd()
        hsSgnf.Draw("nostack E1")
        legHtg.Draw('same')
        self.hsSgnf = hsSgnf
        self.cSgnf = cSgnf

        cEnergy = ROOT.TCanvas("cEnergy_%s" % self.name, "Counts in fine energy bin of %s" % self.name, 800, 600)
        hsEnergy = ROOT.THStack("hsEnergy_%s" % self.name, "Counts in fine energy bin of %s;log_{10}Energy[MeV];[counts]" % self.name)
        for aHtg in self.aaHtgEnergy:
            for htg in aHtg:
                hsEnergy.Add(htg)
        cEnergy.cd()
        hsEnergy.Draw("nostack E1")
        legHtg.Draw('same')
        self.hsEnergy = hsEnergy
        self.cEnergy = cEnergy

        cLCCount = ROOT.TCanvas("cLCCount_%s" % self.name, "Light curve (excess counts) of %s" % self.name, 800, 600)
        hsLCCount = ROOT.THStack("hsLCCount_%s" % self.name, "Light curve (excess counts) of %s;%s;[counts]" % (self.name, self.aaHtgLCCount[0][0].GetXaxis().GetTitle()))
        for aHtg in self.aaHtgLCCount:
            for htg in aHtg:
                hsLCCount.Add(htg)
        cLCCount.cd()
        hsLCCount.Draw("nostack")
        legHtgFill = ROOT.TLegend(0.7, 0.85, 0.95, 0.55, "Event class", 'NDC')
        for iaHtg in range(len(self.aaHtgLCCount)):
            for iHtg in range(len(self.aaHtgLCCount[iaHtg])):
                legHtgFill.AddEntry(self.aaHtgLCCount[iaHtg][iHtg], self.aaClass[iaHtg][iHtg], 'lpf')
        legHtgFill.Draw('same')
        self.hsLCCount = hsLCCount
        self.cLCCount = cLCCount
        self.legHtgFill = legHtgFill

    def writeObjects(self):
        self.trObj.Write("", ROOT.TObject.kOverwrite)
        for aaHtg in self.lstt_hp_htg:
            for aHtg in aaHtg:
                for htg in aHtg:
                    htg.Write("", ROOT.TObject.kOverwrite)
        for aHtg in self.aaHtgNumOn:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)                    
        for aHtg in self.aaHtgNumOff:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)                    
        for aHtg in self.aaHtgNumSig:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)                    
        for aHtg in self.aaHtgNumBkg:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)                    
        for aHtg in self.aaHtgSBratio:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)                    
        for aHtg in self.aaHtgSgnf:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)                    
        for aHtg in self.aaHtgEnergy:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)
        for aHtg in self.aaHtgLCCount:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)
        for aHtg in self.aaHtgEnergyOn:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)
        for aHtg in self.aaHtgLCCountOn:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)
        for aHtg in self.aaHtgEnergyOff:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)
        for aHtg in self.aaHtgLCCountOff:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)

        self.cNumOn.Write("", ROOT.TObject.kOverwrite)
        self.cNumOff.Write("", ROOT.TObject.kOverwrite)
        self.cNumSig.Write("", ROOT.TObject.kOverwrite)
        self.cNumBkg.Write("", ROOT.TObject.kOverwrite)
        self.cSBratio.Write("", ROOT.TObject.kOverwrite)
        self.cSgnf.Write("", ROOT.TObject.kOverwrite)
        self.cEnergy.Write("", ROOT.TObject.kOverwrite)
        self.cLCCount.Write("", ROOT.TObject.kOverwrite)

        # for FITS
        print "Making", self.name, "'s FITS files."
        for iCat in range(len(self.aClass)):
            print self.aClass[iCat]
            for iCla in range(len(self.aaClass[iCat])):
                print self.aaClass[iCat][iCla]
                tbhdr = fits.Header()
                tbhdr.set("TELESCOP", "GLAST", "name of telescope generating data")
                tbhdr.set("INSTRUME", "LAT", "name of instrument generating data")
                tbhdr.set("OBSERVER", "Peter Michelson", "GLAST/LAT PI")
                tbhdr.set("ORIGIN", "Mitsunari Takahashi", "CalOnly analysis developer")
                tbhdr.set("HDUCLAS1", "EVENTS", "extension contains events")
                tbhdr.set("TIMEUNIT", "s", "units for the time related keywords")
                tbhdr.set("TUNIT1", "MeV", "physical unit of field")
                tbhdr.set("TLIMIN1", "0.", "minimum value")
                tbhdr.set("TLIMAX1", "10000000.", "maximum value")
                tbhdr.set("TUNIT2", "deg", "physical unit of field")
                tbhdr.set("TLIMIN2", "0.", "minimum value")
                tbhdr.set("TLIMAX2", "360.", "maximum value")
                tbhdr.set("TUNIT3", "deg", "physical unit of field")
                tbhdr.set("TLIMIN3", "-90.", "minimum value")
                tbhdr.set("TLIMAX3", "90.", "maximum value")
                tbhdr.set("TUNIT4", "deg", "physical unit of field")
                tbhdr.set("TLIMIN4", "0.", "minimum value")
                tbhdr.set("TLIMAX4", "360.", "maximum value")
                tbhdr.set("TUNIT5", "deg", "physical unit of field")
                tbhdr.set("TLIMIN5", "-90.", "minimum value")
                tbhdr.set("TLIMAX5", "90.", "maximum value")
                tbhdr.set("TUNIT6", "deg", "physical unit of field")
                tbhdr.set("TLIMIN6", "0.", "minimum value")
                tbhdr.set("TLIMAX6", "180.", "maximum value")
                tbhdr.set("TUNIT7", "deg", "physical unit of field")
                tbhdr.set("TLIMIN7", "0.", "minimum value")
                tbhdr.set("TLIMAX7", "180.", "maximum value")
                tbhdr.set("TUNIT8", "s", "physical unit of field")
                tbhdr.set("TLIMIN8", "0.", "minimum value")
                tbhdr.set("TLIMAX8", "10000000000.", "maximum value")
                aCol = []
                aCol.append(fits.Column(name='ENERGY', format='E', array=self.npaaENERGY[iCat][iCla]))
                aCol.append(fits.Column(name='RA', format='E', array=self.npaaRA[iCat][iCla]))
                aCol.append(fits.Column(name='DEC', format='E', array=self.npaaDEC[iCat][iCla]))
                aCol.append(fits.Column(name='L', format='E', array=self.npaaL[iCat][iCla]))
                aCol.append(fits.Column(name='B', format='E', array=self.npaaB[iCat][iCla]))
                aCol.append(fits.Column(name='THETA', format='E', array=self.npaaTHETA[iCat][iCla]))
                aCol.append(fits.Column(name='ZENITH_ANGLE', format='E', array=self.npaaZENITH_ANGLE[iCat][iCla]))
                aCol.append(fits.Column(name='TIME', format='D', array=self.npaaTIME[iCat][iCla]))
                aCol.append(fits.Column(name='EVENT_CLASS', format='32X', array=self.npaaEVENT_CLASS[iCat][iCla]))
                aCol.append(fits.Column(name='EVENT_TYPE', format='32X', array=self.npaaEVENT_TYPE[iCat][iCla]))
                cols = fits.ColDefs(aCol)
                tbhdu = fits.BinTableHDU.from_columns(cols, tbhdr)
                tbhdu.name = "EVENTS"
                strNameFitsFile = self.name + "_" + self.aClass[iCat] + "_" + self.aaClass[iCat][iCla] + ".fits"
                tbhdu.writeto('fits/{0}'.format(strNameFitsFile), clobber=True)

    def makeHealCube(self):
        """Make NumPy arrays for HEALPix plots
"""
        print "Making NumPy arrays for HEALPix plots..."
        self.lstt_narr_hp = []
        nmeter = int(self.NPIX/100+1) # For showing progress
        for (icat, lstcat) in enumerate(self.aaClass):
            self.lstt_narr_hp.append([])
            for (icla, cla) in enumerate(lstcat):
                print cla
                path_np_map = 'HEALPixMap_{0}_{1}_NSIDE{2}.npy'.format(self.name, cla, self.NSIDE)
                if os.path.exists(path_np_map):
                    self.lstt_narr_hp[-1].append(np.load(path_np_map))
                    print 'NumPy array file has been loaded.'
                else:
                    self.lstt_narr_hp[-1].append(np.zeros((self.energyPlot.nBin, self.costhetaPlot.nBin, self.NPIX)))
                    htg = self.lstt_hp_htg[icat][icla]
                    for ienr in range(self.energyPlot.nBin):
                        print '  {0}th energy bin'.format(ienr)
                        for icth in range(self.costhetaPlot.nBin):
                            print '  {0}th cos(theta) bin'.format(icth)
                            for ipix in range(self.NPIX):
                                if ipix in self.tpl_pix_valid:
                                    self.lstt_narr_hp[-1][-1][ienr][icth][ipix] = htg.GetBinContent(ienr+1, icth+1, ipix+1)
                                    if ipix%nmeter is 0:
                                        sys.stdout.write('x')
                                        sys.stdout.flush()
                                else:
                                    self.lstt_narr_hp[-1][-1][ienr][icth][ipix] = hppf.UNSEEN
                                    if ipix%nmeter is 0:
                                        sys.stdout.write('.')
                                        sys.stdout.flush()
                            print ''
                    np.save(path_np_map, self.lstt_narr_hp[-1][-1])

        return self.lstt_narr_hp # In radians


    def setHealCube(self, lstt_cube):
        self.lstt_narr_hp = lstt_cube
        print "HEALPix cube has been set."


    def smear(self, path_king='/disk/gamma/cta/store/takhsm/FermiMVA/Dispersion/AG_dispersion.root'):
        """Smear HELAPix map in ROOT histogram.
"""
        if self.NARR_PIXS_ANGDIST is None:
            print "No pixel distance table!!!"
            return 1
        if self.lstt_narr_hp[1] is None:
            print "No CalOnly map arrays."
            return 1

        print "Smearing is ongoing..."
        FILE_KING = ROOT.TFile(path_king, 'READ')
        TP_HTG_KING = (FILE_KING.Get('htgKingN'), FILE_KING.Get('htgKingS'), FILE_KING.Get('htgKingG'))
        fc_King_annulus = ROOT.TF1("fc_King_annulus", "TMath::Sin(x)*[0]*(1.-1./[2])*pow(1.+(x/[1])**2/2./[2],-[2])/[1]**2", 0, math.pi)
        fc_King = ROOT.TF1("fc_King", "[0]*(1.-1./[2])*pow(1.+(x/[1])**2/2./[2],-[2])/2./TMath::Pi()/[1]**2", 0, math.pi)

        sa_pix = hppf.nside2pixarea(self.NSIDE) # Solid angle of a pixel [sr]
        lst_narr_CalOnly_smr = []
        for (icla,cla) in enumerate(self.aaClass[1]):
            print " ", cla
            #htg3D = self.lstt_hp_htg[1][icla]
            #lst_hp_htgCalOnly_smr.append(htg3D.Clone('{0}_smeared'.format(htg3D.GetName())))
            lst_narr_CalOnly_smr.append(np.zeros((self.energyPlot.nBin, self.costhetaPlot.nBin, self.NPIX)))
            for ienr in range(self.energyPlot.nBin):
                print '  {0}th energy bin'.format(ienr)
                kxbin = TP_HTG_KING[0].GetXaxis().FindBin(self.energyPlot.getBinCenter(ienr))
                for icth in range( htg3D.GetYaxis().GetNbins()):
                    print '    {0}th cos(theta) bin'.format(icth)
                    kybin = TP_HTG_KING[0].GetYaxis().FindBin(htg3D.GetYaxis().GetBinCenter(icth))
                    if kxbin>0 and kybin>0:
                        for ipar in range(3): # Setting the parameters of King function
                            # PSF
                            #print kxbin, kybin
                            par_value = TP_HTG_KING[ipar].GetBinContent(kxbin, kybin)
                            #print '    Parameter No.{0}:'.format(ipar), par_value
                            fc_King_annulus.FixParameter(ipar, par_value)
                            fc_King.FixParameter(ipar, par_value)
                        factor_norm = 1.0/fc_King_annulus.Integral(0, math.pi)
                       #print "Normalization factor:", factor_norm
                        for (i, ipix) in enumerate(self.tpl_pix_valid):
                            cnt = self.lstt_narr_hp[1][icla][ienr][icth][ipix]#htg3D.GetBinContent(ienr, icth, ipix+1)
                            if cnt>0:
                            #sys.stdout.write('.')
                                for (j, jpix) in enumerate(self.tpl_pix_valid):
                                    angdist = self.NARR_PIXS_ANGDIST[i+1][j]
                                    lst_narr_CalOnly_smr[icla][ienr][icth][jpix] += cnt*fc_King.Eval(angdist)*factor_norm*sa_pix

        self.lstt_narr_hp[1] = lst_narr_CalOnly_smr
        print self.lstt_narr_hp[1]
        return lst_narr_CalOnly_smr
        

    def drawHealPlot(self):
        """Draw and save HEALPix plots.
"""
        if not self.lstt_hp_htg==None:
            print "Making HEALPix plots..."
            for (icat,lstcat) in enumerate(self.aaClass):
                for (icla, cla) in enumerate(lstcat):
                    print " ", cla
                    spcube = self.lstt_hp_htg[icat][icla]
                    cube = np.sum(spcube, axis=1) # Energy vs. cosTheta vs. pixel =>  Energy vs. Pixel
                    for ienr in range(self.energyPlot.nBin+1):
                        if ienr is self.energyPlot.nBin:
                            enrlow = self.energyPlot.getBin(0)[0]
                            enrup = self.energyPlot.getBin(self.energyPlot.nBin-1)[1]
                            hpmap = np.sum(cube[ienr], axis=0) # Energy vs. Pixel => Pixel
                        else:
                            enrlow = self.energyPlot.getBin(ienr)[0]
                            enrup = self.energyPlot.getBin(ienr)[1]
                            hpmap = cube[ienr] # Energy vs. Pixel => Pixel
                        print "    {0:.1f} - {1:.1f} GeV".format(10**(enrlow-3), 10**(enrup-3))
                        if mollweide==False:
                            hp.visufunc.cartview(hpmap, rot=(self.loncntr, self.latcntr, 0), coord=self.coord, lonra=[-self.radius, self.radius], latra=[-self.radius, self.radius], min=0, flip='astro', title="{0} {1} ({2:.1f} - {3:.1f} GeV)".format(cla, self.name, 10**(enrlow-3), 10**(enrup-3)), unit='counts')
                        else:
                            hp.visufunc.mollview(hpmap, rot=(self.loncntr, self.latcntr, 0), coord=self.coord, min=0, flip='astro', title="{0} {1} ({2:.1f} - {3:.1f} GeV)".format(cla, self.name, 10**(enrlow-3), 10**(enrup-3)), unit='counts')
                        plt.savefig("png/{0}_{1}_NSIDE{2}_E{3}-{4}.png".format(self.name, cla, self.NSIDE, int(100*enrlow+0.5), int(100*enrup+0.5)))
                        plt.close()
        else:
            print "Making HEALPix plots is skipped because no map arrays."
                            

class PointSource(Target):
    def __init__(self, strName, raTgt, decTgt, glTgt, gbTgt, zRedshift=0, rAppa=12., rOffMax=[1.7, 12.], rOffMin=[1., 8.], rOnMax=[0.45, 4.0], config = ClassConfig(), eRegion=EnergyLogRegion(3,4.75,0.25), ePlotRegion=EnergyLogRegion(3, 4.75, 0.25), perf=["/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R100_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R30_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R10_perf.root"], tStart=252460800.0, tEnd=504921604.0, nhpside=64):
        Target.__init__(self, strName, config, eRegion, ePlotRegion, tStart, tEnd, nhpside)
        #self.name = strName
        self.perf = perf
        self.glCntr = glTgt
        self.gbCntr = gbTgt
        self.raCntr = raTgt
        self.decCntr = decTgt
        self.zRedshift = zRedshift
        if rAppa>=90:
            print "Too large radius."
        self.radius = rAppa
        self.HEALPLOT_MOLLWEIDE = False
        self.loncntr = self.raCntr
        self.latcntr = self.decCntr
        self.coord = 'E'
        self.TP_DIR_TRUE = (math.pi/2.-radians(self.decCntr), radians(self.raCntr))
       # self.aaClass = config.aaStrSelect
       # self.energySet = self.eRegion
        self.rOffMax = []#rOffMax
        self.rOffMin = []#rOffMin
        self.rOnMax = []#rOnMax
        self.saOn = []
        self.saOff = []
        for ict in range(len(self.aaClass)):
            self.rOffMax.append([])
            self.rOffMin.append([])
            self.rOnMax.append([])
            self.saOn.append([])
            self.saOff.append([])
            for icl in range(len(self.aaClass[ict])):
                print '***', self.aaClass[ict][icl], '***'
                self.rOffMax[ict].append([])
                self.rOffMin[ict].append([])
                self.rOnMax[ict].append([])
                self.saOn[ict].append([])
                self.saOff[ict].append([])
                for ie in range(self.energySet.nBin):
                    if ict==0:
                        self.rOffMax[ict][icl].append(self.radius)
                        self.rOffMin[ict][icl].append(self.perf.getPSF95(ict, icl, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.rOnMax[ict][icl].append(min(self.ulOn, self.perf.getPSF68(ict, icl, self.energySet.aBin[ie]+self.energySet.wBin/2.0)*self.factorOn))
                    elif ict==1:
                        self.rOffMax[ict][icl].append(self.radius)
                        self.rOffMin[ict][icl].append(self.perf.getPSF95(ict, icl, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.rOnMax[ict][icl].append(min(self.ulOn, self.perf.getPSF68(ict, icl, self.energySet.aBin[ie]+self.energySet.wBin/2.0)*self.factorOn))
                        if not (self.rOffMax[ict][icl][-1]<=self.radius and self.rOffMin[ict][icl][-1]<self.rOffMax[ict][icl][-1] and self.rOnMax[ict][icl][-1]<=self.rOffMin[ict][icl][-1] and self.rOnMax[ict][icl][-1]>0):
                            print "Bad region setup!!"
                            sys.exit(1)
                    print 'ON region:0.0 -', self.rOnMax[ict][icl][-1], 'deg'
                    print 'OFF region:', self.rOffMin[ict][icl][-1], '-', self.rOffMax[ict][icl][-1], 'deg'
                    print 'cf. PSF68:', self.perf.getPSF68(ict, icl, self.energySet.aBin[ie]+self.energySet.wBin/2.0), 'deg', 'PSF95:', self.perf.getPSF95(ict, icl, self.energySet.aBin[ie]+self.energySet.wBin/2.0), 'deg'
                    self.saOn[ict][icl].append( 2.0 * math.pi * ( cos(radians(0.0)) - cos(radians(self.rOnMax[ict][icl][-1])) ) )
                    self.saOff[ict][icl].append( 2.0 * math.pi * ( cos(radians(self.rOffMin[ict][icl][-1])) - cos(radians(self.rOffMax[ict][icl][-1])) ) )
            
#        for ic in range(max(len(rOffMax), len(rOffMin), len(rOnMax))):
#            if not (self.rOffMax[ic]<=self.radius and self.rOffMin[ic]<rOffMax[ic] and self.rOnMax[ic]<=self.rOffMin[ic] and self.rOnMax[ic]>0):
 #               print "Bad region setup!!"
#        self.saOn = []
#        self.saOff = []
#        for hS in range(len(self.aaClass)):
#            self.saOn.append( 2.0 * math.pi * ( cos(radians(0.0)) - cos(radians(self.rOnMax[hS])) ) )
 #           self.saOff.append( 2.0 * math.pi * ( cos(radians(self.rOffMin[hS])) - cos(radians(self.rOffMax[hS])) ) )
        print "Solid angle of ON region:", self.saOn
        print "Solid angle of OFF region:", self.saOff

        aaaGrpMap = []
        aaGreHighEnergy = []
        aaaHtgMap = []
        aaaHtgTheta = []
        aaaHtgMapFolded = []
        aaHtgMapAll = []
        aaHtgMapAllFolded = []
        hwMap = 15
        aaaHtgPsf = []
        for iS in range(len(self.aaClass)):
            aaaGrpMap.append([])
            aaGreHighEnergy.append([])
            aaaHtgMap.append([])
            aaaHtgMapFolded.append([])
            aaHtgMapAll.append([])
            aaHtgMapAllFolded.append([])
            aaaHtgPsf.append([])
            aaaHtgTheta.append([])
#            aaaHtgFolded.append([])
            for jS in range(len(self.aaClass[iS])):
                aaaGrpMap[-1].append([])
                aaaHtgMap[-1].append([])
                aaaHtgMapFolded[-1].append([])
                #aaHtgMapAll[-1].append([])
                #aaHtgMapAllFolded[-1].append([])
                aaaHtgPsf[-1].append([])
                aaaHtgTheta[-1].append([])
#                aaaHtgFolded[-1].append([])
                aaGreHighEnergy[-1].append(ROOT.TGraphAsymmErrors())#TGraphErrors())
                aaGreHighEnergy[-1][-1].SetName("greHighEnergy%s_%s_%s" % (self.name, iS, jS))
                if (2-iS)*jS!=2:
                    aaGreHighEnergy[-1][-1].SetMarkerColor(kGray+(2-iS)*jS*2)
                else:
                    aaGreHighEnergy[-1][-1].SetMarkerColor(kBlack)
                aaGreHighEnergy[-1][-1].SetMarkerStyle(20+iS)
                aaGreHighEnergy[-1][-1].SetMarkerSize(1.0-0.2*iS)
                for iE in range(self.energyPlot.nBin):
                    aaaGrpMap[-1][-1].append(ROOT.TGraphPolar())
                    aaaGrpMap[-1][-1][-1].SetName("grpMap%s_%s_%s_%s" % (self.name, iS, jS, iE))
                    #aaaGrpMap[-1][-1][-1].SetTitle("{0} map around {1} in {2:.1f} - {3:.1f} GeV".format(self.aaClass[iS][jS], self.name, 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)))
                    aaaGrpMap[-1][-1][-1].SetMaxRadial(self.radius)
                    aaaGrpMap[-1][-1][-1].SetMarkerColor(pColor.akColor(iE))
                    aaaGrpMap[-1][-1][-1].SetMarkerStyle(7)
                    if abs(self.decCntr)<60:
                        xLowMap = math.floor(-hwMap/cos(abs(radians(self.decCntr))))
                        xUpMap = math.ceil(hwMap/cos(abs(radians(self.decCntr))))
                        #xLowMap = math.floor(self.raCntr-hwMap/cos(abs(radians(self.decCntr))))
                        #xUpMap = math.ceil(self.raCntr+hwMap/cos(abs(radians(self.decCntr))))
                        nxMap = int(xUpMap - xLowMap)*2
                        yLowMap = self.decCntr - hwMap
                        yUpMap = self.decCntr + hwMap
                        nyMap = int(yUpMap - yLowMap)*2
                        #print "Map histogram:", nxMap, xLowMap, xUpMap, nyMap, yLowMap, yUpMap
                        aaaHtgMap[-1][-1].append(ROOT.TH2D("hMapCel%s_%s_%s_%s" % (self.name, iS, jS, iE), "{0} count map around {1} in {2:.1f} - {3:.1f} GeV;-(RA-{4})[#circ];DEC[#circ];[counts]".format(self.aaClass[iS][jS], self.name, 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3), self.raCntr), nxMap, xLowMap, xUpMap, nyMap, yLowMap, yUpMap))
                        aaaHtgMapFolded[-1][-1].append(aaaHtgMap[-1][-1][-1].Clone("hMapCelFolded%s_%s_%s_%s" % (self.name, iS, jS, iE)))
                    else:
                        xLowMap = math.floor(-hwMap/cos(abs(radians(self.gbCntr))))
                        xUpMap = math.ceil(hwMap/cos(abs(radians(self.gbCntr))))
                        #xLowMap = math.floor(self.glCntr-hwMap/cos(abs(radians(self.gbCntr))))
                        #xUpMap = math.ceil(self.glCntr+hwMap/cos(abs(radians(self.gbCntr))))
                        nxMap = int(xUpMap - xLowMap)*2
                        yLowMap = self.gbCntr - hwMap
                        yUpMap = self.gbCntr + hwMap
                        nyMap = int(yUpMap - yLowMap)*2
                        aaaHtgMap[-1][-1].append(ROOT.TH2D("hMapGal%s_%s_%s_%s" % (self.name, iS, jS, iE), "{0} count map around {1} in {2:.1f} - {3:.1f} GeV;-(Galactic longitude - {4}) [#circ];Galactic latitude [#circ];[counts]".format(self.aaClass[iS][jS], self.name, 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3), self.glCntr), nxMap, xLowMap, xUpMap, nyMap, yLowMap, yUpMap))
                        aaaHtgMapFolded[-1][-1].append(aaaHtgMap[-1][-1][-1].Clone("hMapGalFolded%s_%s_%s_%s" % (self.name, iS, jS, iE)))
                    aaaHtgPsf[-1][-1].append(ROOT.TH2D(aaaHtgMap[-1][-1][-1].Clone("hPsf%s_%s_%s_%s" % (self.name, iS, jS, iE))))
                    #aaaHtgMap[-1][-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                    #aaaHtgMap[-1][-1][-1].SetFillColor(pColor.akColor(jS+(1-iS)*jS))
                    #aaaHtgMap[-1][-1][-1].SetLineWidth(3-iS)
                    #aaaHtgMap[-1][-1][-1].SetLineStyle(2-iS)
                    aaaHtgTheta[-1][-1].append(ROOT.TH1D("hTheta%s_%s_%s_%s" % (self.name, iS, jS, iE), "{0} #theta^2 plot of {1} in {2:.1f} - {3:.1f} GeV;#theta^2 [#circ];[counts]".format(self.aaClass[iS][jS], self.name, 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), int(rAppa**2)*10, 0, int(rAppa**2)))
                    aaaHtgTheta[-1][-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                    aaaHtgTheta[-1][-1][-1].SetLineWidth(3-iS)
                    aaaHtgTheta[-1][-1][-1].SetLineStyle(2-iS)
                # else:
                #         xLowMap = math.floor(-hwMap/cos(abs(radians(self.gbCntr))))
                #         xUpMap = math.ceil(hwMap/cos(abs(radians(self.gbCntr))))
                #         #xLowMap = math.floor(self.glCntr-hwMap/cos(abs(radians(self.gbCntr))))
                #         #xUpMap = math.ceil(self.glCntr+hwMap/cos(abs(radians(self.gbCntr))))
                #         nxMap = int(xUpMap - xLowMap)*2
                #         yLowMap = self.gbCntr - hwMap
                #         yUpMap = self.gbCntr + hwMap
                #         nyMap = int(yUpMap - yLowMap)*2
                #         aaaHtgMap[-1][-1].append(ROOT.TH2D("hMapGal%s_%s_%s_%s" % (self.name, iS, jS, iE), "{0} count map around {1} in {2:.1f} - {3:.1f} GeV;-(Galactic longitude - {4}) [#circ];Galactic latitude [#circ];[counts]".format(self.aaClass[iS][jS], self.name, 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3), self.glCntr), nxMap, xLowMap, xUpMap, nyMap, yLowMap, yUpMap))
                #         aaaHtgMapFolded[-1][-1].append(aaaHtgMap[-1][-1][-1].Clone("hMapGalFolded%s_%s_%s_%s" % (self.name, iS, jS, iE)))
                #     aaaHtgPsf[-1][-1].append(ROOT.TH2D(aaaHtgMap[-1][-1][-1].Clone("hPsf%s_%s_%s_%s" % (self.name, iS, jS, iE))))
                #     #aaaHtgMap[-1][-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                #     #aaaHtgMap[-1][-1][-1].SetFillColor(pColor.akColor(jS+(1-iS)*jS))
                #     #aaaHtgMap[-1][-1][-1].SetLineWidth(3-iS)
                #     #aaaHtgMap[-1][-1][-1].SetLineStyle(2-iS)
                #     aaaHtgTheta[-1][-1].append(ROOT.TH1D("hTheta%s_%s_%s_%s" % (self.name, iS, jS, iE), "{0} #theta^2 plot of {1} in {2:.1f} - {3:.1f} GeV;#theta^2 [#circ];[counts]".format(self.aaClass[iS][jS], self.name, 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), int(rAppa**2)*10, 0, int(rAppa**2)))
                #     aaaHtgTheta[-1][-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                #     aaaHtgTheta[-1][-1][-1].SetLineWidth(3-iS)
                #     aaaHtgTheta[-1][-1][-1].SetLineStyle(2-iS)
                if abs(self.decCntr)<60:
                    aaHtgMapAll[-1].append(ROOT.TH2D("hMapAllCel%s_%s_%s" % (self.name, iS, jS), "{0} count map around {1} in {2:.1f} - {3:.1f} GeV;-(RA-{4})[#circ];DEC[#circ];[counts]".format(self.aaClass[iS][jS], self.name, 10**(self.energyPlot.getRegionLowEdge()-3), 10**(self.energyPlot.getRegionUpEdge()-3), self.raCntr), nxMap, xLowMap, xUpMap, nyMap, yLowMap, yUpMap))
                    aaHtgMapAllFolded[-1].append(aaHtgMapAll[-1][-1].Clone("hMapAllFoldedCel%s_%s_%s" % (self.name, iS, jS)))
                else:
                    aaHtgMapAll[-1].append(ROOT.TH2D("hMapAllGall%s_%s_%s" % (self.name, iS, jS), "{0} count map around {1} in {2:.1f} - {3:.1f} GeV;-(Galactic longitude -{4})[#circ];Galactic latitude[#circ];[counts]".format(self.aaClass[iS][jS], self.name, 10**(self.energyPlot.getRegionLowEdge()-3), 10**(self.energyPlot.getRegionUpEdge()-3), self.glCntr), nxMap, xLowMap, xUpMap, nyMap, yLowMap, yUpMap))
                    aaHtgMapAllFolded[-1].append(aaHtgMapAll[-1][-1].Clone("hMapAllFoldedGal%s_%s_%s" % (self.name, iS, jS)))
#                    aaaHtgFolded[-1][-1].append(ROOT.TH2D("hFolded%s_%s_%s_%s" % (self.name, iS, jS, iE), "{0} PSF folded map around {1} in {2:.1f} - {3:.1f} GeV;RA[#circ];DEC[#circ];[counts]".format(self.aaClass[iS][jS], self.name, 10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)), 12,))
        self.aaaGrpMap = aaaGrpMap
        self.aaGreHighEnergy = aaGreHighEnergy
        self.aaaHtgMap = aaaHtgMap
        self.aaaHtgMapFolded = aaaHtgMapFolded
        self.aaHtgMapAll = aaHtgMapAll
        self.aaHtgMapAllFolded = aaHtgMapAllFolded
        self.aaaHtgPsf = aaaHtgPsf
        self.aaaHtgTheta = aaaHtgTheta        
                        
    def fill(self, raEvent, decEvent, glEvent, gbEvent, eEvent, ctEvent, clEvent, zEvent, tEvent, cthEvent):
        #Target.fill(self, lEvent, bEvent, eEvent, ctEvent, clEvent)
        self.flagOn[0] = kFALSE
        self.flagOff[0] = kFALSE
        binE = self.energySet.findBin(eEvent)
        plotE = self.energyPlot.findBin(eEvent)
        raRad = math.radians(raEvent)
        decRad = math.radians(decEvent)
        glRad = math.radians(glEvent)
        gbRad = math.radians(gbEvent)
        vecTgt = np.array([cos(radians(self.decCntr))*cos(radians(self.raCntr)), cos(radians(self.decCntr))*sin(radians(self.raCntr)), sin(radians(self.decCntr))])
        vecEve = np.array([cos(decRad)*cos(raRad), cos(decRad)*sin(raRad), sin(decRad)])
        radTheta = acos(np.dot(vecTgt, vecEve))
        if tEvent>=self.tStart and tEvent<= self.tEnd and binE>=0 and binE<self.energySet.nBin and zEvent<self.zCut[ctEvent-1]:
            if degrees(radTheta)<=self.radius:
                vecNorth = np.array([cos(radians(self.decCntr+90))*cos(radians(self.raCntr)), cos(radians(self.decCntr+90))*sin(radians(self.raCntr)), sin(radians(self.decCntr+90))])
                if cos(radTheta)!=0:
                    vecProj = vecEve/cos(radTheta) - vecTgt #[vecEve[0]-vecTgt[0]*cos(radTheta[0]), vecEve[1]-vecTgt[1]*cos(radTheta[1]), vecEve[2]-vecTgt[2]*cos(radTheta[2])]
                else:
                    print "Zero division!!"
                vecCross = np.cross(vecProj, vecNorth)
                radPhi = acos(np.dot(vecProj, vecNorth) / (np.linalg.norm(vecProj)*np.linalg.norm(vecNorth)))
                if np.dot(vecCross, vecNorth)<0:
                    radPhi = -radPhi +  2.*math.pi

                for clEventPlus in range(clEvent-int(ctEvent==1 and clEvent==3)):
                    self.lstt_hp_htg[ctEvent-1][clEventPlus].Fill(eEvent, cthEvent, hppf.ang2pix(self.NSIDE, math.pi/2.-math.radians(decEvent), math.radians(raEvent))+0.5)
                    self.aaaHtgTheta[ctEvent-1][clEventPlus][plotE].Fill(math.degrees(radTheta)**2)
                    self.aaaGrpMap[ctEvent-1][clEventPlus][plotE].SetPoint(self.aaaGrpMap[ctEvent-1][clEventPlus][plotE].GetN(), radPhi, math.degrees(radTheta))
                    if abs(self.decCntr)<60:
                        if raEvent<self.raCntr+self.aaaHtgMap[ctEvent-1][clEventPlus][plotE].GetXaxis().GetBinLowEdge(1):
                            raEventCor = (360+raEvent)
                        elif raEvent>self.raCntr+self.aaaHtgMap[ctEvent-1][clEventPlus][plotE].GetXaxis().GetBinUpEdge(self.aaaHtgMap[ctEvent-1][clEventPlus][plotE].GetNbinsX()):
                            raEventCor = (360-raEvent)
                        else:
                            raEventCor = raEvent
                        self.aaaHtgMap[ctEvent-1][clEventPlus][plotE].Fill(-(raEventCor-self.raCntr), decEvent)
                    else:
                        if glEvent<self.aaaHtgMap[ctEvent-1][clEventPlus][plotE].GetXaxis().GetBinLowEdge(1):
                            glEventCor = (360+glEvent)
                        elif glEvent>self.aaaHtgMap[ctEvent-1][clEventPlus][plotE].GetXaxis().GetBinUpEdge(self.aaaHtgMap[ctEvent-1][clEventPlus][plotE].GetNbinsX()):
                            glEventCor = (360-glEvent)
                        else:
                            glEventCor = glEvent
                        self.aaaHtgMap[ctEvent-1][clEventPlus][plotE].Fill(-(glEventCor-self.glCntr), gbEvent)
                if math.degrees(radTheta) < self.rOffMax[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]:
                    # for FITS
                    self.npaaENERGY[ctEvent-1][clEventPlus] = np.append(self.npaaENERGY[ctEvent-1][clEventPlus], eEvent)
                    self.npaaRA[ctEvent-1][clEventPlus] = np.append(self.npaaRA[ctEvent-1][clEventPlus], raEvent)
                    self.npaaDEC[ctEvent-1][clEventPlus] = np.append(self.npaaDEC[ctEvent-1][clEventPlus], decEvent)
                    self.npaaL[ctEvent-1][clEventPlus] = np.append(self.npaaL[ctEvent-1][clEventPlus], glEvent)
                    self.npaaB[ctEvent-1][clEventPlus] = np.append(self.npaaB[ctEvent-1][clEventPlus], gbEvent)
                    self.npaaTHETA[ctEvent-1][clEventPlus] = np.append(self.npaaTHETA[ctEvent-1][clEventPlus], cthEvent)
                    self.npaaZENITH_ANGLE[ctEvent-1][clEventPlus] = np.append(self.npaaZENITH_ANGLE[ctEvent-1][clEventPlus], zEvent)
                    self.npaaTIME[ctEvent-1][clEventPlus] = np.append(self.npaaTIME[ctEvent-1][clEventPlus], tEvent)
                    if ctEvent==1:
                        self.npaaEVENT_CLASS[ctEvent-1][clEventPlus] = np.append(self.npaaEVENT_CLASS[ctEvent-1][clEventPlus], 128*2**(clEvent-1)/2)
                    elif ctEvent==2:
                        self.npaaEVENT_CLASS[ctEvent-1][clEventPlus] = np.append(self.npaaEVENT_CLASS[ctEvent-1][clEventPlus], 4096*2**(clEvent-1))
                    else:
                        print "Error!!! Event type is neither 1 or 2"
                    if cthEvent<0.7:
                        self.npaaEVENT_TYPE[ctEvent-1][clEventPlus] = np.append(self.npaaEVENT_TYPE[ctEvent-1][clEventPlus], 1)
                    else:
                        self.npaaEVENT_TYPE[ctEvent-1][clEventPlus] = np.append(self.npaaEVENT_TYPE[ctEvent-1][clEventPlus], 2)

                    if math.degrees(radTheta) < self.rOnMax[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][plotE]:
                        self.flagOn[0] = kTRUE
                        for clEventPlus in range(clEvent-int(ctEvent==1 and clEvent==3)):
                            self.aaHtgNumOn[ctEvent-1][clEventPlus].Fill(eEvent)
                        self.aaHtgEnergyOn[ctEvent-1][clEventPlus].Fill(eEvent)# (10**(eEvent-3))**2)
                        self.aaHtgLCCountOn[ctEvent-1][clEventPlus].Fill(ConvertMetToMjd(tEvent)-ConvertMetToMjd(self.tStart))
                    #self.aaHtgLCCountOn[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].Fill(tEvent-self.tStart)
                        self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].SetPoint(self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].GetN(), ConvertMetToMjd(tEvent)-ConvertMetToMjd(self.tStart), eEvent)
                        if ctEvent==1:
                            self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].SetPointEYhigh(self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].GetN()-1, abs(math.log10(1+self.perf.getEdisp68(0, clEvent-int(clEvent==3)-1, eEvent))))
                            self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].SetPointEYlow(self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].GetN()-1, abs(math.log10(1-self.perf.getEdisp68(0, clEvent-int(clEvent==3)-1, eEvent))))
                        elif ctEvent==2:
                            self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].SetPointEYhigh(self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].GetN()-1, abs(math.log10(1+self.perf.getEdisp68_cth(1, clEvent-1, eEvent, cthEvent))))
                            self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].SetPointEYlow(self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].GetN()-1, abs(math.log10(1-self.perf.getEdisp68_cth(1, clEvent-1, eEvent, cthEvent))))
                    elif math.degrees(radTheta) >= self.rOffMin[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]: 
                        self.flagOff[0] = kTRUE
                        for clEventPlus in range(clEvent-int(ctEvent==1 and clEvent==3)):
                            self.aaHtgNumOff[ctEvent-1][clEventPlus].Fill(eEvent)
                            self.aaHtgEnergyOff[ctEvent-1][clEventPlus].Fill(eEvent)
                            self.aaHtgLCCountOff[ctEvent-1][clEventPlus].Fill(ConvertMetToMjd(tEvent)-ConvertMetToMjd(self.tStart))
        self.trObj.Fill()

    def calc(self):
        Target.calc(self)
        fPSF = ROOT.TF1("fPSF", "1. / (2*TMath::Pi()*[0]**2) * TMath::Exp(-(x/[0])**2/2.0)")
        for tE in range(self.energyPlot.nBin):
            print "Energy bin:", self.energyPlot.aBin[tE]
            for kAS in range(len(self.aClass)):
                for kSS in range(len(self.aaClass[kAS])):
                    print self.aaClass[kAS][kSS]
                    abin = self.energyPlot.getBin(tE)
                    slow = self.energySet.findBin(abin[0])
                    sup = self.energySet.findBin(abin[1])
                    degPSF = max(self.perf.integralPSF68(kAS, kSS, slow, sup), self.perf.integralPSF95(kAS, kSS, slow, sup)/2.0)
                    #degPSF = max(self.perf.integralPSF68(kAS, kSS, self.energySet.findBin(self.energyPlot.edgeLow(tE)), self.energySet.findBin(self.energyPlot.edgeUp(tE))), self.perf.integralPSF95(kAS, kSS, self.energySet.findBin(self.energyPlot.edgeLow(tE)), self.energySet.findBin(self.energyPlot.edgeUp(tE)))/2.0)
                    #degPSF = max(self.perf.getPSF68(kAS, kSS, self.energySet.getBinCenter(tE)), self.perf.getPSF95(kAS, kSS, self.energySet.getBinCenter(tE))/2.0)
                    fPSF.FixParameter(0, degPSF)
                    print "PSF:", degPSF
                    xCntp = (self.aaaHtgPsf[kAS][kSS][tE].GetXaxis().GetBinLowEdge(1)+self.aaaHtgPsf[kAS][kSS][tE].GetXaxis().GetBinUpEdge(self.aaaHtgPsf[kAS][kSS][tE].GetNbinsX()))/2.0
                    yCntp = (self.aaaHtgPsf[kAS][kSS][tE].GetYaxis().GetBinLowEdge(1)+self.aaaHtgPsf[kAS][kSS][tE].GetYaxis().GetBinUpEdge(self.aaaHtgPsf[kAS][kSS][tE].GetNbinsY()))/2.0
                    for pY in range(self.aaaHtgPsf[kAS][kSS][tE].GetNbinsY()):
                        yDirp = self.aaaHtgPsf[kAS][kSS][tE].GetYaxis().GetBinCenter(pY+1)
                        for pX in range(self.aaaHtgPsf[kAS][kSS][tE].GetNbinsX()):
                            xDirp = self.aaaHtgMapFolded[kAS][kSS][tE].GetXaxis().GetBinCenter(pX+1)
                            degDistp = pCommon.anglePointsDegToDeg(xCntp, yCntp, xDirp, yDirp)
                            densp = fPSF.Eval(degDistp) * abs(cos(radians(90-abs(self.aaaHtgPsf[kAS][kSS][tE].GetYaxis().GetBinUpEdge(pY+1))))-cos(radians(90-abs(self.aaaHtgPsf[kAS][kSS][tE].GetYaxis().GetBinLowEdge(pY+1))))*(self.aaaHtgPsf[kAS][kSS][tE].GetXaxis().GetBinUpEdge(pX+1)-self.aaaHtgMapFolded[kAS][kSS][tE].GetXaxis().GetBinLowEdge(pX+1)))*180./math.pi
                            self.aaaHtgPsf[kAS][kSS][tE].Fill(xDirp, yDirp, densp)
                    itglPsf = self.aaaHtgPsf[kAS][kSS][tE].Integral()
                    # for lB in range(self.aaaHtgMap[kAS][kSS][tE].GetNbinsY()):
                    #    bPoint = self.aaaHtgMap[kAS][kSS][tE].GetYaxis().GetBinCenter(lB+1)
                    #    for lL in range(self.aaaHtgMap[kAS][kSS][tE].GetNbinsX()):
                    #        lPoint = self.aaaHtgMap[kAS][kSS][tE].GetXaxis().GetBinCenter(lL+1)
                    #        nCount = self.aaaHtgMap[kAS][kSS][tE].GetBinContent(lL+1, lB+1)
                    #        if nCount>0:
                    #            for mB in range(self.aaaHtgMapFolded[kAS][kSS][tE].GetNbinsY()):
                    #                bDire = self.aaaHtgMapFolded[kAS][kSS][tE].GetYaxis().GetBinCenter(mB+1)
                    #                for mL in range(self.aaaHtgMapFolded[kAS][kSS][tE].GetNbinsX()):
                    #                    lDire = self.aaaHtgMapFolded[kAS][kSS][tE].GetXaxis().GetBinCenter(mL+1)
                    #                    degDist = pCommon.anglePointsDegToDeg(lPoint, bPoint, lDire, bDire)
                    #                    dense = fPSF.Eval(degDist) * abs(cos(radians(90-abs(self.aaaHtgMapFolded[kAS][kSS][tE].GetYaxis().GetBinUpEdge(mB+1))))-cos(radians(90-abs(self.aaaHtgMapFolded[kAS][kSS][tE].GetYaxis().GetBinLowEdge(mB+1))))*(self.aaaHtgMapFolded[kAS][kSS][tE].GetXaxis().GetBinUpEdge(mL+1)-self.aaaHtgMapFolded[kAS][kSS][tE].GetXaxis().GetBinLowEdge(mL+1)))*180./math.pi / itglPsf
                    #                    self.aaaHtgMapFolded[kAS][kSS][tE].Fill(lDire, bDire, nCount * dense)
        for kAS in range(len(self.aClass)):
            for kSS in range(len(self.aaClass[kAS])):
                for tE in range(self.energyPlot.nBin):
                    self.aaHtgMapAll[kAS][kSS].Add(self.aaaHtgMap[kAS][kSS][tE])
                    self.aaHtgMapAllFolded[kAS][kSS].Add(self.aaaHtgMapFolded[kAS][kSS][tE])


    def setPixelDistanceTable(self, deg_margin):
        """Set a table of distance between HEALPix pixels of nside and a tuple of valid pixels.
"""
        PATH_NPFILE = 'PixelsAngularDistance_{0}_NSIDE{1}.npy'.format(self.name, self.NSIDE)
        if os.path.exists(PATH_NPFILE):
            self.NARR_PIXS_ANGDIST = np.load(PATH_NPFILE)
            print 'NumPy array file has been loaded.'

        else:
            print 'Creating a table of pixel distance...'
            rad_range = radians(self.radius + deg_margin)
            lst_pix_valid = []
            dct_angdist_src = {}

          # Make a list of pixels within ROI (with your margin)
            for ipix in range(self.NPIX):
                tp_pix = hppf.pix2ang(self.NSIDE, ipix)
                angdist = hp.rotator.angdist(tp_pix, self.TP_DIR_TRUE) # Angular distance between the source and each pixel in radians
                if angdist<rad_range:
                    lst_pix_valid.append(ipix)
                    dct_angdist_src[ipix] = angdist

            self.tpl_pix_valid = tuple(lst_pix_valid)
            self.npix_valid = len(self.tpl_pix_valid)
            print self.npix_valid, 'pixels are considered.'
            narr_dist = np.zeros((self.npix_valid+1, self.npix_valid))
           # Make a table of distance between the pixels in the list
            for (h, hpix) in enumerate(self.tpl_pix_valid):
                narr_dist[0][h] = hpix
            for (j, jpix) in enumerate(self.tpl_pix_valid):
                for (k, kpix) in enumerate(self.tpl_pix_valid):
                    narr_dist[j+1][k] = hp.rotator.angdist(hppf.pix2ang(self.NSIDE, jpix), hppf.pix2ang(self.NSIDE, kpix))
            self.NARR_PIXS_ANGDIST = narr_dist
            np.save(PATH_NPFILE, self.NARR_PIXS_ANGDIST)
                        

    def draw(self):
        # if not self.lstt_hp_htg==None:
        #     for (icat,lstcat) in enumerate(self.aaClass):
        #         for (icla, cla) in enumerate(lstcat):
        #             htg3D = self.lstt_hp_htg[icat][icla]
        #             print htg3D.GetName()
        #             hpcube = Healcube(self.name, cla, htg3D, lon_cntr=self.raCntr, lat_cntr=self.decCntr, char_coord='E', deg_radius=self.radius)
        #             hpcube.setmap()
        #             hpcube.draw()
                    
                    # for binE in range(self.energyPlot.nBin):
                    #     hp.visufunc.cartview(lstt_hp_htg[binE][ctEvent][clEvent], rot=(self.raCntr, self.decCntr, 0), coord='E', lonra=[-self.radius, self.radius], latra=[-self.radius, self.radius], min=0, flip='astro', title="{0} {1} ({2:.1f} - {3:.1f} GeV)".format(self.aaClass[ctEvent][clEvent], self.name, 10**(self.energyPlot.getBin(binE)[0]-3), 10**(self.energyPlot.getBin(binE)[1]-3)), unit='counts')
                    #     plt.savefig("{0}_{1}_E{2}-{3}.png".format(self.name, self.aaClass[ctEvent][clEvent], int(100*self.energyPlot.getBin(binE)[0]+0.5), int(100*self.energyPlot.getBin(binE)[1]+0.5)))
                    #     plt.clf()
        self.fRadiusOnMax = []
        self.fRadiusOffMin = []
        self.fRadiusOffMax = []
        for ctEvent in range(len(self.aaClass)):
            self.fRadiusOnMax.append([])
            self.fRadiusOffMin.append([])
            self.fRadiusOffMax.append([])
            for clEvent in range(len(self.aaClass[ctEvent])):
                self.fRadiusOnMax[ctEvent].append([])
                self.fRadiusOffMin[ctEvent].append([])
                self.fRadiusOffMax[ctEvent].append([])
                for binE in range(self.energyPlot.nBin):
                    if abs(self.decCntr)<60:
                        cntrMap = self.decCntr
                    else:
                        cntrMap = self.gbCntr
                    self.fRadiusOnMax[ctEvent][clEvent].append(ROOT.TF2("fRadiusOnMax%s_%s_%s" % (ctEvent, clEvent, binE), "x**2+(y-%s)**2 - %s**2" % (cntrMap, self.rOnMax[ctEvent][clEvent][binE]), -self.radius, self.radius, -self.radius, self.radius))
                    self.fRadiusOffMin[ctEvent][clEvent].append(ROOT.TF2("fRadiusOffMin%s_%s_%s" % (ctEvent, clEvent, binE), "x**2+(y-%s)**2 - %s**2" % (cntrMap, self.rOffMin[ctEvent][clEvent][binE]), -self.radius, self.radius, -self.radius, self.radius))
                    self.fRadiusOffMax[ctEvent][clEvent].append(ROOT.TF2("fRadiusOffMax%s_%s_%s" % (ctEvent, clEvent, binE), "x**2+(y-%s)**2 - %s**2" % (cntrMap, self.rOffMax[ctEvent][clEvent][binE]), -self.radius, self.radius, -self.radius, self.radius))
                    self.fRadiusOnMax[-1][-1][-1].SetMinimum(0)
                    self.fRadiusOnMax[-1][-1][-1].SetMaximum(0)
                    self.fRadiusOnMax[-1][-1][-1].SetLineWidth(1)
                    self.fRadiusOnMax[-1][-1][-1].SetLineColor(kGray)
                    self.fRadiusOffMin[-1][-1][-1].SetMinimum(0)
                    self.fRadiusOffMin[-1][-1][-1].SetMaximum(0)
                    self.fRadiusOffMin[-1][-1][-1].SetLineWidth(1)
                    self.fRadiusOffMin[-1][-1][-1].SetLineColor(kGray)
                    self.fRadiusOffMax[-1][-1][-1].SetMinimum(0)
                    self.fRadiusOffMax[-1][-1][-1].SetMaximum(0)
                    self.fRadiusOffMax[-1][-1][-1].SetLineWidth(1)
                    self.fRadiusOffMax[-1][-1][-1].SetLineColor(kGray)

        # self.aaGrpRadiusDummy = []
        # self.aGrpRadiusOnMax = []
        # self.aGrpRadiusOffMin = []
        # self.aGrpRadiusOffMax = []
        # nPoint = 100
        # for ctEvent in range(len(self.aaClass)):
        #     self.aaGrpRadiusDummy.append([])
        #     for clEvent in range(len(self.aaClass[ctEvent])):
        #         self.aaGrpRadiusDummy[-1].append(ROOT.TGraphPolar(0))
        #         self.aaGrpRadiusDummy[-1][-1].SetName("grpRadiusDummy%s%s" % (ctEvent, clEvent))
        #         self.aaGrpRadiusDummy[-1][-1].SetTitle("%s map around %s" % (self.aaClass[ctEvent][clEvent], self.name))
        #     self.aGrpRadiusOnMax.append(ROOT.TGraphPolar(nPoint))
        #     self.aGrpRadiusOnMax[-1].SetName("grpRadiusOnMax%s" % self.aClass[ctEvent])
        #     self.aGrpRadiusOnMax[-1].SetTitle("Maximum radius of the %s ON region" % self.aClass[ctEvent])
        #     self.aGrpRadiusOnMax[-1].SetLineWidth(1)
        #     self.aGrpRadiusOnMax[-1].SetLineColor(kWhite)
        #     self.aGrpRadiusOnMax[-1].SetMarkerStyle(1)
        #     self.aGrpRadiusOnMax[-1].SetMarkerColor(kWhite)
        #     self.aGrpRadiusOnMax[-1].SetFillStyle(1001)
        #     self.aGrpRadiusOnMax[-1].SetFillColor(kWhite)
        #     for iPoint in range(nPoint):
        #         self.aGrpRadiusOnMax[-1].SetPoint(iPoint, 2.*math.pi*iPoint/nPoint, self.rOnMax[ctEvent][clEvent][binE])
        #     self.aGrpRadiusOffMax.append(ROOT.TGraphPolar(nPoint))
        #     self.aGrpRadiusOffMax[-1].SetName("grpRadiusOffMax%s" % self.aClass[ctEvent])
        #     self.aGrpRadiusOffMax[-1].SetTitle("Maximum radius of the %s OFF region" % self.aClass[ctEvent])
        #     self.aGrpRadiusOffMax[-1].SetLineWidth(1)
        #     self.aGrpRadiusOffMax[-1].SetLineColor(kGray)
        #     self.aGrpRadiusOffMax[-1].SetMarkerStyle(1)
        #     self.aGrpRadiusOffMax[-1].SetMarkerColor(kGray)
        #     self.aGrpRadiusOffMax[-1].SetFillStyle(1001)
        #     self.aGrpRadiusOffMax[-1].SetFillColor(kGray)
        #     for iPoint in range(nPoint):
        #         self.aGrpRadiusOffMax[-1].SetPoint(iPoint, 2*math.pi*iPoint/nPoint, self.rOffMax[ctEvent][clEvent][binE])
        #     self.aGrpRadiusOffMin.append(ROOT.TGraphPolar(nPoint))
        #     self.aGrpRadiusOffMin[-1].SetName("grpRadiusOffMin%s" % self.aClass[ctEvent])
        #     self.aGrpRadiusOffMin[-1].SetTitle("Minimum radius of the %s OFF region" % self.aClass[ctEvent])
        #     self.aGrpRadiusOffMin[-1].SetLineWidth(1)
        #     self.aGrpRadiusOffMin[-1].SetLineColor(kGray+1)
        #     self.aGrpRadiusOffMin[-1].SetMarkerStyle(1)
        #     self.aGrpRadiusOffMin[-1].SetMarkerColor(kGray+1)
        #     self.aGrpRadiusOffMin[-1].SetFillStyle(1001)
        #     self.aGrpRadiusOffMin[-1].SetFillColor(kGray+1)
        #     for iPoint in range(nPoint):
        #         self.aGrpRadiusOffMin[-1].SetPoint(iPoint, 2*math.pi*iPoint/nPoint, self.rOffMin[ctEvent][clEvent][binE])

        aaTitleMap = []
        for ctEvent in range(len(self.aClass)):
            aaTitleMap.append([])
            for clEvent in range(len(self.aaClass[ctEvent])):
                aaTitleMap[-1].append(ROOT.TPaveText(0.1, 1.05, 0.9, 0.9, "NDC"))
                aaTitleMap[-1][-1].AddText("%s map around %s" % (self.aaClass[ctEvent][clEvent], self.name))
                aaTitleMap[-1][-1].SetFillStyle(1001)
                aaTitleMap[-1][-1].SetFillColor(kWhite)
        self.aaTitleMap = aaTitleMap
        Target.draw(self)
        cMap = ROOT.TCanvas("cMap%s" % self.name, "%s map" % self.name, 900, 600)
        nDown = len(self.aaClass)#aaStrSelect)
        nAcross = 0
        aaMgrMap = []
        aaHtgMapDummy = []
        for aS in self.aaClass:#aaStrSelect:    
            nAcross = max(nAcross, len(aS))
        cMap.Divide(nAcross, nDown)
        legGrpMap = ROOT.TLegend(0.1,0.1,0.9,0.9,"Event energy", 'NDC')
        for iaaG in range(len(self.aaaGrpMap)):
            aaMgrMap.append([])
            aaHtgMapDummy.append([])
            for iaG in range(len(self.aaaGrpMap[iaaG])):
                aaMgrMap[-1].append(ROOT.TMultiGraph("mgrMap%s_%s_%s" % (self.name, iaaG, iaG), "%s map around %s" % (self.aaClass[iaaG][iaG], self.name)))
                aaHtgMapDummy[-1].append(ROOT.TH2D("hMapDummy%s_%s_%s" % (self.name, iaaG, iaG), "{0} map around {1};[#circ];[#circ]".format(self.aaClass[iaaG][iaG], self.name), int(2*self.radius), -self.radius, self.radius, int(2*self.radius), -self.radius, self.radius))
                for iG in range(len(self.aaaGrpMap[iaaG][iaG])):
                    aaMgrMap[-1][-1].Add(self.aaaGrpMap[iaaG][iaG][iG])
                    if iaaG==0 and iaG==0:
                        legGrpMap.AddEntry(self.aaaGrpMap[iaaG][iaG][iG], "{0:.1f} - {1:.1f} GeV".format(10**(self.energyPlot.getBin(iG)[0]-3), 10**(self.energyPlot.getBin(iG)[1]-3)), 'p')
                cMap.cd(1+iaaG*nAcross+iaG)
                cMap.cd(1+iaaG*nAcross+iaG).SetTitle("%s map around %s" % (self.aaClass[iaaG][iaG], self.name))
#                self.aGrpRadiusOffMax[iaaG].Draw('NCF')
#                self.aGrpRadiusOffMin[iaaG].Draw('NCF')
#                self.aGrpRadiusOnMax[iaaG].Draw('NCF')
#                if iaaG==0 and iaG==0:
#                    legGrpMap.AddEntry(self.aGrpRadiusOnMax[iaaG], "ON region", 'f')
#                    legGrpMap.AddEntry(self.aGrpRadiusOffMax[iaaG], "OFF region", 'f')
                aaMgrMap[-1][-1].Draw('NP')
                aaTitleMap[iaaG][iaG].Draw('same')
                gPad.Update()
                for iG in range(len(self.aaaGrpMap[iaaG][iaG])):
                    try:
                        self.aaaGrpMap[iaaG][iaG][iG].GetPolargram().SetRangeRadial(0.0, self.radius)
                        self.aaaGrpMap[iaaG][iaG][iG].GetPolargram().SetNdivPolar(104)
                        self.aaaGrpMap[iaaG][iaG][iG].GetPolargram().SetNdivRadial(101)
                    except Exception:
                        print self.aaaGrpMap[iaaG][iaG][iG].GetName(), "Polargram failed."
        cMap.cd(3)
        legGrpMap.Draw()
        self.aaMgrMap = aaMgrMap
        self.cMap = cMap
        self.legGrpMap = legGrpMap
        del aaMgrMap
        del cMap

        cHighEnergy = ROOT.TCanvas("cHighEnergy%s" % self.name, "%s Very high energy photon plot" % self.name, 900, 600)
        mgrHighEnergy = ROOT.TMultiGraph("mgrHighEnergy%s" % self.name, "Very high energy photons from %s (Z=%s)" % (self.name, self.zRedshift))
        legGreHighEnergy = ROOT.TLegend(0.7, 0.85, 0.95, 0.55,"Event class", 'NDC')
        for iaG in range(len(self.aaGreHighEnergy)):
            for iG in range(len(self.aaGreHighEnergy[iaG])):
                mgrHighEnergy.Add(self.aaGreHighEnergy[iaG][iG])
                legGreHighEnergy.AddEntry(self.aaGreHighEnergy[iaG][iG], self.aaClass[iaG][iG], 'p')
        cHighEnergy.cd()
        mgrHighEnergy.Draw("AP")
        mgrHighEnergy.GetXaxis().SetTitle("MJD - %s [day]" % ConvertMetToMjd(self.tStart))
        mgrHighEnergy.GetYaxis().SetTitle("log_{10}(Energy[MeV])")
        legGreHighEnergy.Draw("same")

        cTheta = ROOT.TCanvas("cTheta%s" % self.name, "%s #theta^{2} plots" % self.name, 1200, 600)
        mDown = int(math.sqrt(self.energyPlot.nBin+1))
        mAcross = int(math.ceil((self.energyPlot.nBin+1) / mDown))
        cTheta.Divide(mAcross, mDown)
        aHstaTheta = []
        legTheta = ROOT.TLegend(0.1,0.1,0.9,0.9,"Event class",'NDC')
        for tE in range(self.energyPlot.nBin):
            aHstaTheta.append(ROOT.THStack("hsTheta%s_%s" % (self.name, tE), "#theta^2 plot of {0} in {1:.1f} - {2:.1f} GeV;#theta[#circ]^2;[counts]".format(self.name, 10**(self.energyPlot.getBin(tE)[0]-3), 10**(self.energyPlot.getBin(tE)[1]-3))))
            for iaaH in range(len(self.aaaHtgTheta)):
                for iaH in range(len(self.aaaHtgTheta[iaaH])):
                    self.aaaHtgTheta[iaaH][iaH][tE].Rebin(40)
                    aHstaTheta[-1].Add(self.aaaHtgTheta[iaaH][iaH][tE])
                    if tE==0:
                        legTheta.AddEntry(self.aaaHtgTheta[iaaH][iaH][tE], self.aaClass[iaaH][iaH], 'l')
            cTheta.cd(tE+1)
            aHstaTheta[-1].Draw("nostack")
        cTheta.cd(mAcross*mDown)
        legTheta.Draw()
        self.aHstaTheta = aHstaTheta
        self.legTheta = legTheta
        self.legGreHighEnergy = legGreHighEnergy
        self.cTheta = cTheta
        self.mgrHighEnergy = mgrHighEnergy
        self.cHighEnergy = cHighEnergy
        del aHstaTheta
        del cTheta
        del mgrHighEnergy
        del cHighEnergy

        aCanSpacial = []
        aaaHtgDummySpacial = []

        #legSpacial = ROOT.TLegend(0.1,0.1,0.9,0.9,"Event class",'NDC')
        starSrc = ROOT.TGraph(1)
        starSrc.SetName("Source")
        starSrc.SetTitle("%s" % self.name)
        if abs(self.decCntr)<60:
            starSrc.SetPoint(0, 0, self.decCntr)
        else:
            starSrc.SetPoint(0, 0, self.gbCntr)
        starSrc.SetFillStyle(0)
        starSrc.SetMarkerSize(2)
        starSrc.SetMarkerColor(kMagenta)
        #starSrc.SetLineColorAlpha(kMagenta, 1)
        #starSrc.SetFillColorAlpha(kMagenta, 1)
        starSrc.SetMarkerStyle(29)
        aaaStarRec = []
        legStar = ROOT.TLegend(0.7, 0.8, 0.9, 0.9, "", 'NDC')
        legStar.AddEntry(starSrc)

        lDown = len(self.aaClass)
        lAcross = 0
        for aS in self.aaClass:    
            lAcross = max(lAcross, len(aS))
        for iE in range(self.energyPlot.nBin + 1):
            if iE==self.energyPlot.nBin:
                aCanSpacial.append(ROOT.TCanvas("cSpacialAll%s" % self.name, "Spacial distribution of {0} in {1:.1f} - {2:.1f} GeV".format(self.name,10**(self.energyPlot.getRegionLowEdge()-3), 10**(self.energyPlot.getRegionUpEdge()-3)), 1200, 800))
            else:
                aCanSpacial.append(ROOT.TCanvas("cSpacial%s_%s" % (self.name,iE), "Spacial distribution of {0} in {1:.1f} - {2:.1f} GeV".format(self.name,10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), 1200, 800))
            aCanSpacial[-1].Divide(lAcross, lDown)
            #aaaHtgDummySpacial.append([])
            aaaStarRec.append([])
            for iAS in range(len(self.aaClass)):
                #aaaHtgDummySpacial[-1].append([])
                aaaStarRec[-1].append([])
                for iSS in range(len(self.aaClass[iAS])):
                    #aaaHtgDummySpacial[-1][-1].append(ROOT.TH2D("hDummySpacial%s_%s_%s_%s" % (self.name, iE, iAS, iSS), "{0} Spacial distribution of {1} in {2:.1f} - {3:.1f} GeV;[#circ];[#circ]".format(self.aaClass[iAS][iSS], self.name, 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), 2*self.radius, -self.radius, self.radius, 2*self.radius, -self.radius, self.radius))
                    aCanSpacial[-1].cd(1+iAS*lAcross+iSS)
                    aCanSpacial[-1].cd(1+iAS*lAcross+iSS).SetGridx(1)
                    aCanSpacial[-1].cd(1+iAS*lAcross+iSS).SetGridy(1)
                    #aaaHtgDummySpacial[-1][-1][-1].Draw()
                    if iE==self.energyPlot.nBin:
                        self.aaHtgMapAll[iAS][iSS].SetStats(kFALSE)
                    else:
                        self.aaaHtgMap[iAS][iSS][iE].SetStats(kFALSE)

                    if abs(self.decCntr)<60:
                        if iE==self.energyPlot.nBin:
                            self.aaHtgMapAll[iAS][iSS].GetXaxis().SetRangeUser(-14./cos(abs(radians(self.decCntr))), 14./cos(abs(radians(self.decCntr))))
                            self.aaHtgMapAll[iAS][iSS].GetYaxis().SetRangeUser(self.decCntr-14, self.decCntr+14)
                        else:
                            self.aaaHtgMap[iAS][iSS][iE].GetXaxis().SetRangeUser(-14./cos(abs(radians(self.decCntr))), 14./cos(abs(radians(self.decCntr))))
                            self.aaaHtgMap[iAS][iSS][iE].GetYaxis().SetRangeUser(self.decCntr-14, self.decCntr+14)
                        aaaStarRec[-1][-1].append(ROOT.TGraph(1))
                        aaaStarRec[-1][-1][-1].SetName("starRec%s_%s_%s" %(iE,iAS,iSS))
                        aaaStarRec[-1][-1][-1].SetTitle("Maximum")
                    else:
                        if iE==self.energyPlot.nBin:
                            self.aaHtgMapAll[iAS][iSS].GetXaxis().SetRangeUser(-14./cos(abs(radians(self.gbCntr))), 14./cos(abs(radians(self.gbCntr))))
                            self.aaHtgMapAll[iAS][iSS].GetYaxis().SetRangeUser(self.gbCntr-14, self.gbCntr+14)
                        else:
                            self.aaaHtgMap[iAS][iSS][iE].GetXaxis().SetRangeUser(-14./cos(abs(radians(self.gbCntr))), 14./cos(abs(radians(self.gbCntr))))
                            self.aaaHtgMap[iAS][iSS][iE].GetYaxis().SetRangeUser(self.gbCntr-14, self.gbCntr+14)
                        aaaStarRec[-1][-1].append(ROOT.TGraph(1))
                        aaaStarRec[-1][-1][-1].SetName("starRec%s_%s_%s" %(iE,iAS,iSS))
                        aaaStarRec[-1][-1][-1].SetTitle("Maximum")
                    aaaStarRec[-1][-1][-1].SetFillStyle(0)
                    aaaStarRec[-1][-1][-1].SetMarkerSize(1)
                    aaaStarRec[-1][-1][-1].SetMarkerColor(kBlue)
                    #aaaStarRec[-1][-1][-1].SetLineColorAlpha(kBlue, 1)
                    #aaaStarRec[-1][-1][-1].SetFillColorAlpha(kBlue, 1)
                    aaaStarRec[-1][-1][-1].SetMarkerStyle(21)
                    locmax = ROOT.Long()
                    locmay = ROOT.Long()
                    locmaz = ROOT.Long()
                    if iE==self.energyPlot.nBin:
                        self.aaHtgMapAll[iAS][iSS].GetMaximumBin(locmax, locmay, locmaz)
                        aaaStarRec[-1][-1][-1].SetPoint(0, self.aaHtgMapAll[iAS][iSS].GetXaxis().GetBinCenter(locmax), self.aaHtgMapAll[iAS][iSS].GetYaxis().GetBinCenter(locmay))
                    else:
                        self.aaaHtgMapFolded[iAS][iSS][iE].GetMaximumBin(locmax, locmay, locmaz)
                        aaaStarRec[-1][-1][-1].SetPoint(0, self.aaaHtgMap[iAS][iSS][iE].GetXaxis().GetBinCenter(locmax), self.aaaHtgMap[iAS][iSS][iE].GetYaxis().GetBinCenter(locmay))
                    #print "Maximum point: (", locmax, locmay, ")"
                    if iE==0 and iAS==0 and iSS==0:
                        legStar.AddEntry(aaaStarRec[-1][-1][-1])
                    if iE==self.energyPlot.nBin:
                        self.aaHtgMapAll[iAS][iSS].Draw('aitoff')
                        self.aaHtgMapAll[iAS][iSS].SetMaximum(self.aaHtgMapAll[iAS][0].GetMaximum())
                    else:
                        self.aaaHtgMap[iAS][iSS][iE].Draw('aitoff')
                        self.aaaHtgMap[iAS][iSS][iE].SetMaximum(self.aaaHtgMap[iAS][0][iE].GetMaximum())
                    starSrc.Draw('P same')
                    aaaStarRec[iE][iAS][iSS].Draw('P same')
                    legStar.Draw('same')
                    #self.aaaHtgMap[iAS][iSS][iE].Draw('POL COLZ SAMES')
                    if not iE==self.energyPlot.nBin:
                        self.fRadiusOnMax[iAS][iSS][iE].Draw('same')
                    #self.fRadiusOffMin[iAS][iSS][iE].Draw('same')
                    #self.fRadiusOffMax[iAS][iSS][iE].Draw('same')
                    gPad.Update()
                    #palette = self.aaaHtgMapFolded[iAS][iSS][iE].GetListOfFunctions().FindObject("palette")
                    #try:
                    #    palette.SetX1NDC(0.89)
                    #    palette.SetX2NDC(0.93)
                    #except Exception:
                    #    print "Setting of", palette, "failed."
        self.aCanSpacial = aCanSpacial
        self.aaaHtgDummySpacial = aaaHtgDummySpacial
        self.starSrc = starSrc
        self.aaaStarRec = aaaStarRec
        self.legStar = legStar

    def writeObjects(self):
        Target.writeObjects(self)
        self.cMap.Write("", ROOT.TObject.kOverwrite)
        self.cHighEnergy.Write("", ROOT.TObject.kOverwrite)
        self.cTheta.Write("", ROOT.TObject.kOverwrite)
        self.starSrc.Write("", ROOT.TObject.kOverwrite)
        for aaHtg in self.aaaHtgMap:
            for aHtg in aaHtg:
                for htg in aHtg:
                    htg.Write("", ROOT.TObject.kOverwrite)
        for aHtg in self.aaHtgMapAll:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)
        for aHtg in self.aaHtgMapAllFolded:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)
        for aaHtg in self.aaaHtgPsf:
            for aHtg in aaHtg:
                for htg in aHtg:
                    htg.Write("", ROOT.TObject.kOverwrite)
        for aaHtg in self.aaaHtgMapFolded:
            for aHtg in aaHtg:
                for htg in aHtg:
                    htg.Write("", ROOT.TObject.kOverwrite)
        for aaHtg in self.aaaHtgTheta:
            for aHtg in aaHtg:
                for htg in aHtg:
                    htg.Write("", ROOT.TObject.kOverwrite)
        for aaStar in self.aaaStarRec:
            for aStar in aaStar:
                for star in aStar:
                    star.Write("", ROOT.TObject.kOverwrite)
        for aaGrp in self.aaaGrpMap:
            for aGrp in aaGrp:
                for grp in aGrp:
                    grp.Write("", ROOT.TObject.kOverwrite)
        for aGr in self.aaGreHighEnergy:
            for gr in aGr:
                gr.Write("", ROOT.TObject.kOverwrite)

        #for iE in range(self.energyPlot.nBin):
            #self.aCanSpacial[iE].Write("", ROOT.TObject.kOverwrite)
        for cS in self.aCanSpacial:
            cS.Write("", ROOT.TObject.kOverwrite)
        for aaFc in self.fRadiusOnMax:
            for aFc in aaFc:
                for fc in aFc:
                    fc.Write("", ROOT.TObject.kOverwrite)                

class EarthLimb(Target):
    def __init__(self, strName, zOff1Min=[90, 90], zOff1Max=[100, 100],zOff2Min=[120, 120], zOff2Max=[130, 130], zOnMin=[111.1002, 108.1002], zOnMax=[112.9545, 115.9545], config = ClassConfig(), eRegion=EnergyLogRegion(3,4.75,0.25), ePlotRegion=EnergyLogRegion(3, 4.75, 0.25), perf=["/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R100_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R30_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R10_perf.root"]):
        Target.__init__(self, strName, config, eRegion, ePlotRegion)
        self.perf = perf
        self.zOnMin = []
        self.zOnMax = []
        self.zOff1Min = []
        self.zOff1Max = []
        self.zOff2Min = []
        self.zOff2Max = []
        self.saOn = []
        self.saOff = []
        for hS in range(len(self.aaClass)):
            self.zOnMin.append([])
            self.zOnMax.append([])
            self.zOff1Min.append([])
            self.zOff1Max.append([])
            self.zOff2Min.append([])
            self.zOff2Max.append([])
            self.saOn.append([])
            self.saOff.append([])
            for hSS in range(len(self.aaClass[hS])):
                self.zOnMin[hS].append([])
                self.zOnMax[hS].append([])
                self.zOff1Min[hS].append([])
                self.zOff1Max[hS].append([])
                self.zOff2Min[hS].append([])
                self.zOff2Max[hS].append([])
                self.saOn[hS].append([])
                self.saOff[hS].append([])
                for ie in range(self.energySet.nBin):
                    if hS==0:
                        self.zOnMin[hS][hSS].append(111.10)
                        self.zOnMax[hS][hSS].append(112.95)
                        self.zOff1Min[hS][hSS].append(108.66)
                        self.zOff1Max[hS][hSS].append(109.57)
                        self.zOff2Min[hS][hSS].append(114.52)
                        self.zOff2Max[hS][hSS].append(115.47)
                    elif hS==1:
                        self.zOnMin[hS][hSS].append(111.10+0.1-self.perf.getPSF68(hS, hSS, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.zOnMax[hS][hSS].append(112.95-0.1+self.perf.getPSF68(hS, hSS, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.zOff1Min[hS][hSS].append(95)
                        self.zOff1Max[hS][hSS].append(111.10-self.perf.getPSF95(hS, hSS, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.zOff2Min[hS][hSS].append(112.95+self.perf.getPSF95(hS, hSS, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.zOff2Max[hS][hSS].append(125)
                    self.saOn[hS][hSS].append( 2.0 * math.pi * ( cos(radians(self.zOnMin[hS][hSS][-1])) - cos(radians(self.zOnMax[hS][hSS][-1])) ) )
                    self.saOff[hS][hSS].append( 2.0 * math.pi * ( cos(radians(self.zOff1Min[hS][hSS][-1])) - cos(radians(self.zOff1Max[hS][hSS][-1])) ) + 2.0 * math.pi * ( cos(radians(self.zOff2Min[hS][hSS][-1])) - cos(radians(self.zOff2Max[hS][hSS][-1])) ))
        print "Solid angle of ON region:", self.saOn
        print "Solid angle of OFF region:", self.saOff
        aaHtgZenithTheta = []

    def fill(self, raEvent, decEvent, lEvent, bEvent, eEvent, ctEvent, clEvent, zEvent, tEvent, cthEvent):
        self.flagOn[0] = kFALSE
        self.flagOff[0] = kFALSE
        binE = self.energySet.findBin(eEvent)
        plotE = self.energyPlot.findBin(eEvent)
        if binE>=0 and binE<self.energySet.nBin and tEvent>=self.tStart and tEvent<= self.tEnd:
            zOff1Min = self.zOff1Min[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]
            zOff1Max = self.zOff1Max[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]
            zOff2Min = self.zOff2Min[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]
            zOff2Max = self.zOff2Max[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]
            zOnMax = self.zOnMax[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]
            zOnMin = self.zOnMin[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]

            if zEvent>=zOff1Min:
                for clEventPlus in range(clEvent-int(ctEvent==1 and clEvent==3)):
                    # for FITS
                    self.npaaENERGY[ctEvent-1][clEventPlus] = np.append(self.npaaENERGY[ctEvent-1][clEventPlus], eEvent)
                    self.npaaRA[ctEvent-1][clEventPlus] = np.append(self.npaaRA[ctEvent-1][clEventPlus], raEvent)
                    self.npaaDEC[ctEvent-1][clEventPlus] = np.append(self.npaaDEC[ctEvent-1][clEventPlus], decEvent)
                    self.npaaL[ctEvent-1][clEventPlus] = np.append(self.npaaL[ctEvent-1][clEventPlus], lEvent)
                    self.npaaB[ctEvent-1][clEventPlus] = np.append(self.npaaB[ctEvent-1][clEventPlus], bEvent)
                    self.npaaTHETA[ctEvent-1][clEventPlus] = np.append(self.npaaTHETA[ctEvent-1][clEventPlus], cthEvent)
                    self.npaaZENITH_ANGLE[ctEvent-1][clEventPlus] = np.append(self.npaaZENITH_ANGLE[ctEvent-1][clEventPlus], zEvent)
                    self.npaaTIME[ctEvent-1][clEventPlus] = np.append(self.npaaTIME[ctEvent-1][clEventPlus], tEvent)
                    if ctEvent==1:
                        self.npaaEVENT_CLASS[ctEvent-1][clEventPlus] = np.append(self.npaaEVENT_CLASS[ctEvent-1][clEventPlus], 128*2**(clEvent-1)/2)
                    elif ctEvent==2:
                        self.npaaEVENT_CLASS[ctEvent-1][clEventPlus] = np.append(self.npaaEVENT_CLASS[ctEvent-1][clEventPlus], 4096*2**(clEvent-1))
                    if cthEvent<0.7:
                        self.npaaEVENT_TYPE[ctEvent-1][clEventPlus] = np.append(self.npaaEVENT_TYPE[ctEvent-1][clEventPlus], 1)
                    else:
                        self.npaaEVENT_TYPE[ctEvent-1][clEventPlus] = np.append(self.npaaEVENT_TYPE[ctEvent-1][clEventPlus], 2)

                    if zEvent>=zOnMin and zEvent<zOnMax:
                        self.flagOn[0] = kTRUE
                        self.aaHtgNumOn[ctEvent-1][clEventPlus].Fill(eEvent)
                        self.aaHtgEnergyOn[ctEvent-1][clEventPlus].Fill(eEvent)
                        self.aaHtgLCCountOn[ctEvent-1][clEventPlus].Fill(ConvertMetToMjd(tEvent)-ConvertMetToMjd(self.tStart))
                    elif (zEvent>=zOff1Min and zEvent<zOff1Max) or (zEvent>=zOff2Min and zEvent<zOff2Max):
                        self.flagOff[0] = kTRUE
                        self.aaHtgNumOff[ctEvent-1][clEventPlus].Fill(eEvent)
                        self.aaHtgEnergyOff[ctEvent-1][clEventPlus].Fill(eEvent)
                        self.aaHtgLCCountOff[ctEvent-1][clEventPlus].Fill(ConvertMetToMjd(tEvent)-ConvertMetToMjd(self.tStart))
        self.trObj.Fill()


class GalacticRidge(Target):
    def __init__(self, strName, bOffMin=[50., 50.], lOffMin=[90., 90.], lOffMax=[-90., -90.], bOnMax=[1.5, 3.0], lOnMin=[-50., -51.5], lOnMax=[40., 41.5], config = ClassConfig(), eRegion=EnergyLogRegion(3,4.75,0.25), ePlotRegion=EnergyLogRegion(3, 4.75, 0.25), perf=["/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R100_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R30_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R10_perf.root"], tStart=252460800.0, tEnd=504921604.0):
        Target.__init__(self, strName, config, eRegion, ePlotRegion, tStart, tEnd)
        self.perf = perf
        self.HEALPLOT_MOLLWEIDE = True
        self.loncntr = 0.0
        self.latcntr = 0.0
        self.coord = ('E', 'G')
        self.saOn = []
        self.saOff = []
        self.bOffMin = []#bOffMin
        self.lOffMin = []#lOffMin
        self.lOffMax = []#lOffMax
        self.bOnMax = []#bOnMax
        self.lOnMin = []#lOnMin
        self.lOnMax = []#lOnMax
        for hS in range(len(self.aaClass)):
            self.bOffMin.append([])
            self.lOffMin.append([])
            self.lOffMax.append([])
            self.bOnMax.append([])
            self.lOnMin.append([])
            self.lOnMax.append([])
            self.saOn.append([])
            self.saOff.append([])
            for hSS in range(len(self.aaClass[hS])):
                self.bOffMin[hS].append([])
                self.lOffMin[hS].append([])
                self.lOffMax[hS].append([])
                self.bOnMax[hS].append([])
                self.lOnMin[hS].append([])
                self.lOnMax[hS].append([])
                self.saOn[hS].append([])
                self.saOff[hS].append([])
                for ie in range(self.energySet.nBin):
                    if hS==0:
                        self.bOffMin[hS][hSS].append(50.)
                        self.lOffMin[hS][hSS].append(90.)
                        self.lOffMax[hS][hSS].append(-90.)
                        self.bOnMax[hS][hSS].append(1.5)
                        self.lOnMin[hS][hSS].append(-50.)
                        self.lOnMax[hS][hSS].append(40.)
                    elif hS==1:
                        self.bOffMin[hS][hSS].append(self.bOffMin[0][0][ie])
                        self.lOffMin[hS][hSS].append(self.lOffMin[0][0][ie])
                        self.lOffMax[hS][hSS].append(self.lOffMax[0][0][ie])
                        #print hSS, (self.energySet.aBin[ie]+self.energySet.wBin/2.0)
                        self.bOnMax[hS][hSS].append(self.bOnMax[0][0][ie]-0.1+self.perf.getPSF68(hS, hSS, (self.energySet.aBin[ie]+self.energySet.wBin/2.0)) )
                        self.lOnMin[hS][hSS].append(self.lOnMin[0][0][ie]+0.1-self.perf.getPSF68(hS, hSS, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.lOnMax[hS][hSS].append(self.lOnMax[0][0][ie]-0.1+self.perf.getPSF68(hS, hSS, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                    self.saOn[hS][hSS].append( 2.0 * 2.0 * math.pi * ( cos(radians(90.0-self.bOnMax[hS][hSS][-1])) - cos(radians(90.0)) ) * (self.lOnMax[hS][hSS][-1]-self.lOnMin[hS][hSS][-1])/360. )
                    self.saOff[hS][hSS].append( 2.0 * 2.0 * math.pi * ( cos(radians(0.0)) - cos(radians(90.0-self.bOffMin[hS][hSS][-1])) ) * (self.lOffMax[hS][hSS][-1]+180.+180.-self.lOffMin[hS][hSS][-1])/360. )
#            self.saOn.append( 2.0 * 2.0 * math.pi * ( cos(radians(90.0-self.bOnMax[hS])) - cos(radians(90.0)) ) * (self.lOnMax[hS]-self.lOnMin[hS])/360. )
#            self.saOff.append( 2.0 * 2.0 * math.pi * ( cos(radians(0.0)) - cos(radians(90.0-self.bOffMin[hS])) ) * (self.lOffMax[hS]+180.+180.-self.lOffMin[hS])/360. )
        print "Solid angle of ON region:", self.saOn
        print "Solid angle of OFF region:", self.saOff
        aaaHtgLatitude = []
        aaaHtgMapAllSky = []
        aaaHtgMapAllCel = []
        aaaHtgFoldedMapAllSky = []
        aaaHtgFoldedMapAllCel = []
        for iS in range(len(self.aaClass)):
            aaaHtgLatitude.append([])
            aaaHtgMapAllSky.append([])
            aaaHtgMapAllCel.append([])
            aaaHtgFoldedMapAllSky.append([])
            aaaHtgFoldedMapAllCel.append([])
            for jS in range(len(self.aaClass[iS])):
                aaaHtgLatitude[-1].append([])
                aaaHtgMapAllSky[-1].append([])
                aaaHtgFoldedMapAllSky[-1].append([])
                aaaHtgMapAllCel[-1].append([])
                aaaHtgFoldedMapAllCel[-1].append([])
                for iE in range(self.energyPlot.nBin):
                    aaaHtgLatitude[-1][-1].append(ROOT.TH1D("hLatitude%s_%s_%s" % (iS,jS,iE), "{0} Galactic latitude distribution in {1:.1f} - {2:.1f} GeV;sin(b);counts [counts]".format(self.aaClass[iS][jS], 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), 200, -1, 1))
                    aaaHtgLatitude[-1][-1][-1]
                    aaaHtgLatitude[-1][-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                    aaaHtgLatitude[-1][-1][-1].SetLineWidth(3-iS)
                    aaaHtgLatitude[-1][-1][-1].SetLineStyle(2-iS)
                        
                    aaaHtgMapAllSky[-1][-1].append(ROOT.TH2D("hMapAllSky%s_%s_%s" % (iS,jS,iE), "{0} all sky map in {1:.1f} - {2:.1f} GeV;- Galactic longitude [#circ];Galactic latitude [#circ]".format(self.aaClass[iS][jS], 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), 400, -200, 200, 180, -90, 90))
                    #aaaHtgMapAllSky[-1][-1].append(ROOT.TH2D("hMapAllSky%s_%s_%s" % (iS,jS,iE), "{0} all sky map in {1:.1f} - {2:.1f} GeV;- Galactic longitude [#circ];Galactic latitude [#circ]".format(self.aaClass[iS][jS], 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), 360, -180, 180, 180, -90, 90))
                    aaaHtgFoldedMapAllSky[-1][-1].append(ROOT.TH2D("hFoldedMapAllSky%s_%s_%s" % (iS,jS,iE), "{0} folded all sky map in {1:.1f} - {2:.1f} GeV;- Galactic longitude [#circ];Galactic latitude [#circ]".format(self.aaClass[iS][jS], 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), 360, -180, 180, 140, -70, 70))
                    aaaHtgMapAllCel[-1][-1].append(ROOT.TH2D("hMapAllCel%s_%s_%s" % (iS,jS,iE), "{0} all celestial map in {1:.1f} - {2:.1f} GeV;- RA [#circ];DEC [#circ]".format(self.aaClass[iS][jS], 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), 400, -200, 200, 140, -70, 70))
                    #aaaHtgMapAllCel[-1][-1].append(ROOT.TH2D("hMapAllCel%s_%s_%s" % (iS,jS,iE), "{0} all celestial map in {1:.1f} - {2:.1f} GeV;- RA [#circ];DEC [#circ]".format(self.aaClass[iS][jS], 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), 360, -180, 180, 140, -70, 70))
                    aaaHtgFoldedMapAllCel[-1][-1].append(aaaHtgMapAllCel[-1][-1][-1].Clone("hFoldedMapAllCel%s_%s_%s" % (iS, jS, iE)))
                        #ROOT.TH2D("hFolded%s_%s_%s_%s" % (self.name, iS, jS, iE), "{0} PSF folded map around {1} in {2:.1f} - {3:.1f} GeV;RA[#circ];DEC[#circ];[counts]".format(self.aaClass[iS][jS], self.name, 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), 12,))
        self.aaaHtgLatitude = aaaHtgLatitude
        self.aaaHtgMapAllSky = aaaHtgMapAllSky
        self.aaaHtgFoldedMapAllSky = aaaHtgFoldedMapAllSky
        self.aaaHtgMapAllCel = aaaHtgMapAllCel
        self.aaaHtgFoldedMapAllCel = aaaHtgFoldedMapAllCel

        aaGrPerformance = []
        for tE in range(self.energyPlot.nBin):
            aaGrPerformance.append([])
            for kAS in range(len(self.aClass)):
                aaGrPerformance[-1].append(ROOT.TGraphErrors(len(self.aaClass[kAS])))
                aaGrPerformance[-1][-1].SetName("grPerformance_%s_%s_%s" % (self.name, tE, kAS))
                aaGrPerformance[-1][-1].SetTitle("{0} performance plot of {1} in {2:.1f} - {3:.1f} GeV".format(self.aClass[kAS], self.name, 10**(self.energyPlot.getBin(tE)[0]-3), 10**(self.energyPlot.getBin(tE)[1]-3)))
                aaGrPerformance[-1][-1].GetXaxis().SetTitle("Number of signal events")
                aaGrPerformance[-1][-1].GetYaxis().SetTitle("Number of residual OFF events")
                aaGrPerformance[-1][-1].SetMarkerStyle(25-kAS*5)
                aaGrPerformance[-1][-1].SetLineStyle(2-kAS)
                aaGrPerformance[-1][-1].SetLineWidth(2)
        self.aaGrPerformance = aaGrPerformance
        #print self.aaGrPerformance
        print "Construction finished."

    def fill(self, raEvent, decEvent, lEvent, bEvent, eEvent, ctEvent, clEvent, zEvent, tEvent, cthEvent):
        self.flagOn[0] = kFALSE
        self.flagOff[0] = kFALSE
        clEvent = int(clEvent)
        binE = self.energySet.findBin(eEvent)
        plotE = self.energyPlot.findBin(eEvent)
        lEventCor = lEvent
        raEventCor = raEvent
        if binE>=0 and binE<self.energySet.nBin and tEvent>=self.tStart and tEvent<= self.tEnd:
            if lEvent>180:
                lEventCor = -360+lEvent
            if raEvent>180:
                raEventCor = -360+raEvent
            bOffMin = self.bOffMin[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]
            lOffMin = self.lOffMin[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]
            lOffMax = self.lOffMax[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]
            bOnMax = self.bOnMax[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]
            lOnMin = self.lOnMin[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]
            lOnMax = self.lOnMax[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]
            zCut = self.zCut[ctEvent-1]
            if zEvent<zCut:
            #if binE>=0 and zEvent<zCut:
                for clEventPlus in range(clEvent-int(ctEvent==1 and clEvent==3)):
                    self.aaaHtgMapAllSky[ctEvent-1][clEventPlus][plotE].Fill(-lEventCor, bEvent)
                    self.aaaHtgMapAllCel[ctEvent-1][clEventPlus][plotE].Fill(-raEventCor, decEvent)
                    self.aaaHtgLatitude[ctEvent-1][clEventPlus][plotE].Fill(sin(radians(bEvent)))
                    # for FITS
                    self.npaaENERGY[ctEvent-1][clEventPlus] = np.append(self.npaaENERGY[ctEvent-1][clEventPlus], eEvent)
                    self.npaaRA[ctEvent-1][clEventPlus] = np.append(self.npaaRA[ctEvent-1][clEventPlus], raEvent)
                    self.npaaDEC[ctEvent-1][clEventPlus] = np.append(self.npaaDEC[ctEvent-1][clEventPlus], decEvent)
                    self.npaaL[ctEvent-1][clEventPlus] = np.append(self.npaaL[ctEvent-1][clEventPlus], lEvent)
                    self.npaaB[ctEvent-1][clEventPlus] = np.append(self.npaaB[ctEvent-1][clEventPlus], bEvent)
                    self.npaaTHETA[ctEvent-1][clEventPlus] = np.append(self.npaaTHETA[ctEvent-1][clEventPlus], cthEvent)
                    self.npaaZENITH_ANGLE[ctEvent-1][clEventPlus] = np.append(self.npaaZENITH_ANGLE[ctEvent-1][clEventPlus], zEvent)
                    self.npaaTIME[ctEvent-1][clEventPlus] = np.append(self.npaaTIME[ctEvent-1][clEventPlus], tEvent)
                    if ctEvent==1:
                        self.npaaEVENT_CLASS[ctEvent-1][clEventPlus] = np.append(self.npaaEVENT_CLASS[ctEvent-1][clEventPlus], 128*2**(clEvent-1)/2)
                    elif ctEvent==2:
                        self.npaaEVENT_CLASS[ctEvent-1][clEventPlus] = np.append(self.npaaEVENT_CLASS[ctEvent-1][clEventPlus], 4096*2**(clEvent-1))
                    if cthEvent<0.7:
                        self.npaaEVENT_TYPE[ctEvent-1][clEventPlus] = np.append(self.npaaEVENT_TYPE[ctEvent-1][clEventPlus], 1)
                    else:
                        self.npaaEVENT_TYPE[ctEvent-1][clEventPlus] = np.append(self.npaaEVENT_TYPE[ctEvent-1][clEventPlus], 2)

                    if lEventCor>lOnMin and lEventCor<lOnMax and bEvent>-bOnMax and bEvent<bOnMax:
                        self.flagOn[0] = kTRUE
                        self.aaHtgNumOn[ctEvent-1][clEventPlus].Fill(eEvent)
                        self.aaHtgEnergyOn[ctEvent-1][clEventPlus].Fill(eEvent)#, (10**(eEvent-3))**2)
                        self.aaHtgLCCountOn[ctEvent-1][clEventPlus].Fill(ConvertMetToMjd(tEvent)-ConvertMetToMjd(self.tStart))
                    elif (lEventCor<lOffMax and (bEvent>bOffMin or bEvent<-bOffMin)) or (lEventCor>lOffMin and (bEvent>bOffMin or bEvent<-bOffMin)):
                        self.flagOff[0] = kTRUE
                        self.aaHtgNumOff[ctEvent-1][clEventPlus].Fill(eEvent)
                        self.aaHtgEnergyOff[ctEvent-1][clEventPlus].Fill(eEvent)#, (10**(eEvent-3))**2)
                        self.aaHtgLCCountOff[ctEvent-1][clEventPlus].Fill(ConvertMetToMjd(tEvent)-ConvertMetToMjd(self.tStart))
#    def __init__(self, strName, bOffMin=[50., 50.], lOffMin=[90., 90.], lOffMax=[-90., -90.], bOnMax=[1.5, 5.5], lOnMin=[-50., -54], lOnMax=[40., 44.], config = ClassConfig(), self.energySet=EnergyLogRegion(3,4.75,0.25)):
        self.trObj.Fill()

    def calc(self):
        Target.calc(self)
        #fPSF = ROOT.TF1("fPSF", "1. / (2*TMath::Pi())**0.5 / [0] * TMath::Exp(-(x/[0])**2/2.0)")
        for tE in range(self.energyPlot.nBin):
            for kAS in range(len(self.aClass)):
                for kSS in range(len(self.aaClass[kAS])):
                    #if kAS==0:
                    #    degPSF = 0.1
                    #elif kAS==1:
                    #    degPSF = max(self.perf.getPSF68(kSS, self.energySet.getBinCenter(tE)), self.perf.getPSF95(kSS, self.energySet.getBinCenter(tE))/2.0)
                    #fPSF.FixParameter(0, radians(degPSF))
                    self.aaGrPerformance[tE][kAS].SetPoint(kSS, self.aaHtgNumSig[kAS][kSS].GetBinContent(tE+1), self.aaHtgNumOff[kAS][kSS].GetBinContent(tE+1))
                    self.aaGrPerformance[tE][kAS].SetPointError(kSS, self.aaHtgNumSig[kAS][kSS].GetBinError(tE+1), self.aaHtgNumOff[kAS][kSS].GetBinError(tE+1))
        #print self.aaGrPerformance
        print "Calculation finished."


    def draw(self):
        Target.draw(self)
        if not self.lstt_hp_htg==None:
            for (icat,lstcat) in enumerate(self.aaClass):
                for (icla, cla) in enumerate(lstcat):
                    htg3D = self.lstt_hp_htg[icat][icla]
                    hpcube = Healcube('AllSky', cla, htg3D, lon_cntr=0, lat_cntr=0, char_coord=('E','G'))
                    hpcube.setmap()
                    hpcube.draw(mollweide=True)

#         if not lstt_hp_htg==None:
#             for binE in range(self.energyPlot.nBin):
#                 for ctEvent in range(len(self.aaClass)):
#                     for clEvent in range(len(self.aaClass[ctEvent])):
#                         hp.visufunc.mollview(lstt_hp_htg[binE][ctEvent][clEvent], min=0, flip='astro', title="{0} all sky map ({1:.1f} - {2:.1f} GeV)".format(self.aaClass[ctEvent][clEvent], 10**(self.energyPlot.getBin(binE)[0]-3), 10**(self.energyPlot.getBin(binE)[1]-3)), coord=('E','G'), unit='counts')
#                         plt.savefig("{0}_{1}_E{2}-{3}.png".format(self.name
#                         plt.clf()
# , self.aaClass[ctEvent][clEvent], int(100*self.energyPlot.getBin(binE)[0]+0.5), int(100*self.energyPlot.getBin(binE)[1]+0.5)))

        aCanAllSky = []
        nDown = len(self.aaClass)
        nAcross = 0
        aaaHtgMapAllSky = []
        for aS in self.aaClass:    
            nAcross = max(nAcross, len(aS))
        cLatitude = ROOT.TCanvas("cLatitude", "Galactic latitude", 1200, 600)
        mDown = int(math.sqrt(self.energyPlot.nBin+1))
        mAcross = int(math.ceil((self.energyPlot.nBin+1) / mDown))
        cLatitude.Divide(mAcross, mDown)
        aHstaLatitude = []
        legLatitude = ROOT.TLegend(0.1,0.1,0.9,0.9,"Event class",'NDC')

        for iE in range(self.energyPlot.nBin):
            aCanAllSky.append(ROOT.TCanvas("cAllSky%s" % iE, "All sky map in {0:.1f} - {1:.1f} GeV".format(10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), 1200, 400))
            aCanAllSky[-1].Divide(nAcross, nDown)
            aaaHtgMapAllSky.append([])
            aHstaLatitude.append(ROOT.THStack("hsLatitude%s" % iE, "Galactic latitude in {0:.1f} - {1:.1f} GeV;sin(b);counts [counts]".format(10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3))))
            for iAS in range(len(self.aaClass)):
                aaaHtgMapAllSky[-1].append([])
                for iSS in range(len(self.aaClass[iAS])):
                    aCanAllSky[-1].cd(1+iAS*nAcross+iSS)
                    aCanAllSky[-1].cd(1+iAS*nAcross+iSS).SetGridx(0)
                    aCanAllSky[-1].cd(1+iAS*nAcross+iSS).SetGridy(0)
                    self.aaaHtgMapAllSky[iAS][iSS][iE].GetXaxis().SetRangeUser(-180, 180)
                    self.aaaHtgMapAllSky[iAS][iSS][iE].Draw('aitoff Z')
                   # ch.Draw("Latitude:-(Longitude<180?Longitude:(-360+Longitude))>>%s" % self.aaaHtgMapAllSky[iAS][iSS][iE].GetName(), "Energy>=%s && Energy<%s && ZenithAngle<100 && Category==%s && Class>=%s" % (self.energySet.getBin(iE)[0], self.energySet.getBin(iE)[1], iAS+1, iSS+1+int(iAS==0 and iSS==1)), 'aitoff')
                    self.aaaHtgMapAllSky[iAS][iSS][iE].SetMinimum(0)
                    self.aaaHtgMapAllSky[iAS][iSS][iE].SetMaximum(self.aaaHtgMapAllSky[0][0][iE].GetMaximum())
                    aHstaLatitude[iE].Add(self.aaaHtgLatitude[iAS][iSS][iE])
                    if iE==0:
                        legLatitude.AddEntry(self.aaaHtgLatitude[iAS][iSS][iE], self.aaClass[iAS][iSS], 'l')
            cLatitude.cd(iE+1)
            aHstaLatitude[iE].Draw('nostack')
        # for iAS in range(len(self.aaClass)):
        #     for iSS in range(len(self.aaClass[iAS])):
        #         hsEnergy.Add(self.aaHtgEnergy[iAS][iSS])
      #  cEnergy.cd()
      #  hsEnergy.Draw('nostack')
        cLatitude.cd(4)
        legLatitude.Draw()
        self.aCanAllSky = aCanAllSky
        self.legLatitude = legLatitude
        self.cLatitude = cLatitude
      #  self.cEnergy = cEnergy

        aMgrPerformance = []
        cPerformance = ROOT.TCanvas("cPerformance_%s" % self.name, "Performance plots of %s" % self.name, 1200, 800)# in {1:.1f} - {2:.1f} GeV".format(self.name, 10**(self.energySet.getBin(tE)[0]-3), 10**(self.energySet.getBin(tE)[1]-3)), 600, 400)
        lDown = int(math.sqrt(self.energyPlot.nBin+1))
        lAcross = int(math.ceil((self.energyPlot.nBin+1) / lDown))
        cPerformance.Divide(lAcross, lDown)
        legPerformance = ROOT.TLegend(0.1,0.1,0.9,0.9,"Event category",'NDC')
        for tE in range(self.energyPlot.nBin):
            aMgrPerformance.append(ROOT.TMultiGraph("mgrPerformance_%s_%s" % (self.name, tE), "Performance plots with {0} in {1:.1f} - {2:.1f} GeV".format(self.name, 10**(self.energyPlot.getBin(tE)[0]-3), 10**(self.energyPlot.getBin(tE)[1]-3))))
            for kAS in range(len(self.aClass)):
                aMgrPerformance[tE].Add(self.aaGrPerformance[tE][kAS])
                if tE==0:
                    legPerformance.AddEntry(self.aaGrPerformance[tE][kAS], self.aClass[kAS], 'lp')
            cPerformance.cd(tE+1)
            aMgrPerformance[tE].Draw("APL")
            aMgrPerformance[tE].GetXaxis().SetTitle("Number of signal events [counts]")
            aMgrPerformance[tE].GetYaxis().SetTitle("Number of OFF events [counts]")
        cPerformance.cd(4)
        legPerformance.Draw()
        self.aMgrPerformance = aMgrPerformance
        self.cPerformance = cPerformance
        self.aHstaLatitude = aHstaLatitude
        self.legPerformance = legPerformance
        #print self.aaGrPerformance
        print "Drawing finished."

    def writeObjects(self):
        Target.writeObjects(self)
        self.cLatitude.Write("", ROOT.TObject.kOverwrite)
        self.cPerformance.Write("", ROOT.TObject.kOverwrite)
        for iE in range(self.energyPlot.nBin):
            self.aCanAllSky[iE].Write("", ROOT.TObject.kOverwrite)
            for iAS in range(len(self.aaClass)):
                for iSS in range(len(self.aaClass[iAS])):
                    self.aaaHtgMapAllSky[iAS][iSS][iE].Write("", ROOT.TObject.kOverwrite)
                    self.aaaHtgFoldedMapAllSky[iAS][iSS][iE].Write("", ROOT.TObject.kOverwrite)
                    self.aaaHtgMapAllCel[iAS][iSS][iE].Write("", ROOT.TObject.kOverwrite)
                    self.aaaHtgFoldedMapAllCel[iAS][iSS][iE].Write("", ROOT.TObject.kOverwrite)
                    self.aaaHtgLatitude[iAS][iSS][iE].Write("", ROOT.TObject.kOverwrite)
                try:
                    #print "Writing", self.aaGrPerformance[iE][iAS].GetName()
                    self.aaGrPerformance[iE][iAS].Write()
                except Exception:
                    print "Writing failed."
                #else:
                    #print "Writing suceeded."
        print "Writing finished."

class InnerGalaxy(Target):
    def __init__(self, strName, radius=41., config = ClassConfig(), eRegion=EnergyLogRegion(3,4.75,0.25), ePlotRegion=EnergyLogRegion(3, 4.75, 0.25), perf=["/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R100_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R30_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R10_perf.root"], tStart=252460800.0, tEnd=504921604.0):
        Target.__init__(self, strName, config, eRegion, ePlotRegion, tStart, tEnd)
        self.perf = perf
        self.saOn = []
        self.saOff = []
        self.bOffMin = []
        self.lOffMin = []
        self.lOffMax = []
        self.rOnMax = []
        self.lMaskMin = []
        self.bMaskMax = []
        for hS in range(len(self.aaClass)):
            self.bOffMin.append([])
            self.lOffMin.append([])
            self.lOffMax.append([])
            self.rOnMax.append([])
            self.lMaskMin.append([])
            self.bMaskMax.append([])
            self.saOn.append([])
            self.saOff.append([])
            for hSS in range(len(self.aaClass[hS])):
                self.bOffMin[hS].append([])
                self.lOffMin[hS].append([])
                self.lOffMax[hS].append([])
                self.rOnMax[hS].append([])
                self.lMaskMin[hS].append([])
                self.bMaskMax[hS].append([])
                self.saOn[hS].append([])
                self.saOff[hS].append([])
                for ie in range(self.energySet.nBin):
                    if hS==0:
                        self.bOffMin[hS][hSS].append(50.)
                        self.lOffMin[hS][hSS].append(90.)
                        self.lOffMax[hS][hSS].append(-90.)
                        self.rOnMax[hS][hSS].append(radius)
                        self.lMaskMin[hS][hSS].append(5.)
                        self.bMaskMax[hS][hSS].append(6.)
                    elif hS==1:
                        self.bOffMin[hS][hSS].append(self.bOffMin[0][0][ie])
                        self.lOffMin[hS][hSS].append(self.lOffMin[0][0][ie])
                        self.lOffMax[hS][hSS].append(self.lOffMax[0][0][ie])
                        self.rOnMax[hS][hSS].append(radius)
                        self.lMaskMin[hS][hSS].append(0.)
                        self.bMaskMax[hS][hSS].append(self.bMaskMax[0][0][ie]+self.perf.getPSF68(hS, hSS, self.energySet.getBinCenter(ie)))
                    self.saOn[hS][hSS].append( 2.0*math.pi*(cos(0.0)-cos(radians(self.rOnMax[hS][hSS][-1]))) - (2.0 * 2.0 * math.pi * ( cos(radians(90.0-self.bMaskMax[hS][hSS][-1])) - cos(radians(90.0)) ) * (self.rOnMax[hS][hSS][-1]-self.lMaskMin[hS][hSS][-1])*2./360.) )
                    self.saOff[hS][hSS].append( 2.0 * 2.0 * math.pi * ( cos(radians(0.0)) - cos(radians(90.0-self.bOffMin[hS][hSS][-1])) ) * (self.lOffMax[hS][hSS][-1]+180.+180.-self.lOffMin[hS][hSS][-1])/360. )
        print "Solid angle of ON region:", self.saOn
        print "Solid angle of OFF region:", self.saOff
        aaaHtgMap = []
        for iS in range(len(self.aaClass)):
            aaaHtgMap.append([])
            for jS in range(len(self.aaClass[iS])):
                aaaHtgMap[-1].append([])
                for iE in range(self.energyPlot.nBin):
                    aaaHtgMap[-1][-1].append(ROOT.TH2D("hMap%s_%s_%s_%s" % (self.name, iS,jS,iE), "{0} all sky map in {1:.1f} - {2:.1f} GeV;- Galactic longitude [#circ];Galactic latitude [#circ]".format(self.aaClass[iS][jS], 10**(self.energyPlot.getBin(iE)[0]-3), 10**(self.energyPlot.getBin(iE)[1]-3)), 360, -180, 180, 180, -90, 90))
        self.aaaHtgMap = aaaHtgMap
        print "Construction finished."

    def fill(self, raEvent, decEvent, lEvent, bEvent, eEvent, ctEvent, clEvent, zEvent, tEvent):
        clEvent = int(clEvent)
        binE = self.energySet.findBin(eEvent)
        plotE = self.energyPlot.findBin(eEvent)
        if binE>=0 and binE<self.energySet.nBin and tEvent>=self.tStart and tEvent<= self.tEnd:
            if lEvent>180:
                lEvent = -360+lEvent
            zCut = self.zCut[ctEvent-1]
            if zEvent<zCut:
                for clEventPlus in range(clEvent-int(ctEvent==1 and clEvent==3)):
                    bOffMin = self.bOffMin[ctEvent-1][clEventPlus][binE]
                    lOffMin = self.lOffMin[ctEvent-1][clEventPlus][binE]
                    lOffMax = self.lOffMax[ctEvent-1][clEventPlus][binE]
                    rOnMax = self.rOnMax[ctEvent-1][clEventPlus][binE]
                    lMaskMin = self.lMaskMin[ctEvent-1][clEventPlus][binE]
                    bMaskMax = self.bMaskMax[ctEvent-1][clEventPlus][binE]
                    if pCommon.anglePointsDegToDeg(0., 0., lEvent, bEvent)<rOnMax and (abs(bEvent)>=bMaskMax or abs(lEvent)<lMaskMin):
                        self.flagOn[0] = kTRUE
                        self.aaaHtgMap[ctEvent-1][clEventPlus][plotE].Fill(-lEvent, bEvent)
                        self.aaHtgNumOn[ctEvent-1][clEventPlus].Fill(eEvent)
                        self.aaHtgEnergyOn[ctEvent-1][clEventPlus].Fill(eEvent)#, (10**(eEvent-3))**2)
                        self.aaHtgLCCountOn[ctEvent-1][clEventPlus].Fill(ConvertMetToMjd(tEvent)-ConvertMetToMjd(self.tStart))
                    elif (lEvent<lOffMax and (bEvent>bOffMin or bEvent<-bOffMin)) or (lEvent>lOffMin and (bEvent>bOffMin or bEvent<-bOffMin)):
                        self.flagOff[0] = kTRUE
                        self.aaaHtgMap[ctEvent-1][clEventPlus][plotE].Fill(-lEvent, bEvent)
                        self.aaHtgNumOff[ctEvent-1][clEventPlus].Fill(eEvent)
                        self.aaHtgEnergyOff[ctEvent-1][clEventPlus].Fill(eEvent)#, (10**(eEvent-3))**2)
                        self.aaHtgLCCountOff[ctEvent-1][clEventPlus].Fill(ConvertMetToMjd(tEvent)-ConvertMetToMjd(self.tStart))
        self.trObj.Fill()

    def calc(self):
        Target.calc(self)
        aHtgData = []
        aHtgData.append(self.aaHtgEnergy[0][1].Clone("hEnergy_CalTkr"))
        aHtgData[-1].SetTitle("Regular")
        aHtgData[-1].SetLineColor(kBlack)
        aHtgData[-1].SetMarkerColor(kBlack)
        aHtgData.append(self.aaHtgEnergy[1][1].Clone("hEnergy_CalOnly"))
        aHtgData[-1].SetTitle("CalOnly")
        aHtgData[-1].SetLineColor(kBlue)
        aHtgData[-1].SetMarkerColor(kBlue)
        aHtgData.append(aHtgData[0].Clone("hEnergy_Combined"))
        aHtgData[-1].Reset()
        aHtgData[-1].SetTitle("Combined")
        aHtgData[-1].SetLineColor(kRed)
        aHtgData[-1].SetMarkerColor(kRed)
        aHtgData[-1].Add(aHtgData[0], 1)
        aHtgData[-1].Add(aHtgData[1], 1)
        aFcExpo = []
        aHtgDevi = []
        for htg in aHtgData:
            print htg.GetName()
            htg.SetLineStyle(1)
            htg.Fit("expo")
            aFcExpo.append(htg.GetFunction("expo"))
            aFcExpo[-1].SetLineWidth(1)
            aHtgDevi.append(htg.Clone("%sDevi" % htg.GetName()))
            aHtgDevi[-1].Reset()
            for ix in range(1, aHtgDevi[-1].GetNbinsX()+1):
                vModel = aFcExpo[-1].Integral(aHtgDevi[-1].GetBinLowEdge(ix), aHtgDevi[-1].GetBinLowEdge(ix+1))/aHtgDevi[-1].GetBinWidth(ix)
                aHtgDevi[-1].SetBinContent(ix, (htg.GetBinContent(ix)-vModel)/vModel)
                aHtgDevi[-1].SetBinError(ix, htg.GetBinError(ix)/vModel)
        aFcExpo[0].SetLineColor(kBlack)
        aFcExpo[1].SetLineColor(kBlue)
        aFcExpo[2].SetLineColor(kRed)
        self.aHtgEnergyBoth = aHtgData
        self.aHtgEnergyDevi = aHtgDevi

    def draw(self):
        Target.draw(self)
        self.hsEnergyBoth = ROOT.THStack("hsEnergyBoth", "Counts in fine energy bin of %s;log_{10}Energy[MeV];[counts]" % self.name)
        for hData in self.aHtgEnergyBoth:
            self.hsEnergyBoth.Add(hData)
        self.hsEnergyDevi = ROOT.THStack("hsEnergyDevi", "Deviation from power low fit of %s;log_{10}Energy[MeV];(Counts-PL)/PL" % self.name)
        for hDevi in self.aHtgEnergyDevi:
            self.hsEnergyDevi.Add(hDevi)
        cEnergyBoth = ROOT.TCanvas("cEnergyBoth%s" % self.name, "Counts in fine energy bin of %s" % self.name, 800, 600)
        cEnergyBoth.Divide(1, 2)
        cEnergyBoth.cd(1)
        self.hsEnergyBoth.Draw('E1 nostack')
        cEnergyBoth.cd(2)
        self.hsEnergyDevi.Draw('E1 nostack')
        self.cEnergyBoth = cEnergyBoth

    def writeObjects(self):
        Target.writeObjects(self)
        # for iE in range(self.energyPlot.nBin):
        #     for iAS in range(len(self.aaClass)):
        #         for iSS in range(len(self.aaClass[iAS])):
        #             self.aaaHtgMap[iAS][iSS][iE].Write("", ROOT.TObject.kOverwrite)
        for aaHtg in self.aaaHtgMap:
            for aHtg in aaHtg:
                for htg in aHtg:
                    htg.Write("", ROOT.TObject.kOverwrite)
        for htg in self.aHtgEnergyBoth:
            htg.Write("", ROOT.TObject.kOverwrite)
        for htg in self.aHtgEnergyDevi:
            htg.Write("", ROOT.TObject.kOverwrite)
        self.cEnergyBoth.Write("", ROOT.TObject.kOverwrite)
        print "Writing finished."
