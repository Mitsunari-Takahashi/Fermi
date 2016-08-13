import ROOT
import numpy
from pAnalysisConfig import *
import pColor
import pCommon
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
from ROOT import kWhite, kBlack, kGray, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE

class Target:
    def __init__(self, strName, config = ClassConfig(), eRegion=EnergyLogRegion(3,4.75,0.25), tStart=252460800.0, tEnd=378691200.0):
        self.name = strName
        self.aaClass = config.aaStrSelect
        self.aClass = config.aStrSelect
        self.energySet = eRegion
        self.saOn = [1.,1.]
        self.saOff = [1.,1.]
        self.zCut = config.zCut
        self.tStart = tStart
        self.tEnd = tEnd
        self.tDuration = tEnd - tStart
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
        for iS in range(len(self.aaClass)):
            aaHtgNumOn.append([])
            aaHtgNumOff.append([])
            aaHtgNumSig.append([])
            aaHtgNumBkg.append([])
            aaHtgSBratio.append([])
            aaHtgSgnf.append([])
            aaHtgEnergy.append([])
            aaHtgLCCount.append([])
            for jS in range(len(self.aaClass[iS])):
                aaHtgNumOn[-1].append(ROOT.TH1D("hNumOn%s_%s_%s" % (self.name,iS,jS), "Number of %s ON events of %s;log_{10}Energy[MeV];[counts]" % (self.aaClass[iS][jS], self.name), self.energySet.nBin, self.energySet.edgeLow, self.energySet.edgeUp))
                aaHtgNumOn[-1][-1].SetFillStyle(0)
                aaHtgNumOn[-1][-1].SetLineWidth(3-iS)
                aaHtgNumOn[-1][-1].SetLineStyle(2-iS)
                aaHtgNumOn[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumOn[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgNumOn[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgNumOn[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumOff[-1].append(ROOT.TH1D("hNumOff%s_%s_%s" % (self.name,iS,jS), "Number of %s OFF events of %s;log_{10}Energy[MeV];[counts]" % (self.aaClass[iS][jS], self.name), self.energySet.nBin, self.energySet.edgeLow, self.energySet.edgeUp))
                aaHtgNumOff[-1][-1].SetFillStyle(0)
                aaHtgNumOff[-1][-1].SetLineWidth(3-iS)
                aaHtgNumOff[-1][-1].SetLineStyle(2-iS)
                aaHtgNumOff[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumOff[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgNumOff[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgNumOff[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumSig[-1].append(ROOT.TH1D("hNumSig%s_%s_%s" % (self.name,iS,jS), "Number of %s Signal events of %s;log_{10}Energy[MeV];[counts]" % (self.aaClass[iS][jS], self.name), self.energySet.nBin, self.energySet.edgeLow, self.energySet.edgeUp))
                aaHtgNumSig[-1][-1].SetFillStyle(0)
                aaHtgNumSig[-1][-1].SetLineWidth(3-iS)
                aaHtgNumSig[-1][-1].SetLineStyle(2-iS)
                aaHtgNumSig[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumSig[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgNumSig[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgNumSig[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumBkg[-1].append(ROOT.TH1D("hNumBkg%s_%s_%s" % (self.name,iS,jS), "Number of %s Background events of %s;log_{10}Energy[MeV];[counts]" % (self.aaClass[iS][jS], self.name), self.energySet.nBin, self.energySet.edgeLow, self.energySet.edgeUp))
                aaHtgNumBkg[-1][-1].SetFillStyle(0)
                aaHtgNumBkg[-1][-1].SetLineWidth(3-iS)
                aaHtgNumBkg[-1][-1].SetLineStyle(2-iS)
                aaHtgNumBkg[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgNumBkg[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgNumBkg[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgNumBkg[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgSBratio[-1].append(ROOT.TH1D("hSBratio%s_%s_%s" % (self.name,iS,jS), "%s Signal/Background ratio of %s;log_{10}Energy[MeV];Ratio" % (self.aaClass[iS][jS], self.name), self.energySet.nBin, self.energySet.edgeLow, self.energySet.edgeUp))
                aaHtgSBratio[-1][-1].SetFillStyle(0)
                aaHtgSBratio[-1][-1].SetLineWidth(3-iS)
                aaHtgSBratio[-1][-1].SetLineStyle(2-iS)
                aaHtgSBratio[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgSBratio[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgSBratio[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgSBratio[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgSgnf[-1].append(ROOT.TH1D("hSignificance%s_%s_%s" % (self.name,iS,jS), "%s Significance of %s;log_{10}Energy[MeV];[#sigma]" % (self.aaClass[iS][jS], self.name), self.energySet.nBin, self.energySet.edgeLow, self.energySet.edgeUp))
                aaHtgSgnf[-1][-1].SetFillStyle(0)
                aaHtgSgnf[-1][-1].SetLineWidth(3-iS)
                aaHtgSgnf[-1][-1].SetLineStyle(2-iS)
                aaHtgSgnf[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgSgnf[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgSgnf[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgSgnf[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgEnergy[-1].append(ROOT.TH1D("hEnergy%s_%s_%s" % (self.name,iS,jS), "%s Energy plot of %s;log_{10}Energy[MeV];NE^{2}" % (self.aaClass[iS][jS], self.name), self.energySet.nBin*5, self.energySet.edgeLow, self.energySet.edgeUp))
                aaHtgEnergy[-1][-1].SetFillStyle(0)
                aaHtgEnergy[-1][-1].SetLineWidth(3-iS)
                aaHtgEnergy[-1][-1].SetLineStyle(2-iS)
                aaHtgEnergy[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgEnergy[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgEnergy[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgEnergy[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
                # Light curve in counts
                if self.tDuration >= 12 * pCommon.JulYrInSec: # Year bin
                    hLCCountBinWidth = pCommon.JulYrInSec
                elif self.tDuration >= 12 * pCommon.JulYrInSec/12.: # Month bin
                    hLCCountBinWidth = pCommon.JulYrInSec
                elif self.tDuration>=12 * 24*60*60: # Day bin
                    hLCCountBinWidth = 24.*60.*60.
                else: # Hour bin
                    hLCCountBinWidth = 60.*60.
                hLCCountNumBin = self.tDuration / hLCCountBinWidth                    
                aaHtgLCCount[-1].append(ROOT.TH1D("hLCCount%s_%s_%s" % (self.name,iS,jS), "%s Light curve (in counts) of %s;After MET%s [sec];[counts]" % (self.aaClass[iS][jS], self.name, self.tStart), hLCCountNumBin, 0, self.tDuration))
                aaHtgLCCount[-1][-1].SetFillStyle(3002-2001*iS)
                aaHtgLCCount[-1][-1].SetLineWidth(3-iS)
                aaHtgLCCount[-1][-1].SetLineStyle(2-iS)
                aaHtgLCCount[-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgLCCount[-1][-1].SetFillColor(pColor.akColor(jS+(1-iS)*jS))
                aaHtgLCCount[-1][-1].SetMarkerStyle(pColor.aakMarkerStyle(iS,jS))
                aaHtgLCCount[-1][-1].SetMarkerSize(pColor.aakMarkerSize(iS,jS))
                aaHtgLCCount[-1][-1].SetMarkerColor(pColor.akColor(jS+(1-iS)*jS))
        self.aaHtgNumOn = aaHtgNumOn
        self.aaHtgNumOff = aaHtgNumOff
        self.aaHtgNumSig = aaHtgNumSig
        self.aaHtgNumBkg = aaHtgNumBkg
        self.aaHtgSBratio = aaHtgSBratio
        self.aaHtgSgnf = aaHtgSgnf
        self.aaHtgEnergy = aaHtgEnergy
        self.aaHtgLCCount = aaHtgLCCount

    def calc(self):
        for iS in range(len(self.aaClass)):
            for jS in range(len(self.aaClass[iS])):
                self.aaHtgNumSig[iS][jS].Reset()
                self.aaHtgNumBkg[iS][jS].Reset()
                self.aaHtgSBratio[iS][jS].Reset()
                self.aaHtgSgnf[iS][jS].Reset()
                #self.aaHtgEnergy[iS][jS].Reset()

        for iS in range(len(self.aaClass)):
            print "==========", self.aaClass[iS], "=========="
            for jS in range(len(self.aaClass[iS])):
                print "----------", self.aaClass[iS][jS], "----------"
                for iE in range(self.energySet.nBin):
                    print self.energySet.getBin(iE)
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
                    #sgnf = math.sqrt( 2 * (nOn*math.log((saOff/saOn+1)*nOn/(nOn+nOff)) + nOff*math.log((1+saOn/saOff)*nOff/(nOn+nOff)) ) )
                    #print "Significance:", sgnf
                    #self.aaHtgSgnf[iS][jS].SetBinContent(iE+1, sgnf)
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

        cNumSig = ROOT.TCanvas("cNumSig_%s" % self.name, "Number of Signal events of %s" % self.name, 800, 600)
        hsNumSig = ROOT.THStack("hsNumSig_%s" % self.name, "Number of Signal events of %s;log_{10}Energy[MeV];[counts]" % self.name)
        for aHtg in self.aaHtgNumSig:
            for htg in aHtg:
                hsNumSig.Add(htg)
        cNumSig.cd()
        hsNumSig.Draw("nostack E1")
        legHtg.Draw('same')
        self.hsNumSig = hsNumSig
        self.cNumSig = cNumSig

        cNumBkg = ROOT.TCanvas("cNumBkg_%s" % self.name, "Number of Background events of %s" % self.name, 800, 600)
        hsNumBkg = ROOT.THStack("hsNumBkg_%s" % self.name, "Number of Background events of %s;log_{10}Energy[MeV];[counts]" % self.name)
        for aHtg in self.aaHtgNumBkg:
            for htg in aHtg:
                hsNumBkg.Add(htg)
        cNumBkg.cd()
        hsNumBkg.Draw("nostack E1")
        legHtg.Draw('same')
        self.hsNumBkg = hsNumBkg
        self.cNumBkg = cNumBkg

        cSBratio = ROOT.TCanvas("cSBratio_%s" % self.name, "Signal/Background ratio of %s" % self.name, 800, 600)
        hsSBratio = ROOT.THStack("hsSBratio_%s" % self.name, "Signal/Background ratio of %s;log_{10}Energy[MeV];Ratio" % self.name)
        for aHtg in self.aaHtgSBratio:
            for htg in aHtg:
                hsSBratio.Add(htg)
        cSBratio.cd()
        hsSBratio.Draw("nostack E1")
        legHtg.Draw('same')
        self.hsSBratio = hsSBratio
        self.cSBratio = cSBratio

        # cSgnf = ROOT.TCanvas("cSgnf", "Significance of %s" % self.name, 800, 600)
        # hSgnfDummy = ROOT.TH1D("hSgnfDummy", "Significance of %s;log_{10}Energy[MeV];[#sigma]" % self.name, self.energySet.nBin, self.energySet.edgeLow, self.energySet.edgeUp)
        # hSgnfDummy.GetYaxis().SetRangeUser(0., self.aaHtgSgnf[0][1].GetMaximum()*1.1)
        # cSgnf.cd()
        # hSgnfDummy.Draw()
        # for aHtg in self.aaHtgSgnf:
        #     for htg in aHtg:
        #         htg.Draw("sames")
        # self.cSgnf = cSgnf

        cEnergy = ROOT.TCanvas("cEnergy_%s" % self.name, "Energy plot of %s" % self.name, 800, 600)
        hsEnergy = ROOT.THStack("hsEnergy_%s" % self.name, "Energy plot of %s;log_{10}Energy[MeV];NE^{2}" % self.name)
        for aHtg in self.aaHtgEnergy:
            for htg in aHtg:
                hsEnergy.Add(htg)
        cEnergy.cd()
        hsEnergy.Draw("nostack E1")
        legHtg.Draw('same')
        self.hsEnergy = hsEnergy
        self.cEnergy = cEnergy

        cLCCount = ROOT.TCanvas("cLCCount_%s" % self.name, "Light curve in counts of %s" % self.name, 800, 600)
        hsLCCount = ROOT.THStack("hsLCCount_%s" % self.name, "Light curve in counts of %s;%s;[counts]" % (self.name, self.aaHtgLCCount[0][0].GetXaxis().GetTitle()))
        for aHtg in self.aaHtgLCCount:
            for htg in aHtg:
                hsLCCount.Add(htg)
        cLCCount.cd()
        hsLCCount.Draw()
        legHtgFill = ROOT.TLegend(0.7, 0.85, 0.95, 0.55, "Event class", 'NDC')
        for iaHtg in range(len(self.aaHtgLCCount)):
            for iHtg in range(len(self.aaHtgLCCount[iaHtg])):
                legHtgFill.AddEntry(self.aaHtgLCCount[iaHtg][iHtg], self.aaClass[iaHtg][iHtg], 'lpf')
        legHtgFill.Draw('same')
        self.hsLCCount = hsLCCount
        self.cLCCount = cLCCount
        self.legHtgFill = legHtgFill

    def writeObjects(self):
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
        # for aHtg in self.aaHtgSgnf:
        #     for htg in aHtg:
        #         htg.Write("", ROOT.TObject.kOverwrite)                    
        for aHtg in self.aaHtgEnergy:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)
        for aHtg in self.aaHtgLCCount:
            for htg in aHtg:
                htg.Write("", ROOT.TObject.kOverwrite)

        self.cNumOn.Write("", ROOT.TObject.kOverwrite)
        self.cNumOff.Write("", ROOT.TObject.kOverwrite)
        self.cNumSig.Write("", ROOT.TObject.kOverwrite)
        self.cNumBkg.Write("", ROOT.TObject.kOverwrite)
        self.cSBratio.Write("", ROOT.TObject.kOverwrite)
       # self.cSgnf.Write("", ROOT.TObject.kOverwrite)
        self.cEnergy.Write("", ROOT.TObject.kOverwrite)
        self.cLCCount.Write("", ROOT.TObject.kOverwrite)
        
class PointSource(Target):
    def __init__(self, strName, raTgt, decTgt, zRedshift=0, rAppa=12., rOffMax=[1.7, 12.], rOffMin=[1., 8.], rOnMax=[0.45, 4.0], config = ClassConfig(), eRegion=EnergyLogRegion(3,4.75,0.25), perf=["/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R100_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R30_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R10_perf.root"], tStart=252460800.0, tEnd=378691200.0):
        Target.__init__(self, strName, config, eRegion, tStart, tEnd)
        #self.name = strName
        self.perf = perf
#        self.lCntr = lTgt
#        self.bCntr = bTgt
        self.raCntr = raTgt
        self.decCntr = decTgt
        self.zRedshift = zRedshift
        if rAppa>=90:
            print "Too large radius."
        self.radius = rAppa
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
                self.rOffMax[ict].append([])
                self.rOffMin[ict].append([])
                self.rOnMax[ict].append([])
                self.saOn[ict].append([])
                self.saOff[ict].append([])
                for ie in range(self.energySet.nBin):
                    if ict==0:
                        self.rOffMax[ict][icl].append(self.radius)
                        self.rOffMin[ict][icl].append(0.5)
                        self.rOnMax[ict][icl].append(0.1)
                    elif ict==1:
                        self.rOffMax[ict][icl].append(self.radius)
                        self.rOffMin[ict][icl].append(self.perf.getPSF95(icl, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.rOnMax[ict][icl].append(self.perf.getPSF68(icl, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        if not (self.rOffMax[ict][icl][-1]<=self.radius and self.rOffMin[ict][icl][-1]<self.rOffMax[ict][icl][-1] and self.rOnMax[ict][icl][-1]<=self.rOffMin[ict][icl][-1] and self.rOnMax[ict][icl][-1]>0):
                            print "Bad region setup!!"
                            sys.exit(-1)
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
        aaaHtgFolded = []
        for iS in range(len(self.aaClass)):
            aaaGrpMap.append([])
            aaGreHighEnergy.append([])
            aaaHtgMap.append([])
            aaaHtgTheta.append([])
#            aaaHtgFolded.append([])
            for jS in range(len(self.aaClass[iS])):
                aaaGrpMap[-1].append([])
                aaaHtgMap[-1].append([])
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
                for iE in range(self.energySet.nBin):
                    aaaGrpMap[-1][-1].append(ROOT.TGraphPolar())
                    aaaGrpMap[-1][-1][-1].SetName("grpMap%s_%s_%s_%s" % (self.name, iS, jS, iE))
                    #aaaGrpMap[-1][-1][-1].SetTitle("{0} map around {1} in {2:.1f} - {3:.1f} GeV".format(self.aaClass[iS][jS], self.name, 10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)))
                    aaaGrpMap[-1][-1][-1].SetMaxRadial(self.radius)
                    aaaGrpMap[-1][-1][-1].SetMarkerColor(pColor.akColor(iE))
                    aaaGrpMap[-1][-1][-1].SetMarkerStyle(7)
                    aaaHtgMap[-1][-1].append(ROOT.TH2D("hMap%s_%s_%s_%s" % (self.name, iS, jS, iE), "{0} spacial distribution of {1} in {2:.1f} - {3:.1f} GeV;;;[counts/sr]".format(self.aaClass[iS][jS], self.name, 10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)), 12, 0., 2.*math.pi, int(rAppa), 0, int(rAppa)))
                    aaaHtgMap[-1][-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                    aaaHtgMap[-1][-1][-1].SetFillColor(pColor.akColor(jS+(1-iS)*jS))
                    aaaHtgMap[-1][-1][-1].SetLineWidth(3-iS)
                    aaaHtgMap[-1][-1][-1].SetLineStyle(2-iS)
                    aaaHtgTheta[-1][-1].append(ROOT.TH1D("hTheta%s_%s_%s_%s" % (self.name, iS, jS, iE), "{0} #theta^2 plot of {1} in {2:.1f} - {3:.1f} GeV;#theta^2 [#circ];[counts]".format(self.aaClass[iS][jS], self.name, 10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)), int(rAppa**2)*10, 0, int(rAppa**2)))
                    aaaHtgTheta[-1][-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                    aaaHtgTheta[-1][-1][-1].SetLineWidth(3-iS)
                    aaaHtgTheta[-1][-1][-1].SetLineStyle(2-iS)
#                    aaaHtgFolded[-1][-1].append(ROOT.TH2D("hFolded%s_%s_%s_%s" % (self.name, iS, jS, iE), "{0} PSF folded map around {1} in {2:.1f} - {3:.1f} GeV;RA[#circ];DEC[#circ];[counts]".format(self.aaClass[iS][jS], self.name, 10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)), 12,))
        self.aaaGrpMap = aaaGrpMap
        self.aaGreHighEnergy = aaGreHighEnergy
        self.aaaHtgMap = aaaHtgMap
        self.aaaHtgTheta = aaaHtgTheta

    def fill(self, lEvent, bEvent, eEvent, ctEvent, clEvent, zEvent, tEvent, cthEvent):
        #Target.fill(self, lEvent, bEvent, eEvent, ctEvent, clEvent)
        clEvent = int(clEvent)
        binE = self.energySet.findBin(eEvent)
        lRad = math.radians(lEvent)
        bRad = math.radians(bEvent)
        vecTgt = numpy.array([cos(radians(self.decCntr))*cos(radians(self.raCntr)), cos(radians(self.decCntr))*sin(radians(self.raCntr)), sin(radians(self.decCntr))])
        vecEve = numpy.array([cos(bRad)*cos(lRad), cos(bRad)*sin(lRad), sin(bRad)])
        radTheta = acos(numpy.dot(vecTgt, vecEve))

        if tEvent>=self.tStart and tEvent<= self.tEnd and (degrees(radTheta)<=self.radius) and binE>=0 and binE<self.energySet.nBin and zEvent<self.zCut[ctEvent-1]:
            vecNorth = numpy.array([cos(radians(self.decCntr+90))*cos(radians(self.raCntr)), cos(radians(self.decCntr+90))*sin(radians(self.raCntr)), sin(radians(self.decCntr+90))])
            if cos(radTheta)!=0:
                vecProj = vecEve/cos(radTheta) - vecTgt #[vecEve[0]-vecTgt[0]*cos(radTheta[0]), vecEve[1]-vecTgt[1]*cos(radTheta[1]), vecEve[2]-vecTgt[2]*cos(radTheta[2])]
            else:
                print "Zero division!!"
            vecCross = numpy.cross(vecProj, vecNorth)
            radPhi = acos(numpy.dot(vecProj, vecNorth) / (numpy.linalg.norm(vecProj)*numpy.linalg.norm(vecNorth)))
            if numpy.dot(vecCross, vecNorth)<0:
                radPhi = -radPhi +  2.*math.pi
             #   radPhi = 2.*math.pi+radPhi
            #print "Phi angle :", radPhi, "rad"
            #print ""
            #print self.aaGreHighEnergy
            if math.degrees(radTheta) < self.rOnMax[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]:
                #print "(", ctEvent-1, clEvent-1, ")"
                self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].SetPoint(self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].GetN(), tEvent-self.tStart, eEvent)
                if ctEvent==2:
                    #self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].SetPointError(self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].GetN()-1, 0, self.perf.getEdisp68_cth(clEvent-1, eEvent, cthEvent))
                    self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].SetPointEYhigh(self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].GetN()-1, abs(math.log10(1+self.perf.getEdisp68_cth(clEvent-1, eEvent, cthEvent))))
                    self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].SetPointEYlow(self.aaGreHighEnergy[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].GetN()-1, abs(math.log10(1-self.perf.getEdisp68_cth(clEvent-1, eEvent, cthEvent))))
                self.aaHtgLCCount[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1].Fill(tEvent-self.tStart)
            for clEventPlus in range(clEvent-int(ctEvent==1 and clEvent==3)):
                if binE<self.energySet.nBin:
                    self.aaaHtgTheta[ctEvent-1][clEventPlus][binE].Fill(math.degrees(radTheta)**2)
                #if binE>=0 and binE<self.energySet.nBin:
                    self.aaaGrpMap[ctEvent-1][clEventPlus][binE].SetPoint(self.aaaGrpMap[ctEvent-1][clEventPlus][binE].GetN(), radPhi, math.degrees(radTheta))
                    self.aaaHtgMap[ctEvent-1][clEventPlus][binE].Fill(radPhi, math.degrees(radTheta), 1./(sin(radTheta)/self.aaaHtgMap[ctEvent-1][clEventPlus][binE].GetNbinsX()*2.*math.pi/self.aaaHtgMap[ctEvent-1][clEventPlus][binE].GetNbinsY()) )
                if math.degrees(radTheta) < self.rOnMax[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]:
                #for clEventPlus in range(clEvent-int(ctEvent==1 and clEvent==3)):
                    #self.aaGreHighEnergy[ctEvent-1][clEventPlus].SetPoint(self.aaGreHighEnergy[ctEvent-1][clEventPlus].GetN(), tEvent, eEvent)
                    self.aaHtgNumOn[ctEvent-1][clEventPlus].Fill(eEvent)
                    self.aaHtgEnergy[ctEvent-1][clEventPlus].Fill(eEvent, (10**(eEvent-3))**2)
                elif math.degrees(radTheta) < self.rOffMax[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE] and math.degrees(radTheta) >= self.rOffMin[ctEvent-1][clEvent-int(ctEvent==1 and clEvent==3)-1][binE]:
                #for clEventPlus in range(clEvent-int(ctEvent==1 and clEvent==3)):
                    self.aaHtgNumOff[ctEvent-1][clEventPlus].Fill(eEvent)

        # for aaHtgMap in self.aaaHtgMap:
        #     for aHtgMap in aaHtgMap:
        #         for hMap in aHtgMap:
        #             for iphi in range(hMap.GetNbinsX()):
        #                 for itheta in range(hMap.GetNbinsY()):
        #                     hMap.SetBinContent(iphi+1, itheta+1, hMap.GetBinContent(iphi+1, itheta+1)/cos(radians(hMap.GetYaxis().GetBinCenter(itheta+1))))

    def draw(self):
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
                for binE in range(self.energySet.nBin):
                    self.fRadiusOnMax[ctEvent][clEvent].append(ROOT.TF2("fRadiusOnMax%s_%s_%s" % (ctEvent, clEvent, binE), "x**2+y**2 - %s**2" % self.rOnMax[ctEvent][clEvent][binE], -self.radius, self.radius, -self.radius, self.radius))
                    self.fRadiusOnMax[-1][-1][-1].SetMinimum(0)
                    self.fRadiusOnMax[-1][-1][-1].SetMaximum(0)
                    self.fRadiusOnMax[-1][-1][-1].SetLineWidth(1)
                    self.fRadiusOnMax[-1][-1][-1].SetLineColor(kGray)
                    self.fRadiusOffMin[ctEvent][clEvent].append(ROOT.TF2("fRadiusOffMin%s_%s_%s" % (ctEvent, clEvent, binE), "x**2+y**2 - %s**2" % self.rOffMin[ctEvent][clEvent][binE], -self.radius, self.radius, -self.radius, self.radius))
                    self.fRadiusOffMin[-1][-1][-1].SetMinimum(0)
                    self.fRadiusOffMin[-1][-1][-1].SetMaximum(0)
                    self.fRadiusOffMin[-1][-1][-1].SetLineWidth(1)
                    self.fRadiusOffMin[-1][-1][-1].SetLineColor(kGray)
                    self.fRadiusOffMax[ctEvent][clEvent].append(ROOT.TF2("fRadiusOffMax%s_%s_%s" % (ctEvent, clEvent, binE), "x**2+y**2 - %s**2" % self.rOffMax[ctEvent][clEvent][binE], -self.radius, self.radius, -self.radius, self.radius))
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
                aaHtgMapDummy[-1].append(ROOT.TH2D("hMapDummy%s_%s_%s" % (self.name, iaaG, iaG), "{0} map around {1};[#circ];[#circ]".format(self.aaClass[iaaG][iaG], self.name), 2*self.radius, -self.radius, self.radius, 2*self.radius, -self.radius, self.radius))
                for iG in range(len(self.aaaGrpMap[iaaG][iaG])):
                    aaMgrMap[-1][-1].Add(self.aaaGrpMap[iaaG][iaG][iG])
                    if iaaG==0 and iaG==0:
                        legGrpMap.AddEntry(self.aaaGrpMap[iaaG][iaG][iG], "{0:.1f} - {1:.1f} GeV".format(10**(self.energySet.getBin(iG)[0]-3), 10**(self.energySet.getBin(iG)[1]-3)), 'p')
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
        mgrHighEnergy.GetXaxis().SetTitle("Time after MET%s [sec]" % self.tStart)
        mgrHighEnergy.GetYaxis().SetTitle("log_{10}(Energy[MeV])")
        legGreHighEnergy.Draw("same")

        cTheta = ROOT.TCanvas("cTheta%s" % self.name, "%s #theta^{2} plots" % self.name, 1200, 600)
        mDown = int(math.sqrt(self.energySet.nBin+1))
        mAcross = math.ceil((self.energySet.nBin+1) / mDown)
        cTheta.Divide(mAcross, mDown)
        aHstaTheta = []
        legTheta = ROOT.TLegend(0.1,0.1,0.9,0.9,"Event class",'NDC')
        for tE in range(self.energySet.nBin):
            aHstaTheta.append(ROOT.THStack("hsTheta%s_%s" % (self.name, tE), "#theta^2 plot of {0} in {1:.1f} - {2:.1f} GeV;#theta[#circ]^2;[counts]".format(self.name, 10**(self.energySet.getBin(tE)[0]-3), 10**(self.energySet.getBin(tE)[1]-3))))
            for iaaH in range(len(self.aaaHtgTheta)):
                for iaH in range(len(self.aaaHtgTheta[iaaH])):
                    self.aaaHtgTheta[iaaH][iaH][tE].Rebin(30)
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
        lDown = len(self.aaClass)
        lAcross = 0
        for aS in self.aaClass:    
            lAcross = max(lAcross, len(aS))
        for iE in range(self.energySet.nBin):
            aCanSpacial.append(ROOT.TCanvas("cSpacial%s_%s" % (self.name,iE), "Spacial distribution of {0} in {1:.1f} - {2:.1f} GeV".format(self.name,10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)), 1200, 800))
            aCanSpacial[-1].Divide(lAcross, lDown)
            aaaHtgDummySpacial.append([])
            for iAS in range(len(self.aaClass)):
                aaaHtgDummySpacial[-1].append([])
                for iSS in range(len(self.aaClass[iAS])):
                    aaaHtgDummySpacial[-1][-1].append(ROOT.TH2D("hDummySpacial%s_%s_%s_%s" % (self.name, iE, iAS, iSS), "{0} Spacial distribution of {1} in {2:.1f} - {3:.1f} GeV;[#circ];[#circ]".format(self.aaClass[iAS][iSS], self.name, 10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)), 2*self.radius, -self.radius, self.radius, 2*self.radius, -self.radius, self.radius))
                    aCanSpacial[-1].cd(1+iAS*lAcross+iSS)
                    aCanSpacial[-1].cd(1+iAS*lAcross+iSS).SetGridx(0)
                    aCanSpacial[-1].cd(1+iAS*lAcross+iSS).SetGridy(0)
                    aaaHtgDummySpacial[-1][-1][-1].SetStats(kFALSE)
                    aaaHtgDummySpacial[-1][-1][-1].Draw()
                    #self.aaaHtgMap[iAS][iSS][iE].SetMinimum(0)
                    #self.aaaHtgMap[iAS][iSS][iE].SetMaximum(self.aaaHtgMap[0][0][iE].GetMaximum())
                    self.aaaHtgMap[iAS][iSS][iE].SetStats(kFALSE)
                    self.aaaHtgMap[iAS][iSS][iE].Draw('POL COLZ SAMES')
                    self.fRadiusOnMax[iAS][iSS][iE].Draw('same')
                    self.fRadiusOffMin[iAS][iSS][iE].Draw('same')
                    self.fRadiusOffMax[iAS][iSS][iE].Draw('same')
                    gPad.Update()
                    palette = self.aaaHtgMap[iAS][iSS][iE].GetListOfFunctions().FindObject("palette")
                    try:
                        palette.SetX1NDC(0.89)
                        palette.SetX2NDC(0.93)
                    except Exception:
                        print "Setting of", palette, "failed."
        self.aCanSpacial = aCanSpacial
        self.aaaHtgDummySpacial = aaaHtgDummySpacial

    def writeObjects(self):
        Target.writeObjects(self)
        self.cMap.Write("", ROOT.TObject.kOverwrite)
        self.cHighEnergy.Write("", ROOT.TObject.kOverwrite)
        self.cTheta.Write("", ROOT.TObject.kOverwrite)
        for aaHtg in self.aaaHtgMap:
            for aHtg in aaHtg:
                for htg in aHtg:
                    #print htg.GetName()
                    htg.Write("", ROOT.TObject.kOverwrite)
        for aaHtg in self.aaaHtgTheta:
            for aHtg in aaHtg:
                for htg in aHtg:
                    #print htg.GetName()
                    htg.Write("", ROOT.TObject.kOverwrite)
        for aaGrp in self.aaaGrpMap:
            for aGrp in aaGrp:
                for grp in aGrp:
                    #print grp.GetName()
                    grp.Write("", ROOT.TObject.kOverwrite)

        for aGr in self.aaGreHighEnergy:
            for gr in aGr:
                gr.Write("", ROOT.TObject.kOverwrite)

        for iE in range(self.energySet.nBin):
            self.aCanSpacial[iE].Write("", ROOT.TObject.kOverwrite)

class EarthLimb(Target):
    def __init__(self, strName, zOff1Min=[1.7, 12.], zOff1Max=[1., 8.],zOff2Min=[1.7, 12.], zOff2Max=[1., 8.], zOnMin=[], zOnMax=[0.45, 4.0], config = ClassConfig(), eRegion=EnergyLogRegion(3,4.75,0.25), perf=["/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R100_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R30_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R10_perf.root"]):
        Target.__init__(self, strName, config, eRegion)
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
            self.saOn.append([])
            for hSS in range(len(self.aaClass[hS])):
                self.zOnMin[hS].append([])
                self.zOnMax[hS].append([])
                self.zOff1Min[hS].append([])
                self.zOff1Max[hS].append([])
                self.zOff2Min[hS].append([])
                self.zOff2Max[hS].append([])
                self.saOn[hS].append([])
                self.saOn[hSS].append([])
                for ie in range(self.energySet.nBin):
                    if hS==0:
                        self.zOnMin[hS][hSS].append(111.10)
                        self.zOnMax[hS][hSS].append(112.95)
                        self.zOff1Min[hS][hSS].append(108.66)
                        self.zOff1Max[hS][hSS].append(109.57)
                        self.zOff2Min[hS][hSS].append(114.52)
                        self.zOff2Max[hS][hSS].append(115.47)
                    elif hS==1:
                        self.zOnMin[hS][hSS].append(111.10+0.1-self.perf.getPSF68(hSS, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.zOnMax[hS][hSS].append(112.95-0.1+self.perf.getPSF68(hSS, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.zOff1Min[hS][hSS].append(105)
                        self.zOff1Max[hS][hSS].append(111.10-self.perf.getPSF95(hSS, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.zOff2Min[hS][hSS].append(112.95+self.perf.getPSF95(hSS, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.zOff2Max[hS][hSS].append(120)
                    self.saOn[hS][hSS].append( 2.0 * math.pi * ( cos(radians(self.zOnMin[hs][hSS][-1])) - cos(radians(self.zOnMax[hS][hSS][-1])) ) )
                    self.saOff[hS][hSS].append( 2.0 * math.pi * ( cos(radians(self.zOffMin1[hS][hSS][-1])) - cos(radians(self.zOffMax1[hS][hSS][-1])) ) + 2.0 * math.pi * ( cos(radians(self.zOffMin2[hS][hSS][-1])) - cos(radians(self.zOffMax2[hS][hSS][-1])) ))
        print "Solid angle of ON region:", self.saOn
        print "Solid angle of OFF region:", self.saOff
        aaHtgZenithTheta = []

class GalacticRidge(Target):
    def __init__(self, strName, bOffMin=[50., 50.], lOffMin=[90., 90.], lOffMax=[-90., -90.], bOnMax=[1.5, 3.0], lOnMin=[-50., -51.5], lOnMax=[40., 41.5], config = ClassConfig(), eRegion=EnergyLogRegion(3,4.75,0.25), perf=["/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R100_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R30_perf.root", "/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r9p9_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_WP8CalOnlyLogEnergy_CalOnly_R10_perf.root"], tStart=252460800.0, tEnd=378691200.0):
        Target.__init__(self, strName, config, eRegion, tStart, tEnd)
        self.perf = perf
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
                        self.bOnMax[hS][hSS].append(self.bOnMax[0][0][ie]-0.1+self.perf.getPSF68(hSS, (self.energySet.aBin[ie]+self.energySet.wBin/2.0)) )
                        self.lOnMin[hS][hSS].append(self.lOnMin[0][0][ie]+0.1-self.perf.getPSF68(hSS, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
                        self.lOnMax[hS][hSS].append(self.lOnMax[0][0][ie]-0.1+self.perf.getPSF68(hSS, self.energySet.aBin[ie]+self.energySet.wBin/2.0))
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
                for iE in range(self.energySet.nBin):
                    aaaHtgLatitude[-1][-1].append(ROOT.TH1D("hLatitude%s_%s_%s" % (iS,jS,iE), "{0} Galactic latitude distribution in {1:.1f} - {2:.1f} GeV;sin(b);counts [counts]".format(self.aaClass[iS][jS], 10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)), 200, -1, 1))
                    aaaHtgLatitude[-1][-1][-1]
                    aaaHtgLatitude[-1][-1][-1].SetLineColor(pColor.akColor(jS+(1-iS)*jS))
                    aaaHtgLatitude[-1][-1][-1].SetLineWidth(3-iS)
                    aaaHtgLatitude[-1][-1][-1].SetLineStyle(2-iS)
                        
                    aaaHtgMapAllSky[-1][-1].append(ROOT.TH2D("hMapAllSky%s_%s_%s" % (iS,jS,iE), "{0} all sky map in {1:.1f} - {2:.1f} GeV;- Galactic longitude [#circ];Galactic latitude [#circ]".format(self.aaClass[iS][jS], 10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)), 360, -180, 180, 180, -90, 90))
                    aaaHtgFoldedMapAllSky[-1][-1].append(ROOT.TH2D("hFoldedMapAllSky%s_%s_%s" % (iS,jS,iE), "{0} folded all sky map in {1:.1f} - {2:.1f} GeV;- Galactic longitude [#circ];Galactic latitude [#circ]".format(self.aaClass[iS][jS], 10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)), 360, -180, 180, 140, -70, 70))
                    #aaaHtgFoldedMapAllSky[-1][-1].append(aaaHtgMapAllSky[-1][-1][-1].Clone("hFoldedMapAllSky%s_%s_%s" % (iS, jS, iE)))
                    aaaHtgMapAllCel[-1][-1].append(ROOT.TH2D("hMapAllCel%s_%s_%s" % (iS,jS,iE), "{0} all celestial map in {1:.1f} - {2:.1f} GeV;- RA [#circ];DEC [#circ]".format(self.aaClass[iS][jS], 10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)), 360, -180, 180, 140, -70, 70))
                    aaaHtgFoldedMapAllCel[-1][-1].append(aaaHtgMapAllSCel[-1][-1][-1].Clone("hFoldedMapAllCel%s_%s_%s" % (iS, jS, iE)))
                        #ROOT.TH2D("hFolded%s_%s_%s_%s" % (self.name, iS, jS, iE), "{0} PSF folded map around {1} in {2:.1f} - {3:.1f} GeV;RA[#circ];DEC[#circ];[counts]".format(self.aaClass[iS][jS], self.name, 10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)), 12,))
        self.aaaHtgLatitude = aaaHtgLatitude
        self.aaaHtgMapAllSky = aaaHtgMapAllSky
        self.aaaHtgFoldedMapAllSky = aaaHtgFoldedMapAllSky
        self.aaaHtgMapAllCel = aaaHtgMapAllCel
        self.aaaHtgFoldedMapAllCel = aaaHtgFoldedMapAllCel

        aaGrPerformance = []
        for tE in range(self.energySet.nBin):
            aaGrPerformance.append([])
            for kAS in range(len(self.aClass)):
                aaGrPerformance[-1].append(ROOT.TGraphErrors(len(self.aaClass[kAS])))
                aaGrPerformance[-1][-1].SetName("grPerformance_%s_%s_%s" % (self.name, tE, kAS))
                aaGrPerformance[-1][-1].SetTitle("{0} performance plot of {1} in {2:.1f} - {3:.1f} GeV".format(self.aClass[kAS], self.name, 10**(self.energySet.getBin(tE)[0]-3), 10**(self.energySet.getBin(tE)[1]-3)))
                aaGrPerformance[-1][-1].GetXaxis().SetTitle("Number of signal events")
                aaGrPerformance[-1][-1].GetYaxis().SetTitle("Number of residual OFF events")
                aaGrPerformance[-1][-1].SetMarkerStyle(25-kAS*5)
                aaGrPerformance[-1][-1].SetLineStyle(2-kAS)
                aaGrPerformance[-1][-1].SetLineWidth(2)
        self.aaGrPerformance = aaGrPerformance
        #print self.aaGrPerformance
        print "Construction finished."

    def fill(self, lEvent, bEvent, raEvent, decEvent, eEvent, ctEvent, clEvent, zEvent, tEvent):
        clEvent = int(clEvent)
        binE = self.energySet.findBin(eEvent)
        if binE>=0 and binE<self.energySet.nBin and tEvent>=self.tStart and tEvent<= self.tEnd:
            if lEvent>180:
                lEvent = -360+lEvent
            if raEvent>180:
                raEvent = -360+raEvent
            #lRad = math.radians(lEvent)
            #bRad = math.radians(bEvent)
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
                    #if binE<self.energySet.nBin:
                    self.aaaHtgMapAllSky[ctEvent-1][clEventPlus][binE].Fill(-lEvent, bEvent)
                    self.aaaHtgMapAllCel[ctEvent-1][clEventPlus][binE].Fill(-raEvent, decEvent)
                    self.aaaHtgLatitude[ctEvent-1][clEventPlus][binE].Fill(sin(radians(bEvent)))
                    if lEvent>lOnMin and lEvent<lOnMax and bEvent>-bOnMax and bEvent<bOnMax:
                        self.aaHtgNumOn[ctEvent-1][clEventPlus].Fill(eEvent)
                        self.aaHtgEnergy[ctEvent-1][clEventPlus].Fill(eEvent, (10**(eEvent-3))**2)
                        self.aaHtgLCCount[ctEvent-1][clEventPlus].Fill(tEvent-self.tStart)
                    elif (lEvent<lOffMax and (bEvent>bOffMin or bEvent<-bOffMin)) or (lEvent>lOffMin and (bEvent>bOffMin or bEvent<-bOffMin)):
                        self.aaHtgNumOff[ctEvent-1][clEventPlus].Fill(eEvent)
#    def __init__(self, strName, bOffMin=[50., 50.], lOffMin=[90., 90.], lOffMax=[-90., -90.], bOnMax=[1.5, 5.5], lOnMin=[-50., -54], lOnMax=[40., 44.], config = ClassConfig(), self.energySet=EnergyLogRegion(3,4.75,0.25)):

    def calc(self):
        Target.calc(self)
        fPSF = ROOT.TF1("fPSF", "1. / (2*TMath::Pi())**0.5 / [0] * TMath::Exp(-(x/[0])**2/2.0)")
        for tE in range(self.energySet.nBin):
            for kAS in range(len(self.aClass)):
                for kSS in range(len(self.aaClass[kAS])):
                    # if kAS==0:
                    #     degPSF = 0.1
                    # elif kAS==1:
                    #     degPSF = max(self.perf.getPSF68(kSS, self.energySet.getBinCenter(tE)), self.perf.getPSF95(kSS, self.energySet.getBinCenter(tE))/2.0)
                    # fPSF.FixParameter(0, radians(degPSF))
                    self.aaGrPerformance[tE][kAS].SetPoint(kSS, self.aaHtgNumSig[kAS][kSS].GetBinContent(tE+1), self.aaHtgNumOff[kAS][kSS].GetBinContent(tE+1))
                    self.aaGrPerformance[tE][kAS].SetPointError(kSS, self.aaHtgNumSig[kAS][kSS].GetBinError(tE+1), self.aaHtgNumOff[kAS][kSS].GetBinError(tE+1))
                    #for lB in range(self.aaaHtgMapAllSky[kAS][kSS][tE].GetNbinsY()):
                    #         bPoint = self.aaaHtgMapAllSky[kAS][kSS][tE].GetYais().GetBinCenter(lB+1)
                    #         bRoiLow = bPoint-degPSF*4
                    #         if bRoiLow<self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetYaxis().GetBinLowEdge(1):
                    #             bRoiLow = self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetYaxis().GetBinUpEdge(self.aaaHtgFoldedMapAllSky[kAS][kSS][tE].GetNbinsX()) + lRoiLow
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
        #print self.aaGrPerformance
        print "Calculation finished."

    def draw(self):
        Target.draw(self)
        aCanAllSky = []
        nDown = len(self.aaClass)
        nAcross = 0
        aaaHtgMapAllSky = []
        for aS in self.aaClass:    
            nAcross = max(nAcross, len(aS))
        cLatitude = ROOT.TCanvas("cLatitude", "Galactic latitude", 1200, 600)
        mDown = int(math.sqrt(self.energySet.nBin+1))
        mAcross = math.ceil((self.energySet.nBin+1) / mDown)
        cLatitude.Divide(mAcross, mDown)
        aHstaLatitude = []
        legLatitude = ROOT.TLegend(0.1,0.1,0.9,0.9,"Event class",'NDC')

        for iE in range(self.energySet.nBin):
            aCanAllSky.append(ROOT.TCanvas("cAllSky%s" % iE, "All sky map in {0:.1f} - {1:.1f} GeV".format(10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3)), 1200, 400))
            aCanAllSky[-1].Divide(nAcross, nDown)
            aaaHtgMapAllSky.append([])
            aHstaLatitude.append(ROOT.THStack("hsLatitude%s" % iE, "Galactic latitude in {0:.1f} - {1:.1f} GeV;sin(b);counts [counts]".format(10**(self.energySet.getBin(iE)[0]-3), 10**(self.energySet.getBin(iE)[1]-3))))
            for iAS in range(len(self.aaClass)):
                aaaHtgMapAllSky[-1].append([])
                for iSS in range(len(self.aaClass[iAS])):
                    aCanAllSky[-1].cd(1+iAS*nAcross+iSS)
                    aCanAllSky[-1].cd(1+iAS*nAcross+iSS).SetGridx(0)
                    aCanAllSky[-1].cd(1+iAS*nAcross+iSS).SetGridy(0)
                    self.aaaHtgMapAllSky[iAS][iSS][iE].Draw('aitoff')
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
        lDown = int(math.sqrt(self.energySet.nBin+1))
        lAcross = math.ceil((self.energySet.nBin+1) / lDown)
        cPerformance.Divide(lAcross, lDown)
        legPerformance = ROOT.TLegend(0.1,0.1,0.9,0.9,"Event category",'NDC')
        for tE in range(self.energySet.nBin):
            aMgrPerformance.append(ROOT.TMultiGraph("mgrPerformance_%s_%s" % (self.name, tE), "Performance plots with {0} in {1:.1f} - {2:.1f} GeV".format(self.name, 10**(self.energySet.getBin(tE)[0]-3), 10**(self.energySet.getBin(tE)[1]-3))))
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
        for iE in range(self.energySet.nBin):
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
