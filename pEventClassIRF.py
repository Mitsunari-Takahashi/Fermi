#!/usr/bin/env python

import sys
import ROOT
import ROOT
from ROOT import TTree
from ROOT import TCut
import numpy
import yaml
from array import array
import math
#from pCutBDT import cutBDT
from pAnalysisConfig import *

class EventClassIRF:
    def __init__(self, nameClass, aCut, aHtgIRF, nameBDT, energyRegion):
        self.energyRegion = energyRegion
        self.aCutValue = aCutBDT
        self.aTCut = []
        self.htgAcceptance = aHtgIRF[0]
        # PSF
        self.htgPSF68 = aHtgIRF[2]
        fPSFScale = ROOT.TF1("fPSFScale", "sqrt(pow([0]*pow(pow(x, 10)/100, -[2]), 2) + pow([1], 2))", self.energyRegion.edgeLow, self.energyRegion.edgeUp)
        fPSFScale->SetParameter(2, 0.8)
        self.htgPSF68.ProjectionX().Fit(fPSFScale)
        self.fPSFScale = fPSFScale

        # Energy dispersion
        self.htgEdisp68 = aHtgIRF[1]
        fEdispScale = ROOT.TF2("fEdipsScale", "[0]*pow(x,2) + [1]*pow(y,2) + [2]*x + [3]*y + [4]*x*y + [5]", self.energyRegion.edgeLow, self.energyRegion.edgeUp, 0.0, 1.0)
        self.htgEdisp68.Fit(fEdispScale)
        self.fEdispScale = fEdispScale

        for iE in range(self.energyRegion.nBin):
            self.aTCut.append(ROOT.TCut(""))
#            self.aTCut.append(ROOT.TCut(Form("(%s>=%d)", nameBDT, self.aCutValue[iE])))
            
    def addCut(aCutIn):
        if self.energyRegion.nBin == len(aCutIn):
            self.cTCut = aCutIn
        else:
            print "Array length are mismatched!!!"
            self.cTCut = aCutIn

    def makeTable(self, trAG):
        htgPSF = ROOT.TH1D("htgPSF", "htgPSH", 1000, 0.0, 100)
        htgEdisp = ROOT.TH1D("htgEdisp", "htgEdisp", 1100, -1.0, 10)
        self.aaPSF = []
        self.aaEdisp = []
        for iE in range(self.energyRegion.nBin):
            self.aaPSF.append([])
            self.aaEdisp.append([])
            for iC in range(10):
                # PSF
                strPSFDivider = self.fPSFScale.GetExpFormula()
                strPSFDivider.RelaceAll("[0]", Form("%f", fPSFScale->GetParameter(0)))
                strPSFDivider.RelaceAll("[1]", Form("%f", fPSFScale->GetParameter(1)))
                strPSFDivider.RelaceAll("(x", "(log10(WP8CalOnlyEnergy)")
                strPSFDivider.RelaceAll("*x", "*log10(WP8CalOnlyEnergy)")
                trAG.Draw(Form("TMath.ACos(-Cal1MomXDir*McXDir-Cal1MomYDir*McYDir-Cal1MomZDir*McZDir)*TMath.RadToDeg()/(%s)>>%s", strPSFDivider, htgPSF.GetName()), self.aTCut[iE])
                self.aaPSF[-1].append(htgPSF.GetMean())

                # Edisp
                strEdispDivider = self.fEdispScale.GetExpFormula()
                strEdispDivider.RelaceAll("[0]", Form("%f", fEdispScale->GetParameter(0)))
                strEdispDivider.RelaceAll("[1]", Form("%f", fEdispScale->GetParameter(1)))
                strEdispDivider.RelaceAll("[2]", Form("%f", fEdispScale->GetParameter(2)))
                strEdispDivider.RelaceAll("[3]", Form("%f", fEdispScale->GetParameter(3)))
                strEdispDivider.RelaceAll("[4]", Form("%f", fEdispScale->GetParameter(4)))
                strEdispDivider.RelaceAll("[5]", Form("%f", fEdispScale->GetParameter(5)))
                strEdispDivider.RelaceAll("(x", "(log10(WP8CalOnlyEnergy)")
                strEdispDivider.RelaceAll("(y", "((-Cal1MomZDir)")
                strEdispDivider.RelaceAll("*x", "*log10(WP8CalOnlyEnergy)")
                strEdispDivider.RelaceAll("*y", "*(-Cal1MomZDir)")
                trAG.Draw(Form("((WP8CalOnlyEnergy-McEnergy)/McEnergy)/(%s)>>%s", strEdispDivider, htgEdisp.GetName()), self.aTCut[iE])
                self.aaEdisp[-1].append(htgEdisp.GetMean())

    def writeFITS(self):
        
