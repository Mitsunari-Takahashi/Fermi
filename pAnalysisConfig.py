from ROOT import TH2D
from ROOT import TH2F
from ROOT import TAxis
import ROOT
#from array import array
import numpy as np
class EnergyLogRegion:
    def __init__(self, numEnergyLogBin, energyLogStart, energyLogWidth=0.25):
        self.nBin = int(numEnergyLogBin)
        self.edgeLow = float(energyLogStart)
        self.wBin = float(energyLogWidth)
        self.edgeUp = self.edgeLow + self.nBin * self.wBin
        self.aBin = [self.edgeLow]
        for iBin in range(self.nBin):
            self.aBin.append(self.aBin[-1]+self.wBin)
    def getRegionLowEdge(self):
        return self.edgeLow
    def getRegionUpEdge(self):
        return self.edgeUp
    def printR(self):
        print self.aBin
    def getBin(self, jBin):
        return self.aBin[jBin:jBin+2]
    def getBinCenter(self, jBin):
        return (self.aBin[jBin]+self.aBin[jBin+1])/2.0
    def findBin(self, energyLog):
        return int((energyLog-self.edgeLow)/self.wBin)
    def findBinForce(self, energyLog):
        return min(max(0, int((energyLog-self.edgeLow)/self.wBin)), nBin-1)

class EnergyLogRegionPlot:
    def __init__(self, numEnergyLogBin, energyLogStart, energyLogWidth=0.25):
        self.nBin = int(numEnergyLogBin)
        self.edgeLow = float(energyLogStart)
        self.wBin = float(energyLogWidth)
        self.edgeUp = self.edgeLow + self.nBin * self.wBin
        self.aBin = [self.edgeLow]
        for iBin in range(self.nBin):
            self.aBin.append(self.aBin[-1]+self.wBin)
    def printR(self):
        print self.aBin
    def getBin(self, jBin):
        return self.aBin[jBin:jBin+2]
    def getBinCenter(self, jBin):
        return (self.aBin[jBin]+self.aBin[jBin+1])/2.0
    def findBin(self, energyLog):
        return int((energyLog-self.edgeLow)/self.wBin)
    def findBinForce(self, energyLog):
        return min(max(0, int((energyLog-self.edgeLow)/self.wBin)), nBin-1)

class ClassConfig:
    def __init__(self, eventCategory="Both", cutLevel=[10, 2, 1], binStart=1, zCut=[100.0, 100.0]):
        self.aCutEGB = cutLevel
        self.aStrSelectCalOnly = []
        for vCutEGB in self.aCutEGB:
            self.aStrSelectCalOnly.append("CalOnly_R%s0" % vCutEGB)
        self.aStrSelectCalTkr = ["P8R1_TRANSIENT_R100", "P8R1_SOURCE"]
        self.zCut = []
        if eventCategory == "Both":
            self.aStrSelect = ["CalTkr", "CalOnly"]
            self.aaStrSelect = [self.aStrSelectCalTkr, self.aStrSelectCalOnly]
            self.zCut = zCut
        elif eventCategory == "CalTkr":
            self.aStrSelect = ["CalTkr"]
            self.aaStrSelect = [self.aStrSelectCalTkr]
            self.zCut.append(zCut[0])
        elif eventCategory == "CalOnly":
            self.aStrSelect = ["CalOnly"]
            self.aaStrSelect = [self.aStrSelectCalOnly]
            self.zCut.append(zCut[1])
        else:
            print 'Please input "Both", "CalTkr" or "CalOnly" as the 1st arg.'
        self.nStartBin = binStart

class CutBDT:
    def __init__(self, resultBDT, aCutEGB=[10, 3, 1]):
        self.aCutEGB = aCutEGB
        if resultBDT[-5:] == ".root":
            fileRoc = ROOT.TFile(resultBDT, 'READ')
            print fileRoc.GetName(), " was opened."
            print "Cutting at", aCutEGB, "xEGB level."
            h2Sig = fileRoc.Get("sig_acc")
            h2Bkg = fileRoc.Get("bkg_rate")
            h2Count = fileRoc.Get("bkg_counts_cum_hist")
            h2Egb = fileRoc.Get("egb_rate")

            self.aaValCutBDT = [{'edgeLow':h2Sig.GetXaxis().GetBinLowEdge(1), 'widthBin':h2Sig.GetXaxis().GetBinWidth(1), 'numBin':h2Sig.GetXaxis().GetNbins()}]
            self.aaNumGamCut = []
            self.aaNumRemCut = []
            self.aaNumCountCut = []
            for ie in range(self.aaValCutBDT[0]['numBin']):
                ebLow = self.aaValCutBDT[0]['widthBin']*ie+self.aaValCutBDT[0]['edgeLow']
                ebUp = self.aaValCutBDT[0]['widthBin']*(ie+1)+self.aaValCutBDT[0]['edgeLow']
                print '=========', "{0:.1f}".format(10.**(ebLow-3.)),' - ', "{0:.1f}".format(10.**(ebUp-3.)), ' GeV ========='
                self.aaValCutBDT.append([])
                abFound = []
                self.aaNumGamCut.append([])
                self.aaNumRemCut.append([])
                self.aaNumCountCut.append([])
                for jc in range(len(aCutEGB)):
                    abFound.append(False)
                h1Sig = h2Sig.ProjectionY("h1SigAcc", ie+1, ie+1)
                h1Bkg = h2Bkg.ProjectionY("h1BkgAcc", ie+1, ie+1)
                h1Count = h2Count.ProjectionY("h1BkgCount", ie+1, ie+1)
                h1Egb = h2Egb.ProjectionY("h1EgbRate", ie+1, ie+1)
                for ib in range(1, h1Sig.GetNbinsX()+1):
                    nEgb = h1Egb.GetBinContent(ib)
                    nSig = h1Sig.GetBinContent(ib)
                    nRem = h1Bkg.GetBinContent(ib)
                    nCount = h1Count.GetBinContent(ib)
                    for ic in range(len(aCutEGB)):
                        if (abFound[ic]==False and nRem<=aCutEGB[ic]*nEgb):
                            abFound[ic] = True
                            self.aaValCutBDT[ie+1].append(h1Sig.GetBinCenter(ib))
                            self.aaNumGamCut[-1].append(nSig)
                            self.aaNumRemCut[-1].append(nRem)
                            self.aaNumCountCut[-1].append(nCount)
                            print "--------- Background level:", aCutEGB[ic], "x EGB ---------"
                            print "Cut value:", self.aaValCutBDT[ie+1][-1]
                            print "Acceptance:", self.aaNumGamCut[-1][-1], "[m^2 sr]"
                            print "Background rate:", self.aaNumRemCut[-1][-1], "[MeV sr^-1 s^-1]"
                            print "Background count:", self.aaNumCountCut[-1][-1], "[events]"
        elif resultBDT[-2:] == ".P":
            fileIn = open(resultBDT, 'r')
            self.aaValCutBDT = pickle.load(fileIn)
            print self.aaValCutBDT
            fileIn.close()
                
    def write(self, nameFileOut):
        fileOut = open(nameFileOut, 'w')
        pickle.dump(self.aaValCutBDT, fileOut)
        print "BDT cut values have been  written in", nameFileOut
        fileOut.close()

class CutBDT_2D:
    def __init__(self, aCutBDT, aCosZ, strSuffix=""):
        arCosZ = np.array(aCosZ, dtype=np.float)
        fileOut = ROOT.TFile("hBDT_CutValues_{0}.root".format(strSuffix), 'RECREATE')
        self.aCutEGB = aCutBDT[0].aCutEGB
        self.aaValCutBDT = []
        self.aHtgCutBDT = []
        self.aHtgGamCut = []
        self.aHtgRemCut = []
        self.aHtgCountCut = []
        if len(aCutBDT)==len(aCosZ)-1:
            self.aaValCutBDT.append(aCutBDT[0].aaValCutBDT[0])
            self.aaValCutBDT[0]['aEdgeLowY'] = aCosZ[0:-1]
            self.aaValCutBDT[0]['aEdgeUpY'] = aCosZ[1:]
            self.aaValCutBDT[0]['numBinY'] = len(aCosZ)-1
            for ie in range(self.aaValCutBDT[0]['numBin']):
                self.aaValCutBDT.append([])
                for iEGB in range(len(self.aCutEGB)):
                    self.aaValCutBDT[-1].append([])
                    for ith in range(self.aaValCutBDT[0]['numBinY']):
#                        print self.aaValCutBDT
#                        print self.aaValCutBDT[-1][-1]
#                        print ie+1, iEGB
                        self.aaValCutBDT[-1][-1].append(aCutBDT[ith].aaValCutBDT[ie+1][iEGB])
            print ""
            print "##############################################"
            print "#############   Combined cut   ###############"
            print "##############################################"
            print ""
            for iEGB in range(len(self.aCutEGB)):
                self.aHtgCutBDT.append(ROOT.TH2D("htgCutBDT{0}".format(iEGB), "BDT cut value for CalOnly_R{0}".format(int(self.aCutEGB[iEGB]*10)), self.aaValCutBDT[0]['numBin'], self.aaValCutBDT[0]['edgeLow'], self.aaValCutBDT[0]['edgeLow']+self.aaValCutBDT[0]['numBin']*self.aaValCutBDT[0]['widthBin'], self.aaValCutBDT[0]['numBinY'], arCosZ))
                self.aHtgGamCut.append(ROOT.TH2D("htgGamCut{0}".format(iEGB), "Gamma-ray acceptance value of CalOnly_R{0} [m^2 sr];logWP8CalOnlyEnergy;Cal1MomZDir".format(int(self.aCutEGB[iEGB]*10)), self.aaValCutBDT[0]['numBin'], self.aaValCutBDT[0]['edgeLow'], self.aaValCutBDT[0]['edgeLow']+self.aaValCutBDT[0]['numBin']*self.aaValCutBDT[0]['widthBin'], self.aaValCutBDT[0]['numBinY'], arCosZ))
                self.aHtgRemCut.append(ROOT.TH2D("htgRemCut{0}".format(iEGB), "Residual background rate of CalOnly_R{0} [MeV sr^-1 s^-1];logWP8CalOnlyEnergy;Cal1MomZDir".format(int(self.aCutEGB[iEGB]*10)), self.aaValCutBDT[0]['numBin'], self.aaValCutBDT[0]['edgeLow'], self.aaValCutBDT[0]['edgeLow']+self.aaValCutBDT[0]['numBin']*self.aaValCutBDT[0]['widthBin'], self.aaValCutBDT[0]['numBinY'], arCosZ))
                self.aHtgCountCut.append(ROOT.TH2D("htgCountCut{0}".format(iEGB), "Residual background counts of CalOnly_R{0} [events];logWP8CalOnlyEnergy;Cal1MomZDir".format(int(self.aCutEGB[iEGB]*10)), self.aaValCutBDT[0]['numBin'], self.aaValCutBDT[0]['edgeLow'], self.aaValCutBDT[0]['edgeLow']+self.aaValCutBDT[0]['numBin']*self.aaValCutBDT[0]['widthBin'], self.aaValCutBDT[0]['numBinY'], arCosZ))
                for ie in range(self.aaValCutBDT[0]['numBin']):
                    print '=========', "{0:.1f}".format(10.**(self.aHtgGamCut[0].GetXaxis().GetBinLowEdge(ie+1)-3.)),' - ', "{0:.1f}".format(10.**(self.aHtgGamCut[0].GetXaxis().GetBinUpEdge(ie+1)-3.)), ' GeV ========='
                    for ith in range(self.aaValCutBDT[0]['numBinY']):
                        self.aHtgCutBDT[-1].SetBinContent(ie+1, ith+1, aCutBDT[ith].aaValCutBDT[ie+1][iEGB])
                        self.aHtgGamCut[-1].SetBinContent(ie+1, ith+1, aCutBDT[ith].aaNumGamCut[ie][iEGB])
                        self.aHtgRemCut[-1].SetBinContent(ie+1, ith+1, aCutBDT[ith].aaNumRemCut[ie][iEGB])
                        self.aHtgCountCut[-1].SetBinContent(ie+1, ith+1, aCutBDT[ith].aaNumCountCut[ie][iEGB])
                    print "--------- Background level:", self.aCutEGB[iEGB], "x EGB ---------"
                    print "Cut value:", self.aaValCutBDT[ie+1][iEGB]
                    print "Acceptance:", self.aHtgGamCut[iEGB].Integral(ie+1, ie+1, 1, self.aHtgGamCut[iEGB].GetNbinsY()), "[m^2 sr]"
                    print "Background rate:", self.aHtgRemCut[iEGB].Integral(ie+1, ie+1, 1, self.aHtgRemCut[iEGB].GetNbinsY()), "[MeV sr^-1 s^-1]"
                    print "Background count:", self.aHtgCountCut[iEGB].Integral(ie+1, ie+1, 1, self.aHtgCountCut[iEGB].GetNbinsY()), "[events]"
                fileOut.cd()
                self.aHtgCutBDT[-1].Write()
                self.aHtgGamCut[-1].Write()
                self.aHtgRemCut[-1].Write()
                self.aHtgCountCut[-1].Write()
    def getHtgCutBDT(self):
        print self.aHtgGamCut
        return self.aHtgCutBDT

class CutPerformanceHtg:
#    def __init__(self, aPathFilePerf, classCut, rEnergy):
    def __init__(self, aPathFilePerf):
        self.aRootFile = []
        self.aHtgAccept = []
        self.aHtgEdisp68 = []
        self.aHtgEdisp95 = []
        self.aHtgPSF68 = []
        self.aHtgPSF95 = []
        self.aHtgEdisp68_cth = []
        self.aHtgEdisp95_cth = []
        self.aHtgPSF68_cth = []
        self.aHtgPSF95_cth = []
        self.aFcPSF_cth = []
        for fcl in aPathFilePerf:
            self.aRootFile.append([])
            self.aHtgAccept.append([])
            self.aHtgEdisp68.append([])
            self.aHtgEdisp95.append([])
            self.aHtgPSF68.append([])
            self.aHtgPSF95.append([])
            self.aHtgEdisp68_cth.append([])
            self.aHtgEdisp95_cth.append([])
            self.aHtgPSF68_cth.append([])
            self.aHtgPSF95_cth.append([])
            self.aFcPSF_cth.append([])
            for gcl in fcl:
                self.aRootFile[-1].append(ROOT.TFile(gcl, 'READ'))
                self.aHtgAccept[-1].append(self.aRootFile[-1][-1].Get("acc_hist"))
                self.aHtgEdisp68[-1].append(self.aRootFile[-1][-1].Get("edisp_q68_hist"))
                self.aHtgEdisp95[-1].append(self.aRootFile[-1][-1].Get("edisp_q95_hist"))
                self.aHtgPSF68[-1].append(self.aRootFile[-1][-1].Get("psf_q68_hist"))
                self.aHtgPSF95[-1].append(self.aRootFile[-1][-1].Get("psf_q95_hist"))
                self.aHtgEdisp68_cth[-1].append(self.aRootFile[-1][-1].Get("edisp_cth_q68_hist"))
                self.aHtgEdisp95_cth[-1].append(self.aRootFile[-1][-1].Get("edisp_cth_q95_hist"))
                self.aHtgPSF68_cth[-1].append(self.aRootFile[-1][-1].Get("psf_cth_q68_hist"))
                self.aHtgPSF95_cth[-1].append(self.aRootFile[-1][-1].Get("psf_cth_q95_hist"))
            #self.aFcPSF_cth.append(ROOT.TF1("fcPSF_%s" % fcl.GetName()[:-5], "1. / (2*TMath::Pi())**0.5 / [0] * TMath::Exp(-(x/[0])**2/2.0)"))
            #self.aFcPSF_cth[-1].SetParameter(0, 3.0)
            #self.aFcPSF_cth[-1].Fit
            #self.aRootFile[-1].Close()

    def getAcceptance(self, iCat, iCut, energy):
        return self.aHtgAccept[iCat][iCut].GetBinContent(self.aHtgAccept[iCat][iCut].FindBin(energy))
    def getEdisp68(self, iCat, iCut, energy):
        return self.aHtgEdisp68[iCat][iCut].GetBinContent(self.aHtgEdisp68[iCat][iCut].FindBin(energy))
    def integralEdisp68(self, iCat, iCut, iebin0, iebin1):
        if not iebin0<iebin1:
            print "Integral window is mismatched."
            sys.exit(1)
        wsum = 0.0
        asum = 0.0
        for jBin in range(iebin0, iebin1+1):
            wsum = wsum + self.aHtgEdisp68[iCat][iCut].GetBinContent(jBin)*self.aHtgAccept[iCat][iCut].GetBinContent(jBin)
            asum = asum + self.aHtgAccept[iCat][iCut].GetBinContent(jBin)
        return wsum/asum
    def getEdisp95(self, iCat, iCut, energy):
        return self.aHtgEdisp95[iCat][iCut].GetBinContent(self.aHtgEdisp95[iCat][iCut].FindBin(energy))
    def integralEdisp95(self, iCat, iCut, iebin0, iebin1):
        if not iebin0<iebin1:
            print "Integral window is mismatched."
            sys.exit(1)
        wsum = 0.0
        asum = 0.0
        for jBin in range(iebin0, iebin1+1):
            wsum = wsum + self.aHtgEdisp95[iCat][iCut].GetBinContent(jBin)*self.aHtgAccept[iCat][iCut].GetBinContent(jBin)
            asum = asum + self.aHtgAccept[iCat][iCut].GetBinContent(jBin)
        return wsum/asum
    def getPSF68(self, iCat, iCut, energy):
        return self.aHtgPSF68[iCat][iCut].GetBinContent(self.aHtgPSF68[iCat][iCut].FindBin(energy))
    def integralPSF68(self, iCat, iCut, iebin0, iebin1):
        if not iebin0<iebin1:
            print "Integral window is mismatched."
            sys.exit(1)
        wsum = 0.0
        asum = 0.0
        for jBin in range(iebin0, iebin1+1):
            wsum = wsum + self.aHtgPSF68[iCat][iCut].GetBinContent(jBin)*self.aHtgAccept[iCat][iCut].GetBinContent(jBin)
            asum = asum + self.aHtgAccept[iCat][iCut].GetBinContent(jBin)
        return wsum/asum
    def getPSF95(self, iCat, iCut, energy):
        return self.aHtgPSF95[iCat][iCut].GetBinContent(self.aHtgPSF95[iCat][iCut].FindBin(energy))
    def integralPSF95(self, iCat, iCut, iebin0, iebin1):
        if not iebin0<iebin1:
            print "Integral window is mismatched."
            sys.exit(1)
        wsum = 0.0
        asum = 0.0
        for jBin in range(iebin0, iebin1+1):
            wsum = wsum + self.aHtgPSF95[iCat][iCut].GetBinContent(jBin)*self.aHtgAccept[iCat][iCut].GetBinContent(jBin)
            asum = asum + self.aHtgAccept[iCat][iCut].GetBinContent(jBin)
        return wsum/asum
    def getEdisp68_cth(self, iCat, iCut, energy, cth):
        return self.aHtgEdisp68_cth[iCat][iCut].GetBinContent(self.aHtgEdisp68_cth[iCat][iCut].FindBin(energy, cth))
    def getEdisp95_cth(self, iCat, iCut, energy, cth):
        return self.aHtgEdisp95_cth[iCat][iCut].GetBinContent(self.aHtgEdisp95_cth[iCat][iCut].FindBin(energy, cth))
    def getPSF68_cth(self, iCat, iCut, energy, cth):
        return self.aHtgPSF68_cth[iCat][iCut].GetBinContent(self.aHtgPSF68_cth[iCat][iCut].FindBin(energy, cth))
    def getPSF95_cth(self, iCat, iCut, energy, cth):
        return self.aHtgPSF95_cth[iCat][iCut].GetBinContent(self.aHtgPSF95_cth[iCat][iCut].FindBin(energy, cth))
