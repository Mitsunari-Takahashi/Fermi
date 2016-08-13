#!/usr/bin/env python
"""For study of CalOnly direction dispersion
The histogram of the mis-reconstructed direction is fitted with single King fucntion.
"""
import sys
import math
from math import log10, acos, atan2, pi, degrees, sin, cos, asin
import numpy as np
from array import array
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TTree
from array import array
from pColor import *
ROOT.gROOT.SetBatch() 

par = sys.argv
print par
strPathFileIn = "/nfs/farm/g/glast/u/mtakahas/data/MC/AG200909_62_2016Jun.root"
#strPathFileFr = "/u/gl/mtakahas/work/data/MC/AG200909_62_2016Jun_S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_BDTG1000D06.root"
#strCut="( FswGamState==0 && (TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy  >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000)) ) && Cal1RawEnergySum>=20000 ) && ( (log10(WP8CalOnlyEnergy)>=4.35&&log10(WP8CalOnlyEnergy)<4.55&&S16RAWE20ZDIR010ZCS000catTwoZDIR050_15_BDTG1000D06Log>=2.75875) || (log10(WP8CalOnlyEnergy)>=4.55&&log10(WP8CalOnlyEnergy)<4.75&&S16RAWE20ZDIR010ZCS000catTwoZDIR050_15_BDTG1000D06Log>=2.56125) || (log10(WP8CalOnlyEnergy)>=4.75&&log10(WP8CalOnlyEnergy)<4.95&&S16RAWE20ZDIR010ZCS000catTwoZDIR050_15_BDTG1000D06Log>=2.38875) || (log10(WP8CalOnlyEnergy)>=4.95&&log10(WP8CalOnlyEnergy)<5.15&&S16RAWE20ZDIR010ZCS000catTwoZDIR050_15_BDTG1000D06Log>=2.44875) || (log10(WP8CalOnlyEnergy)>=5.15&&log10(WP8CalOnlyEnergy)<5.35&&S16RAWE20ZDIR010ZCS000catTwoZDIR050_15_BDTG1000D06Log>=2.71875) || (log10(WP8CalOnlyEnergy)>=5.35&&log10(WP8CalOnlyEnergy)<5.55&&S16RAWE20ZDIR010ZCS000catTwoZDIR050_15_BDTG1000D06Log>=2.79375) || (log10(WP8CalOnlyEnergy)>=5.55&&log10(WP8CalOnlyEnergy)<5.75&&S16RAWE20ZDIR010ZCS000catTwoZDIR050_15_BDTG1000D06Log>=2.91625) )"
strCut="( (TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy  >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000)) ) && Cal1RawEnergySum>=20000 && FswGamState==0 )"

fileIn = ROOT.TFile(strPathFileIn, 'READ')
trAG  = fileIn.Get('MeritTuple')
#if len(strPathFileFr)>0:
#    trAG.AddFriend("frTemp = MeritTuple", strPathFileFr)

fileOut = ROOT.TFile("AG_dispersion.root", "RECREATE")
trOut = TTree('trOut','trKingParameters')
e = np.zeros(1, dtype=float)
th = np.zeros(1, dtype=float)
apar = []
for ip in range(6):
    apar.append(np.zeros(1, dtype=float))
trOut.Branch('Energy',e,'Energy/D')
trOut.Branch('Theta',th,'Theta/D')
trOut.Branch('nCore',apar[0],'nCore/D')
trOut.Branch('sCore',apar[1],'sCore/D')
trOut.Branch('gCore',apar[2],'gCore/D')
fileOut.cd()

nBinE = 7
nBinCth = 8
lAxE = 4.35
uAxE = 5.75
lAxCth = 0.2
uAxCth = 1.0
hPSF = ROOT.TH3D("hPSF", "AG MC", nBinE, lAxE, uAxE, nBinCth, lAxCth, uAxCth, 60, 0, 0.300)
trAG.Draw("TMath::ACos(-McXDir*Cal1MomXDir-McYDir*Cal1MomYDir-McZDir*Cal1MomZDir):-McZDir:McLogEnergy>>{0}".format(hPSF.GetName()), strCut, "goff")
print hPSF.GetEntries(), "events"

aCanPSF = []
#cPSF.Divide(4,2)
aHsPsfProj = []
hsPsfProjAve = ROOT.THStack("hsPsfProjAve", "Averaged over #theta")
aHtgPsfProj = []
aFcKing = []
aFcKingAve = []
aHtgParKing = []
aHtgParKing.append(ROOT.TH2D('htgKingN', 'Normalization of King function;McLogEnergy;-McZDir', nBinE, lAxE, uAxE, nBinCth, lAxCth, uAxCth))
aHtgParKing.append(ROOT.TH2D('htgKingS', 'Sigma of King function;McLogEnergy;-McZDir', nBinE, lAxE, uAxE, nBinCth, lAxCth, uAxCth))
aHtgParKing.append(ROOT.TH2D('htgKingG', 'Gamma of King function;McLogEnergy;-McZDir', nBinE, lAxE, uAxE, nBinCth, lAxCth, uAxCth))
aHtgParKingAve = []
aHtgParKingAve.append(ROOT.TH1D('htgKingAveN', 'Normalization averaged over #theta;McLogEnergy', nBinE, lAxE, uAxE))
aHtgParKingAve.append(ROOT.TH1D('htgKingAveS', 'Sigma averaged over #theta;McLogEnergy', nBinE, lAxE, uAxE))
aHtgParKingAve.append(ROOT.TH1D('htgKingAveG', 'Gamma averaged over #theta;McLogEnergy', nBinE, lAxE, uAxE))

prmKing = [5e-05, 1.6e-02, 2.3e+02]

for iE in range(nBinE):
    print "================================"
    print "{0:.1f} - {1:.1f} GeV".format(pow(10, hPSF.GetXaxis().GetBinLowEdge(iE+1)-3), pow(10, hPSF.GetXaxis().GetBinLowEdge(iE+2)-3))
    aHsPsfProj.append(ROOT.THStack("hsPsfProj{0}".format(iE+1), "PSF dispersion in {0:.1f} - {1:.1f} GeV".format(pow(10, hPSF.GetXaxis().GetBinLowEdge(iE+1)-3), pow(10, hPSF.GetXaxis().GetBinLowEdge(iE+2)-3))))
    e[0] = hPSF.GetXaxis().GetBinCenter(iE+1)
    aCanPSF.append(ROOT.TCanvas("cPSF{0}".format(iE+1), "PSF dispersion in {0:.1f} - {1:.1f} GeV".format(pow(10, hPSF.GetXaxis().GetBinLowEdge(iE+1)-3), pow(10, hPSF.GetXaxis().GetBinLowEdge(iE+2)-3)), 1200, 600))
    aHtgPsfProj.append([])
    aFcKing.append([])
    for iTh in range(nBinCth+1):
        if iTh==0:
            th[0] = 0
            print "Cos(theta): Averaged."
            aHtgPsfProj[-1].append(hPSF.ProjectionZ("hPsfProj_{0}_{1}".format(iE+1, iTh), iE+1, iE+1, 1, hPSF.GetNbinsY()))
            hsPsfProjAve.Add(aHtgPsfProj[-1][-1])
        else:
            th[0] = hPSF.GetYaxis().GetBinCenter(iTh)
            print "Cos(theta):", hPSF.GetYaxis().GetBinLowEdge(iTh), hPSF.GetYaxis().GetBinLowEdge(iTh+1)
            aHtgPsfProj[-1].append(hPSF.ProjectionZ("hPsfProj_{0}_{1}".format(iE+1, iTh), iE+1, iE+1, iTh, iTh))
        aFcKing[-1].append(ROOT.TF1("fcKing{0}_{1}".format(iE+1, iTh), "TMath::Sin(x)*([0]*(1.-1./[2])*pow(1.+(x/[1])**2/2./[2], -[2])/2./TMath::Pi()/[1]**2)", 0, 1.0))
        print aHtgPsfProj[-1][-1].GetEntries(), "entries."

        for iPrm in range(len(prmKing)):
            aFcKing[-1][-1].SetParameter(iPrm, prmKing[iPrm])
        aFcKing[-1][-1].SetParLimits(0, 0, 1e-1)
        aFcKing[-1][-1].SetParLimits(1, 0, 1e-1)
        if aHtgPsfProj[-1][-1].Integral()>0:
            aHtgPsfProj[-1][-1].Sumw2()
            aHtgPsfProj[-1][-1].Scale(1./aHtgPsfProj[-1][-1].Integral())
            aHtgPsfProj[-1][-1].Write()
            aHtgPsfProj[-1][-1].SetLineColor(akColor(iTh))
            aHtgPsfProj[-1][-1].SetMarkerColor(akColor(iTh))
            aHtgPsfProj[-1][-1].SetLineWidth(3)
            aHtgPsfProj[-1][-1].SetFillStyle(0)
            aFcKing[-1][-1].SetLineColor(akColor(iTh))
            aFcKing[-1][-1].SetLineStyle(1)
            
            aHtgPsfProj[-1][-1].Fit(aFcKing[-1][-1], "", "", 0, aHtgPsfProj[-1][-1].GetXaxis().GetBinUpEdge(aHtgPsfProj[-1][-1].GetNbinsX()))
            for jp in range(len(prmKing)):
                apar[jp][0] = aFcKing[-1][-1].GetParameter(jp)
                print "Parameter", jp, ":", apar[jp][0]
                if iTh>0:
                    aHtgParKing[jp].SetBinContent(iE+1, iTh, apar[jp][0])
                else:
                    aHtgParKingAve[jp].SetBinContent(iE+1, apar[jp][0])
            trOut.Fill()
            aHsPsfProj[-1].Add(aHtgPsfProj[-1][-1])
            aFcKing[-1][-1].Write()
    aHsPsfProj[-1].Write()
    #aHtgPsfProjAve[-1].Write()
    aCanPSF[-1].cd()
    aHsPsfProj[-1].Draw("nostack")
    aCanPSF[-1].Write()
hPSF.Write()
#cPSF.Write()
trOut.Write()
for hpk in aHtgParKing:
    hpk.Write()
for hpk in aHtgParKingAve:
    hpk.Write()
aFcAve = []
cPsfAve = ROOT.TCanvas("cPsfAve", "PSF dispersion averaged over #theta", 1200, 600)
print "================================"
for iE in range(hPSF.GetNbinsX()):
    print "{0:.1f} - {1:.1f} GeV".format(pow(10, hPSF.GetXaxis().GetBinLowEdge(iE+1)-3), pow(10, hPSF.GetXaxis().GetBinLowEdge(iE+2)-3))
    aHtgPsfProj[iE][0].SetLineColor(akColor(iE))
    aHtgPsfProj[iE][0].SetMarkerColor(akColor(iE))
    aFcAve.append(aHtgPsfProj[iE][0].GetFunction(aFcKing[iE][0].GetName()))
    aFcAve[-1].SetLineColor(akColor(iE))
    aFcAve[-1].SetLineStyle(1)
    cPsfAve.cd()
    hsPsfProjAve.Draw("nostack")
    aFcAve[-1].Draw("same")
    iStep=1
    rad68 = -1
    rad95 = -1
    while iStep<1000 and (rad68==-1 or rad95==-1):
        rCut = float(iStep)/1000.
        limUp = hPSF.GetXaxis().GetBinUpEdge(hPSF.GetNbinsX())
        ratio = aFcAve[-1].Integral(0, rCut*limUp) / aFcAve[-1].Integral(0, limUp)
        if ratio>0.68 and rad68==-1:
            rad68 = rCut*limUp
            print "PSF 68%:", rad68*180./pi, "deg"
        if ratio>0.95 and rad95==-1:
            rad95 = rCut*limUp
            print "PSF 95%:", rad95*180./pi, "deg"
        iStep = iStep+1
    aFcAve[-1].Write()
    cPsfAve.Write()
hsPsfProjAve.Write()
