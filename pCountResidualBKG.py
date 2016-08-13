#!/usr/bin/env python
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

pathFileBKG = '/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_Combined.root'
fileBKG = ROOT.TFile(pathFileBKG, 'READ')
trBKG = fileBKG.Get('MeritTuple')
print trBKG.GetName(), 'is found.'

pathFriendBDT = '/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_Combined_S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060.root'
nameBDT = 'S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060'
pathFileOut = '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/{0}_Residual.root'.format(nameBDT)
fileOut = ROOT.TFile(pathFileOut, 'UPDATE')

fileRaw = ROOT.TFile('/nfs/farm/g/glast/u/mtakahas/data/MC/hRawBKG.root', 'READ')
hRaw = fileRaw.Get('hRawBKG')

dictPathFileAcc = { 'CalOnlyR100':'/u/gl/mtakahas/work/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R100_perf.root', 'CalOnlyR30':'/u/gl/mtakahas/work/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R30_perf.root', 'CalOnlyR10':'/u/gl/mtakahas/work/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R10_perf.root' }

trBKG.AddFriend("{0} = MeritTuple".format(nameBDT), pathFriendBDT)
trBKG.SetAlias("S18ZDIR020catTwoZDIR060Log", "-log10(1.0-(1.0+{0})/2.0)".format(nameBDT))
dictStrCutBDT = { 'CalOnlyR100':'((Cal1MomZDir>=0.2&&Cal1MomZDir<0.6)&&((log10(WP8CalOnlyEnergy)>=4.35&&log10(WP8CalOnlyEnergy)<4.55&&S18ZDIR020catTwoZDIR060Log>=1.42875) || (log10(WP8CalOnlyEnergy)>=4.55&&log10(WP8CalOnlyEnergy)<4.75&&S18ZDIR020catTwoZDIR060Log>=1.52375) || (log10(WP8CalOnlyEnergy)>=4.75&&log10(WP8CalOnlyEnergy)<4.95&&S18ZDIR020catTwoZDIR060Log>=1.72125) || (log10(WP8CalOnlyEnergy)>=4.95&&log10(WP8CalOnlyEnergy)<5.15&&S18ZDIR020catTwoZDIR060Log>=1.97875) || (log10(WP8CalOnlyEnergy)>=5.15&&log10(WP8CalOnlyEnergy)<5.35&&S18ZDIR020catTwoZDIR060Log>=2.12875) || (log10(WP8CalOnlyEnergy)>=5.35&&log10(WP8CalOnlyEnergy)<5.55&&S18ZDIR020catTwoZDIR060Log>=2.25875) || (log10(WP8CalOnlyEnergy)>=5.55&&log10(WP8CalOnlyEnergy)<5.75&&S18ZDIR020catTwoZDIR060Log>=2.36125))) || ((Cal1MomZDir>=0.6)&&((log10(WP8CalOnlyEnergy)>=4.35&&log10(WP8CalOnlyEnergy)<4.55&&S18ZDIR020catTwoZDIR060Log>=1.12625) || (log10(WP8CalOnlyEnergy)>=4.55&&log10(WP8CalOnlyEnergy)<4.75&&S18ZDIR020catTwoZDIR060Log>=1.52125) || (log10(WP8CalOnlyEnergy)>=4.75&&log10(WP8CalOnlyEnergy)<4.95&&S18ZDIR020catTwoZDIR060Log>=1.73875) || (log10(WP8CalOnlyEnergy)>=4.95&&log10(WP8CalOnlyEnergy)<5.15&&S18ZDIR020catTwoZDIR060Log>=1.94625) || (log10(WP8CalOnlyEnergy)>=5.15&&log10(WP8CalOnlyEnergy)<5.35&&S18ZDIR020catTwoZDIR060Log>=2.20125) || (log10(WP8CalOnlyEnergy)>=5.35&&log10(WP8CalOnlyEnergy)<5.55&&S18ZDIR020catTwoZDIR060Log>=2.21375) || (log10(WP8CalOnlyEnergy)>=5.55&&log10(WP8CalOnlyEnergy)<5.75&&S18ZDIR020catTwoZDIR060Log>=2.31125)))', 
                  'CalOnlyR30':' ((Cal1MomZDir>=0.2&&Cal1MomZDir<0.6)&&((log10(WP8CalOnlyEnergy)>=4.35&&log10(WP8CalOnlyEnergy)<4.55&&S18ZDIR020catTwoZDIR060Log>=1.97125) || (log10(WP8CalOnlyEnergy)>=4.55&&log10(WP8CalOnlyEnergy)<4.75&&S18ZDIR020catTwoZDIR060Log>=2.11875) || (log10(WP8CalOnlyEnergy)>=4.75&&log10(WP8CalOnlyEnergy)<4.95&&S18ZDIR020catTwoZDIR060Log>=2.26125) || (log10(WP8CalOnlyEnergy)>=4.95&&log10(WP8CalOnlyEnergy)<5.15&&S18ZDIR020catTwoZDIR060Log>=2.52625) || (log10(WP8CalOnlyEnergy)>=5.15&&log10(WP8CalOnlyEnergy)<5.35&&S18ZDIR020catTwoZDIR060Log>=2.64625) || (log10(WP8CalOnlyEnergy)>=5.35&&log10(WP8CalOnlyEnergy)<5.55&&S18ZDIR020catTwoZDIR060Log>=2.45875) || (log10(WP8CalOnlyEnergy)>=5.55&&log10(WP8CalOnlyEnergy)<5.75&&S18ZDIR020catTwoZDIR060Log>=2.78875))) || ((Cal1MomZDir>=0.6)&&((log10(WP8CalOnlyEnergy)>=4.35&&log10(WP8CalOnlyEnergy)<4.55&&S18ZDIR020catTwoZDIR060Log>=1.57375) || (log10(WP8CalOnlyEnergy)>=4.55&&log10(WP8CalOnlyEnergy)<4.75&&S18ZDIR020catTwoZDIR060Log>=2.15625) || (log10(WP8CalOnlyEnergy)>=4.75&&log10(WP8CalOnlyEnergy)<4.95&&S18ZDIR020catTwoZDIR060Log>=2.27625) || (log10(WP8CalOnlyEnergy)>=4.95&&log10(WP8CalOnlyEnergy)<5.15&&S18ZDIR020catTwoZDIR060Log>=2.43625) || (log10(WP8CalOnlyEnergy)>=5.15&&log10(WP8CalOnlyEnergy)<5.35&&S18ZDIR020catTwoZDIR060Log>=2.62625) || (log10(WP8CalOnlyEnergy)>=5.35&&log10(WP8CalOnlyEnergy)<5.55&&S18ZDIR020catTwoZDIR060Log>=2.58125) || (log10(WP8CalOnlyEnergy)>=5.55&&log10(WP8CalOnlyEnergy)<5.75&&S18ZDIR020catTwoZDIR060Log>=2.67625)))', 
                  'CalOnlyR10':'((Cal1MomZDir>=0.2&&Cal1MomZDir<0.6)&&((log10(WP8CalOnlyEnergy)>=4.35&&log10(WP8CalOnlyEnergy)<4.55&&S18ZDIR020catTwoZDIR060Log>=2.53875) || (log10(WP8CalOnlyEnergy)>=4.55&&log10(WP8CalOnlyEnergy)<4.75&&S18ZDIR020catTwoZDIR060Log>=2.54375) || (log10(WP8CalOnlyEnergy)>=4.75&&log10(WP8CalOnlyEnergy)<4.95&&S18ZDIR020catTwoZDIR060Log>=2.59375) || (log10(WP8CalOnlyEnergy)>=4.95&&log10(WP8CalOnlyEnergy)<5.15&&S18ZDIR020catTwoZDIR060Log>=3.05625) || (log10(WP8CalOnlyEnergy)>=5.15&&log10(WP8CalOnlyEnergy)<5.35&&S18ZDIR020catTwoZDIR060Log>=2.95875) || (log10(WP8CalOnlyEnergy)>=5.35&&log10(WP8CalOnlyEnergy)<5.55&&S18ZDIR020catTwoZDIR060Log>=2.75125) || (log10(WP8CalOnlyEnergy)>=5.55&&log10(WP8CalOnlyEnergy)<5.75&&S18ZDIR020catTwoZDIR060Log>=3.58375))) || ((Cal1MomZDir>=0.6)&&((log10(WP8CalOnlyEnergy)>=4.35&&log10(WP8CalOnlyEnergy)<4.55&&S18ZDIR020catTwoZDIR060Log>=2.28125) || (log10(WP8CalOnlyEnergy)>=4.55&&log10(WP8CalOnlyEnergy)<4.75&&S18ZDIR020catTwoZDIR060Log>=2.66625) || (log10(WP8CalOnlyEnergy)>=4.75&&log10(WP8CalOnlyEnergy)<4.95&&S18ZDIR020catTwoZDIR060Log>=2.70875) || (log10(WP8CalOnlyEnergy)>=4.95&&log10(WP8CalOnlyEnergy)<5.15&&S18ZDIR020catTwoZDIR060Log>=2.90625) || (log10(WP8CalOnlyEnergy)>=5.15&&log10(WP8CalOnlyEnergy)<5.35&&S18ZDIR020catTwoZDIR060Log>=2.96375) || (log10(WP8CalOnlyEnergy)>=5.35&&log10(WP8CalOnlyEnergy)<5.55&&S18ZDIR020catTwoZDIR060Log>=2.89625) || (log10(WP8CalOnlyEnergy)>=5.55&&log10(WP8CalOnlyEnergy)<5.75&&S18ZDIR020catTwoZDIR060Log>=2.80875)))' }
liNameClass = dictStrCutBDT.keys()
liCutBDT = dictStrCutBDT.values()


for cl in liNameClass:
    print cl
    h = ROOT.TH2D('hResidual_{0}'.format(cl), 'Residual background counts in {0} [events]'.format(cl), 7, 4.35, 5.75, 10, 0.0, 1.0)
    h.GetXaxis().SetTitle('log_{10}WP8CalOnlyEnergy')
    h.GetYaxis().SetTitle('Cal1MomZDir')
    trBKG.Draw("Cal1MomZDir:log10(WP8CalOnlyEnergy)>>{0}".format(h.GetName()), "(McSourceId!=7000 && Cal1MomZCrossSide840>=0.0 && (TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000)) ) && Cal1RawEnergySum>=20000 && Cal1MomNumIterations>0 && Cal1TransRms>=10 && Cal1TransRms<70 && Acd2Cal1VetoSigmaHit>0 && FswGamState==0 && (GltGemSummary&0x20)==0) && ({0})".format(dictStrCutBDT[cl]), "goff")
    fileOut.cd()
    h.Write()

    hRej = h.Clone('hRejection_{0}'.format(cl))
    hRej.SetTitle('Residual / Raw MC background counts in {0}'.format(cl))
    hRej.Divide(hRaw)
    fileOut.cd()
    hRej.Write()

    fileAcc = ROOT.TFile(dictPathFileAcc[cl], 'READ')
    hAcc = fileAcc.Get('acc_cth_hist')
    hDivAcc = h.Clone('hResidualPerAcc_{0}'.format(cl))
    hDivAcc.SetTitle('Residual backgrounds / Accepatnce')
    hDivAcc.Divide(hAcc)
    fileOut.cd()
    hDivAcc.Write()

    hPro = ROOT.TH2D('hResidualPro_{0}'.format(cl), 'Residual proton counts in {0} [events]'.format(cl), 7, 4.35, 5.75, 10, 0.0, 1.0)
    hPro.GetXaxis().SetTitle('log_{10}WP8CalOnlyEnergy')
    hPro.GetYaxis().SetTitle('Cal1MomZDir')
    trBKG.Draw("Cal1MomZDir:log10(WP8CalOnlyEnergy)>>{0}".format(hPro.GetName()), "(McSourceId==1000 && Cal1MomZCrossSide840>=0.0 && (TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000)) ) && Cal1RawEnergySum>=20000 && Cal1MomNumIterations>0 && Cal1TransRms>=10 && Cal1TransRms<70 && Acd2Cal1VetoSigmaHit>0 && FswGamState==0 && (GltGemSummary&0x20)==0) && ({0})".format(dictStrCutBDT[cl]), "goff")
    fileOut.cd()
    hPro.Write()

    hProDivAcc = h.Clone('hResidualProPerAcc_{0}'.format(cl))
    hProDivAcc.SetTitle('Residual protons / Accepatnce in {0}'.format(cl))
    hProDivAcc.Divide(hAcc)
    fileOut.cd()
    hProDivAcc.Write()

    hEle = ROOT.TH2D('hResidualEle_{0}'.format(cl), 'Residual electron counts in {0} [events]'.format(cl), 7, 4.35, 5.75, 10, 0.0, 1.0)
    hEle.GetXaxis().SetTitle('log_{10}WP8CalOnlyEnergy')
    hEle.GetYaxis().SetTitle('Cal1MomZDir')
    trBKG.Draw("Cal1MomZDir:log10(WP8CalOnlyEnergy)>>{0}".format(hEle.GetName()), "(McSourceId==2000 && Cal1MomZCrossSide840>=0.0 && (TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000)) ) && Cal1RawEnergySum>=20000 && Cal1MomNumIterations>0 && Cal1TransRms>=10 && Cal1TransRms<70 && Acd2Cal1VetoSigmaHit>0 && FswGamState==0 && (GltGemSummary&0x20)==0) && ({0})".format(dictStrCutBDT[cl]), "goff")
    fileOut.cd()
    hEle.Write()

    hEleDivAcc = h.Clone('hResidualElePerAcc_{0}'.format(cl))
    hEleDivAcc.SetTitle('Residual electrons / Accepatnce in {0}'.format(cl))
    hEleDivAcc.Divide(hAcc)
    fileOut.cd()
    hEleDivAcc.Write()

    hRatioPro = hPro.Clone('hRatioPro_{0}'.format(cl))
    hRatioPro.SetTitle('Ratio of protons in redisual background of {0}'.format(cl))
    hRatioPro.Divide(h)
    fileOut.cd()
    hRatioPro.Write()

    hRatioEle = hEle.Clone('hRatioEle_{0}'.format(cl))
    hRatioEle.SetTitle('Ratio of electorons in redisual background of {0}'.format(cl))
    hRatioEle.Divide(h)
    fileOut.cd()
    hRatioEle.Write()

    hRatioNonProEle = hPro.Clone('hRatioNonProEle_{0}'.format(cl))
    hRatioNonProEle.SetTitle('Ratio of non-p/e CRs in redisual background of {0}'.format(cl))
    for ix in range(hRatioNonProEle.GetNbinsX()):
        for iy in range(hRatioNonProEle.GetNbinsY()):
            hRatioNonProEle.SetBinContent(ix+1, iy+1, 1.0-(hRatioPro.GetBinContent(ix+1, iy+1)+hRatioEle.GetBinContent(ix+1, iy+1)))
    fileOut.cd()
    hRatioNonProEle.Write()
