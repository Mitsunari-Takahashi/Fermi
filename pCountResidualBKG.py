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

pathFileBKG = '/u/gl/mtakahas/work/data/MC/BKG200909_62MCE2e4_Combined_DownlikeCalOnlyZD020CalZC840_CR.root'  #'/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_Combined.root'
fileBKG = ROOT.TFile(pathFileBKG, 'READ')
trBKG = fileBKG.Get('MeritTuple')
print trBKG.GetName(), 'is found.'

pathFriendBDT = '/u/gl/mtakahas/work/data/MC/BKG200909_62MCE2e4_Combined_DownlikeCalOnlyZD020CalZC840_CR_S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR020060.root'
#'/u/gl/mtakahas/work/data/MC/BKG200909_62MCE2e4_Combined_DownlikeCalOnlyZD020CalZC840_CR_S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR060100.root' #'/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_Combined_S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060.root'
nameBDT = 'S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR020060_BDTG500D06' #'S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060'
pathFileOut = '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S20/S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir/{0}_Residual.root'.format(nameBDT) #'/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/{0}_Residual.root'.format(nameBDT)
fileOut = ROOT.TFile(pathFileOut, 'UPDATE')

fileRaw = ROOT.TFile('/nfs/farm/g/glast/u/mtakahas/data/MC/hRawBKG.root', 'READ')
hRaw = fileRaw.Get('hRawBKG')

dictPathFileAcc = { 'CalOnly_R100': '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S20/S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir/S20UDZDLog_ZDIR020100_E7bins/S20UDZDLog_ZDIR020100_E7bins_CalOnly_R100_perf.root', #'/u/gl/mtakahas/work/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28binx_Cth40bins_CalOnly_R100_perf.root', 
                    'CalOnly_R030': '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S20/S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir/S20UDZDLog_ZDIR020100_E7bins/S20UDZDLog_ZDIR020100_E7bins_CalOnly_R030_perf.root', #'/u/gl/mtakahas/work/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28binx_Cth40bins_CalOnly_R030_perf.root', 
                    'CalOnly_R010': '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S20/S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir/S20UDZDLog_ZDIR020100_E7bins/S20UDZDLog_ZDIR020100_E7bins_CalOnly_R010_perf.root', #'/u/gl/mtakahas/work/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28binx_Cth40bins_CalOnly_R010_perf.root',
                    'CalOnly_R003': '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S20/S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir/S20UDZDLog_ZDIR020100_E7bins/S20UDZDLog_ZDIR020100_E7bins_CalOnly_R003_perf.root'} #'/u/gl/mtakahas/work/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28binx_Cth40bins_CalOnly_R003_perf.root' }

trBKG.AddFriend("{0} = MeritTuple".format(nameBDT), pathFriendBDT)
trBKG.SetAlias("S20UDZDLog_ZDIR020060", "-log10(1.0-(1.0+{0})/2.0)".format(nameBDT))
dictStrCutBDT = { 'CalOnly_R100':'(log10(WP8CalOnlyEnergy)>=4.35&&log10(WP8CalOnlyEnergy)<4.55&&S20UDZDLog_ZDIR020060>=2.11125) || (log10(WP8CalOnlyEnergy)>=4.55&&log10(WP8CalOnlyEnergy)<4.75&&S20UDZDLog_ZDIR020060>=2.10125) || (log10(WP8CalOnlyEnergy)>=4.75&&log10(WP8CalOnlyEnergy)<4.95&&S20UDZDLog_ZDIR020060>=2.12875) || (log10(WP8CalOnlyEnergy)>=4.95&&log10(WP8CalOnlyEnergy)<5.15&&S20UDZDLog_ZDIR020060>=2.23625) || (log10(WP8CalOnlyEnergy)>=5.15&&log10(WP8CalOnlyEnergy)<5.35&&S20UDZDLog_ZDIR020060>=2.26125) || (log10(WP8CalOnlyEnergy)>=5.35&&log10(WP8CalOnlyEnergy)<5.55&&S20UDZDLog_ZDIR020060>=2.32125) || (log10(WP8CalOnlyEnergy)>=5.55&&log10(WP8CalOnlyEnergy)<5.75&&S20UDZDLog_ZDIR020060>=2.51125)', 
                  'CalOnly_R030':'(log10(WP8CalOnlyEnergy)>=4.35&&log10(WP8CalOnlyEnergy)<4.55&&S20UDZDLog_ZDIR020060>=2.67125) || (log10(WP8CalOnlyEnergy)>=4.55&&log10(WP8CalOnlyEnergy)<4.75&&S20UDZDLog_ZDIR020060>=2.51125) || (log10(WP8CalOnlyEnergy)>=4.75&&log10(WP8CalOnlyEnergy)<4.95&&S20UDZDLog_ZDIR020060>=2.64375) || (log10(WP8CalOnlyEnergy)>=4.95&&log10(WP8CalOnlyEnergy)<5.15&&S20UDZDLog_ZDIR020060>=2.52625) || (log10(WP8CalOnlyEnergy)>=5.15&&log10(WP8CalOnlyEnergy)<5.35&&S20UDZDLog_ZDIR020060>=2.42625) || (log10(WP8CalOnlyEnergy)>=5.35&&log10(WP8CalOnlyEnergy)<5.55&&S20UDZDLog_ZDIR020060>=2.47625) || (log10(WP8CalOnlyEnergy)>=5.55&&log10(WP8CalOnlyEnergy)<5.75&&S20UDZDLog_ZDIR020060>=2.95125)', 
                  'CalOnly_R010':'(log10(WP8CalOnlyEnergy)>=4.35&&log10(WP8CalOnlyEnergy)<4.55&&S20UDZDLog_ZDIR020060>=2.85375) || (log10(WP8CalOnlyEnergy)>=4.55&&log10(WP8CalOnlyEnergy)<4.75&&S20UDZDLog_ZDIR020060>=2.74125) || (log10(WP8CalOnlyEnergy)>=4.75&&log10(WP8CalOnlyEnergy)<4.95&&S20UDZDLog_ZDIR020060>=2.73125) || (log10(WP8CalOnlyEnergy)>=4.95&&log10(WP8CalOnlyEnergy)<5.15&&S20UDZDLog_ZDIR020060>=2.67375) || (log10(WP8CalOnlyEnergy)>=5.15&&log10(WP8CalOnlyEnergy)<5.35&&S20UDZDLog_ZDIR020060>=2.70625) || (log10(WP8CalOnlyEnergy)>=5.35&&log10(WP8CalOnlyEnergy)<5.55&&S20UDZDLog_ZDIR020060>=2.65375) || (log10(WP8CalOnlyEnergy)>=5.55&&log10(WP8CalOnlyEnergy)<5.75&&S20UDZDLog_ZDIR020060>=3.22125)',
                  'CalOnly_R003': '(log10(WP8CalOnlyEnergy)>=4.35&&log10(WP8CalOnlyEnergy)<4.55&&S20UDZDLog_ZDIR020060>=2.98375) || (log10(WP8CalOnlyEnergy)>=4.55&&log10(WP8CalOnlyEnergy)<4.75&&S20UDZDLog_ZDIR020060>=2.96875) || (log10(WP8CalOnlyEnergy)>=4.75&&log10(WP8CalOnlyEnergy)<4.95&&S20UDZDLog_ZDIR020060>=2.78875) || (log10(WP8CalOnlyEnergy)>=4.95&&log10(WP8CalOnlyEnergy)<5.15&&S20UDZDLog_ZDIR020060>=2.89125) || (log10(WP8CalOnlyEnergy)>=5.15&&log10(WP8CalOnlyEnergy)<5.35&&S20UDZDLog_ZDIR020060>=3.04625) || (log10(WP8CalOnlyEnergy)>=5.35&&log10(WP8CalOnlyEnergy)<5.55&&S20UDZDLog_ZDIR020060>=2.67375) || (log10(WP8CalOnlyEnergy)>=5.55&&log10(WP8CalOnlyEnergy)<5.75&&S20UDZDLog_ZDIR020060>=3.22125)'}
liNameClass = dictStrCutBDT.keys()
liCutBDT = dictStrCutBDT.values()


for cl in liNameClass:
    print cl
    h = ROOT.TH2D('hResidual_{0}'.format(cl), 'Residual background counts in {0} [events]'.format(cl), 7, 4.35, 5.75, 10, 0.0, 1.0)
    h.GetXaxis().SetTitle('log_{10}WP8CalOnlyEnergy')
    h.GetYaxis().SetTitle('Cal1MomZDir')
    trBKG.Draw("Cal1MomZDir:log10(WP8CalOnlyEnergy)>>{0}".format(h.GetName()), "Cal1MomZDir<0.6 && ({0})".format(dictStrCutBDT[cl]), "goff")
#    trBKG.Draw("Cal1MomZDir:log10(WP8CalOnlyEnergy)>>{0}".format(h.GetName()), "(McSourceId!=7000 && Cal1MomZCrossSide840>=0.0 && (TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000)) ) && Cal1RawEnergySum>=20000 && Cal1MomNumIterations>0 && Cal1TransRms>=10 && Cal1TransRms<70 && Acd2Cal1VetoSigmaHit>0 && FswGamState==0 && (GltGemSummary&0x20)==0) && ({0})".format(dictStrCutBDT[cl]), "goff")
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
    trBKG.Draw("Cal1MomZDir:log10(WP8CalOnlyEnergy)>>{0}".format(hPro.GetName()), "(McSourceId==1000) && Cal1MomZDir<0.6 && ({0})".format(dictStrCutBDT[cl]), "goff")
#    trBKG.Draw("Cal1MomZDir:log10(WP8CalOnlyEnergy)>>{0}".format(hPro.GetName()), "(McSourceId==1000 && Cal1MomZCrossSide840>=0.0 && (TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000)) ) && Cal1RawEnergySum>=20000 && Cal1MomNumIterations>0 && Cal1TransRms>=10 && Cal1TransRms<70 && Acd2Cal1VetoSigmaHit>0 && FswGamState==0 && (GltGemSummary&0x20)==0) && ({0})".format(dictStrCutBDT[cl]), "goff")
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
#    trBKG.Draw("Cal1MomZDir:log10(WP8CalOnlyEnergy)>>{0}".format(hEle.GetName()), "(McSourceId==2000 && Cal1MomZCrossSide840>=0.0 && (TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000)) ) && Cal1RawEnergySum>=20000 && Cal1MomNumIterations>0 && Cal1TransRms>=10 && Cal1TransRms<70 && Acd2Cal1VetoSigmaHit>0 && FswGamState==0 && (GltGemSummary&0x20)==0) && ({0})".format(dictStrCutBDT[cl]), "goff")
    trBKG.Draw("Cal1MomZDir:log10(WP8CalOnlyEnergy)>>{0}".format(hEle.GetName()), "(McSourceId==2000) && Cal1MomZDir<0.6 && ({0})".format(dictStrCutBDT[cl]), "goff")
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
