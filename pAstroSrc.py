#!/usr/bin/env python

import sys
import ROOT
import yaml
sys.path.append("/Users/Mitsunari/FermiMVA_store")
ROOT.gROOT.SetBatch()
#from myUtils.__ROOT__ import *
#from myUtils.effUtils import *
from array import array
import math

ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)
log = open("allSky%s.log" % ,"w")

aCutEnergy = [4.75, 5.00, 5.25, 5.50]
aStrEnergy = []
for ice in range(len(aCutEnergy)-1):
    aStrEnergy.append("%3.1fGeV<=Energy<%3.1fGeV" % (pow(10, aCutEnergy[ice]-3), pow(10, aCutEnergy[ice+1]-3)))
vStepEnergy = [float((4.50+4.75)/2.0), float((4.75+5.00)/2.0), float((5.00+5.25)/2.0), float((5.25+5.50)/2.0)]
#aCutEnergy = ["EvtJointLogEnergy>=4.50 && EvtJointLogEnergy<4.75", "EvtJointLogEnergy>=4.75 && EvtJointLogEnergy<5.00", "EvtJointLogEnergy>=5.00 && EvtJointLogEnergy<5.25", "EvtJointLogEnergy>=5.25 && EvtJointLogEnergy<5.5"]

#-------------------- BDT cut study
fileRoc = ROOT.TFile('v20r09p09_G1haB1_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_1_S10WOrw_1Log_roc.root', "READ")
h2Sig = fileRoc.Get("sig_acc")
h2Bkg = fileRoc.Get("bkg_rate")
h2Count = fileRoc.Get("bkg_counts_cum_hist")
h2Egb = fileRoc.Get("egb_rate")

aCutEgb = [10, 2, 1]
aaValCutBDT = []

for ie in range(len(aCutEnergy)-1):
    print >> log , '=========', aCutEnergy[ie],' - ', aCutEnergy[ie+1], '========='
    aaValCutBDT.append([])
    abFound = []
    aNumGamCut = []
    aNumRemCut = []
    aNumCountCut = []
    for jc in range(len(aCutEgb)):
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
        for ic in range(len(aCutEgb)):
            if (abFound[ic]==False and nRem<=aCutEgb[ic]*nEgb):
                abFound[ic] = True
                aaValCutBDT[ie].append(h1Sig.GetBinCenter(ib))
                aNumGamCut.append(nSig)
                aNumRemCut.append(nRem)
                aNumCountCut.append(nCount)
                print >> log , "--------- Background level:", aCutEgb[ic], "x EGB ---------"
                print >> log , "Cut value:", aaValCutBDT[ie][-1]
                print >> log , "Acceptance:", aNumGamCut[-1], "[m^2 sr]"
                print >> log , "Background rate:", aNumRemCut[-1], "[MeV sr^-1 s^-1]"
                print >> log , "Background count:", aNumCountCut[-1], "[events]"

#-------------------- ON/OFF Definition
raSrc = 166.079 # in degree
decSrc = 38.1947 # in degree

vOnMaxCalTkr = 1.0
vOffMinCalTkr = 5.0
vOffMaxCalTkr = 10.0

vOnMaxCalOnly = 5.25
vOffMinCalOnly = 10.0
vOffMaxCalOnly = 20.0

saOnCalTkr = 2.0 * ROOT.TMath.Pi() * ( ROOT.TMath.Cos(0.0/180.*ROOT.TMath.Pi()) - ROOT.TMath.Cos(vOnMaxCalTkr/180.*ROOT.TMath.Pi()) )
saOffCalTkr = 2.0 * ROOT.TMath.Pi() * ( ROOT.TMath.Cos(vOffMinCalTkr/180.*ROOT.TMath.Pi()) - ROOT.TMath.Cos(vOffMaxCalTkr/180.*ROOT.TMath.Pi()) )
saOnCalOnly = 2.0 * ROOT.TMath.Pi() * ( ROOT.TMath.Cos(0.0/180.*ROOT.TMath.Pi()) - ROOT.TMath.Cos(vOnMaxCalOnly/180.*ROOT.TMath.Pi()) )
saOffCalOnly = 2.0 * ROOT.TMath.Pi() * ( ROOT.TMath.Cos(vOffMinCalOnly/180.*ROOT.TMath.Pi()) - ROOT.TMath.Cos(vOffMaxCalOnly/180.*ROOT.TMath.Pi()) )

#-------------------- Add Friends
#nameFileData = ["/Users/Mitsunari/FermiMVA_store/EarthLimbTracks20141125.root"]
#nameFileFriend = []

#-------------------- Load dataS10V200909_020rawe30zcs000nbep006WWOtrkWbkWOmczWOrw_15_BDTG1000D06S10V200909_020rawe30zcs000nbep006WWOtrkWbkWOmczWOrw_15_BDTG1000D06
chainData = ROOT.TChain("MeritTuple")
chainData.Add("/Users/Mitsunari/FermiMVA_store/Mrk421_2014Dec_2009.root")
chainData.Add("/Users/Mitsunari/FermiMVA_store/Mrk421_2014Dec_2010.root")
chainData.Add("/Users/Mitsunari/FermiMVA_store/Mrk421_2014Dec_2011.root")
chainData.Add("/Users/Mitsunari/FermiMVA_store/Mrk421_2014Dec_2012.root")
chainData.Add("/Users/Mitsunari/FermiMVA_store/Mrk421_2014Dec_2013.root")
#chainData.Add("")

#-------------------- Aliases
nameVarBDT = 'S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_BDTG1000D06'
aStrSelect = [["CalTkrFilter", "P8R1_TRANSIENT_R100", "P8R1_SOURCE"], ["CalOnlyFilter", "CalOnly_10xEGB", "CalOnly_2xEGB", "CalOnly_1xEGB"]]

#strCutCalOnly = "FswGamState==0 && (TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy  >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000)))"
#strCutSkim = "Cal1RawEnergySum>30000 && Cal1MomNumIterations>0 && Cal1TransRms>10 && Cal1TransRms<70 && Acd2Cal1VetoSigmaHit>0 && Cal1MomZDir>0.2 && (1.0+WP8CalOnlyBEP_E_BDT)/2.0>0.06"
aliasSelections = yaml.load(open('/Users/Mitsunari/FermiMVA/python/pass8_event_selections.yaml','r'))
for k,v in aliasSelections.iteritems(): 
    chainData.SetAlias(k,v)
#    print >> log , k, chainData.GetAlias(k)
chainData.SetAlias('CalOnlyFilter', 'FswGamFilter&&(CalTrackAngle>0.3||!TracksFilter)')
chainData.SetAlias('SkimCut', '(TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy  >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000))) && FswGamState==0')
chainData.SetAlias('GamProb', '-log10(1.0-(1.0+%s)/2.0)' % "S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_BDTG1000D06")

chainData.SetAlias('CutOnCalOnly', "(acos(sin(FT1CalDec)*sin(decSrc)+cos(FT1CalDec)*cos(decSrc)*cos(FT1CalRa-raSrc))<%f)" % vOnMaxCalOnly)
chainData.SetAlias('CutOffCalOnly', "(acos(sin(FT1CalDec)*sin(decSrc)+cos(FT1CalDec)*cos(decSrc)*cos(FT1CalRa-raSrc))>=%f) && (acos(sin(FT1CalDec)*sin(decSrc)+cos(FT1CalDec)*cos(decSrc)*cos(FT1CalRa-raSrc))<%f)" % (vOffMinCalOnly, vOffMaxCalOnly))
chainData.SetAlias('CutOnCalTkr', "(acos(sin(FT1Dec)*sin(decSrc)+cos(FT1Dec)*cos(decSrc)*cos(FT1Ra-raSrc))<%f)" % vOnMaxCalTkr)
chainData.SetAlias('CutOffCalTkr', "((acos(sin(FT1Dec)*sin(decSrc)+cos(FT1Dec)*cos(decSrc)*cos(FT1Ra-raSrc)))>=%f && (acos(sin(FT1Dec)*sin(decSrc)+cos(FT1Dec)*cos(decSrc)*cos(FT1Ra-raSrc)))<%f)" % (vOffMinCalTkr, vOffMaxCalTkr))

color = [13, 1, 2, 13, 1, 4, 6]
sFill = [0, 3003, 3003, 0, 0, 0, 0]
wLine = [1, 1, 1, 1, 2, 3, 2]
sLine = [3, 2, 2, 3, 1, 1, 1]

#-------------------- 
nameSrc = "Mrk421"
fileOut = ROOT.TFile("%s_compe0.root" % nameSrc, "RECREATE")
cTheta = ROOT.TCanvas("cTheta", "FT1ZenithTheta", 1200, 800)
cTheta.Divide(2,2)

aHtgDummy = []
aHtgZenith = []
# aHtgMap_CalTkrFilter = []
# aHtgMap_TRANSIENT = []
# aHtgMap_SOURCE = []
# aHtgMap_CalOnlyFilter = []
# aHtgMap_CalOnly10 = []
# aHtgMap_CalOnly2 = []
# aHtgMap_CalOnly1 = []
aaHtgMap = [[[], [], []], [[], [], [], []]]
aaNumOn = [[[], [], []], [[], [], [], []]]
aaNumOff = [[[], [], []], [[], [], [], []]]
aaSB = [[[], [], []], [[], [], [], []]]
aaRej = [[[], [], []], [[], [], [], []]]
aaNumSig = [[[], [], []], [[], [], [], []]]
#aaHtgMap = [aHtgMap_CalTkrFilter, aHtgMap_TRANSIENT, aHtgMap_SOURCE, aHtgMap_CalOnlyFilter, aHtgMap_CalOnly10, aHtgMap_CalOnly2, aHtgMap_CalOnly1]
# aSB_CalTkrFilter = []
# aSB_TRANSIENT = []
# aSB_SOURCE = []
# aSB_CalOnlyFilter = []
# aSB_CalOnly10 = []
# aSB_CalOnly2 = []
# aSB_CalOnly1 = []
# aaSB = [aSB_CalTkrFilter, aSB_TRANSIENT, aSB_SOURCE, aSB_CalOnlyFilter, aSB_CalOnly10, aSB_CalOnly2, aSB_CalOnly1]
# aRej_CalTkrFilter = []
# aRej_TRANSIENT = []
# aRej_SOURCE = []
# aRej_CalOnlyFilter = []
# aRej_CalOnly10 = []
# aRej_CalOnly2 = []
# aRej_CalOnly1 = []
# aaRej = [aRej_CalTkrFilter, aRej_TRANSIENT, aRej_SOURCE, aRej_CalOnlyFilter, aRej_CalOnly10, aRej_CalOnly2, aRej_CalOnly1]
# aNumSig_CalTkrFilter = []
# aNumSig_TRANSIENT = []
# aNumSig_SOURCE = []
# aNumSig_CalOnlyFilter = []
# aNumSig_CalOnly10 = []
# aNumSig_CalOnly2 = []
# aNumSig_CalOnly1 = []
# aaNumSig = [aNumSig_CalTkrFilter, aNumSig_TRANSIENT, aNumSig_SOURCE, aNumSig_CalOnlyFilter, aNumSig_CalOnly10, aNumSig_CalOnly2, aNumSig_CalOnly1]
legTheta = ROOT.TLegend(0.55, 0.4, 0.9, 0.6, "Class")
legTheta.SetTextSize(0.03)
legMap = ROOT.TLegend(0.55, 0.4, 0.9, 0.6, "Class")
legMap.SetTextSize(0.03)

aCanMap = []

print >> log , "##################"
print >> log , "Going to analysis."
print >> log , "##################"

for bE in range(len(aCutEnergy)-1):
    print >> log , "==========", aStrEnergy[bE], "=========="
    aCanMap.append(ROOT.TCanvas("cMap%s" % bE, "Map of %s in %s" % (nameSrc, aStrEnergy[bE]), 1200, 600))
    aCanMap[-1].Divide(4, 2)
    aCanMap[-1].cd(4)
    tex = ROOT.TLatex(0.2, 0.5, "%s" % aStrEnergy[bE])
    tex.Draw()
    for iS in range(len(aStrSelect)):
        for jS in range(len(aStrSelect[iS])):
            print >> log , aStrSelect[iS][jS]
            aaHtgMap[iS][jS].append(ROOT.TH2F("aHtgMap%s_%s_%s" % (bE, iS, jS), "%s;RA [deg];DEC[deg]" % aStrSelect[iS][jS], 40, raSrc-20, raSrc+20, 40, decSrc-20, decSrc+20))
            aaHtgMap[iS][jS][-1].SetMarkerStyle(7)
            aaNumOn[iS][jS].append(0)
            aaNumOff[iS][jS].append(0)
            aaNumSig[iS][jS].append(0)
            aaSB[iS][jS].append(0)
            aaRej[iS][jS].append(-999)
            if iS==0 or jS==0:
                strCutRslt = aStrSelect[iS][jS] + " && EvtJointLogEnergy>=%s && EvtJointLogEnergy<%s" % (aCutEnergy[bE],aCutEnergy[bE+1]) + " && FT1ZenithTheta<100"
                print >> log , iS, jS, strCutRslt
                aCanMap[-1].cd(iS*4+jS+1)
                chainData.Draw('FT1Dec:FT1Ra>>%s' % aaHtgMap[iS][jS][-1].GetName(), '%s' % strCutRslt, 'colz')
            else:
                strCutRslt =  "SkimCut==1 && GamProb>%s" % aaValCutBDT[bE][jS-1] + " && EvtJointLogEnergy>=%s && EvtJointLogEnergy<%s" % (aCutEnergy[bE],aCutEnergy[bE+1]) + " && FT1CalZenithTheta<100"
                print >> log , iS, jS, strCutRslt
                aCanMap[-1].cd(iS*4+jS+1)
                chainData.Draw('FT1CalDec:FT1CalRa>>%s' % aaHtgMap[iS][jS][-1].GetName(), '%s' % strCutRslt, 'colz')
                fileOut.cd()
                aaHtgMap[iS][jS][-1].Write()


        
#     print >> log , "Filling events."
#     for iEve in range(chainData.GetEntries()):
#         chainData.GetEntry(iEve)
#         if iEve%10000==0:
#             sys.stdout.write(".")
#         if chainData.EvtJointLogEnergy>=aCutEnergy[bE] and chainData.EvtJointLogEnergy<aCutEnergy[bE+1]: # Energy cut
#             if chainData.CalTkrFilter: # Regular classes
#                 aaHtgMap[0][0][-1].Fill(chainData.FT1Ra, chainData.FT1Dec)
#                 if chainData.CutOnCalTkr:
#                     aaNumOn[0][0][-1]+=1
#                 elif chainData.CutOffCalTkr:
#                     aaNumOff[0][0][-1]+=1
#                 if chainData.P8R1_TRANSIENT_R100: # TRANSIENT
#                     aaHtgMap[0][1][-1].Fill(chainData.FT1Ra, chainData.FT1Dec)
#                     if chainData.CutOnCalTkr:
#                         aaNumOn[0][1][-1]+=1
#                     elif chainData.CutOffCalTkr:
#                         aaNumOff[0][1][-1]+=1
#                     if chainData.P8R1_SOURCE: # SOURCE
#                         aaHtgMap[0][2][-1].Fill(chainData.FT1Ra, chainData.FT1Dec)
#                         if chainData.CutOnCalTkr:
#                             aaNumOn[0][2][-1]+=1
#                         elif chainData.CutOffCalTkr:
#                             aaNumOff[0][2][-1]+=1
#             elif chainData.CalOnlyFilter: # CalOnly classes
#                 aaHtgMap[1][0][-1].Fill(chainData.FT1CalRa, chainData.FT1CalDec)
#                 if chainData.CutOnCalOnly:
#                     aaNumOn[1][0][-1]+=1
#                 elif chainData.CutOffCalOnly:
#                     aaNumOff[1][0][-1]+=1
#                 if chainData.SkimCut: # Skim cuts
#                     if chainData.GamProb>aaValCutBDT[bE][0]: # 10xEGB level cut
#                         aaHtgMap[1][1][-1].Fill(chainData.FT1CalRa, chainData.FT1CalDec)
#                         if chainData.CutOnCalOnly:
#                             aaNumOn[1][1][-1]+=1
#                         elif chainData.CutOffCalOnly:
#                             aaNumOff[1][1][-1]+=1
#                         if chainData.GamProb>aaValCutBDT[bE][1]: # 2xEGB level cut
#                             aaHtgMap[1][2][-1].Fill(chainData.FT1CalRa, chainData.FT1CalDec)
#                             if chainData.CutOnCalOnly:
#                                 aaNumOn[1][2][-1]+=1
#                             elif chainData.CutOffCalOnly:
#                                 aaNumOff[1][2][-1]+=1
#                             if chainData.GamProb>aaValCutBDT[bE][2]: # 1xEGB level cut
#                                 aaHtgMap[1][3][-1].Fill(chainData.FT1CalRa, chainData.FT1CalDec)
#                                 if chainData.CutOnCalOnly:
#                                     aaNumOn[1][3][-1]+=1
#                                 elif chainData.CutOffCalOnly:
#                                     aaNumOff[1][3][-1]+=1
#     print >> log , ""
        
#     for iCan in range(2):
#         if iCan==0:
#             print >> log , "========= Regular classes ========"
#             saOn = saOnCalTkr
#             saOff = saOffCalTkr
#         else:
#             print >> log , "========= CalOnly classes ========"
#             saOn = saOnCalOnly
#             saOff = saOffCalOnly
#         print >> log , "Solid angle of ON region:", saOn
#         print >> log , "Solid angle of OFF region:", saOff
#         for jCan in range(4):
#             aCanMap[-1].cd(iCan*4+jCan+1)
#             aCanMap[-1].cd(iCan*4+jCan+1).SetLogz()
#             if iCan!=0 or jCan!=3:
#                 print >> log , "---------", aStrSelect[iCan][jCan], "----------"
#                 aaHtgMap[iCan][jCan][-1].Draw("colz")
#                 fileOut.cd()
#                 aaHtgMap[iCan][jCan][-1].Write()
#                 print >> log , "Number of ON events:", aaNumOn[iCan][jCan][-1]
#                 print >> log , "Number of OFF events:", aaNumOff[iCan][jCan][-1]
#                 if aaNumOff[iCan][jCan][-1]!=0:
#                     aaSB[iCan][jCan].append((aaNumOn[iCan][jCan][-1]-aaNumOff[iCan][jCan][-1]*saOn/saOff)/aaNumOn[iCan][jCan][-1]*saOn/saOff)
#                     print >> log , "S/B ratio:", aaSB[iCan][jCan][-1]
#                 aaNumSig[iCan][jCan].append(aaNumOn[iCan][jCan][-1]-aaNumOff[iCan][jCan][-1]*saOn/saOff)
#                 print >> log , "Number of signal events:", aaNumSig[iCan][jCan][-1]
#                 if aaNumOff[iCan][jCan][0]!=0:
#                     aaRej[iCan][jCan].append(aaNumOff[iCan][jCan][-1]/aaNumOff[iCan][jCan][0])
#                 else:
#                     print >> log , "Number of raw background is abnormally zero."
#                 print >> log , "Remaning background fraction:", aaRej[iCan][jCan][-1]

    fileOut.cd()
    aCanMap[-1].Write()

#    cTheta.cd(bE+1)
#    cTheta.cd(bE+1).SetLogy()
#    aHtgDummy.append(ROOT.TH1F("aHtgDummy%s" % bE, "%s;Zenith theta [deg];[events]" % aCutEnergy[bE], 120, 80, 140))
#    aHtgDummy[-1].GetYaxis().SetRangeUser(0.5, 1e4)
#    aHtgDummy[-1].Draw()


#   cTheta.cd(bE+1)
#   legTheta.Draw("same")

# fileOut.cd()
# cTheta.Write()

# cSB = ROOT.TCanvas("cSB", "(S-B)/B ratio", 800, 600)
# cNumSig = ROOT.TCanvas("cNumSig", "Number of signal events", 800, 600)
# #cRej = ROOT.TCanvas("cRej", "Background rejection", 800, 600)
# aGrSB = []
# aGrNumSig = []
# #aGrRej = []
# hDummySB = ROOT.TH1F("hDummySB", "(S-B)/B ratio;log10(Energy[MeV])", 4, 4.5, 5.5)
# hDummySB.GetYaxis().SetRangeUser(0, 100)
# cSB.cd()
# hDummySB.Draw()
# hDummyNumSig = ROOT.TH1F("hDummyNumSig", "Number of signal events;log10(Energy[MeV])", 4, 4.5, 5.5)
# hDummyNumSig.GetYaxis().SetRangeUser(0, 5000)
# cNumSig.cd()
# hDummyNumSig.Draw()
# #hDummyRej = ROOT.TH1F("hDummyRej", "Background rejection", 4, 4.5, 5.5)
# #cRej.cd()
# #hDummyRej.Draw()
# legSB = ROOT.TLegend(0.7, 0.3, 0.9, 0.6, "Earthlimb")
# legSB.SetTextSize(0.03)
# for jS in range(len(aStrSelect)):
#     aGrSB.append(ROOT.TGraph(len(vStepEnergy), array('f', vStepEnergy), array('f', aaSB[jS])))
#     aGrSB[-1].SetName(aStrSelect[jS])
#     aGrSB[-1].SetTitle(aStrSelect[jS])
#     aGrSB[-1].SetLineWidth(2)
#     if jS in {0, 2}:
#         aGrSB[-1].SetLineStyle(3)
#     elif jS in {3}:
#         aGrSB[-1].SetLineStyle(2)
#     else:
#         aGrSB[-1].SetLineStyle(1)
#     aGrSB[-1].SetLineColor((jS>=2)+1)
#     aGrSB[-1].SetMarkerColor((jS>=2)+1)
#     aGrSB[-1].SetMarkerStyle(20)
#     cSB.cd()
#     aGrSB[-1].Draw("LP sames")
#     legSB.AddEntry(aGrSB[-1], aStrSelect[jS], "lp")
#     cSB.cd()
#     legSB.Draw("same")

#     aGrNumSig.append(ROOT.TGraph(len(vStepEnergy), array('f', vStepEnergy), array('f', aaNumSig[jS])))
#     aGrNumSig[-1].SetName(aStrSelect[jS])
#     aGrNumSig[-1].SetTitle(aStrSelect[jS])
#     aGrNumSig[-1].SetLineWidth(2)
#     if jS in {0, 2}:
#         aGrNumSig[-1].SetLineStyle(3)
#     elif jS in {3}:
#         aGrNumSig[-1].SetLineStyle(2)
#     else:
#         aGrNumSig[-1].SetLineStyle(1)
#     aGrNumSig[-1].SetLineColor((jS>=2)+1)
#     aGrNumSig[-1].SetMarkerColor((jS>=2)+1)
#     aGrNumSig[-1].SetMarkerStyle(20)
#     aGrNumSig[-1].SetMarkerStyle(20+5*(jS<2))
#     cNumSig.cd()
#     aGrNumSig[-1].Draw("LP sames")
#     cNumSig.cd()
#     legSB.Draw("same")
#     aGrRej.append(ROOT.TGraph(len(vStepEnergy), vStepEnergy), aRej[jS])
#     aGrRej.append(ROOT.TGraph(len(vStepEnergy), vStepEnergy), aaRej[jS])
#     aGrRej[-1].SetName(aStrSelect[jS])
#     aGrRej[-1].SetTitle(aStrSelect[jS])
#     aGrRej[-1].SetLineWidth(2)
#     aGrRej[-1].SetLineColor(iS+1)
#     cRej.cd()
#    aGrRej[-1].Draw("LP sames")

fileOut.cd()
#cSB.Write()
#cNumSig.Write()
#cRej.Write()
  


