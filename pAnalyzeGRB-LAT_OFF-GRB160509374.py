#!/usr/bin/env python

import sys
import ROOT
from ROOT import TTree
from ROOT import TChain
from ROOT import TH2D
import numpy as np
import yaml
import xml.etree.ElementTree as ET
import datetime
from ctypes import *
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
from pColor import *

ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)

#from pCutBDT import cutBDT
from pAnalysisConfig import *

# ----- Event class setup -----
par = sys.argv
cfg = ClassConfig('Both', [10, 3, 1], 1)
aCutEGB = cfg.aCutEGB
aaStrSelect = cfg.aaStrSelect
aStrSelect = cfg.aStrSelect

#IRF
listPathFilePerf = [['/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_TRANSIENT100_P8R2_TRANSIENT100_perf.root', '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_SOURCE_P8R2_SOURCE_perf.root'], 
                    ['/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R100_perf.root', '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R30_perf.root', '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_CalOnly_R10_perf.root']]
htgPerf = CutPerformanceHtg(listPathFilePerf)

# Target
listTgtGRB = [par[1]]

# Catalogue Data
pathList = "/nfs/farm/g/glast/u/mtakahas/data/catalogue/PublicTableGRBs.xml"
fileList = ET.parse(pathList)
rtXml = fileList.getroot()

# OFF regions
#nOff = 4;
#degOffOffset = 14.0;

print "===================="
# Making all sky map

nFile = (len(par)-5)/2
listFileIn = par[6:6+nFile]
listFileDat = par[6+nFile:6+2*nFile]
print listFileIn

aliasSelections = yaml.load(open('/afs/slac.stanford.edu/u/gl/mtakahas/eventSelect/config/pass8_event_selections.yaml','r'))
chIn = ROOT.TChain('EVENTS')
for fileIn in listFileIn:
    chIn.Add(fileIn)
strSuffixOut = par[4]
aCutPsf = [95, 68]

fileRej = ROOT.TFile(par[5], 'READ')
aHtgRej = []
for icc in range(len(aaStrSelect[1])):
    aHtgRej.append(fileRej.Get('h2Rejection{0}'.format(icc)))

for nameGrb in listTgtGRB:
    fileOut = ROOT.TFile('Plot_GRB{0}_{1}.root'.format(nameGrb, strSuffixOut), 'UPDATE')
    #------ Source data -----
    #indexGrbName = nameFileIn.rindex('GRB') + 3
    #indexGrbNameEnd = indexGrbName + 9
    #nameGrb = nameFileIn[indexGrbName:indexGrbNameEnd]
    for grb in rtXml:
        if grb.findtext("./GRBNAME")==nameGrb: 
            trigger_time = float(grb.findtext("./MET"))
            if grb.findtext("./ERROR") == "--" or grb.findtext("./ERROR") == "":
                if grb.findtext("./LATERROR") == "--" or grb.findtext("./LATERROR") == "":
                    err_rad = 0.
                else:
                    err_rad = float(grb.findtext("./LATERROR"))
            else:
                if grb.findtext("./LATERROR") == "--" or grb.findtext("./LATERROR") == "" or float(grb.findtext("./ERROR"))<=float(grb.findtext("./LATERROR")):
                    err_rad = float(grb.findtext("./ERROR"))                    
                    raSrc = float(grb.findtext("./RA"))
                    decSrc = float(grb.findtext("./DEC"))
                else :
                    err_rad = float(grb.findtext("./LATERROR"))                    
                    raSrc = float(grb.findtext("./LATRA"))
                    decSrc = float(grb.findtext("./LATDEC"))
    print ""
    print "==============="
    print "GRB", nameGrb
    print "==============="
    print "(", raSrc, ",", decSrc, "), Error radius:", err_rad, "Trigger MET:", trigger_time 
    if float(par[3])>0 and float(par[2])<0:
        metStart = trigger_time+float(par[2])
        metStop = trigger_time+float(par[3])
    elif float(par[3])>0 and float(par[2])==0:
        metStart = trigger_time
        metStop = trigger_time+float(par[3])
    elif float(par[2])>0 and float(par[3])>float(par[2]):
        metStart = float(par[2])
        metStop = float(par[3])
    else:
        metStart = chIn.GetMinimum("t")-1
        metStop = chIn.GetMaximum("t")+1

    print 'Analysis time domain: MET', metStart, '-', metStop
    nEvent = chIn.GetEntries('TIME>={0} && TIME<{1}'.format(metStart, metStop))
    print "Total number of events in the time domain:", nEvent
    vecTgt = np.array([cos(radians(decSrc))*cos(radians(raSrc)), cos(radians(decSrc))*sin(radians(raSrc)), sin(radians(decSrc))])
    #distCut = min(7.0, 5.0+err_rad)
    cthTgt = 0.986

    # Background estimation
    print "Estimation of unifomal backgrounds"
    nEvtPrecutPSF68 = [0., 0., 0.]
    nEvtPrecutPSF95 = [0., 0., 0.]
    nEvtPostcutPSF68 = [0., 0., 0.]
    nEvtPostcutPSF95 = [0., 0., 0.]
    dictHtgPrecut = {}
    dictHtgPostcut = {}
    for cpc in aCutPsf:
        dictHtgPrecut[cpc] = ROOT.TH3D('hPrecutPSF{0}'.format(cpc), 'Counted event number before the cut within PSF{0}'.format(cpc), 15, 0, 15, 7, 4.35, 5.75, 180, 0, 180)
        dictHtgPostcut[cpc] = ROOT.TH3D('hPostcutPSF{0}'.format(cpc), 'Expected event number after the cut within PSF{0}'.format(cpc), 15, 0, 15, 7, 4.35, 5.75, 180, 0, 180)
    for pathFileDat in listFileDat:
        fileDat = ROOT.TFile(pathFileDat, 'READ')
        print fileDat.GetName()
        chDat = fileDat.Get('MeritTuple')
        nEvtBll = chDat.GetEntries()
        print chDat.GetName(), "has", nEvtBll, "events."
        for jEvt in range(nEvtBll):
            chDat.GetEntry(jEvt)
            if chDat.EvtElapsedTime>=metStart and chDat.EvtElapsedTime<metStop and chDat.FT1CalZenithTheta<90 and chDat.Cal1MomZDir>=0.2 and chDat.Cal1RawEnergySum>=20000 and chDat.EvtJointLogEnergy>=4.35 and chDat.EvtJointLogEnergy<5.75 and (chDat.TkrNumTracks==0 or (math.log10(max(chDat.CalTrackAngle,1E-4)) > (0.529795)*(chDat.EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((chDat.EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((chDat.EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((chDat.EvtJointLogEnergy-3.000000)/0.916667,3))))*(chDat.EvtJointLogEnergy  >= 3.000000 and chDat.EvtJointLogEnergy <= 5.750000) + (-0.398962)*(chDat.EvtJointLogEnergy >  5.750000))) and chDat.FswGamState == 0:
                vecEvtB = np.array([cos(radians(chDat.FT1CalDec))*cos(radians(chDat.FT1CalRa)), cos(radians(chDat.FT1CalDec))*sin(radians(chDat.FT1CalRa)), sin(radians(chDat.FT1CalDec))])
                radThetaB = acos(np.dot(vecTgt, vecEvtB))
                degDistB = degrees(radThetaB)
                aDictDistCutB = []
                for kcc in range(len(aaStrSelect[1])): # Fill each event to every CLASS bin.
                    aDictDistCutB.append({ 'PSF95': (htgPerf.getPSF95_cth(1, kcc-1, chDat.EvtJointLogEnergy, chDat.Cal1MomZDir) + err_rad), 'PSF68': (htgPerf.getPSF68_cth(1, kcc-1, chDat.EvtJointLogEnergy, chDat.Cal1MomZDir) + err_rad) })
                    if degDistB<aDictDistCutB[kcc]['PSF95']:
                        nEvtPrecutPSF95[kcc] = nEvtPrecutPSF95[kcc] + 1
                        dictHtgPrecut[95].Fill(12+kcc, chDat.EvtJointLogEnergy, chDat.Cal1MomZDir)
                        nEvtPostcutPSF95[kcc] = nEvtPostcutPSF95[kcc]+aHtgRej[kcc].GetBinContent(aHtgRej[kcc].GetXaxis().FindBin(chDat.EvtJointLogEnergy), aHtgRej[kcc].GetYaxis().FindBin(chDat.Cal1MomZDir))
                        dictHtgPostcut[95].Fill( 12+kcc, chDat.EvtJointLogEnergy, chDat.Cal1MomZDir, aHtgRej[kcc].GetBinContent(aHtgRej[kcc].GetXaxis().FindBin(chDat.EvtJointLogEnergy),aHtgRej[kcc].GetYaxis().FindBin(chDat.Cal1MomZDir)) )
                        if degDistB<aDictDistCutB[kcc]['PSF68']:
                            nEvtPrecutPSF68[kcc] = nEvtPrecutPSF68[kcc] + 1
                            dictHtgPrecut[68].Fill(12+kcc, chDat.EvtJointLogEnergy, chDat.Cal1MomZDir)
                            nEvtPostcutPSF68[kcc] = nEvtPostcutPSF68[kcc]+aHtgRej[kcc].GetBinContent(aHtgRej[kcc].GetXaxis().FindBin(chDat.EvtJointLogEnergy), aHtgRej[kcc].GetYaxis().FindBin(chDat.Cal1MomZDir))
                            dictHtgPostcut[68].Fill( 12+kcc, chDat.EvtJointLogEnergy, chDat.Cal1MomZDir, aHtgRej[kcc].GetBinContent(aHtgRej[kcc].GetXaxis().FindBin(chDat.EvtJointLogEnergy),aHtgRej[kcc].GetYaxis().FindBin(chDat.Cal1MomZDir)) )
        for kcc in range(len(aaStrSelect[1])):
            print "Number of events within PSF95 after/before cut", aaStrSelect[1][kcc], ":", nEvtPostcutPSF95[kcc],"/",nEvtPrecutPSF95[kcc]
            print "Number of events within PSF68 after/before cut", aaStrSelect[1][kcc], ":", nEvtPostcutPSF68[kcc],"/",nEvtPrecutPSF68[kcc]
    fileOut.cd()
    for cpc in aCutPsf:
        dictHtgPrecut[cpc].Write()
        dictHtgPostcut[cpc].Write()

    # TTree
    fileOut.cd()
    trGRB = ROOT.TTree("trFriendGRB{0}".format(nameGrb), "Friend TTree for GRB{0}".format(nameGrb))
    cdTimeGRB = c_double()
    trGRB.Branch('TIME_GRB', cdTimeGRB, 'TIME_GRB/D')
    cbFlagPSF68 = c_bool()
    trGRB.Branch('FLAG_PSF68', cbFlagPSF68, 'FLAG_PSF68/O')
    cbFlagPSF95 = c_bool()
    trGRB.Branch('FLAG_PSF95', cbFlagPSF95, 'FLAG_PSF95/O')

    # Plot
    #mgr = ROOT.TMultiGraph("mgr", "Gamma-like events within {0} deg".format(distCut))
    mgr = ROOT.TMultiGraph("mgr", "Gamma-like events around GRB{0}".format(nameGrb))
    greOn = []
    mgrZenith = ROOT.TMultiGraph("mgrZenith", "Zenith angle of events around GRB{0}".format(nameGrb))
    greZenith = []
    for cutPsf in aCutPsf:
        print 'PSF cut:', cutPsf, '%'
        #mgr.append(ROOT.TMultiGraph("mgr", "Gamma-like events within PSF{0} from GRB{1}".format(cutPsf, nameGrb)))
        #mgrZenith.append(ROOT.TMultiGraph("mgrZenithPSF{0}".format(cutPsf), "Zenith angle of events within PSF{0} from GRB{1}".format(cutPsf, nameGrb)))
        greOn.append([])
        greZenith.append([])
        for pC in range(len(aaStrSelect)):
            greOn[-1].append([])
            greZenith[-1].append([])
            for qC in range(len(aaStrSelect[pC])):
                greOn[-1][-1].append(ROOT.TGraphErrors())
                greOn[-1][-1][-1].SetName("greOn_{0}_{1}".format(pC, qC))
                greOn[-1][-1][-1].SetTitle("{0}, PSF{1}%".format(aaStrSelect[pC][qC], cutPsf))
                greOn[-1][-1][-1].SetMarkerStyle(20+int(cutPsf==68))
                greZenith[-1][-1].append(ROOT.TGraphErrors())
                greZenith[-1][-1][-1].SetName("greZenith_{0}_{1}".format(pC, qC))
                greZenith[-1][-1][-1].SetTitle("{0}, PSF{1}%".format(aaStrSelect[pC][qC], cutPsf))
                greZenith[-1][-1][-1].SetMarkerStyle(20+int(cutPsf==68))
                if pC==0:
                    greOn[-1][-1][-1].SetMarkerColor(13-12*qC)
                    greZenith[-1][-1][-1].SetMarkerColor(13-12*qC)
                elif pC==1:
                    greOn[-1][-1][-1].SetMarkerColor(kRed+3*(qC-2))
                    greZenith[-1][-1][-1].SetMarkerColor(kRed+3*(qC-2))
                mgr.Add(greOn[-1][-1][-1])
                mgrZenith.Add(greZenith[-1][-1][-1])  

    timeStart = datetime.datetime.now()
    for iEvent in range(nEvent):
        cbFlagPSF68.value = 0
        cbFlagPSF95.value = 0
        chIn.GetEntry(iEvent)
        ngrb = float(nameGrb)
        grbt = chIn.t - trigger_time
        cdTimeGRB.value = grbt
        vecEvt = np.array([cos(radians(chIn.dec))*cos(radians(chIn.ra)), cos(radians(chIn.dec))*sin(radians(chIn.ra)), sin(radians(chIn.dec))])
        dictDistCut = { 'PSF95': (htgPerf.getPSF95_cth(chIn.c-1, 0*(chIn.s==4 or chIn.s==4096)+1*(chIn.s==128 or chIn.s==8192)+2*(chIn.s==16384), chIn.e, cthTgt) + err_rad), 'PSF68': (htgPerf.getPSF68_cth(chIn.c-1, 0*(chIn.s==4 or chIn.s==4096)+1*(chIn.s==128 or chIn.s==8192)+2*(chIn.s==16384), chIn.e, cthTgt) + err_rad) }
        radTheta = acos(np.dot(vecTgt, vecEvt))
        degDist = degrees(radTheta)
        if degDist<dictDistCut['PSF95'] and chIn.t>=metStart and chIn.t<metStop:
            print ""
            print "== ON photon candidate!!! =="
            cbFlagPSF95.value = 1
            bPSF68=0
            if degDist<dictDistCut['PSF68']:
                bPSF68=1
                cbFlagPSF68.value = 1
            if chIn.c == 1:
                if chIn.s == 4:
                    if chIn.z<90:
                        greOn[bPSF68][0][0].SetPoint(greOn[bPSF68][0][0].GetN(), chIn.t-trigger_time, pow(10, chIn.e-3))
                    greZenith[bPSF68][0][0].SetPoint(greZenith[bPSF68][0][0].GetN(), chIn.t-trigger_time, chIn.z)
                elif chIn.s == 128:
                    if chIn.z<90:
                        greOn[bPSF68][0][1].SetPoint(greOn[bPSF68][0][1].GetN(), chIn.t-trigger_time, pow(10, chIn.e-3))
                    greZenith[bPSF68][0][1].SetPoint(greZenith[bPSF68][0][1].GetN(), chIn.t-trigger_time, chIn.z)
            elif chIn.c == 2:
                if chIn.s == 4096:
                    if chIn.z<90:
                        greOn[bPSF68][1][0].SetPoint(greOn[bPSF68][1][0].GetN(), chIn.t-trigger_time, pow(10, chIn.e-3))
                    greZenith[bPSF68][1][0].SetPoint(greZenith[bPSF68][1][0].GetN(), chIn.t-trigger_time, chIn.z)
                elif chIn.s == 8192:
                    if chIn.z<90:
                        greOn[bPSF68][1][1].SetPoint(greOn[bPSF68][1][1].GetN(), chIn.t-trigger_time, pow(10, chIn.e-3))
                    greZenith[bPSF68][1][1].SetPoint(greZenith[bPSF68][1][1].GetN(), chIn.t-trigger_time, chIn.z)
                elif chIn.s == 16384:
                    if chIn.z<90:
                        greOn[bPSF68][1][2].SetPoint(greOn[bPSF68][1][2].GetN(), chIn.t-trigger_time, pow(10, chIn.e-3))
                    greZenith[bPSF68][1][2].SetPoint(greZenith[bPSF68][1][2].GetN(), chIn.t-trigger_time, chIn.z)

            print "Run ID:", chIn.run
            print "Event ID:", chIn.evid
            print "Event category:", cfg.aStrSelect[chIn.c-1]
            print "Event class:", chIn.s
            print "Time from the trigger:", chIn.t-trigger_time, "s"
            print "Anglular distance:", degDist, "deg"
            print "PSF68:", htgPerf.getPSF68_cth(chIn.c-1, 0*(chIn.s==4 or chIn.s==4096)+1*(chIn.s==128 or chIn.s==8192)+2*(chIn.s==16384), chIn.e, chIn.cth), "deg"
            print "Energy:", pow(10,chIn.e-3), "GeV"
            print "Edisp68:", 100*htgPerf.getEdisp68_cth(chIn.c-1, 0*(chIn.s==4 or chIn.s==4096)+1*(chIn.s==128 or chIn.s==8192)+2*(chIn.s==16384), chIn.e, chIn.cth), "%"
            print "Cos( inclination angle ):", chIn.cth
            print "Zenith angle:", chIn.z, "deg"
        trGRB.Fill()
        if iEvent%(nEvent/20)==0:
            rate = int((iEvent*100.)/nEvent+0.5)
            if rate>0:
                nt = (datetime.datetime.now() - timeStart).seconds * (100.-rate)/rate
                meter = "\r[{0}{1}] Wait {2} hr {3} min".format("=" * rate, ' ' * (100-rate), int(nt/3600), (int(nt)%3600)/60+1)
            else:
                meter = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
            sys.stdout.write(meter)
            sys.stdout.flush()
    cEvent = ROOT.TCanvas("cEvent", "GRB {0} gamma-like events".format(nameGrb))
    cEvent.cd()
    mgr.Draw("AP")
    mgr.GetXaxis().SetTitle("Time [s]")
    mgr.GetYaxis().SetTitle("Energy [GeV]")
    leg = ROOT.TLegend(0.67, 0.5, 0.88, 0.88)
    for pD in range(len(aaStrSelect)):
        for qD in range(len(aaStrSelect[pD])):
            leg.AddEntry(greOn[0][pD][qD], greOn[0][pD][qD].GetTitle(), "p")
            leg.AddEntry(greOn[1][pD][qD], greOn[1][pD][qD].GetTitle(), "p")
#    for rr in range(nOff):
#        leg.AddEntry(greOff[1][0][rr], "OFF{0}".format(rr+1), "p")
    leg.Draw("same")
    cZenith = ROOT.TCanvas("cZenith", "Zenith angle of ON/OFF events")
    cZenith.cd()
    mgrZenith.Draw("AP")
    mgrZenith.GetXaxis().SetTitle("Time [s]")
    mgrZenith.GetYaxis().SetTitle("Zenith angle [deg]")
    legZenith = ROOT.TLegend(0.67, 0.5, 0.88, 0.88)
    for pD in range(len(aaStrSelect)):
        for qD in range(len(aaStrSelect[pD])):
            legZenith.AddEntry(greZenith[0][pD][qD], greZenith[0][pD][qD].GetTitle(), "p")
            legZenith.AddEntry(greZenith[1][pD][qD], greZenith[1][pD][qD].GetTitle(), "p")
#    for rz in range(nOff+1):
#        legZenith.AddEntry(greZenith[rz], greZenith[rz].GetTitle(), "p")

    legZenith.Draw("same")
    print ""
    fileOut.cd()
#    trm.Write()
    cEvent.Write()
    cZenith.Write()
    mgr.Write()
    mgrZenith.Write()
    trGRB.Write()
    print "Finished!"
