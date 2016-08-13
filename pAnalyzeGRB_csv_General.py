#!/usr/bin/env python

import sys
import ROOT
ROOT.gROOT.SetBatch()
from ROOT import TTree, TChain
import numpy as np
import yaml
import datetime
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees
import pandas
from pMETandMJD import *
from pColor import *

ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)

from pAnalysisConfig import *

# ----- Event class setup -----
par = sys.argv
cfg = ClassConfig('Both', [10, 3, 1], 1)
aCutEGB = cfg.aCutEGB
aaStrSelect = cfg.aaStrSelect
nStartBin = cfg.nStartBin

nameFileRoc = "/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZCS000wwoTRKwoMCZDIR00woRW_15_S11D_catTwoZDIR050Log_roc.root" #par[2]
nameVarBDT = "S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_BDTG1000D06"
nameFileSuffix = par[1]
#vSepCut = int(par[2])
cutMVA = CutBDT(nameFileRoc, aCutEGB)
aaValCutBDT = cutMVA.aaValCutBDT[0:]
print aaValCutBDT
nEnergyBin = cutMVA.aaValCutBDT[0]['numBin'] - nStartBin
vEnergyBinWidth = cutMVA.aaValCutBDT[0]['widthBin']
vEnergyLow = cutMVA.aaValCutBDT[0]['edgeLow'] + nStartBin*vEnergyBinWidth
vEnergyUp = vEnergyLow + nEnergyBin*vEnergyBinWidth
aaNumEventClass=[]
for hS in range(len(aaStrSelect)):
    aaNumEventClass.append([])
    for iS in range(len(aaStrSelect[hS])):
        aaNumEventClass[hS].append(0)

#IRF
listPathFilePerf = [['/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_TRANSIENT100_P8R2_TRANSIENT100_perf.root', '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_SOURCE_P8R2_SOURCE_perf.root'], 
                    ['/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_CalOnly_R100_perf.root', '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_CalOnly_R30_perf.root', '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_CalOnly_R10_perf.root']]
htgPerf = CutPerformanceHtg(listPathFilePerf)

# Data
pathList = "/nfs/farm/g/glast/u/mtakahas/data/lists/GBMcatalogue_FluenceLt1Micro_20160514.csv"
csv = pandas.read_csv(pathList)
num_lines = sum(1 for line in open(pathList))

# OFF regions
nOff = 0 #4;
degOffOffset = 14.0
zenithCut = 90;
print "===================="
# Making all sky map
listFileIn = par[2:]
print listFileIn

aliasSelections = yaml.load(open('/afs/slac.stanford.edu/u/gl/mtakahas/eventSelect/config/pass8_event_selections.yaml','r'))
cEvent = ROOT.TCanvas("cEvent", "Gamma-like events")
cZenith = ROOT.TCanvas("cZenith", "Zenith angle of ON/OFF events")
for nameFileIn in listFileIn:
    cEvent.Clear()
    cZenith.Clear()
    print ""
    fileIn = ROOT.TFile(nameFileIn, "READ")
    chainData = fileIn.Get("MeritTuple")
    if chainData == None:
        print "No events."
        continue
    nameFileFriend = nameFileIn.replace(".root", "_" + nameVarBDT + ".root")
    chainData.AddFriend("friendTemp=MeritTuple", nameFileFriend)
    for k,v in aliasSelections.iteritems(): 
        chainData.SetAlias(k,v)

    nameFileOut = nameFileIn[:-5] + "_PHOTON_" + nameVarBDT + nameFileSuffix + ".root"
    fileOut = ROOT.TFile(nameFileOut, 'UPDATE')

    #------ Source data -----
    # nameGrb = "160509374" 
    # for iGrb in range(num_lines-1):
    #     if int(nameGrb) == int(csv.ix[iGrb,'name']):
    #         raSrc = float(csv.ix[iGrb,'ra'])
    #         decSrc = float(csv.ix[iGrb,'dec']) 
    #         trigger_time = ConvertMjdToMet(float(csv.ix[iGrb,'trigger_time']))
    #         err_rad = float(csv.ix[iGrb,'error_radius'])
    # print ""
    # print "==============="
    # print "GRB", nameGrb
    # print "==============="
    # print "(", raSrc, ",", decSrc, "), Error radius:", err_rad, "Trigger MET:", trigger_time 
    nEvent = chainData.GetEntries()
    print "Total number of events:", nEvent
    # Plot
    mgr = ROOT.TMultiGraph("mgr", "Gamma-like events around the GRB")
    greOn = []
    #greOff=[]
    for pC in range(len(aaStrSelect)):
        greOn.append([])
        #greOff.append([])
        for qC in range(len(aaStrSelect[pC])):
            greOn[-1].append(ROOT.TGraphErrors())
            greOn[-1][-1].SetName("greOn_{0}_{1}".format(pC, qC))
            greOn[-1][-1].SetTitle("{0} ON".format(aaStrSelect[pC][qC]))
            if pC==0:
                greOn[-1][-1].SetMarkerColor(13-12*qC)
            elif pC==1:
                greOn[-1][-1].SetMarkerColor(kRed+3*(qC-2))
            greOn[-1][-1].SetMarkerStyle(6)
            mgr.Add(greOn[-1][-1])
            #greOff[-1].append([])
            # for hRegio in range(nOff):
            #     greOff[-1][-1].append(ROOT.TGraphErrors())
            #     greOff[-1][-1][-1].SetName("greOff_{0}_{1}_{2}".format(pC, qC, hRegio+1))
            #     greOff[-1][-1][-1].SetTitle("{0} Off{1} events".format(aaStrSelect[pC][qC], hRegio+1))
            #     if pC==0:
            #         greOff[-1][-1][-1].SetMarkerColor(13-12*qC)
            #     elif pC==1:
            #         greOff[-1][-1][-1].SetMarkerColor(kRed+3*(qC-2))
            #     greOff[-1][-1][-1].SetMarkerStyle(25+hRegio)
            #     mgr.Add(greOff[-1][-1][-1])
    #mgrZenith = ROOT.TMultiGraph("mgrZenith", "Zenith angle within ON/OFF regions")
    grZenith = ROOT.TGraph()
    grZenith.SetName("grZenith")
    grZenith.SetTitle("ON")
    grZenith.SetMarkerStyle(7)

    #------ TTree setting -----
    trm = []
    c = np.zeros(1, dtype=np.int32)
    s = np.zeros(1, dtype=np.int32)
    ty = np.zeros(1, dtype=np.int32)
    #ngrb = np.zeros(1, dtype=int)
    evid = np.zeros(1, dtype=int)
    run = np.zeros(1, dtype=int)
    e = np.zeros(1, dtype=float)
    t = np.zeros(1, dtype=float)
    lt = np.zeros(1, dtype=float)
    ra = np.zeros(1, dtype=float)
    dec = np.zeros(1, dtype=float)
    l = np.zeros(1, dtype=float)
    b = np.zeros(1, dtype=float)
    z = np.zeros(1, dtype=float)
    az = np.zeros(1, dtype=float)
    bep = np.zeros(1, dtype=float)
    p = np.zeros(1, dtype=float)
    ctp = np.zeros(1, dtype=float)
    rawe = np.zeros(1, dtype=float)
    cth = np.zeros(1, dtype=float)
    th = np.zeros(1, dtype=float)
    phi = np.zeros(1, dtype=float)
    #dist = np.zeros(1, dtype=float)
    #grbt = np.zeros(1, dtype=float)
    #flag = np.zeros(1, dtype=int)

    trm = ROOT.TTree("trGammas", "Gamma-like events")
    trm.Branch('Category',c,'c/I') # 1:CalTkr or 2:CalOnly
    trm.Branch('EVENT_CLASS',s,'s/I') # 4: TRANSIENT100, 128: SOURCE, 4096: CalOnly_10xEGB, 8192: CalOnly_3xEGB, 16384: CalOnly_1xEGB
    trm.Branch('EVENT_TYPE',ty,'ty/I') # 1: FRONT, 2: BACK, 4, PSF0, ... , 32: PSF3, 64: EDISP0, ... , 512: EDISP3
    #trm.Branch('GRB_NAME',ngrb,'ngrb/I') #EvtEventId
    trm.Branch('EVENT_ID',evid,'evid/I') #EvtEventId
    trm.Branch('RUN_ID',run,'run/I') #EvtRun
    trm.Branch('ENERGY',e,'e/D') #FT1Energy
    trm.Branch('TIME',t,'t/D') #EvtElapsedTime
    trm.Branch('LIVETIME',lt,'lt/D') #EvtLiveTime
    trm.Branch('RA',ra,'ra/D') #FT1Ra
    trm.Branch('DEC',dec,'dec/D') #FT1Dec
    trm.Branch('L',l,'l/D') #FT1L
    trm.Branch('B',b,'b/D') #FT1B
    trm.Branch('ZENITH_ANGLE',z,'z/D') #FT1ZenithTheta
    trm.Branch('EARTH_AZIMUTH_ANGLE',az,'az/D') #FT1EarthAzimuth
    trm.Branch('WP8CTCalOnlyBestEnergyProb',bep,'bep/D') # (WP8CalOnlyBEPCaseE_myBDT+1.0)/2.0
    trm.Branch('WP8CTCalOnlyProb',p,'p/D') # (S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_BDTG1000D06+1.0)/2.0
    trm.Branch('WP8CTAllProb',ctp,'ctp/D') #WP8CTAllProb
    trm.Branch('CalEnergyRaw',rawe,'rawe/D') #CalEnergyRaw
    trm.Branch('CosTHETA',cth,'cth/D') #
    trm.Branch('THETA',th,'th/D') # FT1Theta or -Cal1MomZDir
    trm.Branch('PHI',phi,'phi/D') # FT1Phi or Cal1MomYDir/Cal1MomXDir
    #trm.Branch('DIST',dist,'dist/D')
    #trm.Branch('GRB_TIME',grbt,'grbt/D') 
    #trm.Branch('FLAG',flag,'flag/I') #flag for this GRB, 0: On, 1,2,3,4,...: Off, -1: Other
    timeStart = datetime.datetime.now()
    for iEvent in range(nEvent):
        chainData.GetEntry(iEvent)
        #flag[0] = -1;
        e[0] = chainData.EvtJointLogEnergy
        rawe[0] = chainData.CalEnergyRaw
        c[0] = 0
        s[0] = 0
        ty[0] = 0
        #ngrb[0] = float(nameGrb)
        evid[0] = chainData.EvtEventId
        run[0] = chainData.EvtRun
        t[0] = chainData.EvtElapsedTime
        #grbt[0] = t[0] - trigger_time
        lt[0] = chainData.EvtLiveTime
        bep[0] = (chainData.WP8CalOnlyBEPCaseE_myBDT+1.0)/2.0
        p[0] = (chainData.S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_BDTG1000D06+1.0)/2.0
        ctp[0] = chainData.WP8CTAllProb
        binEnergy = max(min(nEnergyBin-1, int((e[0]-vEnergyLow)/vEnergyBinWidth * (int(e[0]<vEnergyLow)*(-2)+1)) ), 0)
        if (chainData.TkrNumTracks>0) and (math.log10(max(chainData.CalTrackAngle,1E-4)) <= (0.529795)*(e[0] < 3.000000) + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((e[0]-3.000000)/0.916667,1))+(0.583401)*(pow((e[0]-3.000000)/0.916667,2))+(-0.075555)*(pow((e[0]-3.000000)/0.916667,3))))*(e[0] >= 3.000000 and e[0] <= 5.750000) + (-0.398962)*(e[0] > 5.750000)) and chainData.EvtCalCsIRLn>4 and chainData.WP8CTPSFTail>0.05 and chainData.WP8CTBestEnergyProb>0.1 and chainData.FswGamState == 0: # CalTkr
            c[0] = 1
            z[0] = chainData.FT1ZenithTheta
            az[0] = chainData.FT1EarthAzimuth
            ra[0] = chainData.FT1Ra
            dec[0] = chainData.FT1Dec
            l[0] = chainData.FT1L
            b[0] = chainData.FT1B
            cth[0] = chainData.Cal1MomZDir
            th[0] = chainData.FT1Theta
            phi[0] = chainData.FT1Phi
            if ( -math.log10(1.0-ctp[0]) >= (0.010000)*(e[0] < 1.250000) + ((e[0] <= 1.750000)*((0.010000)*(1.0)+(0.000000)*(math.pow((e[0]-1.250000)/0.500000,1))+(0.018669)*(math.pow((e[0]-1.250000)/0.500000,2)))+((e[0] > 1.750000)*(e[0] <= 2.250000))*((0.028669)*(1.0)+(0.037338)*(math.pow((e[0]-1.750000)/0.500000,1))+(-0.017111)*(math.pow((e[0]-1.750000)/0.500000,2)))+((e[0] > 2.250000)*(e[0] <= 2.750000))*((0.048897)*(1.0)+(0.003117)*(math.pow((e[0]-2.250000)/0.500000,1))+(0.001967)*(math.pow((e[0]-2.250000)/0.500000,2)))+((e[0] > 2.750000)*(e[0] <= 3.250000))*((0.053980)*(1.0)+(0.007050)*(math.pow((e[0]-2.750000)/0.500000,1))+(-0.003525)*(math.pow((e[0]-2.750000)/0.500000,2)))+((e[0] > 3.250000)*(e[0] <= 3.750000))*((0.057505)*(1.0)+(0.000000)*(math.pow((e[0]-3.250000)/0.500000,1))+(0.121963)*(math.pow((e[0]-3.250000)/0.500000,2)))+((e[0] > 3.750000)*(e[0] <= 4.250000))*((0.179468)*(1.0)+(0.243925)*(math.pow((e[0]-3.750000)/0.500000,1))+(0.493075)*(math.pow((e[0]-3.750000)/0.500000,2)))+((e[0] > 4.250000)*(e[0] <= 4.750000))*((0.916468)*(1.0)+(1.230076)*(math.pow((e[0]-4.250000)/0.500000,1))+(-0.501532)*(math.pow((e[0]-4.250000)/0.500000,2)))+(e[0] > 4.750000)*((1.645012)*(1.0)+(0.227011)*(math.pow((e[0]-4.750000)/0.500000,1))+(0.029483)*(math.pow((e[0]-4.750000)/0.500000,2))))*(e[0] >= 1.250000 and e[0] <= 5.750000) + (2.216967)*(e[0] > 5.750000) ): #P8R1_TRANSIENT_R100
                if ( -math.log10(1.0-ctp[0]) >= (0.080914)*(e[0] < 1.250000) + ((e[0] <= 1.750000)*((0.080914)*(1.0)+(0.108897)*(pow((e[0]-1.250000)/0.500000,1))+(0.377870)*(pow((e[0]-1.250000)/0.500000,2)))+((e[0] > 1.750000)*(e[0] <= 2.250000))*((0.567682)*(1.0)+(0.864637)*(pow((e[0]-1.750000)/0.500000,1))+(-0.182318)*(pow((e[0]-1.750000)/0.500000,2)))+((e[0] > 2.250000)*(e[0] <= 2.750000))*((1.250000)*(1.0)+(0.500000)*(pow((e[0]-2.250000)/0.500000,1))+(-0.085000)*(pow((e[0]-2.250000)/0.500000,2)))+((e[0] > 2.750000)*(e[0] <= 3.250000))*((1.665000)*(1.0)+(0.330000)*(pow((e[0]-2.750000)/0.500000,1))+(-0.165000)*(pow((e[0]-2.750000)/0.500000,2)))+((e[0] > 3.250000)*(e[0] <= 3.750000))*((1.830000)*(1.0)+(0.000000)*(pow((e[0]-3.250000)/0.500000,1))+(0.285000)*(pow((e[0]-3.250000)/0.500000,2)))+((e[0] > 3.750000)*(e[0] <= 4.250000))*((2.115000)*(1.0)+(0.570000)*(pow((e[0]-3.750000)/0.500000,1))+(-0.185000)*(pow((e[0]-3.750000)/0.500000,2)))+((e[0] > 4.250000)*(e[0] <= 4.750000))*((2.500000)*(1.0)+(0.200000)*(pow((e[0]-4.250000)/0.500000,1))+(0.100000)*(pow((e[0]-4.250000)/0.500000,2)))+(e[0] > 4.750000)*((2.800000)*(1.0)+(0.400000)*(pow((e[0]-4.750000)/0.500000,1))+(-0.112171)*(pow((e[0]-4.750000)/0.500000,2))))*(e[0] >= 1.250000 and e[0] <= 5.750000) + (3.151318)*(e[0] > 5.750000) ) and ( chainData.WP8CTAllBkProb >= (0.366167)*(e[0] < 1.250000) + ((e[0] <= 1.541667)*((0.366167)*(1.0)+(0.028500)*(pow((e[0]-1.250000)/0.291667,1))+(-0.056500)*(pow((e[0]-1.250000)/0.291667,2))+(0.106667)*(pow((e[0]-1.250000)/0.291667,3)))+((e[0] > 1.541667)*(e[0] <= 1.833333))*((0.444833)*(1.0)+(0.235500)*(pow((e[0]-1.541667)/0.291667,1))+(0.263500)*(pow((e[0]-1.541667)/0.291667,2))+(-0.162667)*(pow((e[0]-1.541667)/0.291667,3)))+((e[0] > 1.833333)*(e[0] <= 2.125000))*((0.781167)*(1.0)+(0.274500)*(pow((e[0]-1.833333)/0.291667,1))+(-0.224500)*(pow((e[0]-1.833333)/0.291667,2))+(0.072667)*(pow((e[0]-1.833333)/0.291667,3)))+(e[0] > 2.125000)*((0.903833)*(1.0)+(0.043500)*(pow((e[0]-2.125000)/0.291667,1))+(-0.006500)*(pow((e[0]-2.125000)/0.291667,2))+(-0.000333)*(pow((e[0]-2.125000)/0.291667,3))))*(e[0] >= 1.250000 and e[0] <= 3.000000) + (0.966833)*(e[0] > 3.000000) ):  #P8R1_SOURCE_AllProbFilter&&P8R1_SOURCE_AllBkProbFilter
                    s[0] = 128#3
                    aaNumEventClass[0][1] = aaNumEventClass[0][1]+1
                else:
                    s[0] = 4#1
                    aaNumEventClass[0][0] = aaNumEventClass[0][0]+1
        elif chainData.Cal1RawEnergySum>=20000 and chainData.Cal1MomZDir>=0.1 and chainData.Cal1MomZCrossSide840>=0.0 and (chainData.WP8CalOnlyBEPCaseE_myBDT+1.0)/2.0>0.06 and (chainData.TkrNumTracks==0 or (math.log10(max(chainData.CalTrackAngle,1E-4)) > (0.529795)*(e[0] < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((e[0]-3.000000)/0.916667,1))+(0.583401)*(pow((e[0]-3.000000)/0.916667,2))+(-0.075555)*(pow((e[0]-3.000000)/0.916667,3))))*(e[0]  >= 3.000000 and e[0] <= 5.750000) + (-0.398962)*(e[0] >  5.750000))) and chainData.Acd2Cal1VetoSigmaHit>0 and chainData.Cal1TransRms>=10 and chainData.Cal1TransRms<70 and chainData.Cal1MomNumIterations>0 and chainData.FswGamState == 0: # CalOnly
            c[0] = 2
            z[0] = chainData.FT1CalZenithTheta
            az[0] = chainData.FT1CalEarthAzimuth
            ra[0] = chainData.FT1CalRa
            dec[0] = chainData.FT1CalDec
            l[0] = chainData.FT1CalL
            b[0] = chainData.FT1CalB
            cth[0] = chainData.Cal1MomZDir
            th[0] = math.degrees(math.acos(chainData.Cal1MomZDir))
            phi[0] = math.degrees(math.atan2(chainData.Cal1MomYDir, chainData.Cal1MomXDir))
            if -math.log10(1.0-p[0])>aaValCutBDT[binEnergy+nStartBin][0]: #CalOnly_R100
                if -math.log10(1.0-p[0])>aaValCutBDT[binEnergy+nStartBin][1]: #CalOnly_R30
                    if -math.log10(1.0-p[0])>aaValCutBDT[binEnergy+nStartBin][2]: #CalOnly_R10
                        s[0]=16384#3
                        aaNumEventClass[1][2] = aaNumEventClass[1][2]+1
                    else:
                        s[0] = 8192#2
                        aaNumEventClass[1][1] = aaNumEventClass[1][1]+1
                else:
                    s[0] = 4096#1
                    aaNumEventClass[1][0] = aaNumEventClass[1][0]+1
                if(e[0]<4.55 or cth[0]<0.6):
                    ty[0] = ty[0]+2 #BACK
                else:
                    ty[0] = ty[0]+1 #FRONT

        vecEvt = np.array([cos(radians(dec[0]))*cos(radians(ra[0])), cos(radians(dec[0]))*sin(radians(ra[0])), sin(radians(dec[0]))])
        grZenith.SetPoint(grZenith.GetN(), t[0], z[0])
        if s[0]>0 and z[0]<zenithCut:
            #print ""
            #print "== ON photon candidate!!! =="
            if c[0] == 1:
                if s[0] == 4:
                    greOn[0][0].SetPoint(greOn[0][0].GetN(), t[0], pow(10, e[0]-3))
                elif s[0] == 128:
                    greOn[0][1].SetPoint(greOn[0][1].GetN(), t[0], pow(10, e[0]-3))
            elif c[0] == 2:
                if s[0] == 4096:
                    greOn[1][0].SetPoint(greOn[1][0].GetN(), t[0], pow(10, e[0]-3))
                elif s[0] == 8192:
                    greOn[1][1].SetPoint(greOn[1][1].GetN(), t[0], pow(10, e[0]-3))
                elif s[0] == 16384:
                    greOn[1][2].SetPoint(greOn[1][2].GetN(), t[0], pow(10, e[0]-3))
            #print "Event No.", iEvent
            #print "Event category:", cfg.aStrSelect[c[0]-1]
            #print "Event class:", s[0]
            #print "Time:", t[0], "s"
            #print "PSF68:", htgPerf.getPSF68_cth(c[0]-1, 0*(s[0]==4 or s[0]==4096)+1*(s[0]==128 or s[0]==8192)+2*(s[0]==16384), e[0], cth[0]), "deg"
            #print "Energy:", pow(10,e[0]-3), "GeV"
            #print "Edisp68:", 100*htgPerf.getEdisp68_cth(c[0]-1, 0*(s[0]==4 or s[0]==4096)+1*(s[0]==128 or s[0]==8192)+2*(s[0]==16384), e[0], cth[0]), "%"
            #print "Cos( inclination angle ):", cth[0]
            #print "Zenith angle:", z[0], "deg"
            #print "Run ID:", run[0]
            #print "Event ID:", evid[0]
            trm.Fill()
    cEvent.SetTitle("Gamma-like events")
    cEvent.cd()
    mgr.Draw("AP")
    mgr.GetXaxis().SetTitle("Time [s]")
    mgr.GetYaxis().SetTitle("Energy [GeV]")
    mgr.GetYaxis().SetRangeUser(10, 10000)
    cEvent.SetLogy()
    leg = ROOT.TLegend(0.67, 0.5, 0.88, 0.88)
    for pD in range(len(aaStrSelect)):
        for qD in range(len(aaStrSelect[pD])):
            leg.AddEntry(greOn[pD][qD], greOn[pD][qD].GetTitle(), "p")
    leg.Draw("same")
    cZenith.cd()
    grZenith.Draw("AP")
    grZenith.GetXaxis().SetTitle("Time [s]")
    grZenith.GetYaxis().SetTitle("Zenith angle [deg]")
    print ""
    fileOut.cd()
    trm.Write()
    cEvent.Write()
    cZenith.Write()
    print "Finished!"
