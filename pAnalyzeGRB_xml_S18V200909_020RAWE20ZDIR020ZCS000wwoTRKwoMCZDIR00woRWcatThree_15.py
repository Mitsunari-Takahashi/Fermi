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
nStartBin = cfg.nStartBin

nameFileRoc = ["/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatThree_15/S18ZDIR020catThreeZDIR060_ZDIR020to060_S18ZDIR020catThreeZDIR060Log_ZDIR020to060_roc.root", 
               "/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatThree_15/S18ZDIR020catThreeZDIR060_ZDIR060to090_S18ZDIR020catThreeZDIR060Log_ZDIR060to090_roc.root",
               "/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatThree_15/S18ZDIR020catThree_ZDIR090to100_S18ZDIR020catThreeZDIR060Log_ZDIR090to100_roc.root"]
nameVarBDT = "S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatThree_15_BDTG500D06_catZDIR060"
nameFileSuffix = "S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatThree_15_BDTG500D06_catZDIR060"

aCutMVA = []
for iRocFile in range(len(nameFileRoc)):
    aCutMVA.append(CutBDT(nameFileRoc[iRocFile], aCutEGB))
aCosZ=[0.2, 0.6, 0.9, 1.0]
cutMVA_comb = CutBDT_2D(aCutMVA, aCosZ, "S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatThree_15_BDTG500D06_catZDIR060")
print cutMVA_comb
aaaValCutBDT = cutMVA_comb.aaValCutBDT
aHtgCutBDT = []
fileHtgBDT = ROOT.TFile('hBDT_CutValues_{0}.root'.format("S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatThree_15_BDTG500D06_catZDIR060"))
for jEGB in range(len(aCutEGB)):
    htgBDT = fileHtgBDT.Get('htgCutBDT{0}'.format(jEGB))
    aHtgCutBDT.append(htgBDT)
print aHtgCutBDT
print aaaValCutBDT
nEnergyBin = aaaValCutBDT[0]['numBin'] - nStartBin
vEnergyBinWidth = aaaValCutBDT[0]['widthBin']
vEnergyLow = aaaValCutBDT[0]['edgeLow'] + nStartBin*vEnergyBinWidth
vEnergyUp = vEnergyLow + nEnergyBin*vEnergyBinWidth
nCosthBin = aaaValCutBDT[0]['numBinY']
aCosthLow = aaaValCutBDT[0]['aEdgeLowY']
aCosthUp = aaaValCutBDT[0]['aEdgeUpY']
aaNumEventClass=[]
#aColor = []
for hS in range(len(aaStrSelect)):
    aaNumEventClass.append([])
    for iS in range(len(aaStrSelect[hS])):
        aaNumEventClass[hS].append(0)

#IRF
listPathFilePerf = [['/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_TRANSIENT100_P8R2_TRANSIENT100_perf.root', 
                     '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S16/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15/S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_P8R2_SOURCE_P8R2_SOURCE_perf.root'], 
                    ['/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatThree_15/S18ZDIR020catThreeZDIR060_CalOnly_R100_perf.root', 
                     '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatThree_15/S18ZDIR020catThreeZDIR060_CalOnly_R30_perf.root',
                     '/nfs/farm/g/glast/u/mtakahas/v20r09p09_G1haB1/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatThree_15/S18ZDIR020catThreeZDIR060_CalOnly_R10_perf.root']]
htgPerf = CutPerformanceHtg(listPathFilePerf)

# Data
pathList = "/nfs/farm/g/glast/u/mtakahas/data/catalogue/PublicTableGRBs.xml"
fileList = ET.parse(pathList)
rtXml = fileList.getroot()

# OFF regions
nOff = 4;
degOffOffset = 14.0;

print "===================="
# Making all sky map
listFileIn = par[1:]
print listFileIn

aliasSelections = yaml.load(open('/afs/slac.stanford.edu/u/gl/mtakahas/eventSelect/config/pass8_event_selections.yaml','r'))
for nameFileIn in listFileIn:
    print ""
    print "========================================================================"
    fileIn = ROOT.TFile(nameFileIn, "READ")
    print fileIn.GetName()
    print "========================================================================"
    chainData = fileIn.Get("MeritTuple")
    nameFileFriend = nameFileIn.replace(".root", "_" + nameFileSuffix + ".root")
    chainData.AddFriend("friendTemp=MeritTuple", nameFileFriend)
    chainData.GetListOfFriends().Print()
    for k,v in aliasSelections.iteritems(): 
        chainData.SetAlias(k,v)

    nameFileOut = nameFileIn[:-5] + "_PHOTON_" + nameVarBDT + nameFileSuffix + ".root"
    fileOut = ROOT.TFile(nameFileOut, 'UPDATE')

    #------ Source data -----
    indexGrbName = nameFileIn.rindex('GRB') + 3
    indexGrbNameEnd = indexGrbName + 9
    nameGrb = nameFileIn[indexGrbName:indexGrbNameEnd]
    for grb in rtXml: #for iGrb in range(trList.GetEntries())
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
    nEvent = chainData.GetEntries()
    print "Total number of events:", nEvent
    #distCut = min(7.0, 5.0+err_rad)

    # Plot
    #mgr = ROOT.TMultiGraph("mgr", "Gamma-like events within {0} deg".format(distCut))
    mgr = ROOT.TMultiGraph("mgr", "Gamma-like events around the GRB")
    greOn = []
    greOff=[]
    for pC in range(len(aaStrSelect)):
        greOn.append([])
        greOff.append([])
        for qC in range(len(aaStrSelect[pC])):
            greOn[-1].append(ROOT.TGraphErrors())
            greOn[-1][-1].SetName("greOn_{0}_{1}".format(pC, qC))
            greOn[-1][-1].SetTitle("{0} ON".format(aaStrSelect[pC][qC]))
            greOn[-1][-1].SetMarkerStyle(20)
            if pC==0:
                greOn[-1][-1].SetMarkerColor(13-12*qC)
            elif pC==1:
                greOn[-1][-1].SetMarkerColor(kRed+3*(qC-2))
            greOn[-1][-1].SetMarkerStyle(20)
            mgr.Add(greOn[-1][-1])
            greOff[-1].append([])
            for hRegio in range(nOff):
                greOff[-1][-1].append(ROOT.TGraphErrors())
                greOff[-1][-1][-1].SetName("greOff_{0}_{1}_{2}".format(pC, qC, hRegio+1))
                greOff[-1][-1][-1].SetTitle("{0} Off{1} events".format(aaStrSelect[pC][qC], hRegio+1))
                if pC==0:
                    greOff[-1][-1][-1].SetMarkerColor(13-12*qC)
                elif pC==1:
                    greOff[-1][-1][-1].SetMarkerColor(kRed+3*(qC-2))
                greOff[-1][-1][-1].SetMarkerStyle(25+hRegio)
                mgr.Add(greOff[-1][-1][-1])
    mgrZenith = ROOT.TMultiGraph("mgrZenith", "Zenith angle within ON/OFF regions")
    grZenith = []
    for gRegio in range(nOff+1):
        grZenith.append(ROOT.TGraph())
        grZenith[-1].SetName("grZenith{0}".format(gRegio))
        if gRegio==0:
            grZenith[0].SetTitle("ON")
        else:
            grZenith[gRegio].SetTitle("OFF{0}".format(gRegio))
        grZenith[gRegio].SetMarkerStyle(7)
        grZenith[gRegio].SetMarkerColor(akColor(gRegio))
        mgrZenith.Add(grZenith[-1])

    #------ TTree setting -----
    trm = []
    c = np.zeros(1, dtype=np.int32)
    s = np.zeros(1, dtype=np.int32)
    ty = np.zeros(1, dtype=np.int32)
    ngrb = np.zeros(1, dtype=int)
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
    dist = np.zeros(1, dtype=float)
    grbt = np.zeros(1, dtype=float)
    flag = np.zeros(1, dtype=int)

#    for iRegio in range(1+nOff):
    #if iRegio==0:
    #trm.append(ROOT.TTree("trGammas", "Gamma-like events"))
    trm = ROOT.TTree("trGammas", "Gamma-like events")
        #else:
         #   trm.append(ROOT.TTree("trGammasOFF{0}".format(iRegio), "Gamma-like events in the OFF region {0}".format(iRegio)))
    trm.Branch('Category',c,'c/I') # 1:CalTkr or 2:CalOnly
    trm.Branch('EVENT_CLASS',s,'s/I') # 4: TRANSIENT100, 128: SOURCE, 4096: CalOnly_10xEGB, 8192: CalOnly_3xEGB, 16384: CalOnly_1xEGB
    trm.Branch('EVENT_TYPE',ty,'ty/I') # 1: FRONT, 2: BACK, 4, PSF0, ... , 32: PSF3, 64: EDISP0, ... , 512: EDISP3
    trm.Branch('GRB_NAME',ngrb,'ngrb/I') #EvtEventId
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
    trm.Branch('DIST',dist,'dist/D')
    trm.Branch('GRB_TIME',grbt,'grbt/D') 
    trm.Branch('FLAG',flag,'flag/I') #flag for this GRB, 0: On, 1,2,3,4,...: Off, -1: Other

    timeStart = datetime.datetime.now()
    #print timeStart

    vecTgt = []
    vecTgt.append(np.array([cos(radians(decSrc))*cos(radians(raSrc)), cos(radians(decSrc))*sin(radians(raSrc)), sin(radians(decSrc))]))
    vecTgt.append(np.array([cos(radians(decSrc-degOffOffset))*cos(radians(raSrc)), cos(radians(decSrc-degOffOffset))*sin(radians(raSrc)), sin(radians(decSrc-degOffOffset))]))
    vecTgt.append(np.array([cos(radians(decSrc))*cos(radians(raSrc-degOffOffset/cos(radians(decSrc)))), cos(radians(decSrc))*sin(radians(raSrc-degOffOffset/cos(radians(decSrc)))), sin(radians(decSrc))]))
    vecTgt.append(np.array([cos(radians(decSrc+degOffOffset))*cos(radians(raSrc)), cos(radians(decSrc+degOffOffset))*sin(radians(raSrc)), sin(radians(decSrc+degOffOffset))]))
    vecTgt.append(np.array([cos(radians(decSrc))*cos(radians(raSrc+degOffOffset/cos(radians(decSrc)))), cos(radians(decSrc))*sin(radians(raSrc+degOffOffset/cos(radians(decSrc)))), sin(radians(decSrc))]))

    for iEvent in range(nEvent):
        chainData.GetEntry(iEvent)
        flag[0] = -1;
        e[0] = chainData.EvtJointLogEnergy
        rawe[0] = chainData.CalEnergyRaw
        c[0] = 0
        s[0] = 0
        ty[0] = 0
        ngrb[0] = float(nameGrb)
        evid[0] = chainData.EvtEventId
        run[0] = chainData.EvtRun
        t[0] = chainData.EvtElapsedTime
        grbt[0] = t[0] - trigger_time
        lt[0] = chainData.EvtLiveTime
        bep[0] = (chainData.WP8CalOnlyBEPCaseE_myBDT+1.0)/2.0
        p[0] = (chainData.S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatThree_15_BDTG500D06_catZDIR060+1.0)/2.0
        ctp[0] = chainData.WP8CTAllProb
#        binEnergy = max(min(nEnergyBin-1, int((e[0]-vEnergyLow)/vEnergyBinWidth * (int(e[0]<vEnergyLow)*(-2)+1)) ), 0)
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
                #trm.Fill()
        elif chainData.Cal1RawEnergySum>=20000 and chainData.Cal1MomZDir>=0.2 and chainData.Cal1MomZCrossSide840>=0.0 and (chainData.TkrNumTracks==0 or (math.log10(max(chainData.CalTrackAngle,1E-4)) > (0.529795)*(e[0] < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((e[0]-3.000000)/0.916667,1))+(0.583401)*(pow((e[0]-3.000000)/0.916667,2))+(-0.075555)*(pow((e[0]-3.000000)/0.916667,3))))*(e[0]  >= 3.000000 and e[0] <= 5.750000) + (-0.398962)*(e[0] >  5.750000))) and chainData.Acd2Cal1VetoSigmaHit>0 and chainData.Cal1TransRms>=10 and chainData.Cal1TransRms<70 and chainData.Cal1MomNumIterations>0 and chainData.FswGamState == 0 and e[0]>=4.35 and e[0]<5.75: # CalOnly
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
            if -math.log10(1.0-p[0])>=aHtgCutBDT[0].GetBinContent(aHtgCutBDT[0].GetXaxis().FindBin(e[0]), aHtgCutBDT[0].GetYaxis().FindBin(cth[0])): #CalOnly_R100
                if -math.log10(1.0-p[0])>aHtgCutBDT[1].GetBinContent(aHtgCutBDT[1].GetXaxis().FindBin(e[0]), aHtgCutBDT[1].GetYaxis().FindBin(cth[0])): #CalOnly_R30
                    if -math.log10(1.0-p[0])>aHtgCutBDT[2].GetBinContent(aHtgCutBDT[2].GetXaxis().FindBin(e[0]), aHtgCutBDT[2].GetYaxis().FindBin(cth[0])): #CalOnly_R10
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
        aDist = []
        distCut = htgPerf.getPSF95_cth(c[0]-1, 0*(s[0]==4 or s[0]==4096)+1*(s[0]==128 or s[0]==8192)+2*(s[0]==16384), e[0], cth[0]) + err_rad
        for iRegio in range(1+nOff):
            radTheta = acos(np.dot(vecTgt[iRegio], vecEvt))
            aDist.append(degrees(radTheta))
            if iRegio==0:
                dist[0] = aDist[0]
            if aDist[iRegio] < distCut:
                grZenith[iRegio].SetPoint(grZenith[iRegio].GetN(), t[0]-trigger_time, z[0])
                #if s[0]>0:
                    #print "============================"
                if iRegio==0:
                    flag[0] = 0
                    if s[0]>0:
                        print ""
                        print "== ON photon candidate!!! =="
                        if c[0] == 1:
                            if s[0] == 4:
                                greOn[0][0].SetPoint(greOn[0][0].GetN(), t[0]-trigger_time, pow(10, e[0]-3))
                            elif s[0] == 128:
                                greOn[0][1].SetPoint(greOn[0][1].GetN(), t[0]-trigger_time, pow(10, e[0]-3))
                        elif c[0] == 2:
                            if s[0] == 4096:
                                greOn[1][0].SetPoint(greOn[1][0].GetN(), t[0]-trigger_time, pow(10, e[0]-3))
                            elif s[0] == 8192:
                                greOn[1][1].SetPoint(greOn[1][1].GetN(), t[0]-trigger_time, pow(10, e[0]-3))
                            elif s[0] == 16384:
                                greOn[1][2].SetPoint(greOn[1][2].GetN(), t[0]-trigger_time, pow(10, e[0]-3))
                elif iRegio>0:
                    flag[0] = iRegio
                    if s[0]>0:
                        #print ""
                        #print "== OFF{0} photon candidate! ==".format(iRegio)
                        if c[0] == 1:
                            if s[0] == 4:
                                greOff[0][0][iRegio-1].SetPoint(greOff[0][0][iRegio-1].GetN(), t[0]-trigger_time, pow(10, e[0]-3))
                            elif s[0] == 128:
                                greOff[0][1][iRegio-1].SetPoint(greOff[0][1][iRegio-1].GetN(), t[0]-trigger_time, pow(10, e[0]-3))
                        elif c[0] == 2:
                            if s[0] == 4096:
                                greOff[1][0][iRegio-1].SetPoint(greOff[1][0][iRegio-1].GetN(), t[0]-trigger_time, pow(10, e[0]-3))
                            elif s[0] == 8192:
                                greOff[1][1][iRegio-1].SetPoint(greOff[1][1][iRegio-1].GetN(), t[0]-trigger_time, pow(10, e[0]-3))
                            elif s[0] == 16384:
                                greOff[1][2][iRegio-1].SetPoint(greOff[1][2][iRegio-1].GetN(), t[0]-trigger_time, pow(10, e[0]-3))
                if s[0]>0 and iRegio==0:
                    print "Event No.", iEvent
                    print "Event category:", cfg.aStrSelect[c[0]-1]
                    print "Event class:", s[0]
                    print "Time from the trigger:", t[0]-trigger_time, "s"
                    print "Anglular distance:", dist[0], "deg"
                    print "PSF68:", htgPerf.getPSF68_cth(c[0]-1, 0*(s[0]==4 or s[0]==4096)+1*(s[0]==128 or s[0]==8192)+2*(s[0]==16384), e[0], cth[0]), "deg"
                    print "Energy:", pow(10,e[0]-3), "GeV"
                    print "Edisp68:", 100*htgPerf.getEdisp68_cth(c[0]-1, 0*(s[0]==4 or s[0]==4096)+1*(s[0]==128 or s[0]==8192)+2*(s[0]==16384), e[0], cth[0]), "%"
                    print "Cos( inclination angle ):", cth[0]
                    print "Zenith angle:", z[0], "deg"
                    print "Run ID:", run[0]
                    print "Event ID:", evid[0]
                trm.Fill()
#            trm.Fill()
        if iEvent%(nEvent/200)==0:
            rate = int((iEvent*100.)/nEvent+0.5)
            if rate>0:
                nt = (datetime.datetime.now() - timeStart).seconds * (100.-rate)/rate
                meter = "\r[{0}{1}] {2} Wait {3} hr {4} min".format("=" * rate, ' ' * (100-rate), aaNumEventClass, int(nt/3600), (int(nt)%3600)/60+1)
            else:
                meter = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
            sys.stdout.write(meter)
            sys.stdout.flush()
    cEvent = ROOT.TCanvas("cEvent", "GRB {0} gamma-like events within {1} deg".format(nameGrb, distCut))
    cEvent.cd()
    mgr.Draw("AP")
    mgr.GetXaxis().SetTitle("Time [s]")
    mgr.GetYaxis().SetTitle("Energy [GeV]")
    leg = ROOT.TLegend(0.67, 0.5, 0.88, 0.88)
    for pD in range(len(aaStrSelect)):
        for qD in range(len(aaStrSelect[pD])):
            leg.AddEntry(greOn[pD][qD], greOn[pD][qD].GetTitle(), "p")
    for rr in range(nOff):
        leg.AddEntry(greOff[1][0][rr], "OFF{0}".format(rr+1), "p")
    leg.Draw("same")
    cZenith = ROOT.TCanvas("cZenith", "Zenith angle of ON/OFF events")
    cZenith.cd()
    mgrZenith.Draw("AP")
    mgrZenith.GetXaxis().SetTitle("Time [s]")
    mgrZenith.GetYaxis().SetTitle("Zenith angle [deg]")
    legZenith = ROOT.TLegend(0.67, 0.5, 0.88, 0.88)
    for rz in range(nOff+1):
        legZenith.AddEntry(grZenith[rz], grZenith[rz].GetTitle(), "p")
    legZenith.Draw("same")
    print ""
    fileOut.cd()
    trm.Write()
    cEvent.Write()
    cZenith.Write()
    print "Finished!"
