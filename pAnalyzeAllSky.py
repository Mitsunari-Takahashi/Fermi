#!/usr/bin/env python

import sys
import ROOT
from ROOT import TTree
from ROOT import TChain
import numpy as np
import yaml
import xml.etree.ElementTree as ET
import datetime
#sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
#sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()
from array import array
import math

ROOT.gStyle.SetPadGridX(True)
ROOT.gStyle.SetPadGridY(True)
ROOT.gStyle.SetPadTickX(True)
ROOT.gStyle.SetPadTickY(True)

from ctypes import *
#from pCutBDT import cutBDT
from pAnalysisConfig import *

# ----- Event class setup -----
par = sys.argv
cfg = ClassConfig('Both', [10, 3, 1], 1)
aCutEGB = cfg.aCutEGB
aaStrSelect = cfg.aaStrSelect
nStartBin = cfg.nStartBin
#nameVarBDT = "S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060"
nameVarBDT = par[1]
nameFileOutSuffix = par[2]
fileRoc = ROOT.TFile("hBDT_CutValues_{0}.root".format(nameVarBDT), 'READ')
aHtgCutBDT = []
for hCutEGB in range(len(aCutEGB)):
    aHtgCutBDT.append(fileRoc.Get("htgCutBDT{0}".format(hCutEGB)))

print "===================="
# Making all sky map
listFileIn = par[3:]
print listFileIn
chainData = ROOT.TChain("MeritTuple")
for nameFileIn in listFileIn:
    chainData.Add(nameFileIn)

chainFriend = ROOT.TChain("MeritTuple")
for nameFileIn in listFileIn:
    nameFileFriend = nameFileIn.replace(".root", "_" + nameVarBDT + ".root")
    print nameFileFriend
    chainFriend.Add(nameFileFriend)
chainData.AddFriend(chainFriend, "friendTemp=MeritTuple")
print "Friend: ", chainData.GetListOfFriends().FindObject("friendTemp=MeritTuple").Print()
print "---------------------"

aliasSelections = yaml.load(open('/afs/slac.stanford.edu/u/gl/mtakahas/eventSelect/config/pass8_event_selections.yaml','r'))
for k,v in aliasSelections.iteritems(): 
    chainData.SetAlias(k,v)

#nEnergyBin = cutMVA.aaValCutBDT[0]['numBin'] - nStartBin
nEnergyBin = aHtgCutBDT[0].GetNbinsX() - nStartBin
#vEnergyBinWidth = cutMVA.aaValCutBDT[0]['widthBin']
vEnergyBinWidth = aHtgCutBDT[0].GetXaxis().GetBinWidth(1)
#vEnergyLow = cutMVA.aaValCutBDT[0]['edgeLow'] + nStartBin*vEnergyBinWidth
vEnergyLow = aHtgCutBDT[0].GetXaxis().GetBinLowEdge(1) + nStartBin*vEnergyBinWidth
vEnergyUp = vEnergyLow + nEnergyBin*vEnergyBinWidth
#vEnergyUp = aHtgCutBDT[0].GetXaxis().GetBinUpEdge(aHtgCutBDT[0].GetNbinsX())

nameFileOut = "trAllSkyMap_S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_BDTG1000D06_LPA_P301.root"
nameFileOut = "trAllSkyMap_" + nameVarBDT + "_" + nameFileOutSuffix + ".root"
fileOut = ROOT.TFile(nameFileOut, 'UPDATE')
aaNumEventClass=[]
for hS in range(len(aaStrSelect)):
    aaNumEventClass.append([])
    for iS in range(len(aaStrSelect[hS])):
        aaNumEventClass[hS].append(0)

#------ Input TTree setting -----
bdtR = c_float() 
chainData.SetBranchAddress(nameVarBDT,bdtR)

#bdtR = np.zeros(1, dtype=float)
#bdtR = array(nameVarBDT,[0.0])
#chainData.SetBranchAddress(nameVarBDT, bdtR) # S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_BDTG1000D06
#chainData.SetBranchAddress('{0}'.format(nameVarBDT), bdtR) # S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_BDTG1000D06

#------ Output TTree setting -----
trm = ROOT.TTree("EVENTS", "Gamma-like events")
c = np.zeros(1, dtype=np.int32)
trm.Branch('Category',c,'c/I') # 1:CalTkr or 2:CalOnly
s = np.zeros(1, dtype=np.int32)
trm.Branch('EVENT_CLASS',s,'s/I') # 4: TRANSIENT100, 128: SOURCE, 4096: CalOnly_10xEGB, 8192: CalOnly_3xEGB, 16384: CalOnly_1xEGB
ty = np.zeros(1, dtype=np.int32)
trm.Branch('EVENT_TYPE',ty,'ty/I') # 1: FRONT, 2: BACK, 4, PSF0, ... , 32: PSF3, 64: EDISP0, ... , 512: EDISP3
evid = np.zeros(1, dtype=int)
trm.Branch('EVENT_ID',evid,'evid/I') #EvtEventId
run = np.zeros(1, dtype=int)
trm.Branch('RUN_ID',run,'run/I') #EvtRun
e = np.zeros(1, dtype=float)
trm.Branch('ENERGY',e,'e/D') #FT1Energy
t = np.zeros(1, dtype=float)
trm.Branch('TIME',t,'t/D') #EvtElapsedTime
lt = np.zeros(1, dtype=float)
trm.Branch('LIVETIME',lt,'lt/D') #EvtLiveTime
ra = np.zeros(1, dtype=float)
trm.Branch('RA',ra,'ra/D') #FT1Ra
dec = np.zeros(1, dtype=float)
trm.Branch('DEC',dec,'dec/D') #FT1Dec
l = np.zeros(1, dtype=float)
trm.Branch('L',l,'l/D') #FT1L
b = np.zeros(1, dtype=float)
trm.Branch('B',b,'b/D') #FT1B
z = np.zeros(1, dtype=float)
trm.Branch('ZENITH_ANGLE',z,'z/D') #FT1ZenithTheta
az = np.zeros(1, dtype=float)
trm.Branch('EARTH_AZIMUTH_ANGLE',az,'az/D') #FT1EarthAzimuth
bep = np.zeros(1, dtype=float)
trm.Branch('WP8CTCalOnlyBestEnergyProb',bep,'bep/D') # (WP8CalOnlyBEPCaseE_myBDT+1.0)/2.0
p = np.zeros(1, dtype=float)
trm.Branch('WP8CTCalOnlyProb',p,'p/D') # (S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_BDTG1000D06+1.0)/2.0
ctp = np.zeros(1, dtype=float)
trm.Branch('WP8CTAllProb',ctp,'ctp/D') #WP8CTAllProb
rawe = np.zeros(1, dtype=float)
trm.Branch('CalEnergyRaw',rawe,'rawe/D') #CalEnergyRaw
cth = np.zeros(1, dtype=float)
trm.Branch('CosTHETA',cth,'cth/D') #
th = np.zeros(1, dtype=float)
trm.Branch('THETA',th,'th/D') # FT1Theta or -Cal1MomZDir
phi = np.zeros(1, dtype=float)
trm.Branch('PHI',phi,'phi/D') # FT1Phi or Cal1MomYDir/Cal1MomXDir

nEvent = chainData.GetEntries()
print "Total number of events:", nEvent
timeStart = datetime.datetime.now()
print timeStart

for iEvent in range(nEvent):
    chainData.GetEntry(iEvent)
    e[0] = chainData.EvtJointLogEnergy
    rawe[0] = chainData.CalEnergyRaw
    c[0] = 0
    s[0] = 0
    ty[0] = 0
    evid[0] = chainData.EvtEventId
    run[0] = chainData.EvtRun
    t[0] = chainData.EvtElapsedTime
    lt[0] = chainData.EvtLiveTime
    bep[0] = (chainData.WP8CalOnlyBEPCaseE_myBDT+1.0)/2.0
#    p[0] = (chainData.S16V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatTwoZDIR050_15_BDTG1000D06+1.0)/2.0
#    print "bdtR=",bdtR.value
#    print "bdtR=",chainData.S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060
    p[0] = (bdtR.value+1.0)/2.0
#    print "p=",p[0]
    ctp[0] = chainData.WP8CTAllProb
#    binEnergy = max(min(nEnergyBin-1, int((e[0]-vEnergyLow)/vEnergyBinWidth * (int(e[0]<vEnergyLow)*(-2)+1)) ), 0)
    if (chainData.TkrNumTracks>0) and (math.log10(max(chainData.CalTrackAngle,1E-4)) <= (0.529795)*(e[0] < 3.000000) + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((e[0]-3.000000)/0.916667,1))+(0.583401)*(pow((e[0]-3.000000)/0.916667,2))+(-0.075555)*(pow((e[0]-3.000000)/0.916667,3))))*(e[0] >= 3.000000 and e[0] <= 5.750000) + (-0.398962)*(e[0] > 5.750000)) and chainData.EvtCalCsIRLn>4 and chainData.WP8CTPSFTail>0.05 and chainData.WP8CTBestEnergyProb>0.1 and chainData.FswGamState == 0:
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
            trm.Fill()

    elif chainData.Cal1RawEnergySum>=20000 and chainData.Cal1MomZDir>=0.2 and e[0]>=vEnergyLow and e[0]<vEnergyUp and chainData.Cal1MomZCrossSide840>=0.0 and (chainData.TkrNumTracks==0 or (math.log10(max(chainData.CalTrackAngle,1E-4)) > (0.529795)*(e[0] < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((e[0]-3.000000)/0.916667,1))+(0.583401)*(pow((e[0]-3.000000)/0.916667,2))+(-0.075555)*(pow((e[0]-3.000000)/0.916667,3))))*(e[0]  >= 3.000000 and e[0] <= 5.750000) + (-0.398962)*(e[0] >  5.750000))) and chainData.Acd2Cal1VetoSigmaHit>0 and chainData.Cal1TransRms>=10 and chainData.Cal1TransRms<70 and chainData.Cal1MomNumIterations>0 and chainData.FswGamState == 0:
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
#        if -math.log10(1.0-p[0])>aaValCutBDT[binEnergy+nStartBin][0]:
#            if -math.log10(1.0-p[0])>aaValCutBDT[binEnergy+nStartBin][1]:
#                if -math.log10(1.0-p[0])>aaValCutBDT[binEnergy+nStartBin][2]:
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
            trm.Fill()
    if iEvent%(nEvent/200)==0:
        rate = int((iEvent*100.)/nEvent+0.5)
        if rate>0:
            nt = (datetime.datetime.now() - timeStart).seconds * (100.-rate)/rate
            meter = "\r[{0}{1}] {2} Wait {3} hr {4} min".format("=" * rate, ' ' * (100-rate), aaNumEventClass, int(nt/3600), (int(nt)%3600)/60+1)
        else:
            meter = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
        sys.stdout.write(meter)
        sys.stdout.flush()
print ""        
print "Finished!"
fileOut.cd()
trm.Write()
