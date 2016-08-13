#!/usr/bin/env python

import sys
import ROOT
from ROOT import TTree
from ROOT import TChain
import numpy
import yaml
import datetime
sys.path.append("/disk/gamma/cta/store/takhsm/FermiMVA/AllSky")
sys.path.append("/home/takhsm/FermiMVA/python")
ROOT.gROOT.SetBatch()
from array import array
import math

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

nameFileRoc = par[2]
nameVarBDT = "S11V200909_020RAWE30ZCS000wwoTRKwBKwoMCZDIR00woRWdivCALE_15_BDTG1000D06"
nameVarBDT = par[1]
nameFileSuffix = par[3]
#aCutEGB = [10, 2, 1]
#cutMVA = CutBDT("/home/takhsm/FermiMVA/S10/S10V200909_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15/v20r09p09_G1haB1_S10_020rawe30zdir020nbep006WWOtrkWbkWOmczWOrw_15_A_S10WOrw_BDTG1000D06Log_roc.root", aCutEGB)
cutMVA = CutBDT(nameFileRoc, aCutEGB)
aaValCutBDT = cutMVA.aaValCutBDT[0:] #cutBDT(nameFileRoc, aCutEGB)
print aaValCutBDT

print "===================="
# Making all sky map
#listFileIn = ["/disk/gamma/cta/store/takhsm/FermiMVA/AllSky/LPA_P301_15Feb_2009_1.root"]
listFileIn = par[4:]
print listFileIn
chainData = ROOT.TChain("MeritTuple")
for nameFileIn in listFileIn:
    chainData.Add(nameFileIn)

chainFriend = ROOT.TChain("MeritTuple")
#chainFriendBEP = ROOT.TChain("MeritTuple")
for nameFileIn in listFileIn:
    nameFileFriend = nameFileIn.replace(".root", "_" + nameVarBDT + ".root")
    #nameFileFriendBEP = nameFileIn.replace(".root", "_WP8CalOnlyBEPCaseE_myBDT.root")
    print nameFileFriend
    chainFriend.Add(nameFileFriend)
    #print nameFileFriendBEP
    #chainFriendBEP.Add(nameFileFriendBEP)

chainData.AddFriend(chainFriend, "friendTemp=MeritTuple")
print "Friend: ", chainData.GetListOfFriends().FindObject("friendTemp=MeritTuple").Print()
print "---------------------"
#chainData.AddFriend(chainFriendBEP, "friendBEP=MeritTuple")
#print "BEP friend: ", chainData.GetListOfFriends().FindObject("friendBEP=MeritTuple").Print()
#print "---------------------"

aliasSelections = yaml.load(open('/home/takhsm/FermiMVA/python/pass8_event_selections.yaml','r'))
for k,v in aliasSelections.iteritems(): 
    chainData.SetAlias(k,v)

#fileRoc = ROOT.TFile(nameFileRoc, 'READ')
#h2Sig = fileRoc.Get('sig_acc')
nEnergyBin = cutMVA.aaValCutBDT[0]['numBin'] - nStartBin #h2Sig.ProjectionX().GetNbinsX()-nStartBin
vEnergyBinWidth = cutMVA.aaValCutBDT[0]['widthBin'] #h2Sig.GetXaxis().GetBinWidth(1)
vEnergyLow = cutMVA.aaValCutBDT[0]['edgeLow'] + nStartBin*vEnergyBinWidth #h2Sig.GetXaxis().GetBinLowEdge(1+nStartBin)
vEnergyUp = vEnergyLow + nEnergyBin*vEnergyBinWidth #h2Sig.GetXaxis().GetXmax()

nameFileOut = "trAllSkyMap_S11V200909_020RAWE30ZCS000wwoTRKwBKwoMCZDIR00woRWdivCALE_15_BDTG1000D06_LPA_P301_15Feb_2009_2.root"
#nameFileOut = par[1]
nameFileOut = "trAllSkyMap_" + nameVarBDT + "_" + nameFileSuffix + ".root"
fileOut = ROOT.TFile(nameFileOut, 'UPDATE')
#aaaGrMap = []
aaNumEventClass=[]
aColor = []

#------ TTree setting -----
trm = ROOT.TTree("trGammas", "Gamma-like events")
c = numpy.zeros(1, dtype=int)
s = numpy.zeros(1, dtype=int)
e = numpy.zeros(1, dtype=float)
t = numpy.zeros(1, dtype=float)
ra = numpy.zeros(1, dtype=float)
dec = numpy.zeros(1, dtype=float)
l = numpy.zeros(1, dtype=float)
b = numpy.zeros(1, dtype=float)
z = numpy.zeros(1, dtype=float)
bep = numpy.zeros(1, dtype=float)
cth = numpy.zeros(1, dtype=float)
trm.Branch('Category',c,'c/I') # 1:CalTkr or 2:CalOnly
trm.Branch('Class',s,'s/I') # 1: TRANSIENT, CalOnly_10xEGB, 2: CalOnly_2xEGB, 3: SOURCE, CalOnly_1xEGB
trm.Branch('Energy',e,'e/D')
trm.Branch('Time',t,'t/D')
trm.Branch('RightAscension',ra,'ra/D')
trm.Branch('Declination',dec,'dec/D')
trm.Branch('Longitude',l,'l/D')
trm.Branch('Latitude',b,'b/D')
trm.Branch('ZenithAngle',z,'z/D')
trm.Branch('BestEnergyProb',bep,'bep/D')
trm.Branch('CosInclinationAngle',cth,'cth/D')

for hS in range(len(aaStrSelect)):
    aaNumEventClass.append([])
    for iS in range(len(aaStrSelect[hS])):
        aaNumEventClass[hS].append(0)

nEvent = chainData.GetEntries()
print "Total number of events:", nEvent
timeStart = datetime.datetime.now()
print timeStart

for iEvent in range(nEvent):
    chainData.GetEntry(iEvent)
    e[0] = chainData.EvtJointLogEnergy
    c[0] = 0
    s[0] = 0
    t[0] = chainData.EvtElapsedTime
    bep[0] = (chainData.WP8CalOnlyBEPCaseE_myBDT+1)/2.0
    binEnergy = max(min(nEnergyBin-1, int((e[0]-vEnergyLow)/vEnergyBinWidth * (int(e[0]<vEnergyLow)*(-2)+1)) ), 0)
    if (chainData.TkrNumTracks>0) and (math.log10(max(chainData.CalTrackAngle,1E-4)) <= (0.529795)*(e[0] < 3.000000) + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((e[0]-3.000000)/0.916667,1))+(0.583401)*(pow((e[0]-3.000000)/0.916667,2))+(-0.075555)*(pow((e[0]-3.000000)/0.916667,3))))*(e[0] >= 3.000000 and e[0] <= 5.750000) + (-0.398962)*(e[0] > 5.750000)) and chainData.EvtCalCsIRLn>4 and chainData.WP8CTPSFTail>0.05 and chainData.WP8CTBestEnergyProb>0.1 and chainData.FswGamState == 0:
        c[0] = 1
        z[0] = chainData.FT1ZenithTheta
        ra[0] = chainData.FT1Ra
        dec[0] = chainData.FT1Dec
        l[0] = chainData.FT1L
        b[0] = chainData.FT1B
        cth[0] = chainData.Cal1MomZDir
        if ( -math.log10(1.0-chainData.WP8CTAllProb) >= (0.010000)*(e[0] < 1.250000) + ((e[0] <= 1.750000)*((0.010000)*(1.0)+(0.000000)*(math.pow((e[0]-1.250000)/0.500000,1))+(0.018669)*(math.pow((e[0]-1.250000)/0.500000,2)))+((e[0] > 1.750000)*(e[0] <= 2.250000))*((0.028669)*(1.0)+(0.037338)*(math.pow((e[0]-1.750000)/0.500000,1))+(-0.017111)*(math.pow((e[0]-1.750000)/0.500000,2)))+((e[0] > 2.250000)*(e[0] <= 2.750000))*((0.048897)*(1.0)+(0.003117)*(math.pow((e[0]-2.250000)/0.500000,1))+(0.001967)*(math.pow((e[0]-2.250000)/0.500000,2)))+((e[0] > 2.750000)*(e[0] <= 3.250000))*((0.053980)*(1.0)+(0.007050)*(math.pow((e[0]-2.750000)/0.500000,1))+(-0.003525)*(math.pow((e[0]-2.750000)/0.500000,2)))+((e[0] > 3.250000)*(e[0] <= 3.750000))*((0.057505)*(1.0)+(0.000000)*(math.pow((e[0]-3.250000)/0.500000,1))+(0.121963)*(math.pow((e[0]-3.250000)/0.500000,2)))+((e[0] > 3.750000)*(e[0] <= 4.250000))*((0.179468)*(1.0)+(0.243925)*(math.pow((e[0]-3.750000)/0.500000,1))+(0.493075)*(math.pow((e[0]-3.750000)/0.500000,2)))+((e[0] > 4.250000)*(e[0] <= 4.750000))*((0.916468)*(1.0)+(1.230076)*(math.pow((e[0]-4.250000)/0.500000,1))+(-0.501532)*(math.pow((e[0]-4.250000)/0.500000,2)))+(e[0] > 4.750000)*((1.645012)*(1.0)+(0.227011)*(math.pow((e[0]-4.750000)/0.500000,1))+(0.029483)*(math.pow((e[0]-4.750000)/0.500000,2))))*(e[0] >= 1.250000 and e[0] <= 5.750000) + (2.216967)*(e[0] > 5.750000) ): #P8R1_TRANSIENT_R100
            if ( -math.log10(1.0-chainData.WP8CTAllProb) >= (0.080914)*(e[0] < 1.250000) + ((e[0] <= 1.750000)*((0.080914)*(1.0)+(0.108897)*(pow((e[0]-1.250000)/0.500000,1))+(0.377870)*(pow((e[0]-1.250000)/0.500000,2)))+((e[0] > 1.750000)*(e[0] <= 2.250000))*((0.567682)*(1.0)+(0.864637)*(pow((e[0]-1.750000)/0.500000,1))+(-0.182318)*(pow((e[0]-1.750000)/0.500000,2)))+((e[0] > 2.250000)*(e[0] <= 2.750000))*((1.250000)*(1.0)+(0.500000)*(pow((e[0]-2.250000)/0.500000,1))+(-0.085000)*(pow((e[0]-2.250000)/0.500000,2)))+((e[0] > 2.750000)*(e[0] <= 3.250000))*((1.665000)*(1.0)+(0.330000)*(pow((e[0]-2.750000)/0.500000,1))+(-0.165000)*(pow((e[0]-2.750000)/0.500000,2)))+((e[0] > 3.250000)*(e[0] <= 3.750000))*((1.830000)*(1.0)+(0.000000)*(pow((e[0]-3.250000)/0.500000,1))+(0.285000)*(pow((e[0]-3.250000)/0.500000,2)))+((e[0] > 3.750000)*(e[0] <= 4.250000))*((2.115000)*(1.0)+(0.570000)*(pow((e[0]-3.750000)/0.500000,1))+(-0.185000)*(pow((e[0]-3.750000)/0.500000,2)))+((e[0] > 4.250000)*(e[0] <= 4.750000))*((2.500000)*(1.0)+(0.200000)*(pow((e[0]-4.250000)/0.500000,1))+(0.100000)*(pow((e[0]-4.250000)/0.500000,2)))+(e[0] > 4.750000)*((2.800000)*(1.0)+(0.400000)*(pow((e[0]-4.750000)/0.500000,1))+(-0.112171)*(pow((e[0]-4.750000)/0.500000,2))))*(e[0] >= 1.250000 and e[0] <= 5.750000) + (3.151318)*(e[0] > 5.750000) ) and ( chainData.WP8CTAllBkProb >= (0.366167)*(e[0] < 1.250000) + ((e[0] <= 1.541667)*((0.366167)*(1.0)+(0.028500)*(pow((e[0]-1.250000)/0.291667,1))+(-0.056500)*(pow((e[0]-1.250000)/0.291667,2))+(0.106667)*(pow((e[0]-1.250000)/0.291667,3)))+((e[0] > 1.541667)*(e[0] <= 1.833333))*((0.444833)*(1.0)+(0.235500)*(pow((e[0]-1.541667)/0.291667,1))+(0.263500)*(pow((e[0]-1.541667)/0.291667,2))+(-0.162667)*(pow((e[0]-1.541667)/0.291667,3)))+((e[0] > 1.833333)*(e[0] <= 2.125000))*((0.781167)*(1.0)+(0.274500)*(pow((e[0]-1.833333)/0.291667,1))+(-0.224500)*(pow((e[0]-1.833333)/0.291667,2))+(0.072667)*(pow((e[0]-1.833333)/0.291667,3)))+(e[0] > 2.125000)*((0.903833)*(1.0)+(0.043500)*(pow((e[0]-2.125000)/0.291667,1))+(-0.006500)*(pow((e[0]-2.125000)/0.291667,2))+(-0.000333)*(pow((e[0]-2.125000)/0.291667,3))))*(e[0] >= 1.250000 and e[0] <= 3.000000) + (0.966833)*(e[0] > 3.000000) ):  #P8R1_SOURCE_AllProbFilter&&P8R1_SOURCE_AllBkProbFilter
                s[0] = 3
                aaNumEventClass[0][1] = aaNumEventClass[0][1]+1
            else:
                s[0] = 1
                aaNumEventClass[0][0] = aaNumEventClass[0][0]+1
            trm.Fill()

    elif chainData.Cal1RawEnergySum>30000 and chainData.Cal1MomZDir>0.2 and (chainData.WP8CalOnlyBEPCaseE_myBDT+1.0)/2.0>0.06 and (chainData.TkrNumTracks==0 or (math.log10(max(chainData.CalTrackAngle,1E-4)) > (0.529795)*(e[0] < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((e[0]-3.000000)/0.916667,1))+(0.583401)*(pow((e[0]-3.000000)/0.916667,2))+(-0.075555)*(pow((e[0]-3.000000)/0.916667,3))))*(e[0]  >= 3.000000 and e[0] <= 5.750000) + (-0.398962)*(e[0] >  5.750000))) and chainData.Acd2Cal1VetoSigmaHit>0 and chainData.Cal1TransRms>10 and chainData.Cal1TransRms<70 and chainData.Cal1MomNumIterations>0 and chainData.FswGamState == 0:
        c[0] = 2
        z[0] = chainData.FT1CalZenithTheta
        ra[0] = chainData.FT1CalRa
        dec[0] = chainData.FT1CalDec
        l[0] = chainData.FT1CalL
        b[0] = chainData.FT1CalB
        cth[0] = chainData.Cal1MomZDir
        if -math.log10(1.0-(1.0+chainData.S11V200909_020RAWE30ZCS000wwoTRKwBKwoMCZDIR00woRWdivCALE_15_BDTG1000D06)/2.0)>aaValCutBDT[binEnergy+nStartBin][0]:
            if -math.log10(1.0-(1.0+chainData.S11V200909_020RAWE30ZCS000wwoTRKwBKwoMCZDIR00woRWdivCALE_15_BDTG1000D06)/2.0)>aaValCutBDT[binEnergy+nStartBin][1]:
                if -math.log10(1.0-(1.0+chainData.S11V200909_020RAWE30ZCS000wwoTRKwBKwoMCZDIR00woRWdivCALE_15_BDTG1000D06)/2.0)>aaValCutBDT[binEnergy+nStartBin][2]:
                    s[0]=3
                    aaNumEventClass[1][2] = aaNumEventClass[1][2]+1
                else:
                    s[0] = 2
                    aaNumEventClass[1][1] = aaNumEventClass[1][1]+1
            else:
                s[0] = 1
                aaNumEventClass[1][0] = aaNumEventClass[1][0]+1
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
