#!/usr/bin/env python
import os
import sys
import numpy
from math import log10
from array import array
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
from ctypes import *
ROOT.gROOT.SetBatch() 

par = sys.argv
print par
strForce = par[1]
if strForce=='True':
    print "Overwriting is allowed."
else:
    print "Overwriting is NOT allowed."

aPathFiles = par[2:]

#newVar0 = array('f',[0.0])
newVar1 = c_bool() #array('f',[0.0])
newVar2 = c_bool() 

for pathFile in aPathFiles:
    fileIn = ROOT.TFile(pathFile, "READ")
    print fileIn.GetName(), "is opened."
    trDat = fileIn.Get("MeritTuple")
    nEvt = trDat.GetEntries()
    print trDat.GetName(), "has", nEvt, "events."
    pathFileFr = pathFile[0:-5] + "_CalOnlyVar.root"
    if strForce=='False' and os.path.isfile(pathFileFr)==True:
        print pathFileFr, 'already exits.'
        break
    if pathFileFr == pathFile:
        print pathFileFr
        sys.exit(-1)
    fileFr = ROOT.TFile(pathFileFr, "RECREATE")
    print fileFr.GetName(), "is opened."
    trFr = ROOT.TTree("MeritTuple", "Additional glast tuple")
    trFr.Branch('NoTrk', newVar1, 'NoTrk/O')
    trFr.Branch('CalOnly', newVar2, 'CalOnly/O')

    for iEvt in range(nEvt):
        trDat.GetEntry(iEvt)
        newVar1.value = trDat.TkrNumTracks==0 or (log10(max(trDat.CalTrackAngle,1E-4)) > ( (0.529795)*(trDat.EvtJointLogEnergy < 3.000000) + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((trDat.EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((trDat.EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((trDat.EvtJointLogEnergy-3.000000)/0.916667,3))))*(trDat.EvtJointLogEnergy >= 3.000000) * (trDat.EvtJointLogEnergy <= 5.750000) + (-0.398962)*(trDat.EvtJointLogEnergy > 5.750000) ))
        newVar2.value = newVar1.value and trDat.Cal1RawEnergySum>=20000 and trDat.Cal1MomNumIterations>0 and trDat.FswGamState==0
        trFr.Fill();
        if iEvt%(nEvt/200)==0:
            rate = int((iEvt*100.)/nEvt+0.5)
            if rate>0:
                meter = "\r[{0}{1}]".format("=" * rate, ' ' * (100-rate))
                sys.stdout.write(meter)
                sys.stdout.flush()
    print ""
    print nEvt, "events have been filled."
    fileFr.cd()
    trFr.Write()
    fileFr.Close()
    fileIn.Close()
