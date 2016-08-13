#!/usr/bin/env python

import sys
import ROOT
from ROOT import TH2

par = sys.argv

aPathFileIn = par[1:]
nOff = 4
pathFileOut = "ExposureSummedOFF.root"
if len(aPathFileIn)==1:
    pathFileOut = aPathFileIn[0]
    pathFileOut.replace('.root', 'OFF.root')
fileOut = ROOT.TFile(pathFileOut, 'RECREATE')
#fileSed = ROOT.TFile("hSED09-13.root", 'READ')
#hSED = fileSed.Get('hSED')

htgExpOff = ROOT.TH2D("htgExpOff", "Summed OFF exposure", 7, 4.35, 5.75, 180, 0, 180)
for pathFileIn in aPathFileIn:
    fileIn = ROOT.TFile(pathFileIn, 'READ')
    for iOff in range(nOff):
        htgExpTemp = fileIn.Get("htgExp{0}".format(iOff+1))
        htgExpOff.Add(htgExpTemp)
fileOut.cd()
htgExpOff.Write()
