#!/usr/bin/env python

import sys
import ROOT
from ROOT import TH2, TH1

par = sys.argv
pathFileExp = 'ExposureSummedOFF.root'
pathFileOut = 'ExposureWeightedOFF.root'
if len(par)>1:
    pathFileExp = par[1]
    pathFileOut = pathFileExp
    pathFileOut = pathFileOut.replace('.root', 'Weighted.root')
strHtgExpName = 'htgExpOff'
if len(par)>2:
    strHtgExpName = par[2]
fileSed = ROOT.TFile("hSED09-13.root", 'READ')
htgSED = fileSed.Get('hSED')
print htgSED.GetName()
fileExp = ROOT.TFile(pathFileExp, 'READ')
htgExpOffSum = fileExp.Get(strHtgExpName)
print htgExpOffSum.GetName()
htgExpOffWeight = ROOT.TH1D('{0}Weight'.format(htgExpOffSum.GetName()), 'Weighted {0}'.format(htgExpOffSum.GetTitle()), 180, 0, 180)
for iE in range(htgExpOffSum.GetNbinsX()):
    print iE+1
    htgExpProj = htgExpOffSum.ProjectionY("_py{0}".format(iE+1), iE+1, iE+1)
    weight = htgSED.GetBinContent(iE+1)/htgSED.Integral() #/htgExpProj.Integral()
    print "Weight:", weight
    htgExpOffWeight.Add(htgExpProj, weight)
fileOut = ROOT.TFile(pathFileOut, 'RECREATE')
fileOut.cd()
htgExpOffWeight.Write()
