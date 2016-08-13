#!/usr/bin/env python

import sys
from array import array
import numpy as np
import ROOT
ROOT.gROOT.SetBatch()
from ROOT import TTree, TFile, gROOT
par = sys.argv

cat = np.zeros(1, dtype=int)
at = np.zeros(1, dtype=float)
ene = np.zeros(1, dtype=float)
zen = np.zeros(1, dtype=float)
sep = np.zeros(1, dtype=float)
cla = np.zeros(1, dtype=int)
typ = np.zeros(1, dtype=int)
fl = np.zeros(1, dtype=int)

listFileIn = par[1:]
for nameFileIn in listFileIn:
    print "=========="
    print nameFileIn
    fileIn = ROOT.TFile(nameFileIn, 'READ')
    print fileIn.GetName(), "is opened."
    tr = fileIn.Get("trGammas")
    tr.SetBranchAddress("c",cat)
    tr.SetBranchAddress("grbt",at)
    tr.SetBranchAddress("e",ene)
    tr.SetBranchAddress("dist",sep)
    tr.SetBranchAddress("z",zen)
    tr.SetBranchAddress("s",cla)
    tr.SetBranchAddress("ty",typ)
    tr.SetBranchAddress("flag",fl)
    print tr.GetName(), "is found."
    if tr.GetEntries()>0:
        #for evt in tr:
        for iEvt in range(tr.GetEntries()):
            tr.GetEntry(iEvt)
            if fl[0]==0 and at[0]>=0.:
#if tr.FLAG==0 and tr.GRB_TIME>=0.:
                print "----------"
                print "Category:", cat[0]
                print "Arrival time:", at[0]
                print "Energy:", ene[0]
                print "Zenith:", zen[0]
                print "Angular separation:", sep[0]
                print "Event class:", cla[0]
                print "Event type:", typ[0]
#                 print "Category:", evt.Category
#                 print "Arrival time:", evt.GRB_TIME
#                 print "Energy:", evt.ENERGY
#                 print "Zenith:", evt.ZENITH_ANGLE
#                 print "Angular separation:", evt.DIST
#                 print "Event class:", evt.EVENT_CLASS
#                 print "Event type:", evt.EVENT_TYPE
                
