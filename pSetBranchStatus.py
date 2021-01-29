#!/usr/bin/env python

import ROOT
import sys
ROOT.gROOT.SetBatch()
par = sys.argv 
print par

pathFile = par[1]

fileMel = ROOT.TFile(pathFile, "update")
print fileMel.GetName(), "is opened."
trMel = fileMel.Get("MeritTuple")
print trMel.GetName(), "is found."

trMel.SetBranchStatus("*", 0)
trMel.SetBranchStatus("WP8CalOnlyEnergy", 1)
trMel.SetBranchStatus("CalELayer7", 1)
trMel.SetBranchStatus("CalELayer4", 1)
trMel.SetBranchStatus("CalEdgeEnergy", 1)
trMel.SetBranchStatus("CalNewCfpCalSelChiSq", 1)
trMel.SetBranchStatus("Cal1TransRms", 1)
trMel.SetBranchStatus("CalNewCfpCalTmax", 1)
trMel.SetBranchStatus("CalBkHalfRatio", 1)
trMel.SetBranchStatus("Acd2Cal1Energy15", 1)
trMel.SetBranchStatus("Acd2VetoCount", 1)
trMel.SetBranchStatus("Acd2Cal1VetoSigmaHit", 1)
trMel.SetBranchStatus("CalTrSizeCalT95", 1)
trMel.SetBranchStatus("Cal1FitChiSquare", 1)
trMel.SetBranchStatus("Cal1MomNumCoreXtalsFract", 1)
trMel.SetBranchStatus("Acd2TileEnergy", 1)
trMel.SetBranchStatus("Cal1MomZDir", 1)
trMel.SetBranchStatus("GltGemSummary", 1)
trMel.SetBranchStatus("McSourceId", 1)
trMel.SetBranchStatus("EvtJointEnergy", 1)

trMel.Write()
fileMel.Close()
