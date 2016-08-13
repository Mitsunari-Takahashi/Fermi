#!/usr/bin/env python

import ROOT
import numpy as np
import sys
from pAnalysisConfig import *
ROOT.gROOT.SetBatch()
par = sys.argv
print par
strSuffix = par[1]
fileOut = ROOT.TFile(Form("McPsf{0}", strSuffix), 'RECREATE')
fileAG = ROOT.TFile("/nfs/farm/g/glast/u/mtakahas/data/MC/AG200909_62_2015Jan.root", 'READ')
trAG = fileAG.Get("MeritTuple")
if len(par)>2:
    strCut = par[2]
#pathFileRoc = par[1]
#fileRoc = ROOT.TFile(pathFilePerf, 'READ')
#pathFilePerf = par[2]
#filePerf = ROOT.TFile(pathFilePerf, 'READ')

#cutBase = ROOT.TCut(Form("", ))
#aCutEnergy = []
#aCutClass = []
#aCutClass.append(ROOT.TCut("CalOnlyR10", Form()))
#aCutClass.append(ROOT.TCut("CalOnlyR30", Form()))
#aCutClass.append(ROOT.TCut("CalOnlyR100", Form()))
#aCutType = ["FRONT", "BACK", "PSF0" "PSF1", "PSF2", "PSF3", "EDISP0", "EDISP1", "EDISP2", "EDISP3"]

#for cutClass in aCutClass:
    #trAG.Draw("ROOT.TMath.ACos(-Cal1MomXDir*McXDir-Cal1MomYDir*McYDir-Cal1MomZDir*McZDir) * ROOT.TMath.RadToDeg()>>(1800, 0.0, 180.0)", strCut)
vecMc = np.array([McXDir, McYDir, McZDir])
vecCal1Mom = np.array([Cal1MomXDir, Cal1MomYDir, Cal1MomZDir])
vecDiff = vecCal1Mom - vecMc
normVecDiff = np.linalg.norm(vecDiff)
h3PsfTotal = ROOT.TH3D("h3PsfTotal", "Direction misreconstruction;log_{10}E[MeV];cos(#theta);[deg]"; 7, 4.35, 5.75, 9, 0.1, 1.0, 180, 0, 18)
trAG.Draw("ROOT.TMath.ASin(normVecDiff/2.)*2. * ROOT.TMath.RadToDeg():-McZDir:McLogEnergy>>{0}".format(h3PsfTotal.GetName()), strCut)
fileOut.cd()
h3PsfTotal.Write()
h2PsfTotal_projTh = h3PsfTotal.Projection
