#!/usr/bin/env python
import sys
import numpy
from math import log10
from array import array
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
from array import array
ROOT.gROOT.SetBatch() 

par = sys.argv
print par
aPathFiles = par[1:]

#newVar0 = array('f',[0.0])
newVar1 = array('f',[0.0])
newVar2 = array('f',[0.0])
newVar3 = array('i',[0])
newVar4 = array('f',[0.0])
newVar5 = array('f',[0.0])
newVar6 = array('f',[0.0])
newVar7 = array('f',[0.0])
#newVar0Cor = array('f',[0.0])
newVar1Cor = array('f',[0.0])
newVar2Cor = array('f',[0.0])
newVar3Cor = array('i',[0])
newVar4Cor = array('f',[0.0])
newVar5Cor = array('f',[0.0])
newVar6Cor = array('f',[0.0])
numLayerIn = -1
calZTop =  -47.395
cellVertPitch = 21.35

for pathFile in aPathFiles:
    fileIn = ROOT.TFile(pathFile, "READ")
    print fileIn.GetName(), "is opened."
    trDat = fileIn.Get("MeritTuple")
    nEvt = trDat.GetEntries()
    print trDat.GetName(), "has", nEvt, "events."
    pathFileFr = pathFile[0:-5] + "_newVariables.root"
#    pathFileFr.replace(".root", "_newVariables.root")
    if pathFileFr == pathFile:
        print pathFileFr
        sys.exit(-1)
    fileFr = ROOT.TFile(pathFileFr, "RECREATE")
    print fileFr.GetName(), "is opened."
    trFr = ROOT.TTree("MeritTuple", "Additional glast tuple")
    #trFr.Branch('Cal1MomParticleInLayer', newVar0, 'Cal1MomParticleInLayer/I')
    trFr.Branch('Cal1MomZCrossSide840', newVar1, 'Cal1MomZCrossSide840/F')
    trFr.Branch('Cal1MomZCrossSide714', newVar2, 'Cal1MomZCrossSide714/F')
    trFr.Branch('Cal1MomZCrossSide840Cor', newVar1Cor, 'Cal1MomZCrossSide840Cor/F')
    trFr.Branch('Cal1MomZCrossSide714Cor', newVar2Cor, 'Cal1MomZCrossSide714Cor/F')
    trFr.Branch('Cal1MomParticleInLayer', newVar3, 'Cal1MomParticleInLayer/I')
    trFr.Branch('CalELayerCorrInitialRatioLog', newVar4, 'CalELayerCorrInitialRatioLog/F')
    trFr.Branch('CalELayer34afterInitialRatioLog', newVar5, 'CalELayer34afterInitialRatioLog/F')
    trFr.Branch('CalELayerInitialRawRatioLog', newVar6, 'CalELayerInitialRawRatioLog/F')
    trFr.Branch('Cal1MomParticleInLayerCor', newVar3Cor, 'Cal1MomParticleInLayerCor/I')
    trFr.Branch('CalELayerCorrInitialRatioCorLog', newVar4Cor, 'CalELayerCorrInitialRatioCorLog/F')
    trFr.Branch('CalELayer34afterInitialRatioCorLog', newVar5Cor, 'CalELayer34afterInitialRatioCorLog/F')
    trFr.Branch('CalELayerInitialRawRatioCorLog', newVar6Cor, 'CalELayerInitialRawRatioCorLog/F')
    trFr.Branch('Acd2TileEnergyRatioLog', newVar7, 'Acd2TileEnergyRatioLog/F')

    for iEvt in range(nEvt):
        trDat.GetEntry(iEvt)
        if trDat.Cal1MomXDir!=0 and trDat.Cal1MomYDir!=0:
            newVar1[0] = (trDat.Cal1MomZCntr+min((840-((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXCntr)/(((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXDir), (840-((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYCntr)/(((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYDir))*trDat.Cal1MomZDir) 
            newVar2[0] = (trDat.Cal1MomZCntr+min((714-((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXCntr)/(((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXDir), (714-((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYCntr)/(((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYDir))*trDat.Cal1MomZDir)
            newVar1Cor[0] = (trDat.Cal1MomZCntrCor+min((840-((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXCntrCor)/(((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXDir), (840-((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYCntrCor)/(((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYDir))*trDat.Cal1MomZDir) 
            newVar2Cor[0] = (trDat.Cal1MomZCntrCor+min((714-((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXCntrCor)/(((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXDir), (714-((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYCntrCor)/(((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYDir))*trDat.Cal1MomZDir)
        elif trDat.Cal1MomXDir!=0:
            newVar1[0] = (trDat.Cal1MomZCntr+(840-((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXCntr)/(((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXDir)*trDat.Cal1MomZDir)
            newVar2[0] = (trDat.Cal1MomZCntr+(714-((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXCntr)/(((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXDir)*trDat.Cal1MomZDir)
            newVar1Cor[0] = (trDat.Cal1MomZCntrCor+(840-((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXCntrCor)/(((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXDir)*trDat.Cal1MomZDir)
            newVar2Cor[0] = (trDat.Cal1MomZCntrCor+(714-((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXCntrCor)/(((trDat.Cal1MomXDir>0)*2-1)*trDat.Cal1MomXDir)*trDat.Cal1MomZDir)
        elif trDat.Cal1MomYDir!=0:
            newVar1[0] = (trDat.Cal1MomZCntr+(840-((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYCntr)/(((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYDir)*trDat.Cal1MomZDir)
            newVar2[0] = (trDat.Cal1MomZCntr+(714-((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYCntr)/(((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYDir)*trDat.Cal1MomZDir)
            newVar1Cor[0] = (trDat.Cal1MomZCntrCor+(840-((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYCntrCor)/(((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYDir)*trDat.Cal1MomZDir)
            newVar2Cor[0] = (trDat.Cal1MomZCntrCor+(714-((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYCntrCor)/(((trDat.Cal1MomYDir>0)*2-1)*trDat.Cal1MomYDir)*trDat.Cal1MomZDir)
        else:
            newVar1[0] = 1000000
            newVar2[0] = 1000000
            newVar1Cor[0] = 1000000
            newVar2Cor[0] = 1000000
            print ""
            print "Perpendicular event!"
        numLayerIn = min(7.0, max(0.0, -(newVar2[0]-calZTop)/cellVertPitch))
        newVar3[0] = int(numLayerIn)
        newVar4[0] = log10(max(1E-5,trDat.CalEnergyCorr)) - ( ( numLayerIn<1 )*log10(max(1E-5, trDat.CalELayer0)) + ( numLayerIn>=1 )*( numLayerIn<2 )*log10(max(1E-5, trDat.CalELayer1)) + ( numLayerIn>=2 )*( numLayerIn<3 )*log10(max(1E-5, trDat.CalELayer2)) + ( numLayerIn>=3 )*( numLayerIn<4 )*log10(max(1E-5, trDat.CalELayer3)) + ( numLayerIn>=4 )*( numLayerIn<5 )*log10(max(1E-5, trDat.CalELayer4)) + ( numLayerIn>=5 )*( numLayerIn<6 )*log10(max(1E-5, trDat.CalELayer5)) + ( numLayerIn>=6 )*( numLayerIn<7 )*log10(max(1E-5, trDat.CalELayer6)) )

        newVar5[0] = ( numLayerIn<1 ) * ( log10(max(1E-5, trDat.CalELayer4+trDat.CalELayer3))-log10(max(1E-5, trDat.CalELayer0)) ) + ( numLayerIn>=1 ) * ( numLayerIn<2 ) * ( log10(max(1E-5, trDat.CalELayer5+trDat.CalELayer4))-log10(max(1E-5, trDat.CalELayer1)) ) + ( numLayerIn>=2 ) * ( numLayerIn<3 ) * ( log10(max(1E-5, trDat.CalELayer6+trDat.CalELayer5))-log10(max(1E-5, trDat.CalELayer2)) ) + ( numLayerIn>=3 ) * ( log10(max(1E-5, trDat.CalELayer7+trDat.CalELayer6))-log10(max(1E-5, trDat.CalELayer3)) )

        newVar6[0] = log10(max(1E-5,trDat.CalEnergyRaw)) - ( ( numLayerIn<1 )*log10(max(1E-5, trDat.CalELayer0)) + ( numLayerIn>=1 )*( numLayerIn<2 )*log10(max(1E-5, trDat.CalELayer1)) + ( numLayerIn>=2 )*( numLayerIn<3 )*log10(max(1E-5, trDat.CalELayer2)) + ( numLayerIn>=3 )*( numLayerIn<4 )*log10(max(1E-5, trDat.CalELayer3)) + ( numLayerIn>=4 )*( numLayerIn<5 )*log10(max(1E-5, trDat.CalELayer4)) + ( numLayerIn>=5 )*( numLayerIn<6 )*log10(max(1E-5, trDat.CalELayer5)) + ( numLayerIn>=6 )*( numLayerIn<7 )*log10(max(1E-5, trDat.CalELayer6)) )

        numLayerInCor = min(7.0, max(0.0, -(newVar2Cor[0]-calZTop)/cellVertPitch))
        newVar3Cor[0] = int(numLayerInCor)

        newVar4Cor[0] = log10(max(1E-5,trDat.CalEnergyCorr)) - ( ( numLayerInCor<1 )*log10(max(1E-5, trDat.CalELayer0)) + ( numLayerInCor>=1 )*( numLayerInCor<2 )*log10(max(1E-5, trDat.CalELayer1)) + ( numLayerInCor>=2 )*( numLayerInCor<3 )*log10(max(1E-5, trDat.CalELayer2)) + ( numLayerInCor>=3 )*( numLayerInCor<4 )*log10(max(1E-5, trDat.CalELayer3)) + ( numLayerInCor>=4 )*( numLayerInCor<5 )*log10(max(1E-5, trDat.CalELayer4)) + ( numLayerInCor>=5 )*( numLayerInCor<6 )*log10(max(1E-5, trDat.CalELayer5)) + ( numLayerInCor>=6 )*( numLayerInCor<7 )*log10(max(1E-5, trDat.CalELayer6)) )

        newVar5Cor[0] =  ( numLayerInCor<1 ) * ( log10(max(1E-5, trDat.CalELayer4+trDat.CalELayer3))-log10(max(1E-5, trDat.CalELayer0)) ) + ( numLayerInCor>=1 ) * ( numLayerInCor<2 ) * ( log10(max(1E-5, trDat.CalELayer5+trDat.CalELayer4))-log10(max(1E-5, trDat.CalELayer1)) ) + ( numLayerInCor>=2 ) * ( numLayerInCor<3 ) * ( log10(max(1E-5, trDat.CalELayer6+trDat.CalELayer5))-log10(max(1E-5, trDat.CalELayer2)) ) + ( numLayerInCor>=3 ) * ( log10(max(1E-5, trDat.CalELayer7+trDat.CalELayer6))-log10(max(1E-5, trDat.CalELayer3)) )
        
        newVar6Cor[0] = log10(max(1E-5,trDat.CalEnergyRaw)) - ( ( numLayerInCor<1 )*log10(max(1E-5, trDat.CalELayer0)) + ( numLayerInCor>=1 )*( numLayerInCor<2 )*log10(max(1E-5, trDat.CalELayer1)) + ( numLayerInCor>=2 )*( numLayerInCor<3 )*log10(max(1E-5, trDat.CalELayer2)) + ( numLayerInCor>=3 )*( numLayerInCor<4 )*log10(max(1E-5, trDat.CalELayer3)) + ( numLayerInCor>=4 )*( numLayerInCor<5 )*log10(max(1E-5, trDat.CalELayer4)) + ( numLayerInCor>=5 )*( numLayerInCor<6 )*log10(max(1E-5, trDat.CalELayer5)) + ( numLayerInCor>=6 )*( numLayerInCor<7 )*log10(max(1E-5, trDat.CalELayer6)) )

        newVar7[0] = log10(max(trDat.Acd2TileEnergy/max(10., trDat.EvtJointEnergy) * 100, 1E-5))

        trFr.Fill();
        if iEvt%(nEvt/200)==0:
            rate = int((iEvt*100.)/nEvt+0.5)
            if rate>0:
                meter = "\r[{0}{1}]({2:.2f})".format("=" * rate, ' ' * (100-rate), numLayerIn)
                sys.stdout.write(meter)
                sys.stdout.flush()
    print ""
    print nEvt, "events have been filled."
    fileFr.cd()
    trFr.Write()
    fileFr.Close()
    fileIn.Close()
