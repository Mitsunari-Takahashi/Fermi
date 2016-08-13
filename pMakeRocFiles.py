#!/usr/bin/env python

import sys
import numpy
from array import array
import ROOT
ROOT.gROOT.SetBatch()
import subprocess
par = sys.argv
strSuffix = par[1]
aZDIR = [0.1, 0.5, 0.9, 1.0]
nameFriend = par[1]
nameFriendLog = "{0}Log".format(nameFriend)
dictFile = {"BKG_0":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_0.root",
"BKG_1":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_1.root",
"BKG_2":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_2.root",
"BKG_3":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_3.root",
"BKG_4":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62MCE2e4_4.root",
"BKGp2":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62p2MCE2e4.root",
"BKGp3b":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62p3bMCE2e4.root",
"BKGb8":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62b8MCE2e4.root",
"BKGb9":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62b9MCE2e4.root",
"BKGb10":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62b10MCE2e4.root",
"BKGb11":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62b11MCE2e4.root",
"BKGb12":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62b12MCE2e4.root",
"BKGb13":"/nfs/farm/g/glast/u/mtakahas/data/MC/BKG200909_62b13MCE2e4.root"}
liKeys = dictFile.keys()
liValues = dictFile.values()
for iEl in range(len(dictFile)):
    for iDIR in range(len(aZDIR)-1):
        pathFileTgt = liValues[iEl]
        pathFileOut = 'roc_{0}_ZDIR{1}to{2}.yaml'.format(liKeys[iEl], int(100*aZDIR[iDIR]), int(100*aZDIR[iDIR+1]))
        f = open(pathFileOut,"w")
        strTgt = """mode : 'roc'
datasets : 
  sig : 
    file : '/nfs/farm/g/glast/u/mtakahas/data/MC/AG200909_62_2016Jun.root'
    friends :
    cuts : ' Cal1MomZDir>{0} && Cal1MomZDir<={1} && WP8CTCalOnlyBEPProb_E>=0.06 && Cal1MomZCrossSide840>=0.0 && (TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000)) ) && Cal1RawEnergySum>=20000 && Cal1MomNumIterations>0 && Cal1TransRms>=10 && Cal1TransRms<70 && Acd2Cal1VetoSigmaHit>0 && FswGamState==0 && (GltGemSummary&0x20)==0 '
  bkg :
    file : '{2}'
    friends :
    cuts : ' Cal1MomZDir>{0} && Cal1MomZDir<={1} && CRSample && WP8CTCalOnlyBEPProb_E>=0.06 && Cal1MomZCrossSide840>=0.0 && (TkrNumTracks==0 || (log10(max(CalTrackAngle,1E-4)) > (0.529795)*(EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((EvtJointLogEnergy-3.000000)/0.916667,3))))*(EvtJointLogEnergy >= 3.000000 && EvtJointLogEnergy <= 5.750000) + (-0.398962)*(EvtJointLogEnergy >  5.750000)) ) && Cal1RawEnergySum>=20000 && Cal1MomNumIterations>0 && Cal1TransRms>=10 && Cal1TransRms<70 && Acd2Cal1VetoSigmaHit>0 && FswGamState==0 && (GltGemSummary&0x20)==0 '

  default :
    fraction : 0.75
    friend_dir : 'data'
    aliases : 'pass8'
    friends : null
    alias_set:
        {3} : -log10(1.0-({4}+1.0)/2.0)
# Define a list of variables that will be used for calculating a ROC
# curve in bins of energy. 
rocVarList : # TO BE CHANGED!
   - '{3}'

coreVarList : 
  - 'WP8CTPSFCore' 
  - 'WP8CTPSFTail' 

# String prefix prepended to all output files and plots.
outprefix : '{5}_{6}' # TO BE CHANGED!

# Set the energy binning in log10(E/MeV) : min/max/bin_size 
ebins : '4.35/5.75/0.2'
evar : 'log10WP8CalOnlyEnergy'
""".format(aZDIR[iDIR], aZDIR[iDIR+1], pathFileTgt, nameFriendLog, nameFriend, liKeys[iEl], strSuffix)
        f.write(strTgt)
        f.close()
        subprocess.call(['python', '/u/gl/mtakahas/eventSelect/scripts/CutPerformance.py', '--config={0}'.format(pathFileOut)])
