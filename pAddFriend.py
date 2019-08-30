#!/usr/bin/env python

import ROOT
import sys
ROOT.gROOT.SetBatch()
par = sys.argv 
print par
aPathFile = par[1:]

aSuffixFriend=["U01_V200909_62IRFBK_020RAWE20woRW_09", "newVariables", "WP8CalOnlyBEPCaseE_myBDT", "CalOnlyVar"]
aNameFriend=["UpDown", "newVariables", "WP8CalOnlyBEP_E_BDT", "CalOnlyVar"]
if len(aSuffixFriend)!=len(aNameFriend):
    print "Numbers of the name and the suffix of the friend files are different!!"
    sys.exit(1)
nFriend = len(aSuffixFriend)
for pathFile in aPathFile:
    if pathFile[0] is not "/":
        print "Please input ablosute path!"
        sys.exit(1)
    fileMel = ROOT.TFile(pathFile, "update")
    print fileMel.GetName(), "is opened."
    trMel = fileMel.Get("MeritTuple")
    print trMel.GetName(), "is found."
    #trMel.GetListOfFriends().RemoveAll()
    aPathFriend=[]
    for iFr in range(nFriend):
        aPathFriend.append(pathFile[:-5])
        aPathFriend[-1] = aPathFriend[-1] + "_" + aSuffixFriend[iFr] + ".root"
        print "  Friend:", aPathFriend[-1]
        strAddFriend = aNameFriend[iFr] + " = MeritTuple"
        trMel.AddFriend(strAddFriend, aPathFriend[-1])
    listFr = trMel.GetListOfFriends()
    if listFr is None:
        print "No added friends!"
        sys.exit(1)
    else:
        print listFr.Print
        trMel.SetAlias("WP8CTCalOnlyBEPProb_E", "(WP8CalOnlyBEPCaseE_myBDT+1)/2.0")
        trMel.Write()
    fileMel.Close()
