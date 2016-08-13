#!/usr/bin/env python

import ROOT
import sys
ROOT.gROOT.SetBatch()
par = sys.argv #args # 
print par
aPathFile = par[3:]
nFriend = len(par)-2
print nFriend
aSuffixFriend=par[1]
aNameFriend=par[2]
# if options.bep==True:
#     nFriend = nFriend+1
#     aSuffixFriend.append("")
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
#    for iFr in range(nFriend):
    aPathFriend.append(pathFile[:-5])
    aPathFriend[-1] = aPathFriend[-1] + "_" + aSuffixFriend + ".root"
    print "  Friend:", aPathFriend[-1]
    strAddFriend = aNameFriend + " = MeritTuple"
    trMel.AddFriend(strAddFriend, aPathFriend[-1])
    listFr = trMel.GetListOfFriends()
    if listFr is None:
        print "No added friends!"
        sys.exit(1)
    else:
        print listFr.Print()
        trMel.SetAlias("{0}Log".format(aNameFriend), "-log10(1.0-({0}+1.0)/2.0)".format(aNameFriend))
        listAlias = trMel.GetListOfAliases()
        print listAlias.Print()
        trMel.Write()
    fileMel.Close()
