#!/usr/bin/env python

import ROOT
import sys
#from optparse import OptionParser

# if __name__ == '__main__':
#     p = OptionParser(version="ver:%s" % __version__)
#     p.add_option('--BEP', action='store_true', dest="bep", default=False,
#                  help="Add BEP variable as a friend.")
#     (opts, args) = p.parse_args()
ROOT.gROOT.SetBatch()
par = sys.argv #args # 
print par
aPathFile = par[1:]

aSuffixFriend=["newVariables", "WP8CalOnlyBEPCaseE_myBDT", "CalOnlyVar"]
aNameFriend=["newVariables", "WP8CalOnlyBEP_E_BDT", "CalOnlyVar"]
if len(aSuffixFriend)!=len(aNameFriend):
    print "Numbers of the name and the suffix of the friend files are different!!"
    sys.exit(1)
nFriend = len(aSuffixFriend)
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
