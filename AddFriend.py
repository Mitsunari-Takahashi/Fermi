#!/usr/bin/env python

import ROOT
import sys
import click
import subprocess
from pLsList import ls_list
ROOT.gROOT.SetBatch()


@click.command()
@click.argument('pathfiles', type=str)
@click.option('--remove', '-r', is_flag=True, help="Remove the current friends")
def main(pathfiles, remove):
    li_path_files = ls_list(pathfiles)

    nFriend = 2
    aSuffixFriend=["newVariables", "WP8CalOnlyBEPCaseE_myBDT"]
    aNameFriend=["newVariables", "WP8CalOnlyBEP_E_BDT"]

    for pathFile in li_path_files:
        if pathFile[0] is not "/":
            print "Please input ablosute path!"
            sys.exit(1)
        fileMel = ROOT.TFile(pathFile, "update")
        print fileMel.GetName(), "is opened."
        trMel = fileMel.Get("MeritTuple")
        print trMel.GetName(), "is found."
        if remove==True:
            trMel.GetListOfFriends().RemoveAll()
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


if __name__ == '__main__':
    main()
