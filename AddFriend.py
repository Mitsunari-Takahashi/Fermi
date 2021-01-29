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
@click.option('--friend', '-f', type=(str, str), multiple=True, default=[(None, None)], help="(Input suffix, Friend name)")
def main(pathfiles, remove, friend):
    li_path_files = ls_list(pathfiles)

    if friend==[(None, None)]:
        nFriend = 2
        #aSuffixFriend=["newVariables", "WP8CalOnlyBEPCaseE_myBDT"]
        #aNameFriend=["newVariables", "WP8CalOnlyBEP_E_BDT"]
        aSuffixFriend=["newVariables", "U01_V200909_62IRFBK_020RAWE20woRW_09ZDir_MinNode2_D06_MinNode00005"]
        aNameFriend=["newVariables", "UpDownZDir"]
    else:
        nFriend = len(friend)
        aSuffixFriend=[]
        aNameFriend=[]
        for fr in friend:
            aSuffixFriend.append(fr[0])
            aNameFriend.append(fr[1])

    for pathFile in li_path_files:
        if pathFile[0] is not "/":
            print "Please input ablosute path!"
            sys.exit(1)
        fileMel = ROOT.TFile(pathFile, "update")
        print fileMel.GetName(), "is opened."
        trMel = fileMel.Get("MeritTuple")
        print trMel.GetName(), "is found."
        if remove==True and trMel.GetListOfFriends()!=None:
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
            if "WP8CalOnlyBEP_E_BDT" in aNameFriend:
                trMel.SetAlias("WP8CTCalOnlyBEPProb_E", "(WP8CalOnlyBEPCaseE_myBDT+1)/2.0")
            if "UpDownZDir" in aNameFriend:
                trMel.SetAlias("UpDownZDir", "U01_V200909_62IRFBK_020RAWE20woRW_09ZDir_MinNode2_D06_MinNode00005")
            if "S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR060100_BDTG500D06" in aNameFriend:
                trMel.SetAlias("S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR060100_BDTG500D06Log", "-log10(1.-(1.+S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR060100_BDTG500D06)/2.)")
            if "S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR020060_BDTG500D06" in aNameFriend:
                trMel.SetAlias("S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR020060_BDTG500D06Log", "-log10(1.-(1.+S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR020060_BDTG500D06)/2.)")
            if "T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_12Dir_ZDIR020050_BDTG500D06" in aNameFriend:
                trMel.SetAlias("T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_12Dir_ZDIR020050_BDTG500D06Log", "-log10(1.-(1.+T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_12Dir_ZDIR020050_BDTG500D06)/2.)")
            if "T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_12Dir_ZDIR050100_BDTG500D06" in aNameFriend:
                trMel.SetAlias("T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_12Dir_ZDIR050100_BDTG500D06Log", "-log10(1.-(1.+T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_12Dir_ZDIR050100_BDTG500D06)/2.)")

            trMel.Write()
        fileMel.Close()


if __name__ == '__main__':
    main()
