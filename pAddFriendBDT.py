#!/usr/bin/env python

import ROOT
import sys
import os
import os.path
ROOT.gROOT.SetBatch()
par = sys.argv #args # 
print par
aPathFile = par[1:]
#nFriend = len(aPathFile)
#print nFriend
#aSuffixFriend=par[1]
#aNameFriend=par[2]
# if options.bep==True:
#     nFriend = nFriend+1
#     aSuffixFriend.append("")


MAXCTH = 1.0
MINCTH = -1.0
MAXLOGE = 10
MINLOGE = 0

class BDTVariable:
    def __init__(self, name, file_suffix, cth_range=(MINCTH, MAXCTH), loge_range=(MINLOGE, MAXLOGE)):
        self.name = name
        self.file_suffix = file_suffix
        self.cth_min = cth_range[0]
        self.cth_max = cth_range[1]
        self.loge_min = loge_range[0]
        self.loge_max = loge_range[1]


    def get_friendpath(self, filepath):
        return filepath.replace(".root", "_{}.root".format(self.file_suffix))
        #return '{0}_{1}.root'.format(os.path.splitext(filepath), self.file_suffix) #aSuffixFriend)


    def get_str_limit(self):
        list_lim = []
        if self.cth_max<MAXCTH:
            list_lim.append('(Cal1MomZDir<{0}'.format(self.cth_max))
        if self.cth_min>MINCTH:
            list_lim.append('(Cal1MomZDir>={0}'.format(self.cth_min))
        if self.loge_max<MAXLOGE:
            list_lim.append('(log10(WP8CalOnlyEnergy)<{0}'.format(self.loge_max))
        if self.loge_min<MINLOGE:
            list_lim.append('(log10(WP8CalOnlyEnergy)>={0}'.format(self.loge_min))
        str_lim = ' && '.join(list_lim)


        
class BDTVariableSet:
    def __init__(self, name, variables):
        self.name = name
        self.dict_variables = {}
        for v in variables:
            self.dict_variables[v.name] = v


    def add_friends(self, tree_merit, file_origin):
        for n, v in self.dict_variables.items():
            strAddFriend = '='.join([n, "MeritTuple"])
            tree_merit.AddFriend(strAddFriend, v.get_friendpath(file_origin))
        listFr = trMel.GetListOfFriends()
        if listFr is None:
            print "No friend was added!"
            sys.exit(1)
        else:
            listFr.Print()


    def set_aliases(self, tree_merit):
        for n, v in self.dict_variables.items():
            tree_merit.SetAlias("{0}Log".format(n), "-log10(1.0-({0}+1.0)/2.0)".format(n))
        #print trMel.GetListOfAliases().Print()
        trMel.Write()


    #def get_unialias(self):
    #    for n, v in self.dict_variables:
        


v_S20_catTwoZDIR_17_Z1 = BDTVariable(name='S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR020060_BDTG500D06', file_suffix='S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR020060_BDTG500D06', cth_range=(0.2, 0.6))
v_S20_catTwoZDIR_17_Z2 = BDTVariable(name='S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR060100_BDTG500D06', file_suffix='S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR060100_BDTG500D06', cth_range=(0.6, 1.0))
v_T02_catTwoZDIR_11_Z1 = BDTVariable(name='T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_11Dir_ZDIR020050_BDTG500D06', file_suffix='T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_11Dir_ZDIR020050_BDTG500D06', cth_range=(0.2, 0.5))
v_T02_catTwoZDIR_11_Z2 = BDTVariable(name='T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_11Dir_ZDIR050100_BDTG500D06', file_suffix='T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_11Dir_ZDIR050100_BDTG500D06', cth_range=(0.5, 1.0))

vset_S20_catTwoZDIR_17 = BDTVariableSet(name='S20_catTwoZDIR_17', variables=[v_S20_catTwoZDIR_17_Z1, v_S20_catTwoZDIR_17_Z2, v_T02_catTwoZDIR_11_Z1, v_T02_catTwoZDIR_11_Z2])


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

    print vset_S20_catTwoZDIR_17.name
    vset_S20_catTwoZDIR_17.add_friends(trMel, pathFile)
    vset_S20_catTwoZDIR_17.set_aliases(trMel)

    # aPathFriend.append(pathFile[:-5])
    # aPathFriend[-1] = aPathFriend[-1] + "_" + aSuffixFriend + ".root"
    # print "  Friend:", aPathFriend[-1]
    # strAddFriend = aNameFriend + " = MeritTuple"
    # trMel.AddFriend(strAddFriend, aPathFriend[-1])
    # listFr = trMel.GetListOfFriends()
    # if listFr is None:
    #     print "No added friends!"
    #     sys.exit(1)
    # else:
    #     print listFr.Print()
    #     trMel.SetAlias("{0}Log".format(aNameFriend), "-log10(1.0-({0}+1.0)/2.0)".format(aNameFriend))
    #     listAlias = trMel.GetListOfAliases()
    #     print listAlias.Print()
    #     trMel.Write()
    fileMel.Close()
