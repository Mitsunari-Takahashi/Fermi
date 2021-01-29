#!/usr/bin/env python2

import sys
import os
import os.path as path
import subprocess
import time
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TEntryList, TEventList, TTree, TChain, TFile, TList
ROOT.gROOT.SetBatch()
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL

##### Logger #####
logger = getLogger(__name__)
handler = StreamHandler()
loglevel = 'INFO'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)

@click.command()
@click.argument('pathin', type=str)
@click.option('--cut', type=str, default='U01_V200909_62IRFBK_020RAWE20woRW_09_BDTG200D05>=0.0 && EvtEventId%5==0 && CalOnly==1 && Cal1MomZDir>=0.2 && log10(WP8CalOnlyEnergy)>=4.35 && log10(WP8CalOnlyEnergy)<=5.75 && Cal1MomZCrossSide840>=0 && (GltGemSummary&0x20)==0') #"U01_V200909_62IRFBK_020RAWE20woRW_09_BDTG200D05>=0.0 && CalOnly==1 && Cal1MomZDir>=0.2 && log10(WP8CalOnlyEnergy)>=4.35 && log10(WP8CalOnlyEnergy)<=5.75 && McSourceId!=7000 && Cal1MomZCrossSide840>=0 && (GltGemSummary&0x20)==0")
@click.option('--branches', type=str, default='/u/gl/mtakahas/work/data/lists/branchlist_forEval_2019Oct.yaml')
@click.option('--suffix', type=str, default='Skimmed')
@click.option('--friend', '-f', multiple=True, default=[])
@click.option('--skipskim', is_flag=True)
@click.option('--outdir', type=str, default='.')
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
@click.option('--logout', type=str, default='./MeritSkimmerWithFriends.log')
def main(pathin, cut, branches, suffix, friend, skipskim, outdir, loglevel, logout):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)

    filein = ROOT.TFile(pathin, "READ")
    treein  = filein.Get("MeritTuple")
    logger.info('Input tree has {0} events.'.format(treein.GetEntries()))
    pathskimmed = '/'.join([outdir, path.basename(pathin).replace(".root", "_{0}.root".format(suffix))])

    if skipskim==False:
#        mscmd = ['bsub', '-o{0}'.format(logout), '-JskimFr', '-W1000', '-We100', 'python', '/u/gl/mtakahas/eventSelectForTMVA42/scripts/MeritSkimmer.py', '--aliases=pass8', '--selection="{0}"'.format(cut), '--extra_trees=jobinfo', pathin, '--output={0}'.format(pathskimmed), '--branches={0}'.format(branches)]
        mscmd = ['python', '/u/gl/mtakahas/eventSelectForTMVA42/scripts/MeritSkimmer.py', '--aliases=pass8', '--selection="{0}"'.format(cut), '--extra_trees=jobinfo', pathin, '--output={0}'.format(pathskimmed), '--branches={0}'.format(branches)]
        logger.debug(mscmd)
        proc_skim = subprocess.Popen(mscmd)
        #time.sleep(1000*60)
    
    treein.Draw(">>enlist", cut) #, "entrylist")
    enlist = gDirectory.Get("enlist")
    logger.info("Event list: {0}".format(enlist.GetN()))

    dictfriends = {}
    dict_friendsuffix = {'WP8CalOnlyBEP_E_BDT': 'WP8CalOnlyBEPCaseE_myBDT', 
                         'UpDown':'U01_V200909_62IRFBK_020RAWE20woRW_09',
                         'UpDownZDir':'U01_V200909_62IRFBK_020RAWE20woRW_09ZDir_MinNode2_D06_MinNode00005'}
    if len(friend)==0:
        logger.warning('List of friends in {0} is referred.'.format(treein.GetName()))
        rlistfriends = treein.GetListOfFriends()
        for ifr in range(rlistfriends.GetEntries()):
            logger.debug(rlistfriends.At(ifr).GetName())
            dictfriends[rlistfriends.At(ifr).GetName()] = rlistfriends.At(ifr).GetTitle()
    else:
        for fr in friend:
            logger.info(fr)
            if fr in dict_friendsuffix.keys():
                dictfriends[fr] = pathin.replace(".root", "_{0}.root".format(dict_friendsuffix[fr]))
            else:
                dictfriends[fr] = pathin.replace(".root", "_{0}.root".format(fr))

    logger.debug(dictfriends)
    for namefriend, pathfriendin in dictfriends.items():
        logger.info(pathfriendin)
        filefriendin = ROOT.TFile(pathfriendin, "READ")
        treefriendin = filefriendin.Get("MeritTuple")
        logger.info('Friend tree has {0} events.'.format(treefriendin.GetEntries()))
        treefriendin.SetEventList(enlist)
        #treefriendin.SetEntryList(enlist)

        if namefriend in dict_friendsuffix.keys():
            friendsuffix = dict_friendsuffix[namefriend]
        else:
            friendsuffix = namefriend
        pathfriendskimmed = '/'.join([outdir, path.basename(pathskimmed).replace(suffix, '_'.join([suffix,friendsuffix]))]) #pathskimmed.replace(".root", "_{0}.root".format(namefriend))
        logger.debug(pathfriendskimmed)
        filefriendskimmed = ROOT.TFile(pathfriendskimmed, "RECREATE")
        filefriendskimmed.cd()
        treefriendskimmed = treefriendin.CopyTree("")
        logger.info('Skimmed friend tree has {0} events.'.format(treefriendskimmed.GetEntries()))
        treefriendskimmed.Write()
        
        #treeskimmed.AddFriend(treefriendskimmed, namefriend)
    #fileskimmed.cd()
    #treeskimmed.Write()

    #fileskimmed = ROOT.TFile(pathskimmed, "UPDATE")
    #treeskimmed  = fileskimmed.Get("MeritTuple")
    #logger.info('Skimmed tree has {0} events.'.format(treeskimmed.GetEntries()))
    #treeskimmed.GetListOfFriends().RemoveAll()

    if skipskim==False:
        logger.info('Waiting the skimming process...')
        proc_skim.wait()
    

if __name__ == '__main__':
    main()
