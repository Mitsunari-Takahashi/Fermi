#!/usr/bin/env python

import sys
import os
import numpy as np
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi, log10
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TTree, TH2F
ROOT.gROOT.SetBatch()
from pColor import *
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL
from ctypes import *


##### Logger #####
logger = getLogger(__name__)
handler = StreamHandler()
loglevel = 'INFO'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)


def weight_distribution(pathbase_file_out, path_hist_table, name_hist_table, path_file_merit, name_tree_merit='MeritTuple', factor=1.0):
    file_merit = ROOT.TFile(path_file_merit, "READ")
    tree_merit = file_merit.Get(name_tree_merit)
    nEvt = tree_merit.GetEntries()
    print tree_merit.GetName(), "has", nEvt, "events."

    file_out = ROOT.TFile('.'.join([pathbase_file_out, 'root']), "UPDATE")
    print '{0} is opened.'.format(file_out.GetName())

    if path_hist_table==None and name_hist_table==None:
        flag_notable = True
    else:
        flag_notable = False
        file_hist_table = ROOT.TFile(path_hist_table, "READ")
        hist_table = file_hist_table.Get(name_hist_table)
        #hist_table.Scale(1./hist_table.Integral())
        #hist_table.Write()
        
    file_out.cd()
    tree_rannum = ROOT.TTree('MeritTuple', 'weights')
    vweight = c_float()
    str_name_branch = 'weights'
    tree_rannum.Branch(str_name_branch, vweight, str_name_branch+'/F')
    file_out.cd()
    if flag_notable==True:
        for iEvt, evt in enumerate(tree_merit):
            vweight.value = factor
            tree_rannum.Fill()
            if iEvt%(nEvt/1000)==0:
                rate = int((iEvt*100.)/nEvt+0.5)
                if rate>0:
                    meter = "\r[{0}{1}] WEIGHT: {2}".format("=" * rate, ' ' * (100-rate), vweight.value)
                    sys.stdout.write(meter)
                    sys.stdout.flush()
    else:
        for iEvt, evt in enumerate(tree_merit):
            #vweight.value = 1./hist_table.GetBinContent(hist_table.GetXaxis().FindBin(log10(evt.WP8CalOnlyEnergy)), hist_table.GetYaxis().FindBin(evt.Cal1MomZDir)) * factor
            vweight.value = hist_table.GetBinContent(hist_table.GetXaxis().FindBin(log10(evt.WP8CalOnlyEnergy)), hist_table.GetYaxis().FindBin(evt.Cal1MomZDir)) * factor
            #vweight.value = hist_table.GetBinContent(hist_table.GetXaxis().FindBin(evt.McLogEnergy), hist_table.GetYaxis().FindBin(-evt.McZDir)) * factor
            tree_rannum.Fill()
            if iEvt%(nEvt/1000)==0:
                rate = int((iEvt*100.)/nEvt+0.5)
                if rate>0:
                    meter = "\r[{0}{1}] WEIGHT: {2}".format("=" * rate, ' ' * (100-rate), vweight.value)
                    sys.stdout.write(meter)
                    sys.stdout.flush()
    print ""
    print nEvt, "events have been filled."
    tree_rannum.Write()


@click.command()
@click.argument('inpathmerit', type=str)
@click.argument('outpathbase', type=str)
@click.option('--table', type=(str, str), default=(None, None), help='<path> <name> of the table histogram') #, default=('/nfs/farm/g/glast/u/mtakahas/data/MC/MixSamples/ThreeMCSamplesDistribution.root', 'histProtonPrim_weight'))
@click.option('--factor', type=float, default=1.0)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(inpathmerit, outpathbase, table, factor, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)
    weight_distribution(pathbase_file_out=outpathbase, path_hist_table=table[0], name_hist_table=table[1], path_file_merit=inpathmerit, factor=factor)


if __name__ == '__main__':
    main()
