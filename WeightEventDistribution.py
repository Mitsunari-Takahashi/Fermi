#!/usr/bin/env python

import sys
import os
import numpy as np
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TTree
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


def weight_distribution(pathbase_file_out, path_hist_origin, name_hist_origin, path_hist_goal, name_hist_goal, path_file_merit, name_tree_merit='MeritTuple', zdironly=False, factor=1.0):
    file_merit = ROOT.TFile(path_file_merit, "READ")
    tree_merit = file_merit.Get(name_tree_merit)
    nEvt = tree_merit.GetEntries()
    print tree_merit.GetName(), "has", nEvt, "events."
    file_hist_origin = ROOT.TFile(path_hist_origin, "READ")
    hist_origin = file_hist_origin.Get(name_hist_origin)
    hist_origin.RebinY(2)
    file_hist_goal = ROOT.TFile(path_hist_goal, "READ")
    hist_goal = file_hist_goal.Get(name_hist_goal)
    hist_goal.RebinY(2)
    if zdironly==True:
        hist_origin = hist_origin.ProjectionY("{0}_ZDir".format(hist_origin.GetName()), 1, hist_origin.GetNbinsX())
        print hist_origin.GetName()
        for ibin in range(1, hist_origin.GetNbinsX()+1):
            print ibin, hist_origin.GetBinContent(ibin)
        hist_goal = hist_goal.ProjectionY("{0}_ZDir".format(hist_goal.GetName()), 1, hist_goal.GetNbinsX())
        print hist_origin.GetName()
        for ibin in range(1, hist_goal.GetNbinsX()+1):
            print ibin, hist_goal.GetBinContent(ibin)


    hist_ratio = hist_goal.Clone('hist_ratio')
    hist_ratio.SetTitle('Ratio of the original event counts to the goal counts')
    hist_ratio.Divide(hist_origin)
    print hist_ratio.GetName()
    for ibin in range(1, hist_ratio.GetNbinsX()+1):
        print ibin, hist_ratio.GetBinContent(ibin)
    downscale_minimum = hist_ratio.GetMaximum()
    if downscale_minimum<=0:
        logger.critical('Histogram {0} has zero!!!'.format(hist_ratio.GetName()))
    else:
        logger.info('Scaling factor: {0}'.format(downscale_minimum))
    hist_ratio.Scale(factor/downscale_minimum)

    file_out = ROOT.TFile('.'.join([pathbase_file_out, 'root']), "UPDATE")
    file_out.cd()
    can_ratio = ROOT.TCanvas('can_ratio', 'Ratio of the original event counts to the goal counts')
    can_ratio.cd()
    hist_ratio.Draw("colz" if zdironly==False else "")
    can_ratio.SaveAs('.'.join([pathbase_file_out, 'pdf']))

    tree_rannum = ROOT.TTree('weights', 'weights')
    vweight = c_float()
    str_name_branch = 'weights'
    tree_rannum.Branch(str_name_branch, vweight, str_name_branch+'/F')
    file_out.cd()
    for iEvt, evt in enumerate(tree_merit):
        if zdironly==False:
            vweight.value = hist_ratio.GetBinContent(hist_ratio.GetXaxis().FindBin(evt.McLogEnergy), hist_ratio.GetYaxis().FindBin(-evt.McZDir))
        else:
            vweight.value = hist_ratio.GetBinContent(hist_ratio.GetXaxis().FindBin(-evt.McZDir))
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
@click.option('--origin', type=(str, str), help='<path> <name> of the original histogram', default=('/nfs/farm/g/glast/u/mtakahas/data/MC/allPro_tuning/hist_allPro200909_62_p1-7_2016Nov_McZDir_vs_McLogEnergy.root', 'h2_McZDir_vs_McLogEnergy'))
@click.option('--goal', type=(str, str), help='<path> <name> of the goal histogram', default=('/nfs/farm/g/glast/u/mtakahas/data/MC/allPro_tuning/histBKG200909_62MCE2e4_Combined_McZDir_vs_McLogEnergy.root', 'h2_McZDir_vs_McLogEnergy'))
@click.option('--zdironly', is_flag=True, default=False)
@click.option('--factor', type=float, default=1.0)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(inpathmerit, outpathbase, origin, goal, zdironly, factor, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)
    weight_distribution(pathbase_file_out=outpathbase, path_hist_origin=origin[0], name_hist_origin=origin[1], path_hist_goal=goal[0], name_hist_goal=goal[1], path_file_merit=inpathmerit, zdironly=zdironly, factor=factor)


if __name__ == '__main__':
    main()
