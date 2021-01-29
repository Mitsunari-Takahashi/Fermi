#!/usr/bin/env python

import sys
import os.path
import ROOT
from ROOT import TFile, TTree, TChain, TH1, TH2, TH3, TH1F, TH2F, TH3F, TGraph, TGraphErrors, TCanvas, TF1, TF2, kRed, kBlue, kMagenta, kBlack, kWhite
import numpy as np
import commands
import click
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, log, log10, sqrt
import CreateGalOffROC3D

ROOT.gROOT.SetBatch()


def find_point_rbkg(graph, bkg_floor, rbkg):
    d1, d2 = ROOT.Double(0), ROOT.Double(0)
    npoint_threshold = None
    npoint_threshold_error = None
    for ip in range(graph.GetN()):
        graph.GetPoint(ip, d1, d2)
        e2 = graph.GetErrorY(ip)
        #print "{0} +/- {1}".format(d2, e2)
        if (d2-e2)<=bkg_floor*rbkg: #/d1*acc_bkg_floor
            npoint_threshold_error = ip
            break
    for jp in range(ip, graph.GetN()):
        graph.GetPoint(jp, d1, d2)
        if d2<=bkg_floor*rbkg: #/d1*acc_bkg_floor
            npoint_threshold = jp
            break
    #print npoint_threshold_error, npoint_threshold
    return (npoint_threshold, npoint_threshold_error)
    


@click.command()
@click.argument('flightroc', type=str)
@click.option('--outpath', '-o', default=None)
def main(flightroc, outpath):
    # Input
    file_roc = TFile(flightroc, "READ")
    hist_covered = file_roc.Get("hist_covered")
    NSTEP = 1000
    list_rbkg = [3., 10.]
    rbkg_offset = 1.

    if outpath is None:
        outpath = flightroc.replace("ROC", "CUT")
    file_out = TFile(outpath, "RECREATE")
    file_out.cd()
    dict_hist_out = {}
    for rbkg in list_rbkg:
        dict_hist_out[rbkg] = hist_covered.Clone("hist_bdt_cut_R{0:0>3.0f}".format(rbkg*10))
        dict_hist_out[rbkg].SetTitle("Flight data BDT cut R{0:0>3.0f}".format(rbkg*10))
    
    loge_valid_min = dict_hist_out[rbkg].GetXaxis().GetXmax()
    loge_valid_max = 0
    cth_valid_min = dict_hist_out[rbkg].GetYaxis().GetXmax()
    cth_valid_max = 0

    for iloge in range(1, 1+hist_covered.GetXaxis().GetNbins()):
        for jcth in range(1, 1+hist_covered.GetYaxis().GetNbins()):
            str_region = "E{e}_Z{c}".format(e=iloge, c=jcth)
            if hist_covered.GetBinContent(iloge, jcth)>0:
                print "*", str_region

                loge_bin_min = hist_covered.GetXaxis().GetBinLowEdge(iloge)
                if loge_bin_min<loge_valid_min:
                    loge_valid_min = loge_bin_min
                loge_bin_max = hist_covered.GetXaxis().GetBinLowEdge(iloge+1)
                if loge_bin_max>loge_valid_max:
                    loge_valid_max = loge_bin_max
                cth_bin_min = hist_covered.GetYaxis().GetBinLowEdge(jcth)
                if cth_bin_min<cth_valid_min:
                    cth_valid_min = cth_bin_min
                cth_bin_max = hist_covered.GetYaxis().GetBinLowEdge(jcth+1)
                if cth_bin_max>cth_valid_max:
                    cth_valid_max = cth_bin_max

                dir_roc = file_roc.Get(str_region)
                graph_evtdens = dir_roc.Get("graph_evtdens_{0}".format(str_region))
                hist_acc = dir_roc.Get("sig_acc_cth_cut_{0}".format(str_region))
                for rbkg in list_rbkg:
                    print "{0}x BKG floor".format(rbkg)
                    bkg_floor = CreateGalOffROC3D.get_ebl_count_integral(e0=10**hist_covered.GetXaxis().GetBinLowEdge(iloge), e1=10**hist_covered.GetXaxis().GetBinLowEdge(iloge+1))
                    print "BKG floor:", bkg_floor
                    np_cut, np_cut_error = find_point_rbkg(graph_evtdens, bkg_floor=bkg_floor, rbkg=rbkg+rbkg_offset)
                    print np_cut, np_cut_error
                    cut_bdt = hist_acc.GetXaxis().GetBinLowEdge(np_cut)
                    cut_bdt_err = cut_bdt - hist_acc.GetXaxis().GetBinLowEdge(np_cut_error)
                    d1, d2 = ROOT.Double(0), ROOT.Double(0) #ROOT.Double
                    graph_evtdens.GetPoint( np_cut, d1, d2 )
                    print "  BDT cut: {0:1.2E} +/- {1:1.2E} (Acceptance: {2:1.2E}; Background: {3:1.2E})".format(cut_bdt, cut_bdt_err, hist_acc.GetBinContent(np_cut), d2)
                    dict_hist_out[rbkg].SetBinContent(iloge, jcth, cut_bdt)
                    dict_hist_out[rbkg].SetBinError(iloge, jcth, cut_bdt_err)
                print ""

    dict_func_cut = {}
    #dict_func_cut['pol0'] = TF2("pol0", "[0]", loge_valid_min, loge_valid_max, cth_valid_min, cth_valid_max)
    #dict_func_cut['pol1'] = TF2("pol1", "[0]+[1]*x+[2]*y", loge_valid_min, loge_valid_max, cth_valid_min, cth_valid_max)
    #dict_func_cut['pol1'].SetLineColor(kRed)
    #dict_func_cut['pol2'] = TF2("pol2", "[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y", loge_valid_min, loge_valid_max, cth_valid_min, cth_valid_max)
    dict_func_cut['pol3'] = TF2("pol3", "[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y+[6]*x*x*x+[7]*x*x*y+[8]*x*y*y+[9]*y*y*y", loge_valid_min, loge_valid_max, cth_valid_min, cth_valid_max)
    dict_func_cut['pol3'].SetLineColor(kWhite)

    for rbkg in list_rbkg:
        print "* CalOnlyR{0:0>3.0f}".format(rbkg*10)
        dict_hist_out[rbkg].GetXaxis().SetRangeUser(loge_valid_min, loge_valid_max)
        dict_hist_out[rbkg].GetYaxis().SetRangeUser(cth_valid_min, cth_valid_max)
        for func_name, func in dict_func_cut.items():
            print func_name
            dict_hist_out[rbkg].Fit(func, "+", "")
            str_bdt_cut_func = "{0}".format(func.GetExpFormula())
            str_bdt_cut_func = str_bdt_cut_func.replace("*x", "*log10(WP8CalOnlyEnergy)").replace("*y", "*Cal1MomZDir").replace("(x", "(log10(WP8CalOnlyEnergy)").replace("(y", "(Cal1MomZDir")
            for kf in range(0, func.GetNpar()):
                str_bdt_cut_func = str_bdt_cut_func.replace("[{0}]".format(kf), "{0:1.3E}".format(func.GetParameter(kf)))
            print str_bdt_cut_func
        dict_hist_out[rbkg].Write()


if __name__ == '__main__':
    main()
