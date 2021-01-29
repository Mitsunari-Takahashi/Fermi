#!/usr/bin/env python

import sys
import os.path
import ROOT
from ROOT import TFile, TTree, TChain, TH1, TH2, TH3, TH1F, TH2F, TH3F, TGraph, TGraphErrors, TCanvas, TF1, TF2
import numpy as np
import commands
import click
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, log, log10, sqrt


ROOT.gROOT.SetBatch()


def find_first_nonzero(graph, backward=True, stat_factor=2.):
    d1, d2 = ROOT.Double(0), ROOT.Double(0) #ROOT.Double
    e2 = 0.
    for ip in range(graph.GetN()):
        jp =  ip if backward is False else graph.GetN()-ip-1
        graph.GetPoint(jp, d1, d2 )
        e1 = graph.GetErrorX(jp)
        e2 = graph.GetErrorY(jp)
        if d2-e2*stat_factor>0 and d1-e1>0:
            return (jp, d1)


def fit_pol1(graph, func, acc_min, acc_max):
    graph.Fit(func, "Q", "", acc_min, acc_max)
    return func.GetParameter(0), func.GetParameter(1)


def find_point_rbkg(graph, bkg_floor, acc_bkg_floor, rbkg):
    d1, d2 = ROOT.Double(0), ROOT.Double(0)
    npoint_threshold = None
    npoint_threshold_error = None
    for ip in range(graph.GetN()):
        graph.GetPoint(ip, d1, d2)
        e2 = graph.GetErrorY(ip)
        #print "{0} +/- {1}".format(d2, e2)
        if (d2-e2)<=bkg_floor*rbkg/d1*acc_bkg_floor:#*d1/acc_bkg_floor:
            npoint_threshold_error = ip
            break
    for jp in range(ip, graph.GetN()):
        graph.GetPoint(jp, d1, d2)
        if d2<=bkg_floor*rbkg/d1*acc_bkg_floor:#*d1/acc_bkg_floor:
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
    list_rbkg = [1., 10.]

    if outpath is None:
        outpath = flightroc.replace("ROC", "CUT")
    file_out = TFile(outpath, "RECREATE")
    file_out.cd()
    dict_hist_out = {}
    for rbkg in list_rbkg:
        dict_hist_out[rbkg] = hist_covered.Clone("hist_bdt_cut_R{0:0>3.0f}".format(rbkg*10))
        dict_hist_out[rbkg].SetTitle("Flight data BDT cut R{0:0>3.0f}".format(rbkg*10))

    for iloge in range(1, 1+hist_covered.GetXaxis().GetNbins()):
        for jcth in range(1, 1+hist_covered.GetYaxis().GetNbins()):
            str_region = "E{e}_Z{c}".format(e=iloge, c=jcth)
            if hist_covered.GetBinContent(iloge, jcth)>0:
                print "*", str_region
                dir_roc = file_roc.Get(str_region)
                graph_evtdens = dir_roc.Get("graph_evtdens_{0}".format(str_region))
                hist_acc = dir_roc.Get("sig_acc_cth_cut_{0}".format(str_region))
                iacc_min_nonzero, acc_min_nonzero = find_first_nonzero(graph_evtdens, True)
                print "First non-zero acceptance (Point No. {1}): {0:1.2E}".format(acc_min_nonzero, iacc_min_nonzero)#-graph_evtdens.GetErrorX(iacc_min_nonzero))
                acc0, bkg0 = ROOT.Double(0), ROOT.Double(0) #ROOT.Double
                graph_evtdens.GetPoint(0, acc0, bkg0)

                can = TCanvas("can_{0}".format(str_region), graph_evtdens.GetTitle())

                p0, p1 = 0., 0.
                p0_prev, p1_prev = 0., 0.
                acc, acc_prev = 0., 0.
                fpol1 = TF1("fpol1", "[0]*1E-5+[1]*x", 0, acc0)
                fpol1.SetParameter(0, 3)
                fpol1.SetParameter(1, 0.)
                for jacc in range(NSTEP):
                    acc = acc0/2.*float(NSTEP-jacc)/NSTEP
                    fit_pol1(graph_evtdens, fpol1, acc_min_nonzero-graph_evtdens.GetErrorX(iacc_min_nonzero), acc)
                    p0 = (fpol1.Eval(acc_min_nonzero) + fpol1.Eval(acc))/2.
                    p1 = fpol1.GetParameter(1) #fit_pol1(graph_evtdens, fpol1, acc_min_nonzero, acc)
                    e1 = fpol1.GetParError(1)

                    #pmin = fpol1.Eval(acc_min_nonzero)
                    #pmax = fpol1.Eval(acc)

                    if acc<=acc_min_nonzero:
                        break
                    elif abs(p1)<e1:#0:
                        break
                    else:
                        p0_prev, p1_prev = p0, p1
                        acc_prev = acc
                print acc_prev, p0_prev, p1_prev
                print acc, p0, p1

                bkg_floor = p0#_prev
                acc_bkg_floor = acc#_prev
                print "Background floor: {0:1.2E} (Acceptance: {1:1.2E})".format(bkg_floor, acc_bkg_floor)

                can.cd()
                graph_evtdens.Draw("APL")
                graph_evtdens.GetXaxis().SetRangeUser(0, acc0*0.8)
                graph_evtdens.GetYaxis().SetRangeUser(bkg_floor*0.1, bkg_floor*100)
                can.SetLogy()
                can.Write()

                for rbkg in list_rbkg:
                    print "{0}x BKG floor".format(rbkg)
                    np_cut, np_cut_error = find_point_rbkg(graph_evtdens, bkg_floor=bkg_floor, acc_bkg_floor=acc_bkg_floor, rbkg=rbkg)
                    cut_bdt = hist_acc.GetXaxis().GetBinLowEdge(np_cut)
                    cut_bdt_err = cut_bdt - hist_acc.GetXaxis().GetBinLowEdge(np_cut_error)
                    d1, d2 = ROOT.Double(0), ROOT.Double(0) #ROOT.Double
                    graph_evtdens.GetPoint( np_cut, d1, d2 )
                    print "  BDT cut: {0:1.2E} +/- {1:1.2E} (Acceptance: {2:1.2E}; Background: {3:1.2E})".format(cut_bdt, cut_bdt_err, hist_acc.GetBinContent(np_cut), d2)
                    dict_hist_out[rbkg].SetBinContent(iloge, jcth, cut_bdt)
                    dict_hist_out[rbkg].SetBinError(iloge, jcth, cut_bdt_err)
                print ""

    dict_func_cut = {}
    dict_func_cut['pol0'] = TF2("pol0", "[0]", 4.35, 5.75, 0.0, 1.0)
    dict_func_cut['pol1'] = TF2("pol1", "[0]+[1]*x+[2]*y", 4.35, 5.75, 0.0, 1.0)
    dict_func_cut['pol2'] = TF2("pol2", "[0]+[1]*x+[2]*y+[3]*x*x+[4]*x*y+[5]*y*y", 4.35, 5.75, 0.0, 1.0)

    for rbkg in list_rbkg:
        dict_hist_out[rbkg].Fit()
        dict_hist_out[rbkg].Write()


if __name__ == '__main__':
    main()
