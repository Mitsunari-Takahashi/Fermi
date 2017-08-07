#!/usr/bin/env python

import sys
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE
ROOT.gROOT.SetBatch()
from pColor import *
import ReadLTFCatalogueInfo


def plot_grbs_deviation(path_tree, path_table, path_out):
    file_tree = ROOT.TFile(path_tree, "READ")
    tree = file_tree.Get('tr')
    ngrb = tree.GetEntries()/2
    print ngrb, 'GRBs.'
    tb = ReadLTFCatalogueInfo.open_table(1, path_table)
    FLUENCE_CUT = [1.45E-4, 3.70E-5]
    lst_htg_deviation = []
    lst_htg_deviation_signed = []
    lst_gre_deviation_signed_vs_index005GeV = []
    lst_htg_index100GeV = []
    lst_htg_index005GeV = []
    lst_gre_index_100GeV_vs_005GeV = []
    str_cat = ['all categories', 'category 1', 'category 2', 'category 3']
    file_out = ROOT.TFile(path_out, 'RECREATE')

    for icat in (0,1,2,3):
        lst_htg_deviation.append(ROOT.TH1D('htg_deviation_{0}'.format(icat), '{0};Deviation[#sigma],[bursts]'.format(str_cat[icat]), 110, 0, 11))
        lst_htg_deviation_signed.append(ROOT.TH1D('htg_deviation_signed_{0}'.format(icat), 'Signed deviation of counts in 10 GeV - 100 GeV;[#sigma],[bursts]', 160, -11, 5))
        lst_gre_deviation_signed_vs_index005GeV.append(ROOT.TGraphErrors())
        lst_gre_deviation_signed_vs_index005GeV[-1].SetName('gre_deviation_signed_vs_index005GeV_{0}'.format(icat))
        lst_gre_deviation_signed_vs_index005GeV[-1].SetTitle(str_cat[icat])
        lst_gre_deviation_signed_vs_index005GeV[-1].GetYaxis().SetTitle('Deviation')
        lst_gre_deviation_signed_vs_index005GeV[-1].GetXaxis().SetTitle('Index')
        lst_htg_index100GeV.append(ROOT.TH1D('htg_index100GeV_{0}'.format(icat), str_cat[icat]+';Index in 100MeV-100GeV', 80, -8, 0))
        lst_htg_index005GeV.append(ROOT.TH1D('htg_index005GeV_{0}'.format(icat), str_cat[icat]+';Index in 100MeV-5GeV', 80, -8, 0))
        lst_gre_index_100GeV_vs_005GeV.append(ROOT.TGraphErrors())
        lst_gre_index_100GeV_vs_005GeV[-1].SetName('gre_index_100GeV_vs_005GeV_{0}'.format(icat))
        lst_gre_index_100GeV_vs_005GeV[-1].SetTitle(str_cat[icat])
        lst_gre_index_100GeV_vs_005GeV[-1].GetYaxis().SetTitle('Index in 100MeV-5GeV')
        lst_gre_index_100GeV_vs_005GeV[-1].GetXaxis().SetTitle('Index in 100MeV-100GeV')
        
    for i in range(ngrb):
        tree.GetEntry(2*i)
        name = '{0:0>9}'.format(int(tree.name))
        print '=====', name, '====='
        #tb1 = ReadLTFCatalogueInfo.select_one_by_name(int(tree.name), path_table)
        #tb1 = ReadLTFCatalogueInfo.select_one_by_name(path_table, name)
        #fluence_gbm = tb1['FLUENCE']
        ncategory = ReadLTFCatalogueInfo.judge_category_fluence(tb, name, FLUENCE_CUT)
        print 'GBM fluence category:', ncategory
        index005GeV = tree.Index1
        index005GeV_err = tree.Index1_err
        tree.GetEntry(2*i+1)
        ts100GeV = tree.ts
        print 'TS:', ts100GeV
        if ts100GeV<25:
            print ' Skipping.'
            continue
        index100GeV = tree.Index1
        index100GeV_err = tree.Index1_err
        deviation = tree.deviation
        sign_deviation = tree.sign_deviation

        lst_htg_deviation[0].Fill(deviation)
        lst_htg_deviation_signed[0].Fill(sign_deviation*deviation)
        lst_gre_deviation_signed_vs_index005GeV[0].SetPoint(i, index005GeV, sign_deviation*deviation)
        lst_gre_deviation_signed_vs_index005GeV[0].SetPointError(i, index005GeV_err, 0)
        lst_htg_index100GeV[0].Fill(index100GeV)
        lst_htg_index005GeV[0].Fill(index005GeV)
        lst_gre_index_100GeV_vs_005GeV[0].SetPoint(i, index005GeV, index100GeV)
        lst_gre_index_100GeV_vs_005GeV[0].SetPointError(i, index005GeV_err, index100GeV_err)

        lst_htg_deviation[ncategory].Fill(deviation)
        lst_htg_deviation_signed[ncategory].Fill(sign_deviation*deviation)
        lst_gre_deviation_signed_vs_index005GeV[ncategory].SetPoint(i, index005GeV, sign_deviation*deviation)
        lst_gre_deviation_signed_vs_index005GeV[ncategory].SetPointError(i, index005GeV_err, 0)
        lst_htg_index100GeV[ncategory].Fill(index100GeV)
        lst_htg_index005GeV[ncategory].Fill(index005GeV)
        lst_gre_index_100GeV_vs_005GeV[ncategory].SetPoint(i, index005GeV, index100GeV)
        lst_gre_index_100GeV_vs_005GeV[ncategory].SetPointError(i, index005GeV_err, index100GeV_err)

    file_out.cd()
    for icat in (0,1,2,3):    
        lst_htg_deviation[icat].Write()
        lst_htg_deviation_signed[icat].Write()
        lst_gre_deviation_signed_vs_index005GeV[icat].Write()
        lst_htg_index100GeV[icat].Write()
        lst_htg_index005GeV[icat].Write()
        lst_gre_index_100GeV_vs_005GeV[icat].Write()


@click.command()
@click.argument('tree', type=str)
@click.argument('output', type=str)
@click.option('--table', type=str, default='/Users/Mitsunari/FermiAnalysis/catalogue/LAT2CATALOG-v1-LTF.fits')
def main(tree, output, table):
    plot_grbs_deviation(path_tree=tree, path_table=table, path_out=output)


if __name__ == '__main__':
    main()
