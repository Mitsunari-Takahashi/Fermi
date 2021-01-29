#!/usr/bin/env python

import sys
import os.path
import ROOT
from ROOT import gPad, TFile, TTree, TChain, TH1, TH2, TH3, TH1F, TH2F, TH3F, TGraph, TGraphErrors, TCanvas, TF1, TGaxis, kBlack, kRed, kBlue, kMagenta, kViolet, kPink, kGray, kOrange, kAzure
import numpy as np
import commands
import click
from array import array
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees, log, log10, sqrt


ROOT.gROOT.SetBatch()


def get_ebl_count_integral(e0, e1, norm=1.03E-5, slope=-2.41):
    count = norm * (pow(e1, slope+1)-pow(e0, slope+1)) / (pow(10**5, slope+1)-pow(10**2, slope+1)) * 1E4 # m^-2 s^-1 sr^-1
    return count


def create_galoff_roc3D(dir_out, hist_evt, hist_exp, hist_acc, negativecth, logemin, logemax, cthmin, cthmax, hist_bkg=None):

    cthsign = 1 if negativecth is False else -1

    dir_out.cd()
    hist_covered = TH2F('hist_covered', 'hist_covered', hist_acc.GetXaxis().GetNbins(), hist_acc.GetXaxis().GetBinLowEdge(1), hist_acc.GetXaxis().GetBinLowEdge(hist_acc.GetXaxis().GetNbins()+1), hist_acc.GetZaxis().GetNbins(), hist_acc.GetZaxis().GetBinLowEdge(1), hist_acc.GetZaxis().GetBinLowEdge(hist_acc.GetZaxis().GetNbins()+1))
    
    for ix in range(1, 1+hist_exp.GetXaxis().GetNbins()): # Energy
        energy_cen = hist_exp.GetXaxis().GetBinCenter(ix)
        energy_lo = hist_exp.GetXaxis().GetBinLowEdge(ix)
        energy_up = hist_exp.GetXaxis().GetBinLowEdge(ix+1)
        str_label_energy = "{0:1.2f}<log10(Energy)<{1:1.2f}".format(energy_lo, energy_up)
        print '* ' + str_label_energy

        ebl_count = get_ebl_count_integral(10**energy_lo, 10**energy_up)
        print "EGB: {0} m^-2 s^-1 sr^-1".format(ebl_count)
        evtdens_min_ebl = 1.
        evtdens_max_ebl = 100.
        evtdens_min = evtdens_min_ebl*ebl_count
        evtdens_max = evtdens_max_ebl*ebl_count

        if hist_exp.GetXaxis().GetBinCenter<logemin or hist_exp.GetXaxis().GetBinCenter(ix)>logemax:
            print '  Skipped.'
        else:
            for iz in range(1, 1+hist_exp.GetZaxis().GetNbins()): # Inclination
                cth_lo = hist_exp.GetZaxis().GetBinLowEdge(iz)
                cth_up = hist_exp.GetZaxis().GetBinLowEdge(iz+1)
                str_label_cth = "{0:1.2f}<cos#theta<{1:1.2f}".format(cthsign*cth_lo, cthsign*cth_up)
                print '** ' + str_label_cth
                str_label_com = ' '.join([str_label_energy, str_label_cth])
                if hist_exp.GetZaxis().GetBinCenter(iz)<cthmin or hist_exp.GetZaxis().GetBinCenter(iz)>cthmax:
                    print '  Skipped.'
                else:
                    name_dir = 'E{0}_Z{1}'.format(ix, iz)
                    dir_xz = dir_out.mkdir(name_dir, str_label_com)
                    dir_xz.cd()
                    if negativecth is False:
                        hist_evt1D = hist_evt.ProjectionY('{0}_E{1}_Z{2}'.format(hist_evt.GetName(), ix, iz), hist_evt.GetXaxis().FindBin(energy_lo), hist_evt.GetXaxis().FindBin(energy_up), hist_evt.GetZaxis().FindBin(cth_lo), hist_evt.GetZaxis().FindBin(cth_up))
                    else:
                        hist_evt1D = hist_evt.ProjectionY('{0}_E{1}_Z{2}'.format(hist_evt.GetName(), ix, iz), hist_evt.GetXaxis().FindBin(energy_lo), hist_evt.GetXaxis().FindBin(energy_up), hist_evt.GetZaxis().FindBin(-cth_up), hist_evt.GetZaxis().FindBin(-cth_lo))
                    hist_evt1D.SetTitle('{};BDT;Event number [counts]'.format(str_label_com))
                    print 'Event number: {}'.format(hist_evt1D.GetEntries())
                    hist_evt1D.Sumw2()
                    hist_exp1D = hist_exp.ProjectionY('{0}_E{1}_Z{2}'.format(hist_exp.GetName(), ix, iz), ix, ix, iz, iz)
                    hist_exp1D.SetTitle('{};BDT;Exposure [m^2 sr s]'.format(str_label_com))
                    hist_acc1D = hist_acc.ProjectionY('{0}_E{1}_Z{2}'.format(hist_acc.GetName(), ix, iz), ix, ix, iz, iz)
                    hist_acc1D.SetTitle('{};BDT;Acceptance [m^2 sr]'.format(str_label_com))
                    #hist_evt1D_cum = hist_evt1D.GetCumulative(0)
                    #hist_evt1D_cum.SetTitle('{};BDT;Cumulative event number [counts]'.format(str_label_com))
                    # hist_exp1D_cum = hist_exp1D.GetCumulative(0) #TH1F('hist_exp_cum_E{0}_Z{1}'.format(ix, iz), ' '.join([str_label_energy, str_label_cth])+';BDT;Cumulative exposure')
                    # hist_exp1D.SetTitle('{};BDT;Cumulative exposure [m^2 sr s]'.format(str_label_com))
                    # hist_acc1D_cum = hist_acc1D.GetCumulative(0)
                    # hist_acc1D.SetTitle('{};BDT;Cumulative acceptance [m^2 sr]'.format(str_label_com))

                    hist_bkg1D = hist_bkg.ProjectionY('{0}_E{1}_Z{2}'.format(hist_bkg.GetName(), ix, iz), ix, ix, iz, iz)
                    hist_bkg1D.SetTitle('{};BDT;Background rate'.format(str_label_com))
                    hist_bkg1D.SetLineColor(kGray)

                    hist_evt1D.Write()
                    hist_exp1D.Write()
                    hist_acc1D.Write()
                    hist_bkg1D.Write()

                    hist_evt1D_cum = hist_evt1D.GetCumulative(0)
                    hist_evt1D_cum.Write()
                    can_evtbkg = TCanvas('can_evtbkg_E{0}_Z{1}'.format(ix, iz), str_label_com, 800, 600)                    
                    can_evtbkg.cd()
                    hist_evt1D_cum.Draw()
                    hist_bkg1D.Scale(hist_evt1D_cum.GetMaximum()/hist_bkg1D.GetMaximum())
                    hist_bkg1D.Draw("same")
                    can_evtbkg.SetLogy()
                    can_evtbkg.Write()

                    graph_evtrate = TGraphErrors()
                    graph_evtrate.SetName('graph_evtrate_E{0}_Z{1}'.format(ix, iz))
                    graph_evtrate.SetTitle(str_label_com)

                    graph_cut = TGraphErrors()
                    graph_cut.SetName('graph_cut_E{0}_Z{1}'.format(ix, iz))
                    graph_cut.SetTitle(str_label_com)

                    graph_energyrate = TGraph()
                    graph_energyrate.SetName('graph_energyrate_E{0}_Z{1}'.format(ix, iz))
                    graph_energyrate.SetTitle(str_label_com)
                    
                    graph_evtdens = TGraphErrors()
                    graph_evtdens.SetName('graph_evtdens_E{0}_Z{1}'.format(ix, iz))
                    graph_evtdens.SetTitle(str_label_com)

                    graph_energydens = TGraph()
                    graph_energydens.SetName('graph_energydens_E{0}_Z{1}'.format(ix, iz))
                    graph_energydens.SetTitle(str_label_com)

                    can_evtrate = TCanvas('can_evtrate_E{0}_Z{1}'.format(ix, iz), str_label_com, 800, 600)
                    can_energyrate = TCanvas('can_energyrate_E{0}_Z{1}'.format(ix, iz), str_label_com, 800, 600)
                    can_evtdens = TCanvas('can_evtdens_E{0}_Z{1}'.format(ix, iz), str_label_com, 800, 600)
                    can_energydens = TCanvas('can_energydens_E{0}_Z{1}'.format(ix, iz), str_label_com, 800, 600)

                    hist_covered.SetBinContent(ix, iz, 1)
                    flag_zerodiv = False

                    acc_ideal = 0

                    for iy in range(1, 1+hist_exp1D.GetNbinsX()): # BDT
                        bdt_cen = hist_exp1D.GetXaxis().GetBinCenter(iy)
                        bdt_lo = hist_exp1D.GetXaxis().GetBinLowEdge(iy)
                        #bdt_up = hist_exp1D.GetXaxis().GetBinLowEdge(iy+1)
                        exp = hist_exp1D.GetBinContent(iy)
                        exp_err = hist_exp1D.GetBinError(iy)
                        acc = hist_acc1D.GetBinContent(iy)
                        acc_err = hist_acc1D.GetBinError(iy)
                        aeff = acc / (2.*math.pi*(cth_up-cth_lo))
                        aeff_err = acc_err / (2.*math.pi*(cth_up-cth_lo))
                        evt_err = ROOT.Double(0)
                        evt = hist_evt1D.IntegralAndError(hist_evt1D.GetXaxis().FindBin(bdt_lo), hist_evt1D.GetNbinsX()+1, evt_err)

                        if exp>0:
                            #wlvtime = exp / aeff
                            evtenergy = evt*(10**energy_cen)
                            evtrate = evt / exp * aeff
                            evtrate_err = evt_err / exp * aeff
                            energyrate = evtenergy / exp / (10**energy_up/10**energy_lo/math.e) * aeff
                            evtdens = evt / exp
                            if evt>0:
                                evtdens_err = evtdens * sqrt( (evt_err/evt)**2 + (exp_err/exp)**2 )
                            else:
                                evtdens_err = 0
                            energydens = evtenergy / exp / (10**energy_up/10**energy_lo/math.e)
                            #energyrate = evtenergy / wlvtime / (10**energy_up/10**energy_lo/math.e) 
                            #evtrate = evt / exp
                            #energyrate = evtenergy / exp / (10**energy_up/10**energy_lo/math.e)
                        elif evt==0:
                            evtrate = 0
                            evtrate_err = 0
                            energyrate = 0
                            evtdens = 0
                            evtdens_err = 0
                            energydens = 0
                        else:
                            str_zerodiv = '''   BDT: {0}
   Event number: {1}
   Exposure: {2}'''.format(bdt_lo, evt, exp)
                            if flag_zerodiv is False:
                                print str_zerodiv
                                flag_zerodiv = True
                            continue
                        graph_evtrate.SetPoint(iy-1, acc, evtrate)
                        if iy==1:
                            acc_ideal = acc
                        graph_evtrate.SetPointError(iy-1, acc_err, evtrate_err)
                        graph_energyrate.SetPoint(iy-1, acc, energyrate)
                        graph_evtdens.SetPoint(iy-1, acc, evtdens)
                        graph_evtdens.SetPointError(iy-1, acc_err, evtdens_err)
                        graph_energydens.SetPoint(iy-1, acc, energydens)
                        graph_cut.SetPoint(iy-1, acc, bdt_lo)
                        graph_cut.SetPointError(iy-1, acc_err, 0)
                    graph_evtrate.Write()
                    graph_energyrate.Write()
                    graph_evtdens.Write()
                    graph_energydens.Write()
                    graph_cut.Write()
                    
                    can_evtrate.cd()
                    graph_evtrate.Draw("AL")
                    graph_evtrate.GetXaxis().SetRangeUser(0, acc_ideal*0.9)
                    graph_evtrate.GetYaxis().SetRangeUser(1e-6, 0.01)
                    graph_evtrate.GetXaxis().SetTitle("Acceptance [m^{2} sr]")
                    graph_evtrate.GetYaxis().SetTitle("Events [counts sr^{-1} s^{-1}]")
                    can_evtrate.SetLogy()
                    can_evtrate.Write()

                    can_energyrate.cd()
                    graph_energyrate.Draw("AL")
                    graph_energyrate.GetXaxis().SetRangeUser(0, 0.5)
                    graph_energyrate.GetYaxis().SetRangeUser(1e-3, 1e3)
                    graph_energyrate.GetXaxis().SetTitle("Acceptance [m^{2} sr]")
                    graph_energyrate.GetYaxis().SetTitle("E^{2} #times Differencial energy rate [MeV sr^{-1} s^{-1}]")
                    can_energyrate.SetLogy()
                    can_energyrate.Write()

                    can_evtdens.cd()
                    graph_evtdens.Draw("AL")
                    graph_evtdens.GetXaxis().SetRangeUser(0, 0.5)
                    graph_evtdens.GetYaxis().SetRangeUser(evtdens_min, evtdens_max)
                    graph_evtdens.GetXaxis().SetTitle("Acceptance [m^{2} sr]")
                    graph_evtdens.GetYaxis().SetTitle("Events [counts m^{-2} sr^{-1} s^{-1}]")
                    gaxis_ebl = TGaxis(graph_evtdens.GetXaxis().GetXmax(),evtdens_min,graph_evtdens.GetXaxis().GetXmax(),evtdens_max,evtdens_min_ebl,evtdens_max_ebl,510,"GL+")
                    #gaxis_ebl = TGaxis(can_evtdens.GetUxmax(),evtdens_min,can_evtdens.GetUxmax(),evtdens_max,evtdens_min_ebl,evtdens_max_ebl,510,"R-")
                    gaxis_ebl.SetTitle("[EBL level]")
                    gaxis_ebl.Draw()
                    can_evtdens.SetGridy(1)
                    can_evtdens.SetLogy()
                    can_evtdens.Write()

                    can_energydens.cd()
                    graph_energydens.Draw("AL")
                    graph_energydens.GetXaxis().SetRangeUser(0, 0.5)
                    graph_energydens.GetYaxis().SetRangeUser(1, 1e3)
                    graph_energydens.GetXaxis().SetTitle("Acceptance [m^{2} sr]")
                    graph_energydens.GetYaxis().SetTitle("E^{2} #times Differencial energy rate [MeV m^{-2} sr^{-1} s^{-1}]")
                    can_energydens.SetLogy()
                    can_energydens.Write()
    dir_out.cd()
    hist_covered.Write()


@click.command()
@click.argument('eventhist', type=str)
@click.argument('exposure', type=str)
@click.argument('roc', type=str)
@click.option('--outpath', '-o', default=None)
@click.option('--logemin', default=4.35, type=float)
@click.option('--logemax', default=5.75, type=float)
@click.option('--cthmin', default=0.6, type=float)
@click.option('--cthmax', default=1.0, type=float)
@click.option('--negativecth', '-n', is_flag=True)
#@click.option('--bsub', '-s', is_flag=True, default=False)
def main(eventhist, exposure, roc, outpath, logemin, logemax, cthmin, cthmax, negativecth):
    # Input
    file_evt = TFile(eventhist, "READ")
    hist_evt = file_evt.Get("Events")
    file_exp = TFile(exposure, "READ")
    hist_exp = file_exp.Get("exp_cth_cut")
    file_roc = TFile(roc, "READ")
    hist_acc = file_roc.Get("sig_acc_cth_cut")
    hist_bkg = file_roc.Get("bkg_rate_cth_cut")

    # Output
    file_out = TFile(outpath, "RECREATE")
    file_out.cd()

    create_galoff_roc3D(file_out, hist_evt, hist_exp, hist_acc, negativecth, logemin, logemax, cthmin, cthmax, hist_bkg)


if __name__ == '__main__':
    main()
