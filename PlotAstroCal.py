#!/usr/bin/env python

#import sys
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TH1, TDirectoryFile, TDirectory
from ROOT import kWhite, kBlack, kGray, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink
ROOT.gROOT.SetBatch()
import pColor

@click.command()
@click.argument('gamsrc', type=str)
@click.argument('model', type=str)
@click.argument('flight', type=str)
@click.option('--suffix', type=str, default='')
def main(gamsrc, flight, model, suffix):

    gStyle.SetOptStat(kFALSE);
    gStyle.SetPalette(54)
    gStyle.SetTextColor(kRed)
    gStyle.SetPaintTextFormat('+3.2f')

    if suffix != '':
        suffix = '_' + suffix
    STR_OUT = 'Compare_Flight_vs_Model_{0}{1}'.format(gamsrc, suffix)
    FILE_OUT = ROOT.TFile('{0}.root'.format(STR_OUT), 'RECREATE')

    TPL_CLASS = ('CalOnlyR100', 'CalOnlyR30', 'CalOnlyR10')

    FILE_FLIGHT = ROOT.TFile(flight, 'READ')
    DIR_FLIGHT = FILE_FLIGHT.GetDirectory(gamsrc)

    FILE_MODEL = ROOT.TFile(model, 'READ')

    CVS = ROOT.TCanvas('cvs_{0}'.format(gamsrc), 'Comparison of model vs. flight data of {0}'.format(gamsrc), 1600, 1200)
    CVS.Divide(2, 2)
    hs_devi = ROOT.THStack('hs_deviation_{0}'.format(gamsrc), 'Deviation of {0} flight data from model;log_{{10}}Energy [MeV];(N_{{obs}}-N_{{exp}})/N_{{exp}}'.format(gamsrc))
    lst_cvs_cth = []

    lst_htg_model_cth = []
    #lst_htg_model_cth_rebin = []
    lst_htg_model_all = []
    lst_htg_model_on = []
    lst_htg_cth = []
    lst_htg_flight = []
    lst_htg_flight_cth = []
    lst_htg_devi = []
    lst_htg_devi_cth = []
    leg = ROOT.TLegend(0.55, 0.55, 0.875, 0.85, 'Integral count distribution')
    leg_devi = ROOT.TLegend(0.2, 0.65, 0.4, 0.875, 'Class')

    for (iclas, clas) in enumerate(TPL_CLASS):
        print '=====', clas, '====='

        lst_cvs_cth.append(ROOT.TCanvas('cvs_cth_{0}_{1}'.format(gamsrc, clas), 'Comparison of model vs. flight data of {0} {1} depending on energy and cos#theta'.format(clas, gamsrc), 1600, 1200))
        lst_cvs_cth[-1].Divide(2, 2)

        lst_htg_flight_cth.append(DIR_FLIGHT.Get('hNumSig{0}_1_{1}'.format(gamsrc, iclas)))
        #lst_htg_flight_cth[-1].GetZaxis().SetRangeUser(0, lst_htg_flight_cth[-1].GetMaximum())
        lst_htg_flight_cth[-1].SetMarkerColor(kRed)
        lst_cvs_cth[-1].cd(1)
        lst_htg_flight_cth[-1].Draw('COLZ TEXTE')
        lst_htg_flight.append(lst_htg_flight_cth[-1].ProjectionX('{0}_projEnergy'.format(lst_htg_flight_cth[-1].GetName())))
        lst_htg_flight[-1].SetLineColor(kRed)
        lst_htg_flight[-1].SetLineWidth(2)
        lst_htg_flight[-1].SetFillStyle(0)
        lst_htg_flight[-1].SetMarkerColor(kRed)
        lst_htg_flight[-1].SetMarkerStyle(20)
        lst_htg_flight[-1].SetTitle('{0} {1}'.format(clas, gamsrc))
        CVS.cd(iclas+1)
        CVS.cd(iclas+1).SetGridy()
        lst_htg_flight[-1].Draw()
        lst_htg_flight[-1].GetYaxis().SetRangeUser(0, lst_htg_flight[-1].GetMaximum()*1.2)

        lst_htg_model_on.append(FILE_MODEL.Get('htg_sed_model_{0}_ON'.format(clas)))
        lst_htg_model_on[-1].SetLineColor(kWhite)
        lst_htg_model_on[-1].SetFillColor(kBlack)
        lst_htg_model_on[-1].SetFillStyle(3001)
        lst_htg_model_on[-1].SetMarkerColor(kBlack)
        lst_htg_model_on[-1].SetMarkerStyle(25)
        CVS.cd(iclas+1)
        CVS.cd(iclas+1).SetGridy()
        lst_htg_model_on[-1].Draw('same E2')

        lst_htg_model_cth.append(FILE_MODEL.Get('htg_sed_model_{0}'.format(clas)))
        lst_htg_model_cth[-1].SetTitle('Model of {0} {1} [counts];log_{{10}}Energy [MeV];cos#theta'.format(clas, gamsrc))
        lst_cvs_cth[-1].cd(2)
        lst_htg_model_cth[-1].Draw('colz')
        lst_htg_model_all.append(FILE_MODEL.Get('htg_sed_model_{0}_projEnergy'.format(clas)))
        lst_htg_model_all[-1].SetLineColor(kBlack)
        lst_htg_model_all[-1].SetFillColor(kWhite)
        lst_htg_model_all[-1].SetFillStyle(0)
        lst_htg_model_all[-1].SetMarkerColor(kBlack)
        lst_htg_model_all[-1].SetMarkerStyle(23)
        CVS.cd(iclas+1)
        lst_htg_model_all[-1].Draw('same E2')

        lst_htg_devi_cth.append(lst_htg_flight_cth[-1].Clone('htg_deviation_cth_{0}'.format(clas)))
        lst_htg_devi_cth[-1].SetMarkerColor(kRed)
        htg_model_cth_rebin = lst_htg_model_cth[-1].RebinY(lst_htg_model_cth[-1].GetNbinsY()/lst_htg_devi_cth[-1].GetNbinsY(), '{0}_rebinY'.format(lst_htg_model_cth[-1].GetName()))
        lst_htg_devi_cth[-1].Add(htg_model_cth_rebin, -1)
        lst_htg_devi_cth[-1].Divide(htg_model_cth_rebin)
        lst_htg_devi_cth[-1].GetYaxis().SetRangeUser(0.2,1)
        lst_htg_devi_cth[-1].GetZaxis().SetRangeUser(-2,2)
        lst_htg_devi_cth[-1].SetTitle('Fractional deviation of {0} {1} from model'.format(clas, gamsrc))
        lst_cvs_cth[-1].cd(3)
        lst_htg_devi_cth[-1].Draw('COLZ TEXTE')

        lst_htg_devi.append(lst_htg_flight[-1].Clone('htg_deviation_{0}'.format(clas)))
        lst_htg_devi[-1].Add(lst_htg_model_on[-1], -1)
        lst_htg_devi[-1].Divide(lst_htg_model_on[-1])
        lst_htg_devi[-1].SetTitle('{0} {1}'.format(clas, gamsrc))
        lst_htg_devi[-1].SetYTitle('Fractional deviation')
        lst_htg_devi[-1].SetLineColor(pColor.akColor(iclas))
        lst_htg_devi[-1].SetLineWidth(3-iclas)
        lst_htg_devi[-1].SetFillStyle(0)
        lst_htg_devi[-1].SetMarkerColor(pColor.akColor(iclas))
        lst_htg_devi[-1].SetMarkerStyle(pColor.aakMarkerStyle(1, iclas))
        leg_devi.AddEntry(lst_htg_devi[-1], clas, 'lp')
        hs_devi.Add(lst_htg_devi[-1])

        if iclas==0:
            leg.AddEntry(lst_htg_flight[-1], 'Observed excess', 'lp')
            leg.AddEntry(lst_htg_model_on[-1], 'Expected ON', 'lpf')
            leg.AddEntry(lst_htg_model_all[-1], 'Expected all', 'lpf')
        CVS.cd(iclas+1)
        leg.Draw('same')
        FILE_OUT.cd()
        lst_cvs_cth[-1].Write()        
        lst_cvs_cth[-1].SaveAs('{0}_cth_{1}.png'.format(STR_OUT, clas))

    CVS.cd(4)
    CVS.cd(4).SetGridy()
    fc_zero = ROOT.TF1('fc_zero', '0', 4.35, 5.75)
    fc_zero.SetLineColor(kGreen+1)
    fc_zero.SetLineWidth(2)
    hs_devi.Draw('nostack E1')
    fc_zero.Draw('same')
    leg_devi.Draw('same')
    
    FILE_OUT.cd()
    CVS.Write()
    CVS.SaveAs('{0}.png'.format(STR_OUT))
    print 'Done.'

if __name__ == '__main__':
    main()
