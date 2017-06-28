#!/usr/bin/env python
import sys
import click
from math import log10, acos, atan2, pi, degrees, sin, cos, asin, radians
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TTree, TH1
from pColor import *
ROOT.gROOT.SetBatch() 


def PlotPSFandAcceptance(path_perf, emin, emax, path_out):
    file_perf = ROOT.TFile(path_perf, 'READ')
    htg_acc = file_perf.Get('acc_cth_hist')
    htg_psf = file_perf.Get('psf_cth_q68_hist')
    file_out = ROOT.TFile(path_out, 'RECREATE')
    file_out.cd()
    # Acceptance
    pf2d_acc = ROOT.TProfile2D(htg_acc.GetName().replace('_hist', '_prof2d'), htg_acc.GetTitle().replace('_hist', '_prof2d'), htg_acc.GetNbinsX(), htg_acc.GetXaxis().GetXmin(), htg_acc.GetXaxis().GetXmax(), htg_acc.GetNbinsY(), htg_acc.GetYaxis().GetXmin(), htg_acc.GetYaxis().GetXmax())
    for ix in range(1, pf2d_acc.GetNbinsX()+1):
        for iy in range(1, pf2d_acc.GetNbinsY()+1):
            pf2d_acc.Fill(pf2d_acc.GetXaxis().GetBinCenter(ix), pf2d_acc.GetYaxis().GetBinCenter(iy), htg_acc.GetBinContent(ix, iy), pow(10**htg_acc.GetXaxis().GetBinCenter(ix),-1))
    pf2d_acc.Write()
    pf_acc_projCth = pf2d_acc.ProjectionY(pf2d_acc.GetName()+'_projCth', pf2d_acc.GetXaxis().FindBin(emin), pf2d_acc.GetXaxis().FindBin(emax))
    pf_acc_projCth.SetTitle('Acceptance\, weighted\, by\, E^{{-1}}\, for\, {0} \leq log_{{10}}E [MeV] < {1};\cos{{\theta}};[m^2 sr]'.format(emin, emax))
    pf_acc_projCth.Write()
    pf_acc_projCth_cum = pf_acc_projCth.GetCumulative()
    for kx in range(1, 1+pf_acc_projCth_cum.GetNbinsX()):
        pf_acc_projCth_cum.SetBinError(kx,0)
    pf_acc_projCth_cum.Write()
    pf_acc_cum_divRoi = pf_acc_projCth_cum.Clone('pf_acc_cum_divRoi')
    pf_acc_cum_divRoi.SetTitle('Acceptance devided by solid angle of RoI')

    # PSF
    pf2d_psf = ROOT.TProfile2D(htg_psf.GetName().replace('_hist', '_prof2d'), htg_psf.GetTitle().replace('_hist', '_prof2d'), htg_psf.GetNbinsX(), htg_psf.GetXaxis().GetXmin(), htg_psf.GetXaxis().GetXmax(), htg_psf.GetNbinsY(), htg_psf.GetYaxis().GetXmin(), htg_psf.GetYaxis().GetXmax())
    for ix in range(1, pf2d_psf.GetNbinsX()+1):
        for iy in range(1, pf2d_psf.GetNbinsY()+1):
            pf2d_psf.Fill(pf2d_psf.GetXaxis().GetBinCenter(ix), pf2d_psf.GetYaxis().GetBinCenter(iy), htg_psf.GetBinContent(ix, iy), htg_acc.GetBinContent(ix, iy)*pow(10**htg_acc.GetXaxis().GetBinCenter(ix),-1))
    pf2d_psf.Write()
    pf_psf_projCth = pf2d_psf.ProfileY(pf2d_psf.GetName()+'_projCth', pf2d_psf.GetXaxis().FindBin(emin), pf2d_psf.GetXaxis().FindBin(emax))
    pf_psf_projCth.Write()
    htg_psfSq = ROOT.TH1D('htg_cos_psf_q68_cth', '1-\cos{{(PSF68)}} \, weighted\, by\, acceptance\, and\, E^{{-1}}\, for\, {0} \leq log_{{10}}E [MeV] < {1};\cos{{\theta}};Solid angle [sr]'.format(emin, emax), htg_psf.GetNbinsY(), htg_psf.GetYaxis().GetXmin(), htg_psf.GetYaxis().GetXmax())
    for jx in range(1, 1+htg_psfSq.GetNbinsX()):
        htg_psfSq.SetBinContent(jx, (1-cos(radians(pf_psf_projCth.GetBinContent(jx))))*2.*pi)
        #htg_psfSq.SetBinError(jx, 2.*pf_psf_projCth.GetBinContent(jx)*pf_psf_projCth.GetBinError(jx))
    htg_psfSq.Write()
    pf_acc_cum_divRoi.Divide(htg_psfSq)
    for lx in range(1, 1+pf_acc_cum_divRoi.GetNbinsX()):
        pf_acc_cum_divRoi.SetBinError(lx,0)
    pf_acc_cum_divRoi.Write()

    htg_acc_psf = ROOT.TH1D('htg_acc_psf', 'Acceptance weighted by E^{-1} vs. PSF68 cut', 100, 0, 10)
    for ix in range(1, pf2d_acc.GetNbinsX()+1):
        for iy in range(1, pf2d_acc.GetNbinsY()+1):
            htg_acc_psf.Fill(htg_psf.GetBinContent(ix, iy), htg_acc.GetBinContent(ix, iy)*pow(10**htg_acc.GetXaxis().GetBinCenter(ix),-1))
    htg_acc_psf.Write()
    htg_acc_psf_cum = htg_acc_psf.GetCumulative()
    htg_acc_psf_cum.Write()
    htg_acc_psf_cum_divRoi = htg_acc_psf_cum.Clone('{0}_divRoI'.format(htg_acc_psf_cum.GetName()))
    htg_acc_psf_cum_divRoi.SetTitle('Weighted acceptance / solid angle')
    for ibin in range(1, 1+htg_acc_psf_cum_divRoi.GetNbinsX()):
        htg_acc_psf_cum_divRoi.SetBinContent(ibin, htg_acc_psf_cum.GetBinContent(ibin)/((1-cos(radians(htg_acc_psf_cum_divRoi.GetBinCenter(ibin))))*2.*pi))
    htg_acc_psf_cum_divRoi.Write()

    
@click.command()
@click.argument('perf', type=str)
@click.argument('emin', type=float)
@click.argument('emax', type=float)
@click.option('--output', '-o', type=str, default='./FIG_PSF_and_Acceptance.root')
def main(perf, emin, emax, output):
    PlotPSFandAcceptance(perf, emin, emax, output)

    
if __name__ == '__main__':
    main()
