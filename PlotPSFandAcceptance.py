#!/usr/bin/env python
import sys
import click
from math import log10, acos, atan2, pi, degrees, sin, cos, asin
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TTree
from pColor import *
ROOT.gROOT.SetBatch() 


def PlotPSFandAcceptance(path_perf, emin, emax, path_out):
    file_perf = ROOT.TFile(path_perf, 'READ')
    htg_acc = file_perf.Get('acc_cth_hist')
    htg_psf = file_perf.Get('psf_cth_q68_hist')
    file_out = ROOT.TFile(path_out, 'RECREATE')
    file_out.cd()
    # Acceptance
    pf2d_acc = ROOT.TProfile2D(htg_acc.GetName().replace('_hist', '_prof2d'), htg_acc.GetTitle().replace('_hist', '_prof2d'), htg_acc.GetNbinsX(), htg_acc.GetXmin(), htg_acc.GetXmax(), htg_acc.GetNbinsY(), htg_acc.GetYmin(), htg_acc.GetYmax())
    for ix in range(1, pf2d_acc.GetNbinsX()+1):
        for iy in range(1, pf2d_acc.GetNbinsY()+1):
            pf2d_acc.Fill(pf2d_acc.GetBinCenter(ix), pf2d_acc.GetBinCenter(iy), htg_acc.GetBinContent(ix, iy), pow(10**htg_acc.GetBinCenter(ix),-1))
    pf2d_acc.Write()
    pf_acc_projCth = pf2d_acc.ProfileY(pf2d_acc.GetName()+'_projCth', 'Acceptance weighted E^{-1}', pf2d_acc.GetXaxis().FindBin(emin), pf2d_acc.GetXaxis().FindBin(emax))
    pf_acc_projCth.Write()
    
    # PSF
    pf2d_psf = ROOT.TProfile2D(htg_psf.GetName().replace('_hist', '_prof2d'), htg_psf.GetTitle().replace('_hist', '_prof2d'), htg_psf.GetNbinsX(), htg_psf.GetXmin(), htg_psf.GetXmax(), htg_psf.GetNbinsY(), htg_psf.GetYmin(), htg_psf.GetYmax())
    for ix in range(1, pf2d_psf.GetNbinsX()+1):
        for iy in range(1, pf2d_psf.GetNbinsY()+1):
            pf2d_psf.Fill(pf2d_psf.GetBinCenter(ix), pf2d_psf.GetBinCenter(iy), htg_psf.GetBinContent(ix, iy), htg_acc.GetBinContent(ix, iy)*pow(10**htg_acc.GetBinCenter(ix),-1))
    pf2d_psf.Write()
    pf_psf_projCth = pf2d_psf.ProfileY(pf2d_psf.GetName()+'_projCth', 'PSF68 weighted by acceptance and E^{-1}', pf2d_psf.GetXaxis().FindBin(emin), pf2d_psf.GetXaxis().FindBin(emax))
    pf_psf_projCth.Write()
    htg_psfSq = ROOT.TH1D('htg_psf_q68_squared_cth', 'PSF68^{2} weighted by acceptance and E^{-1}', htg_psf.GetNbinsY(), htg_psf.GetYmin(), htg_psf.GetYmax())
    for jx in range(1, 1+htg_psfSq.GetNbinsX()):
        htg_psfSq.SetBinContent(jx, pf_psf_projCth.GetBinContent(jx)**2)
        htg_psfSq.SetBinError(jx, 2.*pf_psf_projCth.GetBinContent(jx)*pf_psf_projCth.GetBinError(jx))
    htg_psfSq.Write()

    

@click.command()
@click.argument('perf', type=str)
@click.argument('emin', type=float)
@click.argument('emax', type=float)
@click.option('--output', '-o', type=str, default='./FIG_PSF_and_Acceptance.root')
def main(perf, emin, emax, output):
    PlotPSFandAcceptance(perf, emin, emax, output)

    
if __name__ == '__main__':
    main()
