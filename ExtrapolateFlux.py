#!/usr/bin/env python

import ROOT
from ROOT import TGraph2D
from ROOT import TH2D
import click
from pEBLatten import EBLmodel
ROOT.gROOT.SetBatch()


def ExtrapolateFlux(energy, dnde, pwlindex, binw=0.1, enr_min=4.25, nbin=16, redshift=0.):
    fc_dnde = ROOT.TF1('fc_dnde_pwl', '[0]*(x/[1])**[2]', 100, 1000000) # x = Energy in MeV (linear)
    fc_dnde.FixParameter(0, dnde)
    fc_dnde.FixParameter(1, energy)
    fc_dnde.FixParameter(2, pwlindex)
    
    LST_BINEDGE_LOW = [ ibin*binw+enr_min for ibin in range(nbin) ]
    lst_flux_itgl = []

    #EBL
    ebl = EBLmodel('Dominguez_etal_2011', '/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/GRB/160509A/CalOnly/Likelihood/10000-177828MeV_770-8233sec/EBL_Dominguez_etal_2011.root')

    for (ibin, edgelow) in enumerate(LST_BINEDGE_LOW):
        low = 10**edgelow
        up = 10**(edgelow+binw)
        if redshift==0.:
            atten=1.
        else:
            e_atten = edgelow+binw/2.0-3.0
            atten = ROOT.TMath.Exp(-ebl.interpolate(e_atten, redshift))
        lst_flux_itgl.append(fc_dnde.Integral(low, up)*atten)
    return lst_flux_itgl


@click.command()
@click.argument('energy', type=float)
@click.argument('dnde', type=float)
@click.option('--pwlindex', type=float, default=-2.)
@click.option('--redshift', '-z', type=float, default=0.)
def main(energy, dnde, pwlindex, redshift):
    """Extrapolates the Power-law function of dN/dE and returns an array of integral flux.
"""
    lst_flux_itgl = ExtrapolateFlux(energy, dnde, pwlindex, redshift)
    print lst_flux_itgl


if __name__ == '__main__':
    main()
