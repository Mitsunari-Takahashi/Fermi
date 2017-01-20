#!/usr/bin/env python

import ROOT
import click
ROOT.gROOT.SetBatch()


def ExtrapolateFlux(energy, dnde, pwlindex):
    fc_dnde = ROOT.TF1('fc_dnde_pwl', '[0]*(x/[1])**(-[2])', 100, 1000000) # x = Energy in MeV (linear)
    fc_dnde.FixParameter(0, dnde)
    fc_dnde.FixParameter(1, energy)
    fc_dnde.FixParameter(2, pwlindex)
    
    BINW = 0.1
    ENR_MIN = 4.25
    NBIN = 16
    LST_BINEDGE_LOW = [ ibin*BINW+ENR_MIN for ibin in range(NBIN) ]
    lst_flux_itgl = []

    for (ibin, edgelow) in enumerate(LST_BINEDGE_LOW):
        lst_flux_itgl.append(fc_dnde.Integral(10**edgelow, 10**(edgelow+BINW)))
    return lst_flux_itgl


@click.command()
@click.argument('energy', type=float)
@click.argument('dnde', type=float)
@click.option('--pwlindex', type=float, default=2.)
def main(energy, dnde, pwlindex):
    """Extrapolates the Power-law function of dN/dE and returns an array of integral flux.
"""
    lst_flux_itgl = ExtrapolateFlux(energy, dnde, pwlindex)
    print lst_flux_itgl


if __name__ == '__main__':
    main()
