#!/usr/bin/env python
"""For showing MVA variables of one event.
"""
import sys
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TTree
ROOT.gROOT.SetBatch() 


DCT_VAR = {'log10CalNewCfpCalSelChiSq':'1*log10(CalNewCfpCalSelChiSq)',
           'log10Cal1TransRms':'log10(Cal1FitChiSquare)',
           'CalNewCfpCalTmax':'CalNewCfpCalTmax',
           'CalBkHalfRatio':'CalBkHalfRatio',
           'CalELayer74RatioLog':'log10(max(-5, CalELayer7))-log10(max(-5, CalELayer4))',
           'Acd2Cal1Energy15Log':'log10(max(Acd2Cal1Energy15,1E-6))',
           'Acd2VetoCountLog':'log10(max(Acd2VetoCount,3E-1))',
           'Acd2Cal1VetoSigmaHitLog':'log10(max(Acd2Cal1VetoSigmaHit,1E-3))',
           'CalELayerCorrInitialRatioLog':'CalELayerCorrInitialRatioLog',
           'CalELayer34afterInitialRatioLog':'CalELayer34afterInitialRatioLog',
           'CalTrSizeCalT95':'CalTrSizeCalT95',
           'log10Cal1FitChiSquare':'log10(max(-5, Cal1FitChiSquare))',
           'Cal1MomNumCoreXtalsFract':'Cal1MomNumCoreXtals/Cal1NumXtals',
           'Acd2TileEnergyRatioLog':'Acd2TileEnergyRatioLog',
           'CalEdgeEnergyLog':'log10(max(-5, CalEdgeEnergy))',           
           }


def ShowEventVariable(pathfilein, nametr, nrunid, nevtid):
    fileIn = ROOT.TFile(pathfilein)
    trIn = fileIn.Get(nametr)
    print trIn.GetName(), 'is found.'
    str_scan = ''
    for key, formula in DCT_VAR.items():
        trIn.SetAlias(key, formula)
        str_scan = str_scan + key + ':'
    str_scan = str_scan[:-1]
    trIn.Scan(str_scan, 'EvtRun=={0} && EvtEventId=={1}'.format(nrunid, nevtid))


@click.command()
@click.argument('pathfilein')
@click.option('--nametr', default='MeritTuple')
@click.argument('nrunid', type=int)
@click.argument('nevtid', type=int)
def main(pathfilein, nametr, nrunid, nevtid):
    ShowEventVariable(pathfilein, nametr, nrunid, nevtid)


if __name__ == '__main__':
    main()
