#!/usr/bin/env python
"""For showing MVA variables of one event.
"""
import sys
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TTree
ROOT.gROOT.SetBatch()
from pVariablesForMVA import DCT_ALIAS, TPL_VAR


def ShowEventVariable(pathfilein, nametr, nrunid, nevtid, comparable):
    fileIn = ROOT.TFile(pathfilein)
    trIn = fileIn.Get(nametr)
    print 'Shown:', TPL_VAR
    print trIn.GetName(), 'is found.'
    str_scan = ''
    for key, formula in DCT_ALIAS.items():
        trIn.SetAlias(key, formula)
    for var in TPL_VAR:
        str_scan = str_scan + var + ':'
    str_scan = str_scan[:-1]
    trIn.Scan(str_scan, 'EvtRun=={0} && EvtEventId=={1}'.format(nrunid, nevtid))


@click.command()
@click.argument('pathfilein')
@click.option('--nametr', default='MeritTuple')
@click.argument('nrunid', type=int)
@click.argument('nevtid', type=int)
@click.option('--comparable', default='')
def main(pathfilein, nametr, nrunid, nevtid, comparable):
    ShowEventVariable(pathfilein, nametr, nrunid, nevtid, comparable)


if __name__ == '__main__':
    main()
