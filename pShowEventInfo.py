#!/usr/bin/env python
"""For showing contents of one event.
"""
import sys
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TTree
ROOT.gROOT.SetBatch() 


@click.command()
@click.argument('pathfilein')
@click.argument('nametr')
@click.argument('nrunid', type=int)
@click.argument('nevtid', type=int)
def main(pathfilein, nametr, nrunid, nevtid):
    fileIn = ROOT.TFile(pathfilein)
    trIn = fileIn.Get(nametr)
    print trIn.GetName(), 'is found.'
    for iEvt in range(trIn.GetEntries()):
        trIn.GetEntry(iEvt)
        if trIn.RUN_ID==nrunid and trIn.EVENT_ID==nevtid:
            trIn.Show(iEvt)


if __name__ == '__main__':
    main()
