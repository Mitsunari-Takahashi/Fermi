#!/usr/bin/env python

import sys
#import numpy
#import math
#from math import cos, sin, tan, acos, asin, atan, radians, degrees, pi
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TRandom2
ROOT.gROOT.SetBatch()


def SimulateBinominalWaitInterval(mu=1., duration=42, samplerate=100, nsim=100, outpath='./waiting_simulations.root'):
    prob=mu/float(samplerate)
    fileout = ROOT.TFile(outpath,'RECREATE')
    dirhtg = fileout.mkdir('Histograms')
    dirhtg.cd()
    hs = ROOT.THStack('hs','{0} times simulations'.format(nsim))
    htgs = []
    for isim in range(nsim):
        gen = TRandom2()
        gen.SetSeed(0);
        htg.append(ROOT.TH1F('htg{0:05d}'.format(isim), '{0}th simulation;log(Waiting);intervals'.format(isim), -3, 2))
        time0 = 0.
        for iinst in range(duration*samlerate):
            suc = gen.Binominal(1, prob)
            if bool(suc)==True:
                time1 = float(iinst)/float(samlerate)
                htg[-1].Fill(log10(time1-time0))
                time0 = time1
        htg[-1].Write()
        hs.Add(htg[-1])
    fileout.cd()
    hs.Write()


@click.command()
@click.option('--mu', '-m', type=float, default=1.)
@click.option('--duration', '-d', type=int, default=42)
@click.option('--samplerate', '-s', type=int, default=100)
@click.option('--nsim', type=int, default=100)
@click.option('--outpath', type=str, default='./waiting_simulations.root')
#@click.option('--language', type=click.Choice(['Japanese', 'English']))
def main(mu, duration, samplerate, nsim, outpath):
    SimulateBinominalWaitInterval(mu, duration, samplerate, nsim, outpath)

    
if __name__ == '__main__':
    main()
