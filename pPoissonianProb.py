#!/usr/bin/env python
"""
For calucuating Possion coincidence probability.
python -m pPoissonianProb <lambda> <threshold>
"""
import sys
import math
import ROOT
ROOT.gROOT.SetBatch()


def calcProb(tp):
    """
Calcurate Poissonian coincidence probability.
pPoissonianProb.calcProb([<lambda>, <threshold>])
"""
    lamb = float(tp[0])
    thre = float(tp[1])
    maxi = thre*100.
#    fPoi = ROOT.TF1("fPoi", "TMath::Poisson(x, [0])", 0, maxi)
    fPoiI = ROOT.TF1("fPoi", "TMath::PoissonI(x, [0])", 0, maxi)
#    fPoi.SetParameter(0, lamb)
    fPoiI.SetParameter(0, lamb)
#    return [fPoi.Integral(thre, maxi)/fPoi.Integral(0, maxi), fPoiI.Integral(thre, maxi)/fPoiI.Integral(0, maxi)]
    return fPoiI.Integral(thre, maxi)/fPoiI.Integral(0, maxi)



if __name__ == '__main__':
    if len(sys.argv)==3:
        result = calcProb(sys.argv[1:])
        print result

