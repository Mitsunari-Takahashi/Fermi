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
@click.argument('energy', type=float)
@click.argument('costheta', type=float)
@click.option('--momcntrx', '-x', type=float, default=None)
@click.option('--momcntry', '-y', type=float, default=None)
@click.option('--tol', '-t', type=float, default=0.1)
@click.option('--maxnum', type=float, default=100)
@click.option('--output', '-o', type=str, default='./similar_events.txt')
def main(pathfilein, energy, costheta, momcntrx, momcntry, tol, maxnum, output):
    fileIn = ROOT.TFile(pathfilein)
    trIn = fileIn.Get('MeritTuple')
    print trIn.GetName(), 'is found.'
    with open(output, 'w') as ou:
        print 'Output file', output, 'is opened.'
        ou.write("# RunID   EvtID   RECON_File   MERIT_File")
        ncount = 0
        for iEvt in range(trIn.GetEntries()):
            trIn.GetEntry(iEvt)
            if energy>=trIn.WP8CalOnlyEnergy*(1.-tol) and energy<trIn.WP8CalOnlyEnergy*(1.+tol):
                if costheta>=trIn.Cal1MomZDir*(1.-0.1) and costheta<trIn.Cal1MomZDir*(1.+0.1):
                    if ( momcntrx is None or (abs(momcntrx)>=abs(trIn.Cal1MomXCntrCor)*(1.-tol) and abs(momcntrx)<abs(trIn.Cal1MomXCntrCor)*(1.+tol)) ) or ( momcntry is None or (abs(momcntry)>=abs(trIn.Cal1MomYCntrCor)*(1.-tol) and abs(momcntry)<abs(trIn.Cal1MomYCntrCor)*(1.+tol)) ) :
                        ou.write("""
{run} {evt} root://glast-rdr//glast/mc/ServiceChallenge/AG-GR-v20r09p09-OVL6p2/recon/AG-GR-v20r09p09-OVL6p2-{run:0>6}-recon.root root://glast-rdr//glast/mc/ServiceChallenge/AG-GR-v20r09p09-OVL6p2/merit/AG-GR-v20r09p09-OVL6p2-{run:0>6}-merit.root""".format(run=trIn.EvtRun, evt=trIn.EvtEventId))
                        ncount += 1
            if ncount>maxnum:
                print 'Number of selected events reached {0}.'.format(ncount)
                break
                    
if __name__ == '__main__':
    main()
