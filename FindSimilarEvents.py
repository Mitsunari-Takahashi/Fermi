#!/usr/bin/env python
"""Finds events similar to your event.
"""
import sys
import math
import click
import healpy as hp
from healpy import pixelfunc as hppf
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TTree
ROOT.gROOT.SetBatch() 
from pFindHEALPix import *


@click.command()
@click.argument('pathfilein')
@click.argument('energy', type=float)
@click.argument('costheta', type=float)
@click.option('--galoff', is_flag=True)
@click.option('--calonly', is_flag=True)
@click.option('--momcntrx', '-x', type=float, default=None)
@click.option('--momcntry', '-y', type=float, default=None)
@click.option('--tol', '-t', type=float, default=0.1)
@click.option('--maxnum', type=float, default=10)
@click.option('--output', '-o', type=str, default='./similar_events.txt')
def main(pathfilein, energy, costheta, momcntrx, momcntry, galoff, calonly, tol, maxnum, output):

    aHpxGalOFF = None
    nhpside = 16
    if galoff is True:
        pathCatalogue = "/nfs/farm/g/glast/u/mtakahas/FermiAnalysis/Catalogues/gll_psch_v13.fit" #3FHL
        aHpxGalOFF = find_galoff_healpxs(nhpside, 0, pathCatalogue)
        print 'Galactic OFF:', aHpxGalOFF

    with open(pathfilein, 'r') as inp:
        print 'Input file', pathfilein, 'is opened.'
        line = inp.readline()

        with open(output, 'w') as ou:
            print 'Output file', output, 'is opened.'
            ou.write("# RunID   EvtID   RECON_File   MERIT_File")

            ncount = 0
            flag_break = False
            while line and flag_break is False:
                line.decode('unicode-escape')
                merit = line[:-1]
                print merit
                fileIn = ROOT.TFile.Open(merit)
                print fileIn.GetName(), 'is opened.'
                trIn = fileIn.Get('MeritTuple')
                print trIn.GetName(), 'is found.'
                for iEvt in range(trIn.GetEntries()):
                    trIn.GetEntry(iEvt)
                    if energy>=trIn.WP8CalOnlyEnergy*(1.-tol) and energy<trIn.WP8CalOnlyEnergy*(1.+tol) and costheta>=trIn.Cal1MomZDir*(1.-0.1) and costheta<trIn.Cal1MomZDir*(1.+0.1) and ( ( momcntrx is None or (abs(momcntrx)>=abs(trIn.Cal1MomXCntrCor)*(1.-tol) and abs(momcntrx)<abs(trIn.Cal1MomXCntrCor)*(1.+tol)) ) or ( momcntry is None or (abs(momcntry)>=abs(trIn.Cal1MomYCntrCor)*(1.-tol) and abs(momcntry)<abs(trIn.Cal1MomYCntrCor)*(1.+tol)) ) ) and trIn.FT1CalZenithTheta<=90.:
                        if calonly is False or (trIn.TkrNumTracks==0 or (math.log10(max(trIn.CalTrackAngle,1E-4)) > (0.529795)*(trIn.EvtJointLogEnergy < 3.000000)  + ((1.0)*((0.529795)*(1.0)+(-1.379791)*(pow((trIn.EvtJointLogEnergy-3.000000)/0.916667,1))+(0.583401)*(pow((trIn.EvtJointLogEnergy-3.000000)/0.916667,2))+(-0.075555)*(pow((trIn.EvtJointLogEnergy-3.000000)/0.916667,3))))*(trIn.EvtJointLogEnergy  >= 3.000000 and trIn.EvtJointLogEnergy <= 5.750000) + (-0.398962)*(trIn.EvtJointLogEnergy >  5.750000)) ):
                            if aHpxGalOFF is not None:
                                nhpx_evt = hppf.ang2pix(nhpside, math.pi/2.-math.radians(trIn.FT1CalDec), math.radians(trIn.FT1CalRa))
                                print nhpx_evt
                                if nhpx_evt in aHpxGalOFF:
                                    ou.write("""
{run} {evt} {reco} {meri}""".format(run=trIn.EvtRun, evt=trIn.EvtEventId, reco=merit.replace('merit', 'recon'), meri=merit))
                                    print ncount
                                    ncount += 1
                            else:
                                ou.write("""
{run} {evt} {reco} {meri}""".format(run=trIn.EvtRun, evt=trIn.EvtEventId, reco=merit.replace('merit', 'recon'), meri=merit))
                                print ncount, 'events found.'
                                ncount += 1
                    if ncount>maxnum:
                        print 'Number of selected events reached {0}.'.format(ncount)
                        flag_break = True
                        break
                line = inp.readline()
            print ''
            print 'Finished.'

                    
if __name__ == '__main__':
    main()
