#!/usr/bin/env python

import sys
import numpy as np
import click
from astropy.io import fits
import astropy.units as u
from math import sqrt
import commands


def get_entries(pathFileEvt, metStart, metStop):
    """Return entry number between metStart and metStop.
"""
    hdulistEVT = fits.open(pathFileEvt)
    tbdataEVT = hdulistEVT[1].data
    aTIME = tbdataEVT.field('TIME')
    n = 0
    for t in range(len(aTIME)):
        if t>=metStart and t<metStop:
            n+=1
    return n


def find_goodstat_periods(pathFileEvt, metStart, metStop, nthreshold, noffset=0): #, ethreshold=100., evtclass=128):
    """Look over spacecraft files and find time intervals which have more than N events.
"""
    periods_goodstat = [metStart]
    nstat = 0
    hdulistEVT = fits.open(pathFileEvt)
    tbdataEVT = hdulistEVT[1].data
    aTIME = tbdataEVT.field('TIME')
    #aENERGY = tbdataEVT.field('ENERGY')
    #aEVTCLASS = tbdataEVT.field('EVENT_CLASS')
    aTIME.sort()
    nEVT = len(aTIME)
    time_prev = 0
    print "  ", pathFileEvt, "(", nEVT, "events )"
    for iEVT in range(nEVT):
        if aTIME[iEVT]<time_prev:
            print 'Odd order!!!'
            return 1
        if aTIME[iEVT]>=metStart and aTIME[iEVT]<metStop: # and aEVTCLASS[iEVT]>=evtclass and aENERGY[iEVT]>ethreshold:
            nstat+=1
            time_fill = sqrt(aTIME[iEVT]*aTIME[iEVT+1])
            if (nstat-noffset)/max(1., pow(sqrt(aTIME[iEVT]*aTIME[iEVT+1])-periods_goodstat[-1], 1./4.)) >= nthreshold and int(time_fill)>=int(periods_goodstat[-1]) and time_fill-periods_goodstat[-1]>=1:
                periods_goodstat.append(time_fill)
                nstat = 0
        time_prev = aTIME[iEVT]
#    if nstat<nthreshold:
    periods_goodstat = periods_goodstat[:-1]
    periods_goodstat.append(metStop)
    return periods_goodstat


@click.command()
@click.argument('evtfile', type=str)
@click.argument('start', type=float)
@click.argument('stop', type=float)
@click.option('--threshold', '-t', type=int, default=10)
#@click.option('--ethreshold', '-e', type=float, default=100.)
#@click.option('--evtclass', '-c', type=int, default=128)
def main(evtfile, start, stop, threshold): #, ethreshold, evtclass):
    gti = find_goodstat_periods(evtfile, start, stop, threshold) #, ethreshold, evtclass)
    print gti


if __name__ == '__main__':
    main()
