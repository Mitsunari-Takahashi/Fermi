#!/usr/bin/env python

import sys
import numpy as np
import click
from astropy.io import fits
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
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


def get_entries_roi(pathFileEvt, tStart, tStop, rlim=0.0, ra=0, dec=0, torigin=0.0):
    """Return entry number between tStart - tStop from torigin, and within rlim degrees from (ra, dec).
"""
    metStart = tStart + torigin
    metStop = tStop + torigin
    hdulistEVT = fits.open(pathFileEvt)
    tbdataEVT = hdulistEVT[1].data
    aTIME = tbdataEVT.field('TIME')
    aRA = tbdataEVT.field('RA')
    aDEC = tbdataEVT.field('DEC')
    nEVT = len(aTIME)
    time_prev = 0
    coords_target = SkyCoord(ra, dec, unit="deg")
    nstat = 0
    for iEVT in range(1, nEVT):
        if aTIME[iEVT]<time_prev:
            print 'Odd order!!!'
            sys.exit(1)
        if aTIME[iEVT]>=metStart and aTIME[iEVT]<metStop:
            coords_event = SkyCoord(aRA[iEVT], aDEC[iEVT], unit="deg")
            sep_ang = coords_target.separation(coords_event)
            sep_deg = float(sep_ang.to_string(unit=u.deg, decimal=True))
            if sep_deg<=rlim:
                nstat+=1
        time_prev = aTIME[iEVT]
    return nstat


def find_goodstat_periods(pathFileEvt, tStart, tStop, nthreshold, rlim=0.0, ra=0, dec=0, torigin=0.0, noffset=0): #, ethreshold=100., evtclass=128):
    """Look over spacecraft files and find time intervals which have more than N events.
"""
    metStart = tStart + torigin
    metStop = tStop + torigin
    periods_goodstat = [[tStart]]
    nstat = 1
    hdulistEVT = fits.open(pathFileEvt)
    tbdataEVT = hdulistEVT[1].data
    aTIME = tbdataEVT.field('TIME')
    aRA = tbdataEVT.field('RA')
    aDEC = tbdataEVT.field('DEC')
    #aTIME.sort()
    nEVT = len(aTIME)
    time_prev = 0
    coords_target = SkyCoord(ra, dec, unit="deg")
    print "  ", pathFileEvt, "(", nEVT, "events )"
    for iEVT in range(1, nEVT):
        if aTIME[iEVT]<time_prev:
            print 'Odd order!!!'
            return 1
        if aTIME[iEVT]>=metStart and aTIME[iEVT]<metStop:
            time_fill = sqrt((aTIME[iEVT-1]-torigin)*(aTIME[iEVT]-torigin))
            if not rlim>0:
                nstat+=1
            else:
                coords_event = SkyCoord(aRA[iEVT], aDEC[iEVT], unit="deg")
                sep_ang = coords_target.separation(coords_event)
                sep_deg = float(sep_ang.to_string(unit=u.deg, decimal=True))
                if sep_deg<=rlim:
                    nstat+=1
                else:
                    continue
            if (nstat-noffset) >= nthreshold and int(time_fill)>=int(periods_goodstat[-1][0]): #and time_fill-periods_goodstat[-1][0]>=1:
                periods_goodstat[-1].append(time_fill)
                periods_goodstat.append([time_fill])
                nstat = 0
        time_prev = aTIME[iEVT]
    if len(periods_goodstat)>1:
        periods_goodstat = periods_goodstat[:-1]
        periods_goodstat[-1][1] = tStop
    else:
        periods_goodstat[-1].append(tStop)
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
