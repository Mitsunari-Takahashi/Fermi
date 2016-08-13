#!/usr/bin/env python
"""
For converting sky coordinates
python -m pConvertSkyCoords <RA in deg> <DEC in deg>
return [L, B]
"""
import sys
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
import ROOT
ROOT.gROOT.SetBatch()


def convert_radec_to_galactic(radec):
    """
Convert a sky position from RADEC to Galactic
pConvertSkyCoords.convert_radec_to_galactic([ra in deg, dec in deg])
return [l in deg, b in deg]
"""
    ra = radec[0]
    dec = radec[1]
    c = SkyCoord(ra, dec, unit="deg")
    l = c.galactic.l.deg
    b = c.galactic.b.deg
    return [l, b]


def convert_galactic_to_radec(gal):
    """
Convert a sky position from Galactic to RADEC
pConvertSkyCoords.convert_radec_to_galactic([l in deg, b in deg])
return [ra in deg, dec in deg]
"""
    l = gal[0]
    b = gal[1]
    c = SkyCoord(l, b, Galactic, unit="deg")
    ra = c.ra.deg
    dec = c.dec.deg
    return [ra, dec]


if __name__ == '__main__':
    if len(sys.argv)==3:
        result = convert_radec_to_galactic(sys.argv[1:])
        print result
    else:
        help(convert_radec_to_galactic)
