import ROOT
import numpy
import sys
import math
from math import cos, sin, tan, acos, asin, atan, radians, degrees

JulYrInSec = 31557600

def anglePointsRadToRad(phi0, theta0, phi1, theta1):
    vec0 = numpy.array([cos(theta0)*cos(phi0), cos(theta0)*sin(phi0), sin(theta0)])
    vec1 = numpy.array([cos(theta1)*cos(phi1), cos(theta1)*sin(phi1), sin(theta1)])
    # if numpy.all(vec0-vec1==0):
    #     return 0
    # elif numpy.all(vec0+vec1==0):
    #     return math.pi
    # else:
    dotp = numpy.dot(vec0, vec1)
    if dotp>=1:
        return 0
    elif dotp<=-1:
        return math.pi
    else:
        return acos(dotp)
#    else:
#        print "Vector0", vec0
#        print "Vector1", vec1
#        print "dotp:", dotp
#        sys.exit(1)

def anglePointsDegToDeg(phi0, theta0, phi1, theta1):
    return degrees(anglePointsRadToRad(radians(phi0), radians(theta0), radians(phi1), radians(theta1)))
