import ROOT
from ROOT import kWhite, kBlack, kGray, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kSpring, kTeal, kAzure, kViolet, kPink

def stepColorRGB(nStep, iStep):
    aColorRGB = []
    for pS in range(nStep):
        aColorRGB.append(ROOT.TColor(3001+pS, 255*(-2.*pS/(nStep-1)+1)*int(2.*pS<nStep-1), 255*(2.*pS/(nStep-1)*int(2.*pS<nStep-1)+(-2.*pS/(nStep-1)+2)*int(2.*pS>=nStep-1)), 255*(2.*pS/(nStep-1)-1)*int(2.*pS>nStep-1)))
    return aColorRGB[iStep].GetNumber()

def akColor(iColor):
    if iColor>=0 and iColor<11:
        akColor = [kBlack, kRed, kBlue, kGreen+3, kMagenta, kCyan+3, kYellow+2, kOrange+2, 28, 39, 46]
        return akColor[iColor]
    else:
        print "No color No.%s is defined." % iColor
        return kWhite

def aakMarkerStyle(iMarker, jMarker):
    if iMarker==0:
        if jMarker==0:
            return 25;
        elif jMarker==1:
            return 5;
        else:
            print "Bad second argument!"
    elif iMarker==1:
        if jMarker==0:
            return 25;
        elif jMarker==1:
            return 20;
        elif jMarker==2:
            return 5;
        else:
            print "Bad second argument!"
    else:
        print "Bad first argument!"

def aakMarkerSize(iMarker, jMarker):
    if iMarker==0:
        if jMarker==0:
            return 1.2;
        elif jMarker==1:
            return 1.5;
        else:
            print "Bad second argument!"
    elif iMarker==1:
        if jMarker==0:
            return 1.2;
        elif jMarker==1:
            return 1.0;
        elif jMarker==2:
            return 1.5;
        else:
            print "Bad second argument!"
    else:
        print "Bad first argument!"
