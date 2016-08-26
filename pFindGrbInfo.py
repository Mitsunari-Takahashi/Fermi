#!/usr/bin/env python

import sys
import xml.etree.ElementTree as ET

def find_grb_info(nameGrb, rtXml):
    for grb in rtXml:
        if grb.findtext("./GRBNAME")==nameGrb: 
            trigger_time = float(grb.findtext("./MET"))
            if grb.findtext("./ERROR") == "--" or grb.findtext("./ERROR") == "":
                if grb.findtext("./LATERROR") == "--" or grb.findtext("./LATERROR") == "":
                    err_rad = 0.
                else:
                    err_rad = float(grb.findtext("./LATERROR"))
            else:
                if grb.findtext("./LATERROR") == "--" or grb.findtext("./LATERROR") == "" or float(grb.findtext("./ERROR"))<=float(grb.findtext("./LATERROR")):
                    err_rad = float(grb.findtext("./ERROR"))                    
                    raSrc = float(grb.findtext("./RA"))
                    decSrc = float(grb.findtext("./DEC"))
                else :
                    err_rad = float(grb.findtext("./LATERROR"))                    
                    raSrc = float(grb.findtext("./LATRA"))
                    decSrc = float(grb.findtext("./LATDEC"))
    dict_src = {"TRIGGER_MET": trigger_time, "ERROR_RADIUS": err_rad, "RA": raSrc, "DEC": decSrc}
    print ""
    print "==============="
    print "GRB", nameGrb
    print "==============="
    print "(", dict_src["RA"], ",", dict_src["DEC"], "), Error radius:", dict_src["ERROR_RADIUS"], "Trigger MET:", dict_src["TRIGGER_MET"]
    return dict_src
