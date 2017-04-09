#!/usr/bin/env python

import ROOT
from ROOT import TGraph2D
import click
from math import log10
ROOT.gROOT.SetBatch()


class EBLmodel:
    def __init__(self, name, filedat):
        self.name = name
        if filedat[-5:] == '.root':
            fileroot = ROOT.TFile(filedat, 'READ')
            self.gr2_tau = fileroot.Get(name)
            print self.gr2_tau.GetName(), 'is found.'
            self.gr2_tau.SetDirectory(0)
        else:
            self.gr2_tau = ROOT.TGraph2D()
            self.gr2_tau.SetName(self.name)
            self.gr2_tau.SetTitle('Optical depth')
            self.gr2_tau.GetXaxis().SetTitle('log_{10}Energy[GeV]')
            self.gr2_tau.GetYaxis().SetTitle('redshift')
            npoint = 0
            with open(filedat, 'r') as model:
                for (ienergy, energy) in enumerate(model):
                    if ienergy==0:
                        lst_str_redshift = energy[16:].split(",")
                        self.nredshift = len(lst_str_redshift)
                        self.redshift = []
                        for zstr in lst_str_redshift:
                            self.redshift.append(float(zstr))
                    else:
                        lst_str_energy = energy.split()
                        e_gam = float(lst_str_energy[0])
                        for iz in range(self.nredshift):
                            if self.redshift[iz]<2.0 and e_gam<1000:
                                tau = float(lst_str_energy[iz+1])
                                print npoint, log10(e_gam), self.redshift[iz], tau
                                self.gr2_tau.SetPoint(npoint, log10(e_gam), self.redshift[iz], tau)
                                npoint+=1
            self.gr2_tau.SaveAs('EBL_{0}.root'.format(self.name))


    def interpolate(self, loge, z):
        """Provide log(Energy) in GeV and z(redshift)
        """
        return self.gr2_tau.Interpolate(loge, z)
                    
                
        
