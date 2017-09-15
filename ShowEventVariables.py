#!/usr/bin/env python
"""For showing MVA variables of one event.
"""
import sys
import click
import ROOT
from ROOT import gROOT, gDirectory, gPad, gSystem, gStyle, kTRUE, kFALSE, TTree
ROOT.gROOT.SetBatch()
from pVariablesForMVA import DCT_ALIAS, TPL_VAR


def ShowEventVariable(pathfilein, nametr, nrunid, nevtid, comparable, pathfileout):
    fileIn = ROOT.TFile(pathfilein)
    trIn = fileIn.Get(nametr)
    print 'Shown:', TPL_VAR
    print trIn.GetName(), 'is found.'
    str_scan = ''
    for key, formula in DCT_ALIAS.items():
        trIn.SetAlias(key, formula)
    for var in TPL_VAR:
        str_scan = str_scan + var + ':'
    str_scan = str_scan[:-1]
    trIn.Scan(str_scan, 'EvtRun=={0} && EvtEventId=={1}'.format(nrunid, nevtid))

    if comparable!='':
        if pathfileout=='':
            pathfileout = comparable.replace('.root', '_run{0:d}_evt{1:d}.root'.format(nrunid, nevtid))
        file_out = ROOT.TFile(pathfileout, 'RECREATE')
        file_out.cd()
        for evt in trIn:
            if evt.EvtRun==nrunid and evt.EvtEventId==nevtid:
                evt_enr = log10(evt.WP8CalOnlyEnergy)
                evt_cth = evt.Cal1MomZDir
                break
        file_htg = ROOT.TFile(comparable, 'READ')
        for var in TPL_VAR:
            print '=====', var, '====='
            hs = ROOT.THStack('hs_{0}'.format(var), var)
            htg3D_gam = file_htg.Get('h3Gam_{0}'.format(var))
            htg1D_gam = htg3D_gam.ProjectionZ('{0}_projZ'.format(htg3D_gam.GetName()), htg3D_gam.GetXaxis().FindBin(evt_enr), htg3D_gam.GetXaxis().FindBin(evt_enr), htg3D_gam.GetYaxis().FindBin(evt_cth), htg3D_gam.GetYaxis().FindBin(evt_cth))
            htg1D_gam.SetTitle('Gammas')
            hs.Add(htg1D_gam)
            htg3D_had = file_htg.Get('h3Had_{0}'.format(var))
            htg1D_had = htg3D_had.ProjectionZ('{0}_projZ'.format(htg3D_had.GetName()), htg3D_had.GetXaxis().FindBin(evt_enr), htg3D_had.GetXaxis().FindBin(evt_enr), htg3D_had.GetYaxis().FindBin(evt_cth), htg3D_had.GetYaxis().FindBin(evt_cth))
            htg1D_had.SetTitle('Hadrons')
            hs.Add(htg1D_had)
            htg3D_lep = file_htg.Get('h3Lep_{0}'.format(var))
            htg1D_lep = htg3D_lep.ProjectionZ('{0}_projZ'.format(htg3D_lep.GetName()), htg3D_lep.GetXaxis().FindBin(evt_enr), htg3D_lep.GetXaxis().FindBin(evt_enr), htg3D_lep.GetYaxis().FindBin(evt_cth), htg3D_lep.GetYaxis().FindBin(evt_cth))
            htg1D_lep.SetTitle('Leptons')
            hs.Add(htg1D_lep)
            max_peak = max(htg1D_gam.GetMaximum(), htg1D_had.GetMaximum(), htg1D_lep.GetMaximum())
            
            htg_evt = htg1D_gam.Clone('htg_run{0:d}_evtid{1:d}'.format(nrunid, nevtid))
            htg_evt.SetTitle('Run {0}, Event {1}'.format(nrunid, nevtid))
            htg_evt.SetLineColor(kBlack)
            htg_evt.SetFillColor(kBlack)
            htg_evt.SetFillStyle(1001)
            trIn.Draw('{0}>>{1}'.format(var, htg_evt.GetName()), 'EvtRun=={0} && EvtEventId=={1}'.format(nrunid, nevtid), 'goff')
            htg_evt.SetMaximum(max_peak)
            hs.Add(htg_evt)
            file_out.cd()
            hs.Write()
            can = ROOT.TCanvas(var, var, 800, 500)
            can.cd()
            hs.Draw('nostack')
            can.BuildLegend()
            can.Write()
            pathplotout = pathfileout.replace('.root', '_{0}.png'.format(var))
            can.SaveAs(pathplotout)


@click.command()
@click.argument('pathfilein')
@click.option('--nametr', default='MeritTuple')
@click.argument('nrunid', type=int)
@click.argument('nevtid', type=int)
@click.option('--reference', '-r', default='')
@click.option('--fileout', '-o', default='')
def main(pathfilein, nametr, nrunid, nevtid, reference, fileout):
    ShowEventVariable(pathfilein, nametr, nrunid, nevtid, reference, fileout)


if __name__ == '__main__':
    main()
