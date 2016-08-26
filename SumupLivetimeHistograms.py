#!/usr/bin/env python

import ROOT
from ROOT import TH1, TH2, TH3
import sys
import click
import subprocess
from pLsList import ls_list
ROOT.gROOT.SetBatch()


@click.command()
@click.argument('pathfiles', type=str)
@click.argument('pathfileout', type=str)
#@click.option('--remove', '-r', is_flag=True, help="Remove the current friends")
def main(pathfiles, pathfileout):
    li_path_files = ls_list(pathfiles)
    li_name_htg = []
    li_htg = []
    file_out = ROOT.TFile(pathfileout, "UPDATE")
    for (ifile, path_file) in enumerate(li_path_files):
        file_htg = ROOT.TFile(path_file)
        print file_htg.GetName(), 'is found.'
        for (ihtg, name_htg) in enumerate(li_name_htg):
            htg = file_htg.Get(name_htg)
            if ifile==0:
                li_htg.append(htg.Clone("{0}_summed".format(htg.GetName())))
            else:
                li_htg[ihtg].Add(htg)
    for htg in li_htg:
        htg.Write()


if __name__ == '__main__':
    main()
