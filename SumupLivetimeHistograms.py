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
#    LI_NAME_HTG = ["htgLt_GRB160509374_0_yx", "htgLt_GRB160509374_0_scaled"]
    LI_NAME_HTG = ["htgLt_GalacticOFF_yx"]
    li_htg = []
    file_out = ROOT.TFile(pathfileout, "UPDATE")
    for path_file in li_path_files:
        file_htg = ROOT.TFile(path_file)
        print file_htg.GetName(), 'is found.'
        for (ihtg, name_htg) in enumerate(LI_NAME_HTG):
            htg = file_htg.Get(name_htg)
            if path_file==li_path_files[0]:
                file_out.cd()
                li_htg.append(htg.Clone("{0}_summed".format(htg.GetName())))
                print li_htg[ihtg].GetName(), "is cloned."
                li_htg[-1].Write()
            else:
                li_htg[ihtg].Add(htg)
    for htg in li_htg:
        file_out.cd()
        htg.Write()
    print "Finished."

if __name__ == '__main__':
    main()
