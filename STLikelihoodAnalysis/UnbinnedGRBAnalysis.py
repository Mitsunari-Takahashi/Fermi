#!/usr/bin/env python

import click
import pLATLikelihoodConfig
#import pyLikelihood
#from UnbinnedAnalysis import *
from STLikelihoodAnalysis import get_module_logger

##### Logger #####
logger = get_module_logger(__name__)


def unbinned_grb_analysis(name, mode, force, suffix, grbcatalogue):
    #self.obs = UnbinnedObs()#self.path_filtered_gti, self.path_ft2, expMap=self.path_exposure, expCube=self.path_livetime, irfs=self.irfs)
    #self.like = UnbinnedAnalysis() #self.obs, self.path_model_xml, optimizer='NewMinuit')
    grb = pLATLikelihoodConfig.GRBTarget(name, grbcatalogue)
    ana = pLATLikelihoodConfig.GRBConfig(grb, mode, emin=100., emax=5623.41, emin_eval=10000., emax_eval=100000.)
    ana.setup(force={'download':False, 'filter':False, 'maketime':False, 'livetime':False, 'exposure':False, 'model_3FGL_sources':False, 'diffuse_responses':False})
    # ana.set_directories()
    # ana.download()
    # ana.filter()
    # ana.maketime()
    # ana.livetime(force)
    # ana.exposure(force)
    # ana.model_3FGL_sources(force)
    # ana.diffuse_responses(force)
    ana.fit(bredo=False)
    ana.plot_countspectra_fitted()
    ana.eval_flux_and_error()


@click.command()
@click.argument('name', type=str)
#@click.argument('dst', nargs=-1)
@click.option('--suffix', type=str, default='')
@click.option('--grbcatalogue', type=str, default=pLATLikelihoodConfig.GRB_CATALOGUE_LTF)
#@click.option('--values', multiple=True)
@click.option('--mode', type=click.Choice(['unified', 'prompt', 'afterglow', 'earlyAG', 'lateAG', 'lightcurve', 'special']))
@click.option('--force', '-f', is_flag=True)
def main(name, mode, force, suffix, grbcatalogue):
    unbinned_grb_analysis(name, mode, force, suffix, grbcatalogue)


if __name__ == '__main__':
    main()
