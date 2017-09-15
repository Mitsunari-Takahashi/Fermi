#!/usr/bin/env python

import sys
import click
import pLATLikelihoodConfig
from STLikelihoodAnalysis import get_module_logger

##### Logger #####
logger = get_module_logger(__name__)


def unbinned_grb_analysis(name, mode, emin, emax, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, masifps):
    if spectraltype=='PowerLaw':
        spectralpars = {'Prefactor':1e-10, 'Index':-2.0, 'Scale':1000.}
    elif spectraltype=='ExpCutoff':
        spectralpars = {'Prefactor':1e-10, 'Index':-2.0, 'Scale':1000., 'Ebreak':10.0, 'P1':10000., 'P2':0, 'P3':0}
    else:
        logger.critical("""{0} is NOT available!!! Use PowerLaw or ExpCutoff.""".format(spectraltype))
        sys.exit(1)
    grb = pLATLikelihoodConfig.GRBTarget(name, grbcatalogue, spectraltype=spectraltype, spectralpars=spectralpars)
    ana = pLATLikelihoodConfig.GRBConfig(grb, mode, emin=emin, emax=emax, deg_roi=roi, psForce=masifps)
    ana.setup(force={'download':False, 'filter':force, 'maketime':force, 'livetime':force, 'exposure':force, 'model_3FGL_sources':True, 'diffuse_responses':force})
    if modelonly==True:
        esl = ana.set_likelihood()
        sys.exit(esl)
    ana.fit(bredo=True)
    ana.plot_countspectra_fitted()
    ana.eval_flux_and_error()


@click.command()
@click.argument('name', type=str)
@click.option('--grbcatalogue', '-c', type=str, default=pLATLikelihoodConfig.GRB_CATALOGUE_LTF)
@click.option('--mode', '-m', type=click.Choice(['unified', 'prompt', 'afterglow', 'earlyAG', 'lateAG', 'lightcurve', 'special']))
@click.option('--emin', type=float, default=100.)
@click.option('--emax', type=float, default=100000.)
@click.option('--roi', type=float, default=12.)
@click.option('--spectraltype', type=click.Choice(['PowerLaw', 'ExpCutoff']))
@click.option('--suffix', '-s', type=str, default='')
@click.option('--force', '-f', is_flag=True)
@click.option('--modelonly', is_flag=True)
@click.option('--masifps', is_flag=True)
@click.option('--refit', '-r', is_flag=True)
@click.option('--outdir', '-o', type=str, default='')
def main(name, mode, emin, emax, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, masifps):
    unbinned_grb_analysis(name, mode, emin, emax, roi, spectraltype, refit, force, suffix, grbcatalogue, modelonly, outdir, masifps)


if __name__ == '__main__':
    main()
