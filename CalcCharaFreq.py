#!/usr/bin/env python

import sys
import os
import numpy as np
import math
from math import sqrt, cos, sin, tan, acos, asin, atan, radians, degrees, pi
from astropy import units as u
from astropy.units import cds
import astropy.constants as constants
import click
from logging import getLogger,StreamHandler,DEBUG,INFO,WARNING,ERROR,CRITICAL

##### Logger #####
logger = getLogger(__name__)
handler = StreamHandler()
loglevel = 'INFO'
handler.setLevel(loglevel)
logger.setLevel(loglevel)
logger.addHandler(handler)


class ShockedFluid:
    def __init__(self, e52, n0, p, epsilonB, epsilonE):
        self.energy = e52 * 1E52 * u.erg
        self.ndensity = n0 * u.cm**-3
        self.p = p
        self.epsilonB = epsilonB
        self.epsilonE = epsilonE


    def nu_c(self, t_day):
        return ( 3.*sqrt(3.) * u.M_e.cgs * constants.e.esu * (cds.c.cgs)**0.5 ) / ( 4. * constants.sigma_T.cgs**2 * (self.energy)**0.5 * self.ndensity * u.M_p.cgs * self.epsilonB**(3./2.) * ((t_day*u.d).to(u.s))**0.5 )

    
    def nu_m(self, t_day):
        return ( 3. * self.energy * self.epsilonB / c**5 / ((t_day*u.d).to(u.s))**3. )**0.5 * self.epsilonE**2 * ((self.p-2.)/(self.p-1.))**2 * (u.M_p.cgs/u.M_e.cgs)**2 * constants.e.esu / u.M_e.cgs / 8. / pi
        

@click.command()
@click.argument('tday', type=float)
@click.option('--e52', type=float, default=1)
@click.option('--n0', type=float, default=1)
@click.option('--pelec', '-p', type=float, default=2.5)
@click.option('--epsilonb', '-b', type=float, default=1E-4)
@click.option('--epsilone', '-e', type=float, default=0.5)
#@click.option('--values', type=(str, int))
#@click.option('--values', multiple=True)
#@click.option('--language', type=click.Choice(['Japanese', 'English']))
#@click.option('--shout', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(tday, e52, n0, pelec, epsilonb, epsilone, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)
    shocked_fluid = ShockedFluid(e52, n0, pelec, epsilonb, epsilone)
    logger.info('v_c = {nuc} at T = {td} days'.format(nuc=shocked_fluid.nu_c(tday).to(u.Hz), td=tday))
    logger.info('v_m = {num} at T = {td} days'.format(num=shocked_fluid.nu_m(tday).to(u.Hz), td=tday))

if __name__ == '__main__':
    main()
