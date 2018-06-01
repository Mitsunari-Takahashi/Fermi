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
    def __init__(self, e52, n0, epsilonB):
        self.energy = e52 * 1E52 * u.erg
        self.ndensity = n0 * u.cm**-3
        self.epsilonB = epsilonB


    def nu_c(self, t_day):
        return ( 3.*sqrt(3.) * u.M_e.cgs * constants.e.esu * (cds.c.cgs)**0.5 ) / ( 4. * constants.sigma_T.cgs**2 * (self.energy)**0.5 * self.ndensity * u.M_p.cgs * self.epsilonB**(3./2.) * ((t_day*u.d).to(u.s))**0.5 )
        #return ( 3.*sqrt(3.) * constants.m_e.si * constants.e.si * (constants.c.si)**0.5 ) / ( 4. * constants.sigma_T.si**2 * (self.energy.si)**0.5 * self.ndensity.si * constants.m_p.si * self.epsilonB**(3./2.) * ((t_day*u.d).to(u.s))**0.5 )
        

@click.command()
@click.argument('tday', type=float)
@click.option('--e52', type=float, default=1)
@click.option('--n0', type=float, default=1)
@click.option('--epsilonb', '-b', type=float, default=1E-4)
#@click.option('--values', type=(str, int))
#@click.option('--values', multiple=True)
#@click.option('--language', type=click.Choice(['Japanese', 'English']))
#@click.option('--shout', is_flag=True)
@click.option('--loglevel', type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'CRITICAL']), default='INFO')
def main(tday, e52, n0, epsilonb, loglevel):
    ##### Logger #####
    handler.setLevel(loglevel)
    logger.setLevel(loglevel)
    logger.addHandler(handler)
    shocked_fluid = ShockedFluid(e52, n0, epsilonb)
    #print shocked_fluid.nu_c(tday)
    #print shocked_fluid.nu_c(tday).cgs
    #print shocked_fluid.nu_c(tday).to(u.Hz)
    logger.info('v_c = {nuc} at T = {td} days'.format(nuc=shocked_fluid.nu_c(tday).to(u.Hz), td=tday))

if __name__ == '__main__':
    main()
