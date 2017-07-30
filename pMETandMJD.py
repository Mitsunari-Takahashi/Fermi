#!/usr/bin/env python

import click

DateOffset = 51910.0
SecPerDay = 86400.

def ConvertMetToMjd(MET):
   #Correction for leap seconds
   nleapsec = 0
   if MET > 157766400: #2006-01-01 00:00:00.000
      nleapsec+=1
   if MET > 252460801.0: #2009-01-01 00:00:00.000
      nleapsec+=1
   if MET > 362793602.0: #2012-07-01 00:00:00.000
      nleapsec+=1
   if MET > 457401603.0: #2015-07-01 00:00:00.000
      nleapsec+=1
   if MET > 504921604.0: #2017-01-01 00:00:00
      nleapsec+=1
   MJD=(MET-nleapsec)/SecPerDay+DateOffset
   return MJD

def ConvertMjdToMet(MJD):
   #Correction for leap seconds
   nleapsec = 0
   if MJD > 53736.0: #2006-01-01 00:00:00.000
      nleapsec+=1
   if MJD > 54832.0: #2009-01-01 00:00:00.000
      nleapsec+=1
   if MJD > 56109.0: #2012-07-01 00:00:00.000
      nleapsec+=1
   if MJD > 57204.0: #2015-07-01 00:00:00.000
      nleapsec+=1
   if MJD > 57754.0: #2017-01-01 00:00:00.000
      nleapsec+=1
   MET = (MJD-DateOffset) * SecPerDay + nleapsec
   return MET

def ConvertMetToFMW(MET):
   DateOffset=51910.0
   WeekOffset=-387
   MJD=ConvertMetToMjd(MET)
   FMW=(MJD-DateOffset)/7.+WeekOffset
   return FMW


@click.command()
@click.argument('source', type=float)
@click.option('--mode', '-m', type=click.Choice(['MJDtoMET', 'METtoMJD']))
def main(source, mode):
   if mode=='MJDtoMET':
      result = ConvertMjdToMet(source)
   elif mode=='METtoMJD':
      result = ConvertMetToMjd(source)  
   print result


if __name__ == '__main__':
    main()
