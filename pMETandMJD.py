def ConvertMetToMjd(MET):
   DateOffset=51910.0
   MJD=MET/3600./24.+DateOffset
   return MJD

def ConvertMjdToMet(MJD):
   DateOffset=51910.0
   MET=(MJD-DateOffset)*24.*3600.
   return MET

def ConvertMetToFMW(MET):
   DateOffset=51910.0
   WeekOffset=-387
   MJD=ConvertMetToMjd(MET)
   FMW=(MJD-DateOffset)/7.+WeekOffset
   return FMW
