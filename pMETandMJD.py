def ConvertMetToMjd(MET):
   DateOffset=51910.0
   MJD=MET/3600./24.+DateOffset
   return MJD

def ConvertMjdToMet(MJD):
   DateOffset=51910.0
   MET=(MJD-DateOffset)*24.*3600.
   return MET
