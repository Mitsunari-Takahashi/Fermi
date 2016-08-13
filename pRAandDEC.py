def ConvertToDeg(RA_hms, DEC_dms):
   DateOffset=51910.0
   RA_deg=RA_hms[0]*15+RA_hms[1]*(15/60.)+RA_hms[2]*(15/3600.)
   if DEC_dms[0]>=0:
      DEC_deg=DEC_dms[0]+DEC_dms[1]/60.+DEC_dms[2]/3600.
   else:
      DEC_deg=DEC_dms[0]-DEC_dms[1]/60.-DEC_dms[2]/3600.
   return [RA_deg, DEC_deg]

#def ConvertToHMS(RA_deg, DEC_deg):
 #  RA_hms=
  # DEC_dms=
   #return [RA_hms, DEC_dms]
