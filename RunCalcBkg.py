#!/usr/bin/env python

import sys
import os.path 
import subprocess
par = sys.argv

grb = par[1]
cls0 = par[2] #R100, R030, R010, R003

if os.path.exists('Plot_GRB{0}_LATT0-10000sec.root'.format(grb))==False:
    subprocess.call(['python', '/home/takhsm/PythonModuleMine/Fermi/AnalyzeGRB-LAT.py', '-s', 'LATT0-10000sec', grb, "/disk/gamma/cta/store/takhsm/FermiMVA/AllSky/events/trAllSkyMap_S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15_catZDIR060_BDTG500D06_catZDIR060_20{0}_???.root".format(grb[:2]), '0', '10000'])
if os.path.exists('Livetime_GRB{0}_T0-10000_CalOnly{1}.root'.format(grb, cls0))==False:
    subprocess.call(['python', '/home/takhsm/PythonModuleMine/Fermi/pLivetimePoint.py', 'CalOnly{0}'.format(cls0), '0', '10000', '0', '0', cls0, grb])

subprocess.call(['python', '/home/takhsm/PythonModuleMine/Fermi/Exposure.py', '--limpsf','2', 'Livetime_GRB{0}_T0-10000_CalOnly{1}.root'.format(grb, cls0), '/disk/gamma/cta/store/takhsm/FermiMVA/MVA/S18/S18V200909_020RAWE20ZDIR020ZCS000wwoTRKwoMCZDIR00woRWcatTwo_15/S18ZDIR020catTwoZDIR060_E28bin_Cth40bins_axisObs_CalOnly_{0}_perf.root'.format(cls0), 'htgLt_0_scaled'])

subprocess.call(['python', '/home/takhsm/PythonModuleMine/Fermi/CalcExpectedBackgroundLocal.py', '--suffix','T0-10000sec_limpsf2deg', '--htgevton','htgEvt_CalOnly_{0}_PSF68'.format(cls0), '--htgevtoff','htgEvt_GalOffCalOnly_{0}'.format(cls0), '--htgexpon','htgExp_0_scaled_limpsf2deg', '--htgexpoff','htgExp_GalacticOFF_yx', '-m','GalacticOFF', '/disk/gamma/cta/store/takhsm/FermiMVA/GRB/{0}/Plot_GRB{0}_LATT0-10000sec.root'.format(grb), '/disk/gamma/cta/store/takhsm/FermiMVA/AllSky/OFF-HighB/Plot_GalOff_MET239557417-501033396.root', 'Exposure_GRB{0}_T0-10000_CalOnly{1}_limpsf2deg.root'.format(grb, cls0), '/disk/gamma/cta/store/takhsm/FermiMVA/AllSky/Livetime_GalacticOff/Exposure_GalacticOff__MET239557417-501033396_sum_CalOnly{0}.root'.format(cls0)])
