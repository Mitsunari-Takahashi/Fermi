#!/usr/bin/env python
"""To be included in scripts for checking MVA variables.
"""
DCT_ALIAS = {'log10CalNewCfpCalSelChiSq':'1*log10(CalNewCfpCalSelChiSq)',
             'log10Cal1TransRms' : '1*log10( Cal1TransRms )',
             'CalELayer74RatioLog':'log10(max(-5, CalELayer7))-log10(max(-5, CalELayer4))',
             'Acd2Cal1Energy15Log':'log10(max(Acd2Cal1Energy15,1E-6))',
             'Acd2VetoCountLog':'log10(max(Acd2VetoCount,3E-1))',
             'Acd2Cal1VetoSigmaHitLog':'log10(max(Acd2Cal1VetoSigmaHit,1E-3))',
             'log10Cal1FitChiSquare':'log10(max(-5, Cal1FitChiSquare))',
             'Cal1MomNumCoreXtalsFract':'Cal1MomNumCoreXtals/Cal1NumXtals',
             'CalEdgeEnergyLog':'log10(max(-5, CalEdgeEnergy))'
             }

#TPL_VAR = ('log10CalNewCfpCalSelChiSq', 'log10Cal1TransRms', 'CalNewCfpCalTmax', 'CalBkHalfRatio', 'CalELayer74RatioLog', 'Acd2Cal1Energy15Log', 'Acd2VetoCountLog', 'Acd2Cal1VetoSigmaHitLog', 'CalELayerCorrInitialRatioLog', 'CalELayer34afterInitialRatioLog', 'CalTrSizeCalT95', 'log10Cal1FitChiSquare', 'Cal1MomNumCoreXtalsFract', 'Acd2TileEnergyRatioLog', 'CalEdgeEnergyLog')
TPL_VAR = ('Acd2TileEnergyRatioLog', 'CalEdgeEnergyLog')
