import ROOT

def cutBDT(nameFileRoc, aCutEgb=[10, 3, 1], nameVar="S18V200909_020RAWE20ZDIR010ZCS000wwoTRKwoMCZDIR00woRWcatThree_15_2_BDTG500D06Log"):
    fileRoc = ROOT.TFile(nameFileRoc, 'READ')
    print fileRoc.GetName(), " was opened."
    print "Cutting at", aCutEgb, "xEGB level."
    h2Sig = fileRoc.Get("sig_acc")
    h2Bkg = fileRoc.Get("bkg_rate")
    h2Count = fileRoc.Get("bkg_counts_cum_hist")
    h2Egb = fileRoc.Get("egb_rate")

    aaValCutBDT = []
    nEnergyBin = h2Sig.ProjectionX().GetNbinsX()

    aStrForYaml = []
    for cc in aCutEgb:
        aStrForYaml.append("CalOnly_{0}xEGB: '".format(cc))

    for ie in range(nEnergyBin):
        ebLow = h2Sig.GetXaxis().GetBinLowEdge(ie+1)
        ebUp = h2Sig.GetXaxis().GetBinLowEdge(ie+2)
        print '=========', ebLow,' - ', ebUp, '========='
        aaValCutBDT.append([])
        abFound = []
        aNumGamCut = []
        aNumRemCut = []
        aNumCountCut = []
        for jc in range(len(aCutEgb)):
            abFound.append(False)
        h1Sig = h2Sig.ProjectionY("h1SigAcc", ie+1, ie+1)
        h1Bkg = h2Bkg.ProjectionY("h1BkgAcc", ie+1, ie+1)
        h1Count = h2Count.ProjectionY("h1BkgCount", ie+1, ie+1)
        h1Egb = h2Egb.ProjectionY("h1EgbRate", ie+1, ie+1)
        for ib in range(1, h1Sig.GetNbinsX()+1):
            nEgb = h1Egb.GetBinContent(ib)
            nSig = h1Sig.GetBinContent(ib)
            nRem = h1Bkg.GetBinContent(ib)
            nCount = h1Count.GetBinContent(ib)
            for ic in range(len(aCutEgb)):
                if (abFound[ic]==False and nRem<=aCutEgb[ic]*nEgb):
                    abFound[ic] = True
                    aaValCutBDT[ie].append(h1Sig.GetBinCenter(ib))
                    aNumGamCut.append(nSig)
                    aNumRemCut.append(nRem)
                    aNumCountCut.append(nCount)
                    print "--------- Background level:", aCutEgb[ic], "x EGB ---------"
                    print "Cut value:", aaValCutBDT[ie][-1]
                    print "Acceptance:", aNumGamCut[-1], "[m^2 sr]"
                    print "Background rate:", aNumRemCut[-1], "[MeV sr^-1 s^-1]"
                    print "Background count:", aNumCountCut[-1], "[events]"
                    if ie != 0:
                        aStrForYaml[ic] = aStrForYaml[ic] + " || "
                    aStrForYaml[ic] = aStrForYaml[ic] + "(log10(WP8CalOnlyEnergy)>={0}&&log10(WP8CalOnlyEnergy)<{1}&&{2}>={3})".format(ebLow, ebUp, nameVar, aaValCutBDT[ie][-1])
    print "For perfromance YAML file:"
    for jc in range(len(aCutEgb)):
        aStrForYaml[jc] = aStrForYaml[jc] + "'"
        print aStrForYaml[jc]
    return aaValCutBDT
