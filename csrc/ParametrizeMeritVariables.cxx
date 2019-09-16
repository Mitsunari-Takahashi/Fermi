#include <stdio.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <sstream>
#include <stdint.h>
#include <getopt.h>
#include "TMath.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2F.h"
#include "TF2.h"
#include "TTree.h"
#include "Variable.h"

using namespace std;


int main(int argc,char *argv[])
{
  string strPathFileIn = (string)argv[1];
  TFile *fileIn = new TFile(strPathFileIn.c_str(), "READ");
  TTree *trDat =  (TTree*)fileIn->Get("MeritTuple");
  const Int_t nEvt = trDat->GetEntries();
  cout << trDat->GetName() << " has " << nEvt <<  " events." << endl;

  string strPathFileFit = argv[2];
  TFile* fileFit = new TFile(strPathFileFit.c_str(), "READ");
  
  Double_t energy;
  Double_t zdir;
  trDat->SetBranchAddress("WP8CalOnlyEnergy", &energy);
  trDat->SetBranchAddress("Cal1MomZDir", &zdir);

  /* VariableCalELayerCorrInitialLog */
  VariableCalELayerCorrInitialLog vCalELayerCorrInitialLog = VariableCalELayerCorrInitialLog("CalELayerCorrInitialLog", "ACD", fileFit, 1E-5);
  trDat->SetBranchAddress("Cal1MomParticleInLayerCor", &vCalELayerCorrInitialLog.cal1MomParticleInLayerCor);
  for(int j=0; j<8; j++)
    {
      trDat->SetBranchAddress(Form("CalELayer%d", j), &vCalELayerCorrInitialLog.calELayer[j]);
    }

  /* VariableFloat */
  VariableFloat vCalBkHalfRatio = VariableFloat("CalBkHalfRatio", "CAL", fileFit);
  trDat->SetBranchAddress("CalBkHalfRatio", &vCalBkHalfRatio.original);

  VariableFloat vCalELayerCorrInitialRatioLog = VariableFloat("CalELayerCorrInitialRatioLog", "ACD", fileFit);
  trDat->SetBranchAddress("CalELayerCorrInitialRatioLog", &vCalELayerCorrInitialRatioLog.original);

  VariableFloat vCalELayer34afterInitialRatioLog = VariableFloat("CalELayer34afterInitialRatioLog", "CAL", fileFit);
  trDat->SetBranchAddress("CalELayer34afterInitialRatioLog", &vCalELayer34afterInitialRatioLog.original);

  VariableFloat vCalTrSizeCalT95 = VariableFloat("CalTrSizeCalT95", "CAL", fileFit);
  trDat->SetBranchAddress("CalTrSizeCalT95", &vCalTrSizeCalT95.original);

  VariableDouble vCal1MomNumCoreXtalsFract = VariableDouble("Cal1MomNumCoreXtalsFract", "CAL", fileFit);
  trDat->SetBranchAddress("Cal1MomNumCoreXtalsFract", &vCal1MomNumCoreXtalsFract.original);

  VariableFloat vAcd2TileEnergyRatioLog = VariableFloat("Acd2TileEnergyRatioLog", "ACD", fileFit);
  trDat->SetBranchAddress("Acd2TileEnergyRatioLog", &vAcd2TileEnergyRatioLog.original);

  VariableFloat vCalNewCfpCalTmax = VariableFloat("CalNewCfpCalTmax", "CAL", fileFit);
  trDat->SetBranchAddress("CalNewCfpCalTmax", &vCalNewCfpCalTmax.original);

  //  VariableFloat vCalELayer0 = VariableFloat("CalELayer0", "ACD", fileFit);
  //&(vCalELayer0.original) = &(vCalELayerCorrInitialLog.calELayer[0]);
  //trDat->SetBranchAddress("CalELayer0", &vCalELayer0.original);

  /* VariableFloatLog10, 1E-5 */
  VariableFloatLog10 vCalEdgeEnergyLog = VariableFloatLog10("CalEdgeEnergyLog", "CAL", fileFit, 1E-5);
  trDat->SetBranchAddress("CalEdgeEnergy", &vCalEdgeEnergyLog.original);

  VariableFloatLog10 vlog10CalELayer0 = VariableFloatLog10("log10CalELayer0", "ACD", fileFit, 1E-5);
  //  vlog10CalELayer0.SetOriginal(vCalELayerCorrInitialLog.calELayer[0]);
  //  &(vlog10CalELayer0.original) = &(vCalELayerCorrInitialLog.calELayer[0]);
  //  trDat->SetBranchAddress("CalELayer0", &vlog10CalELayer0.original);

  VariableFloatLog10 vAcd2TileEnergyLog = VariableFloatLog10("Acd2TileEnergyLog", "ACD", fileFit, 1E-5);
  trDat->SetBranchAddress("Acd2TileEnergy", &vAcd2TileEnergyLog.original);

  VariableFloatLog10 vlog10CalNewCfpCalSelChiSq = VariableFloatLog10("log10CalNewCfpCalSelChiSq", "CAL", fileFit, 1E-5);
  trDat->SetBranchAddress("CalNewCfpCalSelChiSq", &vlog10CalNewCfpCalSelChiSq.original);

  /* VariableFloatLog10, 1E-10 */
  VariableFloatLog10 vlog10Cal1TransRms = VariableFloatLog10("log10Cal1TransRms", "CAL", fileFit, 1E-5);
  trDat->SetBranchAddress("Cal1TransRms", &vlog10Cal1TransRms.original);

  /* VariableFloatLog10, 1E-20 */
  VariableFloatLog10 vlog10Cal1FitChiSquare = VariableFloatLog10("log10Cal1FitChiSquare", "CAL", fileFit, 1E-20);
  trDat->SetBranchAddress("Cal1FitChiSquare", &vlog10Cal1FitChiSquare.original);

  /* VariableFloatLog10, 1E-6 */
  VariableFloatLog10 vAcd2Cal1Energy15Log = VariableFloatLog10("Acd2Cal1Energy15Log", "ACD", fileFit, 1E-6);
  trDat->SetBranchAddress("Acd2Cal1Energy15", &vAcd2Cal1Energy15Log.original);

  /* VariableUIntLog10, 3E-1 */
  VariableUIntLog10 vAcd2VetoCountLog = VariableUIntLog10("Acd2VetoCountLog", "ACD", fileFit, 3E-1);
  trDat->SetBranchAddress("Acd2VetoCount", &vAcd2VetoCountLog.original);

  /* VariableFloatLog10, 1E-3 */
  VariableFloatLog10 vAcd2Cal1VetoSigmaHitLog = VariableFloatLog10("Acd2Cal1VetoSigmaHitLog", "ACD", fileFit, 1E-3);
  trDat->SetBranchAddress("Acd2Cal1VetoSigmaHit", &vAcd2Cal1VetoSigmaHitLog.original);

  /* CalELayer74RatioLog */
  VariableFloatRatioLog10 vCalELayer74RatioLog = VariableFloatRatioLog10("CalELayer74RatioLog", "CAL", fileFit, 1E-5);
  //vCalELayer74RatioLog.SerOriginal(vCalELayerCorrInitialLog.calELayer[7]);
  //vCalELayer74RatioLog.SerOriginal2(vCalELayerCorrInitialLog.calELayer[4]);
  //*(vCalELayer74RatioLog.original) = &(vCalELayerCorrInitialLog.calELayer[7]);
  //*(vCalELayer74RatioLog.original2) = &(vCalELayerCorrInitialLog.calELayer[4]);
  //trDat->SetBranchAddress("CalELayer7", &vCalELayer74RatioLog.original);
  //trDat->SetBranchAddress("CalELayer4", &vCalELayer74RatioLog.original2);

  string strPathFileParam = strPathFileIn.replace((int)strPathFileIn.length()-5, 5, "_parametrized.root");
  TFile* fileParam = new TFile(strPathFileParam.c_str(), "UPDATE");

  TTree* trParam = new TTree("MeritTuple", "Parametrized variables");

  trParam->Branch("CalBkHalfRatio_param", &(vCalBkHalfRatio.parametrized), "CalBkHalfRatio_param/F");
  trParam->Branch("CalELayerCorrInitialRatioLog_param", &(vCalELayerCorrInitialRatioLog.parametrized), "CalELayerCorrInitialRatioLog_param/F");
  trParam->Branch("CalELayer34afterInitialRatioLog_param", &(vCalELayer34afterInitialRatioLog.parametrized), "CalELayer34afterInitialRatioLog_param/F");
  trParam->Branch("CalTrSizeCalT95_param", &(vCalTrSizeCalT95.parametrized), "CalTrSizeCalT95_param/F");
  trParam->Branch("Cal1MomNumCoreXtalsFract_param", &(vCal1MomNumCoreXtalsFract.parametrized), "Cal1MomNumCoreXtalsFract_param/F");
  trParam->Branch("Acd2TileEnergyRatioLog_param", &(vAcd2TileEnergyRatioLog.parametrized), "Acd2TileEnergyRatioLog_param/F");
  trParam->Branch("CalNewCfpCalTmax_param", &(vCalNewCfpCalTmax.parametrized), "CalNewCfpCalTmax_param/F");
  //  trParam->Branch("CalELayer0_param", &(vCalELayer0.parametrized), "CalELayer0_param/F");

  trParam->Branch("CalEdgeEnergyLog_param", &(vCalEdgeEnergyLog.parametrized), "CalEdgeEnergyLog_param/F");
  trParam->Branch("log10CalELayer0_param", &(vlog10CalELayer0.parametrized), "log10CalELayer0_param/F");
  trParam->Branch("Acd2TileEnergyLog_param", &(vAcd2TileEnergyLog.parametrized), "Acd2TileEnergyLog_param/F");
  trParam->Branch("Acd2VetoCountLog_param", &(vAcd2VetoCountLog.parametrized), "Acd2VetoCountLog_param/F");
  trParam->Branch("Acd2Cal1VetoSigmaHitLog_param", &(vAcd2Cal1VetoSigmaHitLog.parametrized), "Acd2Cal1VetoSigmaHitLog_param/F");
  trParam->Branch("log10CalNewCfpCalSelChiSq_param", &(vlog10CalNewCfpCalSelChiSq.parametrized), "log10CalNewCfpCalSelChiSq_param/F");

  trParam->Branch("log10Cal1TransRms_param", &(vlog10Cal1TransRms.parametrized), "log10Cal1TransRms_param/F");
  trParam->Branch("log10Cal1FitChiSquare_param", &(vlog10Cal1FitChiSquare.parametrized), "log10Cal1FitChiSquare_param/F");
  trParam->Branch("CalELayer74RatioLog_param", &(vCalELayer74RatioLog.parametrized), "CalELayer74RatioLog_param/F");
  trParam->Branch("Acd2Cal1Energy15Log_param", &(vAcd2Cal1Energy15Log.parametrized), "Acd2Cal1Energy15Log_param/F");
  trParam->Branch("CalELayerCorrInitialLog_param", &(vCalELayerCorrInitialLog.parametrized), "CalELayerCorrInitialLog_param/F");
  Double_t loge;
  for(int iEvt=0; iEvt<nEvt; iEvt++)
    {
      trDat->GetEntry(iEvt);
      loge = TMath::Log10(energy);

      vCalBkHalfRatio.CalcValue();
      vCalBkHalfRatio.Parametrize(loge, zdir);

      vCalELayerCorrInitialRatioLog.CalcValue();
      vCalELayerCorrInitialRatioLog.Parametrize(loge, zdir);

      vCalELayer34afterInitialRatioLog.CalcValue();			      
      vCalELayer34afterInitialRatioLog.Parametrize(loge, zdir);

      vCalTrSizeCalT95.CalcValue();			      
      vCalTrSizeCalT95.Parametrize(loge, zdir);

      vCal1MomNumCoreXtalsFract.CalcValue();			      
      vCal1MomNumCoreXtalsFract.Parametrize(loge, zdir);

      vAcd2TileEnergyRatioLog.CalcValue();			      
      vAcd2TileEnergyRatioLog.Parametrize(loge, zdir);

      vCalNewCfpCalTmax.CalcValue();			      
      vCalNewCfpCalTmax.Parametrize(loge, zdir);

      // vCalELayer0.original = vCalELayerCorrInitialLog.calELayer[0];
      // vCalELayer0.CalcValue();			      
      // vCalELayer0.Parametrize(loge, zdir);

      vCalEdgeEnergyLog.CalcValue();			      
      vCalEdgeEnergyLog.Parametrize(loge, zdir);

      vlog10CalELayer0.original = vCalELayerCorrInitialLog.calELayer[0];
      vlog10CalELayer0.CalcValue();			      
      vlog10CalELayer0.Parametrize(loge, zdir);

      vAcd2TileEnergyLog.CalcValue();		      
      vAcd2TileEnergyLog.Parametrize(loge, zdir);

      vAcd2VetoCountLog.CalcValue();		      
      vAcd2VetoCountLog.Parametrize(loge, zdir);

      vAcd2Cal1VetoSigmaHitLog.CalcValue();		      
      vAcd2Cal1VetoSigmaHitLog.Parametrize(loge, zdir);

      vlog10CalNewCfpCalSelChiSq.CalcValue();		      
      vlog10CalNewCfpCalSelChiSq.Parametrize(loge, zdir);

      vlog10Cal1TransRms.CalcValue();		      
      vlog10Cal1TransRms.Parametrize(loge, zdir);

      vlog10Cal1FitChiSquare.CalcValue();		      
      vlog10Cal1FitChiSquare.Parametrize(loge, zdir);

      vCalELayer74RatioLog.original = vCalELayerCorrInitialLog.calELayer[7];
      vCalELayer74RatioLog.original2 = vCalELayerCorrInitialLog.calELayer[4];
      vCalELayer74RatioLog.CalcValue();		      
      vCalELayer74RatioLog.Parametrize(loge, zdir);

      vAcd2Cal1Energy15Log.CalcValue();		      
      vAcd2Cal1Energy15Log.Parametrize(loge, zdir);

      vCalELayerCorrInitialLog.CalcValue();		      
      vCalELayerCorrInitialLog.Parametrize(loge, zdir);

      if(iEvt%100000==0)
	{
	  cout << "." << flush;
	  if(iEvt%1000000==0)
	    {
	      cout << Form("\n%d%% have been filled.\n", (trParam->GetEntries()*100)/nEvt) << flush;
	      cout << "Acd2VetoCount: " << vAcd2VetoCountLog.original << endl;
	      vAcd2VetoCountLog.Print();
	    }
	}
	//cout << ".";
      trParam->Fill();
      //
    }
  cout << endl;
  trParam->Write();
  delete trParam;
  fileParam->Close();
  fileFit->Close();
  fileIn->Close();
}
