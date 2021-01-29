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
using namespace TMath;

int main(int argc,char *argv[])
{
  string strPathFileIn = (string)argv[1];
  TFile *fileIn = new TFile(strPathFileIn.c_str(), "READ");
  TTree *trDat =  (TTree*)fileIn->Get("MeritTuple");
  const Int_t nEvt = trDat->GetEntries();
  cout << trDat->GetName() << " has " << nEvt <<  " events." << endl;

  //  string strPathFileFit = argv[2];
  //  TFile* fileFit = new TFile(strPathFileFit.c_str(), "READ");
  Float_t minLimit = 1E-5;
  Double_t maxZCross = 1000000;
  Float_t minDir = Power(10, -10.5);
  
  Double_t WP8CalOnlyEnergy;
  trDat->SetBranchAddress("WP8CalOnlyEnergy", &WP8CalOnlyEnergy);
  Double_t Cal1MomXDir;
  trDat->SetBranchAddress("Cal1MomXDir", &Cal1MomXDir);
  Double_t Cal1MomYDir;
  trDat->SetBranchAddress("Cal1MomYDir", &Cal1MomYDir);
  Double_t Cal1MomZDir;
  trDat->SetBranchAddress("Cal1MomZDir", &Cal1MomZDir);
  Float_t Cal1MomXCntr;
  trDat->SetBranchAddress("Cal1MomXCntr", &Cal1MomXCntr);
  Float_t Cal1MomYCntr;
  trDat->SetBranchAddress("Cal1MomYCntr", &Cal1MomYCntr);
  Float_t Cal1MomZCntr;
  trDat->SetBranchAddress("Cal1MomZCntr", &Cal1MomZCntr);
  Float_t Cal1MomXCntrCor;
  trDat->SetBranchAddress("Cal1MomXCntrCor", &Cal1MomXCntrCor);
  Float_t Cal1MomYCntrCor;
  trDat->SetBranchAddress("Cal1MomYCntrCor", &Cal1MomYCntrCor);
  Float_t Cal1MomZCntrCor;
  trDat->SetBranchAddress("Cal1MomZCntrCor", &Cal1MomZCntrCor);
  Float_t CalELayer[8];
  for(int i=0; i<8; i++)
    trDat->SetBranchAddress(Form("CalELayer%d", i), &CalELayer[i]);
  Float_t CalEnergyCorr;
  trDat->SetBranchAddress("CalEnergyCorr", &CalEnergyCorr);
  //  Float_t CalEnergyRaw;
  //  trDat->SetBranchAddress("CalEnergyRaw", &CalEnergyRaw);
  Float_t Acd2TileEnergy;
  trDat->SetBranchAddress("Acd2TileEnergy", &Acd2TileEnergy);
  Float_t EvtJointEnergy;
  trDat->SetBranchAddress("EvtJointEnergy", &EvtJointEnergy);

  Float_t calZTop = -47.395;
  Float_t cellVertPitch = 21.35;

  string strPathFileNew = strPathFileIn.replace((int)strPathFileIn.length()-5, 5, "_newVariables.root");
  TFile* fileNew = new TFile(strPathFileNew.c_str(), "UPDATE");

  TTree* trNew = new TTree("MeritTuple", "Additional glast tuple");

  Int_t Cal1MomParticleInLayer;
  trNew->Branch("Cal1MomParticleInLayer", &Cal1MomParticleInLayer, "Cal1MomParticleInLayer/I");
  Float_t Cal1MomZCrossSide840;
  trNew->Branch("Cal1MomZCrossSide840", &Cal1MomZCrossSide840, "Cal1MomZCrossSide840/F");
  Float_t Cal1MomZCrossSide714;
  trNew->Branch("Cal1MomZCrossSide714", &Cal1MomZCrossSide714, "Cal1MomZCrossSide714/F");
  // Float_t Cal1MomZCrossSide840Cor;
  // trNew->Branch("Cal1MomZCrossSide840Cor", &Cal1MomZCrossSide840Cor, "Cal1MomZCrossSide840Cor/F");
  // Float_t Cal1MomZCrossSide714Cor;
  // trNew->Branch("Cal1MomZCrossSide714Cor", &Cal1MomZCrossSide714Cor, "Cal1MomZCrossSide714Cor/F");
  Float_t CalELayerCorrInitialRatioLog;
  trNew->Branch("CalELayerCorrInitialRatioLog", &CalELayerCorrInitialRatioLog, "CalELayerCorrInitialRatioLog/F");
  Float_t CalELayer34afterInitialRatioLog;
  trNew->Branch("CalELayer34afterInitialRatioLog", &CalELayer34afterInitialRatioLog, "CalELayer34afterInitialRatioLog/F");
  //  Float_t CalELayerInitialRawRatioLog;
  //  trNew->Branch("CalELayerInitialRawRatioLog", &CalELayerInitialRawRatioLog, "CalELayerInitialRawRatioLog/F");
  //Int_t Cal1MomParticleInLayerCor;
  //trNew->Branch("Cal1MomParticleInLayerCor", &Cal1MomParticleInLayerCor, "Cal1MomParticleInLayerCor/I");
  //Float_t CalELayerCorrInitialRatioCorLog;
  //trNew->Branch("CalELayerCorrInitialRatioCorLog", &CalELayerCorrInitialRatioCorLog, "CalELayerCorrInitialRatioCorLog/F");
  //Float_t CalELayer34afterInitialRatioCorLog;
  //trNew->Branch("CalELayer34afterInitialRatioCorLog", &CalELayer34afterInitialRatioCorLog, "CalELayer34afterInitialRatioCorLog/F");
  //Float_t CalELayerInitialRawRatioCorLog;
  //trNew->Branch("CalELayerInitialRawRatioCorLog", &CalELayerInitialRawRatioCorLog, "CalELayerInitialRawRatioCorLog/F");
  Float_t Acd2TileEnergyRatioLog;
  trNew->Branch("Acd2TileEnergyRatioLog", &Acd2TileEnergyRatioLog, "Acd2TileEnergyRatioLog/F");

  Float_t numLayerIn = -1;
  Float_t numLayerInCor = -1;
  for(int iEvt=0; iEvt<nEvt; iEvt++)
    {
      trDat->GetEntry(iEvt);
      if(Abs(Cal1MomXDir)>minDir && Abs(Cal1MomYDir)>minDir)
	{
	  Cal1MomZCrossSide840 = Min(maxZCross, (Cal1MomZCntr+Min((840-((Cal1MomXDir>0)*2-1)*Cal1MomXCntr)/(((Cal1MomXDir>0)*2-1)*Cal1MomXDir), (840-((Cal1MomYDir>0)*2-1)*Cal1MomYCntr)/(((Cal1MomYDir>0)*2-1)*Cal1MomYDir))*Cal1MomZDir));
	  Cal1MomZCrossSide714 = Min(maxZCross, (Cal1MomZCntr+Min((714-((Cal1MomXDir>0)*2-1)*Cal1MomXCntr)/(((Cal1MomXDir>0)*2-1)*Cal1MomXDir), (714-((Cal1MomYDir>0)*2-1)*Cal1MomYCntr)/(((Cal1MomYDir>0)*2-1)*Cal1MomYDir))*Cal1MomZDir));
	  //Cal1MomZCrossSide840Cor = Min(maxZCross, (Cal1MomZCntrCor+Min((840-((Cal1MomXDir>0)*2-1)*Cal1MomXCntrCor)/(((Cal1MomXDir>0)*2-1)*Cal1MomXDir), (840-((Cal1MomYDir>0)*2-1)*Cal1MomYCntrCor)/(((Cal1MomYDir>0)*2-1)*Cal1MomYDir))*Cal1MomZDir));
	  //Cal1MomZCrossSide714Cor = Min(maxZCross, (Cal1MomZCntrCor+Min((714-((Cal1MomXDir>0)*2-1)*Cal1MomXCntrCor)/(((Cal1MomXDir>0)*2-1)*Cal1MomXDir), (714-((Cal1MomYDir>0)*2-1)*Cal1MomYCntrCor)/(((Cal1MomYDir>0)*2-1)*Cal1MomYDir))*Cal1MomZDir));
	}
      else if(Abs(Cal1MomXDir)>minDir)
	{
	  Cal1MomZCrossSide840 = Min(maxZCross, (Cal1MomZCntr+(840-((Cal1MomXDir>0)*2-1)*Cal1MomXCntr)/(((Cal1MomXDir>0)*2-1)*Cal1MomXDir)*Cal1MomZDir));
	  Cal1MomZCrossSide714 = Min(maxZCross, (Cal1MomZCntr+(714-((Cal1MomXDir>0)*2-1)*Cal1MomXCntr)/(((Cal1MomXDir>0)*2-1)*Cal1MomXDir)*Cal1MomZDir));
	  //Cal1MomZCrossSide840Cor = Min(maxZCross, (Cal1MomZCntrCor+(840-((Cal1MomXDir>0)*2-1)*Cal1MomXCntrCor)/(((Cal1MomXDir>0)*2-1)*Cal1MomXDir)*Cal1MomZDir));
	  //Cal1MomZCrossSide714Cor = Min(maxZCross, (Cal1MomZCntrCor+(714-((Cal1MomXDir>0)*2-1)*Cal1MomXCntrCor)/(((Cal1MomXDir>0)*2-1)*Cal1MomXDir)*Cal1MomZDir));
	}
      else if(Abs(Cal1MomYDir)>minDir)
	{
	  Cal1MomZCrossSide840 = Min(maxZCross, (Cal1MomZCntr+(840-((Cal1MomYDir>0)*2-1)*Cal1MomYCntr)/(((Cal1MomYDir>0)*2-1)*Cal1MomYDir)*Cal1MomZDir));
	  Cal1MomZCrossSide714 = Min(maxZCross, (Cal1MomZCntr+(714-((Cal1MomYDir>0)*2-1)*Cal1MomYCntr)/(((Cal1MomYDir>0)*2-1)*Cal1MomYDir)*Cal1MomZDir));
	  //Cal1MomZCrossSide840Cor = Min(maxZCross, (Cal1MomZCntrCor+(840-((Cal1MomYDir>0)*2-1)*Cal1MomYCntrCor)/(((Cal1MomYDir>0)*2-1)*Cal1MomYDir)*Cal1MomZDir));
	  //Cal1MomZCrossSide714Cor = Min(maxZCross, (Cal1MomZCntrCor+(714-((Cal1MomYDir>0)*2-1)*Cal1MomYCntrCor)/(((Cal1MomYDir>0)*2-1)*Cal1MomYDir)*Cal1MomZDir));
	}
      else
	{
	  Cal1MomZCrossSide840 = maxZCross; //1000000;
	  Cal1MomZCrossSide714 = maxZCross; //1000000;
	  //Cal1MomZCrossSide840Cor = maxZCross; //1000000;
	  //Cal1MomZCrossSide714Cor = maxZCross; //1000000;
	  cout << endl;
	  cout << "Perpendicular event!" << endl;
	}
      numLayerIn = Min(float(7.0), Max(float(0.0), -(Cal1MomZCrossSide714-calZTop)/cellVertPitch));
      Cal1MomParticleInLayer = int(numLayerIn);

      CalELayerCorrInitialRatioLog = (Log10(Max(minLimit,CalEnergyCorr)) - Log10(Max(minLimit, CalELayer[Cal1MomParticleInLayer]))); //Log10(Max(minLimit,CalEnergyCorr)) - ( ( numLayerIn<1 )*Log10(Max(minLimit, CalELayer[0])) + ( numLayerIn>=1 )*( numLayerIn<2 )*Log10(Max(minLimit, CalELayer[1])) + ( numLayerIn>=2 )*( numLayerIn<3 )*log10(max(minLimit, trDat.CalELayer2)) + ( numLayerIn>=3 )*( numLayerIn<4 )*log10(max(minLimit, trDat.CalELayer3)) + ( numLayerIn>=4 )*( numLayerIn<5 )*log10(max(minLimit, trDat.CalELayer4)) + ( numLayerIn>=5 )*( numLayerIn<6 )*log10(max(minLimit, trDat.CalELayer5)) + ( numLayerIn>=6 )*( numLayerIn<7 )*log10(max(minLimit, trDat.CalELayer6)) )

      CalELayer34afterInitialRatioLog = Log10(Max(minLimit, CalELayer[Min(3,Cal1MomParticleInLayer)+4]+CalELayer[Min(3,Cal1MomParticleInLayer)+3]))-Log10(Max(minLimit, CalELayer[Min(3,Cal1MomParticleInLayer)]));
//( numLayerIn<1 ) * ( log10(max(minLimit, trDat.CalELayer4+trDat.CalELayer3))-log10(max(minLimit, trDat.CalELayer0)) ) + ( numLayerIn>=1 ) * ( numLayerIn<2 ) * ( log10(max(minLimit, trDat.CalELayer5+trDat.CalELayer4))-log10(max(minLimit, trDat.CalELayer1)) ) + ( numLayerIn>=2 ) * ( numLayerIn<3 ) * ( log10(max(minLimit, trDat.CalELayer6+trDat.CalELayer5))-log10(max(minLimit, trDat.CalELayer2)) ) + ( numLayerIn>=3 ) * ( log10(max(minLimit, trDat.CalELayer7+trDat.CalELayer6))-log10(max(minLimit, trDat.CalELayer3)) )

      //CalELayerInitialRawRatioLog = ( Log10(Max(minLimit,CalEnergyRaw)) - Log10(Max(minLimit, CalELayer[Cal1MomParticleInLayer]) ) ); //Log10(Max(minLimit,CalEnergyRaw)) - ( ( numLayerIn<1 )*log10(max(minLimit, trDat.CalELayer0)) + ( numLayerIn>=1 )*( numLayerIn<2 )*log10(max(minLimit, trDat.CalELayer1)) + ( numLayerIn>=2 )*( numLayerIn<3 )*log10(max(minLimit, trDat.CalELayer2)) + ( numLayerIn>=3 )*( numLayerIn<4 )*log10(max(minLimit, trDat.CalELayer3)) + ( numLayerIn>=4 )*( numLayerIn<5 )*log10(max(minLimit, trDat.CalELayer4)) + ( numLayerIn>=5 )*( numLayerIn<6 )*log10(max(minLimit, trDat.CalELayer5)) + ( numLayerIn>=6 )*( numLayerIn<7 )*log10(max(minLimit, trDat.CalELayer6)) )

      //numLayerInCor = Min(float(7.0), Max(float(0.0), -(Cal1MomZCrossSide714Cor-calZTop)/cellVertPitch));
      //Cal1MomParticleInLayerCor = int(numLayerInCor);
      //CalELayerCorrInitialRatioCorLog = (Log10(Max(minLimit,CalEnergyCorr)) - Log10(Max(minLimit, CalELayer[Cal1MomParticleInLayerCor])));
      //CalELayer34afterInitialRatioCorLog = Log10(Max(minLimit, CalELayer[Min(3,Cal1MomParticleInLayerCor)+4]+CalELayer[Min(3,Cal1MomParticleInLayerCor)+3])) - Log10(Max(minLimit, CalELayer[Min(3,Cal1MomParticleInLayerCor)]));
      //CalELayerInitialRawRatioCorLog = ( Log10(Max(minLimit,CalEnergyRaw)) - Log10(Max(minLimit, CalELayer[Cal1MomParticleInLayerCor]) ) ); 

      Acd2TileEnergyRatioLog = Log10(Max(Acd2TileEnergy/Max(float(10.), EvtJointEnergy) * 100, minLimit));

      if(iEvt%100000==0)
	{
	  cout << "." << flush;
	  if(iEvt%1000000==0)
	    {
	      cout << Form("\n%d%% have been filled.\n", (trNew->GetEntries()*100)/nEvt) << flush;
	      cout << "Acd2TileEnergyRatioLog: " << Acd2TileEnergyRatioLog << endl;
	    }
	}
      trNew->Fill();
    }
  cout << endl;
  trNew->Write();
  delete trNew;
  fileNew->Close();
  fileIn->Close();
}
