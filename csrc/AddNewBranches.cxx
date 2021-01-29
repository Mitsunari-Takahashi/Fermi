#include <stdio.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <sstream>
#include <stdint.h>
#include <getopt.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"


using namespace std;
using namespace TMath;

int main(int argc,char *argv[])
{
  string strPathFileMain = (string)argv[1];
  TFile *fileMain = new TFile(strPathFileMain.c_str(), "UPDATE");
  TTree *trMain =  (TTree*)fileMain->Get("MeritTuple");
  const Int_t nEvtMain = trMain->GetEntries();
  cout << fileMain->GetName() << " has " << nEvtMain <<  " events." << endl;

  string strPathFileAdd = argv[2];
  TFile* fileAdd = new TFile(strPathFileAdd.c_str(), "READ");
  TTree *trAdd =  (TTree*)fileAdd->Get("MeritTuple");
  const Int_t nEvtAdd = trAdd->GetEntries();
  cout << fileAdd->GetName() << " has " << nEvtAdd <<  " events." << endl;

  if(nEvtMain!=nEvtAdd)
    cout << "The event numbers in the two files do not match to each other!" << endl;

  Int_t nBranchesAdd = argc - 3;
  if(nBranchesAdd<1)
    cout << "Please specify the name of the branches to add." << endl;
  vector<string> strBranchesAdd;
  Float_t values[nBranchesAdd];
  Float_t valuesAdd[nBranchesAdd];
  for(int i=0; i<nBranchesAdd; i++)
    {
      values[i] = -99999999.;
      valuesAdd[i] = -99999999.;
    }

  cout << "Number of branches to add: " << nBranchesAdd << endl;
  TBranch* brNew[nBranchesAdd];
  for(int i=0; i<nBranchesAdd; i++)
    {
      strBranchesAdd.push_back(argv[i+3]);
      cout << strBranchesAdd[i] << endl;
      trAdd->SetBranchAddress(strBranchesAdd[i].c_str(), &(values[i]));
      brNew[i] = (TBranch*)trMain->Branch(strBranchesAdd[i].c_str(), &(valuesAdd[i]), Form("%s/%s", strBranchesAdd[i].c_str(), "F"));
      brNew[i]->Print();
    }
  cout << endl;

  for(int iEvt=0; iEvt<nEvtAdd; iEvt++)
    {
      trAdd->GetEntry(iEvt);
      if(iEvt%100000==0)
	cout << iEvt << "th event" << endl;
      for(int jBr=0; jBr<nBranchesAdd; jBr++)
	{
	  valuesAdd[jBr] = (Float_t)values[jBr];
	  if(iEvt%100000==0)
	    cout << strBranchesAdd[jBr] << "=" << valuesAdd[jBr] << endl;
	  brNew[jBr]->Fill();
	  if(iEvt%10000000==0)
	    brNew[jBr]->Print();
	}
      // if(iEvt%100000==0)
      // 	{
      //      cout << Form("\n%d%% have been filled.\n", (iEvt*100)/nEvtMain) << flush;
      // 	  cout << "." << flush;
      // 	  //if(iEvt%1000000==0)
      // 	}
      //cout << Form("\n%d events have been filled.\n", iEvt) << flush;
      //cout << iEvt << endl;
    }
  cout << endl;
  cout << nEvtAdd << " events have been filled." << endl;
  for(int jBr=0; jBr<nBranchesAdd; jBr++)
    brNew[jBr]->Print();
  fileMain->cd();
  trMain->Write();
  delete trMain;
  delete trAdd;
  fileAdd->Close();
  fileMain->Close();
}
