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


class BDTbranch
{
private:
  std::string name_branch;
  std::string type_value;
  std::string path_file;
  std::string path_rootdir;
  std::string name_tree;
  Float_t cth_min;
  Float_t cth_max;
  Float_t loge_min;
  Float_t loge_max;
  TFile *fileBDT
  TTree *treeBDT;

public:
  BDTbranch(std::string name_branch, std::string path_file, std::string name_tree="MeritTuple", float cth_min=-99999999, float cth_max=99999999, float loge_min=0., float loge_max=10.,std::string  path_rootdir="/");
  ~BDTbranch();
  std::string GetNameBranch();
  std::string GetTypeValue();
  std::string GetPathFile();
  std::string GetPathRootDir();
  std::string GetNameTree();
  Int_t GetEntry(Int_t ievt);
  Int_t GetEntries();
};

BDTbranch::BDTbranch(std::string name_branch, std::string path_file, std::string name_tree="MeritTuple", float cth_min=-99999999, float cth_max=99999999, float loge_min=0., float loge_max=10., std::string  path_rootdir="")
{
  name_branch = name_branch;
  type_value = path_file;
  path_file = name_tree;
  path_rootdir = path_rootdir;
  name_tree = name_tree;
  path_roottree = path_rootdir + "/" + name_tree;
  cth_min = cth_min;
  cth_max = cth_max;
  loge_min = loge_min;
  loge_max = loge_max;
  fileBDT = new TFile(path_file, "READ");
  tree = (TTree*)fileBDT->Get(path_roottree.c_str());
}

BDTbranch::~BDTbranch()
{
  delete fileBDT;
};


std::string BDTbranch::GetNameBranch()
{
  return name_branch;
}

std::string BDTbranch::GeTypeValue()
{
  return type_value;
}

std::string BDTbranch::GetPathFile()
{
  return path_file;
}

std::string BDTbranch::GetPathRootDir()
{
  return path_rootdir;
}

std::string BDTbranch::GetNameTree()
{
  return name_tree;
}

BDTbranch::GetEntry(Int_t ievt)
{
  return tree->GetEntry(ievt);
}

BDTbranch::GetEntries()
{
  return tree->GetEntry();
}


class BDTbranch_float : public BDTbranch
{
private:
  Float_t braddress = 0;

public:
  Int_t GetValue();
};

BDTbranch_float::BDTbranch_float(std::string name_branch, std::string path_file, std::string name_tree="MeritTuple", float cth_min=-99999999, float cth_max=99999999, float loge_min=0., float loge_max=10., std::string  path_rootdir="") //Constructor
{
  BDTbranch:: BDTbranch(name_branch, path_file, name_tree, cth_min, cth_max, loge_min, loge_max, path_rootdir);
  tree->SetBranchAddress(name_branch.c_str(), &braddress);
}

BDTbranch_float::GetValue()
{
  return braddress;
}


int main(int argc,char *argv[])
{
  string strPathFileMain = (string)argv[1];
  TFile *fileMain = new TFile(strPathFileMain.c_str(), "UPDATE");
  TTree *trMain =  (TTree*)fileMain->Get("MeritTuple");
  const Int_t nEvtMain = trMain->GetEntries();
  cout << fileMain->GetName() << " has " << nEvtMain <<  " events." << endl;
  Float_t cal1momzdir=0.;
  Float_t wp8calonlyenergy=0.;
  trMain->SetBranchAddress("Cal1MomZDir", &cal1momzdir);
  trMain->SetBranchAddress("WP8CalOnlyEnergy", &wp8calonlyenergy);

  const Int_t nBDTs = 5;
  const string strNameBDTs[nBDTs] = {"S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR060100_BDTG500D06", "S20_S025B025_CALE435575RAWE20ZDIR020ZCS000UD000catTwoZDIR_17_UpDownZDir_ZDIR020060_BDTG500D06", "T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_11Dir_ZDIR000020_BDTG500D06", "T02V200909_IRFBK_020RAWE20wwoTRKwoRWwUpdown0catTwoZ_11Dir_ZDIR050100_BDTG500D06"};
  const float cthMaxs[nBDTs] = {99999999., 0.6, 0.2, 0.0, -0.2, -0.5};
  const float cthMins[nBDTs] = {0.6, 0.2, 0.0, -0.2, -0.5, -99999999.};
  string strPathFileBDT;
  BDTbranch_float::BDTbranch_float bdt_branches[nBDTs];
  for(int ibdt=0; ibdt<nBDTs; ibdt++)
    {
      strPathFileBDT = strPathFileMain.substr(0, strPathFileMain.length()-5) + "_" + strNameBDTs[ibdt] + ".root";
      bdt_branches[ibdt] = BDTbranch_float:BDTbranch_float(strNameBDTs[ibdt], strPathFileBDT, "MeritTuple", cthMins[ibdt], cMaxs[ibdt]);
      Int_t nEvtBDT = bdt_branches[ibdt]->GetEntries();
      if(nEvtBDT!=nEvtMain)
	{
	  cout << bdt_branches[ibdt]->GetNameBranch() << " has " << nEvtBDT <<  " events!" << endl;
	  cout << "Different from the number of the main tree!!" << endl;
	}
    }

  for(int ievt=0; ievt<NEvtMain; ievt++)
    {
      trMain->GetEntry(ievt);

      

    }

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
