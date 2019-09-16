#include <stdio.h>
//#include <vector>
//#include <fstream>
//#include <sstream>
#include <iostream>
#include <string>
//#include <stdint.h>
#include "TMath.h"
//#include "TH2F.h"
//#include "TF1.h"
#include "TDirectory.h"
//#include "TTree.h"
#include "Variable.h"

//using namespace std;


const Double_t MinValue = -7.;
const Double_t MaxValue = 8.;
const Float_t MIN_LOGENERGY = 4.35;
const Float_t MAX_LOGENERGY = 5.75;
const Float_t MIN_ZDIR = 0.2;
const Float_t MAX_ZDIR = 1.0;

Variable::Variable(std::string name, std::string category, std::string definition)
{
  SetName(name);
  SetDefinition(definition);
  calCategory = false;
  acdCategory = false;
  tkrCategory = false;
  SetCategory(category);
  Print();
}

std::string Variable::GetName()
{
  return name;
}

int Variable::SetName(std::string newname)
{
  name = newname;
}

std::string Variable::GetDefinition()
{
  return definition;
}

int Variable::SetDefinition(std::string newdefinition)
{
  definition = newdefinition;
}

int Variable::SetCategory(std::string category)
{
  if(category=="CAL")
    calCategory = true;
  else if(category=="ACD")
    acdCategory = true;
  else if(category=="TKR")
    tkrCategory = true;
}

Int_t Variable::FindSourceId()
{
  Int_t srcid;
  if(calCategory==true)
    srcid = 1000;
  else if(acdCategory==true)
    srcid = 2000;
  return srcid;
}

void Variable::Print()
{
  std::cout << "Name: " << GetName() << std::endl;
  std::cout << "Definition: " << GetDefinition() << std::endl;
}

VariableFloat::VariableFloat(std::string name, std::string category, TDirectory* rfdir, std::string definition) : Variable(name, category, definition)
{
  value = 0.;
  original = 0.;
  TDirectory* vardir = (TDirectory*)rfdir->GetDirectory(name.c_str());
  f2d_fitSig = (TF2*)vardir->Get("f2D_Parabolic_Gamma");
  if(FindSourceId()==1000)
    f2d_fitBkg = (TF2*)vardir->Get("f2D_Parabolic_Hadron");
  else if(FindSourceId()==2000)
    f2d_fitBkg = (TF2*)vardir->Get("f2D_Parabolic_Lepton");
}

// int VariableFloat::SetOriginal(Float_t &original_ref)
// {
//   original = original_ref;
// }

void VariableFloat::Print()
{
  std::cout << "Name: " << GetName() << ", Definition: " << GetDefinition() << std::endl;
  std::cout << original << " -> " << value << " -> " << parametrized << std::endl;
}

int VariableFloat::CalcValue()
{
  value = original;
}

int VariableFloat::Parametrize(Double_t energy, Double_t zdir)
{
  meanSig = f2d_fitSig->Eval(energy, zdir);
  meanBkg = f2d_fitBkg->Eval(energy, zdir);
  if(meanSig==meanBkg)
    {
      parametrized = -10.;
      return 1;
    }
  a = 1./(meanBkg-meanSig);
  b = -a*meanSig;
  Float_t p = TMath::Min(MaxValue, TMath::Max(MinValue, a*value+b));
  if(std::isnan(p))
    {
      parametrized = -10.;
      return 1;
    }
  else
    {
      parametrized = p;
      return 0;
    }
}


/* VariableDouble */
VariableDouble::VariableDouble(std::string name, std::string category, TDirectory* rfdir, std::string definition) : Variable(name, category, definition)
{
  value = 0.;
  original = 0.;
  TDirectory* vardir = (TDirectory*)rfdir->GetDirectory(name.c_str());
  f2d_fitSig = (TF2*)vardir->Get("f2D_Parabolic_Gamma");
  if(FindSourceId()==1000)
    f2d_fitBkg = (TF2*)vardir->Get("f2D_Parabolic_Hadron");
  else if(FindSourceId()==2000)
    f2d_fitBkg = (TF2*)vardir->Get("f2D_Parabolic_Lepton");
}

void VariableDouble::Print()
{
  std::cout << "Name: " << GetName() << ", Definition: " << GetDefinition() << std::endl;
  std::cout << original << " -> " << value << " -> " << parametrized << std::endl;
}

int VariableDouble::CalcValue()
{
  value = original;
}

int VariableDouble::Parametrize(Double_t energy, Double_t zdir)
{
  meanSig = f2d_fitSig->Eval(energy, zdir);
  meanBkg = f2d_fitBkg->Eval(energy, zdir);
  if(meanSig==meanBkg)
    {
      parametrized = -10.;
      return 1;
    }
  a = 1./(meanBkg-meanSig);
  b = -a*meanSig;
  Float_t p = TMath::Min(MaxValue, TMath::Max(MinValue, a*value+b));
  if(std::isnan(p))
    {
      parametrized = -10.;
      return 1;
    }
  else
    {
      parametrized = p;
      return 0;
    }
}

VariableFloatLog10::VariableFloatLog10(std::string name, std::string category, TDirectory* rfdir, Float_t limit, std::string definition) : VariableFloat(name, category, rfdir, definition)
{
  SetLowLim(limit);
}
 
int VariableFloatLog10::SetLowLim(Float_t limit)
{
  lowlim = limit;
}

Float_t VariableFloatLog10::GetLowLim()
{
  return lowlim;
}

int VariableFloatLog10::CalcValue()
{
  //std::cout << "VariableFloatLog10" << std::endl;
  value = TMath::Log10(TMath::Max(lowlim, original));
}

// int VariableFloatLog10::Parametrize(Double_t energy, Double_t zdir)
// {
//   meanSig = f2d_fitSig->Eval(energy, zdir);
//   meanBkg = f2d_fitBkg->Eval(energy, zdir);
//   a = 1./(meanBkg-meanSig);
//   b = -a*meanSig;
//   parametrized = TMath::Min(MaxValue, TMath::Max(MinValue, a*value+b));
// }

VariableFloatRatioLog10::VariableFloatRatioLog10(std::string name, std::string category, TDirectory* rfdir, Float_t limit, std::string definition) : VariableFloatLog10(name, category, rfdir, limit, definition)
{
  //  SetLowLim(limit);
  original2 = 0.;
}

void VariableFloatRatioLog10::Print()
{
  std::cout << "Name: " << GetName() << ", Definition: " << GetDefinition() << std::endl;
  std::cout << "(" << original << ", " << original2 << ") -> " << value << " -> " << parametrized << std::endl;
}

int VariableFloatRatioLog10::CalcValue()
{
  value = TMath::Log10(TMath::Max(lowlim, original))-TMath::Log10(TMath::Max(lowlim, original2));
}


VariableCalELayerCorrInitialLog::VariableCalELayerCorrInitialLog(std::string name, std::string category, TDirectory* rfdir, Float_t limit, std::string definition) : VariableFloatLog10(name, category, rfdir, limit, definition)
{
  cal1MomParticleInLayerCor = 0;
  for(int i=0; i<8; i++)
    calELayer[i] = 0.;
}

int VariableCalELayerCorrInitialLog::CalcValue()
{
  value = TMath::Log10(TMath::Max(lowlim, calELayer[cal1MomParticleInLayerCor]));
}

/* VariableUInt */
VariableUInt::VariableUInt(std::string name, std::string category, TDirectory* rfdir, std::string definition) : Variable(name, category, definition)
{
  value = 0.;
  original = 0;
  TDirectory* vardir = (TDirectory*)rfdir->GetDirectory(name.c_str());
  f2d_fitSig = (TF2*)vardir->Get("f2D_Parabolic_Gamma");
  if(FindSourceId()==1000)
    f2d_fitBkg = (TF2*)vardir->Get("f2D_Parabolic_Hadron");
  else if(FindSourceId()==2000)
    f2d_fitBkg = (TF2*)vardir->Get("f2D_Parabolic_Lepton");
}

void VariableUInt::Print()
{
  std::cout << "Name: " << GetName() << ", Definition: " << GetDefinition() << std::endl;
  std::cout << original << " -> " << value << " -> " << parametrized << std::endl;
}

int VariableUInt::CalcValue()
{
  value = original;
}

int VariableUInt::Parametrize(Double_t energy, Double_t zdir)
{
  meanSig = f2d_fitSig->Eval(energy, zdir);
  meanBkg = f2d_fitBkg->Eval(energy, zdir);
  if(meanSig==meanBkg)
    {
      parametrized = -10.;
      return 1;
    }
  a = 1./(meanBkg-meanSig);
  b = -a*meanSig;
  Float_t p = TMath::Min(MaxValue, TMath::Max(MinValue, a*value+b));
  if(std::isnan(p))
    {
      parametrized = -10.;
      return 1;
    }
  else
    {
      parametrized = p;
      return 0;
    }
}

/* UIntLog10 */
VariableUIntLog10::VariableUIntLog10(std::string name, std::string category, TDirectory* rfdir, Float_t limit, std::string definition) : VariableUInt(name, category, rfdir, definition)
{
  SetLowLim(limit);
}
 
int VariableUIntLog10::SetLowLim(Float_t limit)
{
  lowlim = limit;
}

Float_t VariableUIntLog10::GetLowLim()
{
  return lowlim;
}

int VariableUIntLog10::CalcValue()
{
  value = TMath::Log10(TMath::Max(lowlim, (Float_t)original));
  //  std::cout << "original=" << original << ", value=" << value << std::endl;
}

// void VariableUIntLog10::Print()
// {
//   std::cout << "Name: " << GetName() << ", Definition: " << GetDefinition() << std::endl;
//   std::cout << original << " -> " << value << " -> " << parametrized << std::endl;
// }
