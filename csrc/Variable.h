#include "TF2.h"

class Variable
{
private:
  std::string name;
  std::string definition;
  Bool_t calCategory;
  Bool_t acdCategory;
  Bool_t tkrCategory;

public:
  Variable(std::string name, std::string category, std::string definition="");
  std::string GetName();
  int SetName(std::string newname);
  std::string GetDefinition();
  int SetDefinition(std::string newdefinition);
  int SetCategory(std::string category);
  Int_t FindSourceId();
  void Print();
};
 

class VariableFloat : public Variable
{
private:
  TF2* f2d_fitSig;
  TF2* f2d_fitBkg;
  Double_t meanSig;
  Double_t meanBkg;
  Double_t a;
  Double_t b;
public:
  VariableFloat(std::string name, std::string category, TDirectory* rfdir, std::string definition="");
  Float_t value;
  Float_t original;
  Float_t parametrized;
  //int SetOriginal(Float_t &original_ref);
  int CalcValue();
  int Parametrize(Double_t energy, Double_t zdir);
  void Print();
  //int SetFitFunction(TDirectory rfdir);
};


class VariableFloatLog10 : public VariableFloat
{
 private:
  Float_t lowlim;
 public:
  //  Float_t value;
  //  Float_t original;
  //  Float_t parametrized;
  VariableFloatLog10(std::string name, std::string category, TDirectory* rfdir, Float_t limit=1E-5, std::string definition="");
  int CalcValue();
  int SetLowLim(Float_t limit);
  Float_t GetLowLim();
};

class VariableFloatRatioLog10 : public VariableFloatLog10
{
 private:
  Float_t lowlim;
 public:
  //  Float_t value;
  //  Float_t original;
  //  Float_t parametrized;
  VariableFloatRatioLog10(std::string name, std::string category, TDirectory* rfdir, Float_t limit=1E-5, std::string definition="");
  void Print();
  int CalcValue();
  int SetLowLim(Float_t limit);
  Float_t GetLowLim();
  Float_t original2;
};

class VariableCalELayerCorrInitialLog : public VariableFloatLog10
{
 private:
  Float_t lowlim;
 public:
  VariableCalELayerCorrInitialLog(std::string name, std::string category, TDirectory* rfdir, Float_t lowlim=1E-5, std::string definition="");
  int CalcValue();
  int SetLowLim(Float_t limit);
  Float_t GetLowLim();
  Int_t cal1MomParticleInLayerCor;
  Float_t calELayer[8];
};

class VariableDouble : public Variable
{
private:
  TF2* f2d_fitSig;
  TF2* f2d_fitBkg;
  Double_t meanSig;
  Double_t meanBkg;
  Double_t a;
  Double_t b;
public:
  VariableDouble(std::string name, std::string category, TDirectory* rfdir, std::string definition="");
  Float_t value;
  Double_t original;
  Float_t parametrized;
  int CalcValue();
  int Parametrize(Double_t energy, Double_t zdir);
  void Print();
  //int SetFitFunction(TDirectory rfdir);
};


class VariableUInt : public Variable
{
private:
  TF2* f2d_fitSig;
  TF2* f2d_fitBkg;
  Double_t meanSig;
  Double_t meanBkg;
  Double_t a;
  Double_t b;
public:
  VariableUInt(std::string name, std::string category, TDirectory* rfdir, std::string definition="");
  Float_t value;
  UInt_t original;
  Float_t parametrized;
  int CalcValue();
  int Parametrize(Double_t energy, Double_t zdir);
  void Print();
  //int SetFitFunction(TDirectory rfdir);
};

class VariableUIntLog10 : public VariableUInt
{
 private:
  Float_t lowlim;
 public:
  //Float_t value;
  //UInt_t original;
  //Float_t parametrized;
  VariableUIntLog10(std::string name, std::string category, TDirectory* rfdir, Float_t limit=1E-5, std::string definition="");
  int CalcValue();
  int SetLowLim(Float_t limit);
  //  void Print();
  //  int Parametrize(Double_t energy, Double_t zdir);
  Float_t GetLowLim();
};
