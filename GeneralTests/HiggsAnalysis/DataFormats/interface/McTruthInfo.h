#ifndef HiggsAnalysis_DataFormats_McTruthInfo_h
#define HiggsAnalysis_DataFormats_McTruthInfo_h

#include <TObject.h>
#include <TLorentzVector.h>
#include <vector>

using namespace std;

#include "HiggsAnalysis/DataFormats/interface/Lepton.h"

//enum flag{one, two};

class McTruthInfo 
{
 public:
  //ClassDef(MyEvent, 1);
  
  McTruthInfo():higgs_(),z1_(),z2_(){;}
  virtual ~McTruthInfo(){;}
  
  void setHiggs(TLorentzVector& input){higgs_ = input;}
  void setZ1(TLorentzVector& input){z1_ = input;}
  void setZ2(TLorentzVector& input){z2_ = input;}

  TLorentzVector& higgs() {return higgs_;}
  TLorentzVector& z1() {return z1_;}
  TLorentzVector& z2() {return z2_;}

  /*
  unsigned short numEs() {return e.size();}
  unsigned short numZs() {return Z.size();}
  unsigned short numHs() {return H.size();}
  
  //HBoson* H1()  {return !H.empty() ? H[0] : new HBoson();}
  
  ZBoson* Z1() {return !H.empty() ? H[0]->Z1() : new ZBoson();}
  ZBoson* Z2() {return !H.empty() ? H[0]->Z2() : new ZBoson();}

  Lepton* e1() {return !H.empty() ? H[0]->Z1()->e1() : new Lepton();}
  Lepton* e2() {return !H.empty() ? H[0]->Z1()->e2() : new Lepton();}
  Lepton* e3() {return !H.empty() ? H[0]->Z2()->e1() : new Lepton();}
  Lepton* e4() {return !H.empty() ? H[0]->Z2()->e2() : new Lepton();}
  */

  TLorentzVector higgs_,z1_,z2_;

  //std::vector<Lepton*> e;
  //std::vector<ZBoson*> Z;
  //std::vector<HBoson*> H;
  

};

#endif
