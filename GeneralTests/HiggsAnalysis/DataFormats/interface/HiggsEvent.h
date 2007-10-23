#ifndef HiggsAnalysis_DataFormats_HiggsEvent_h
#define HiggsAnalysis_DataFormats_HiggsEvent_h

#include <TObject.h>
#include <TLorentzVector.h>
#include <vector>

using namespace std;

#include "HiggsAnalysis/DataFormats/interface/Lepton.h"
#include "HiggsAnalysis/DataFormats/interface/McTruthInfo.h"
//#include "HiggsAnalysis/DataFormats/interface/ZBoson.h"
//#include "HiggsAnalysis/DataFormats/interface/HBoson.h"

//enum flag{one, two};

//class HiggsEvent : public TObject
class HiggsEvent 
{
 public:
  //ClassDef(MyEvent, 1);
  
  HiggsEvent();
  virtual ~HiggsEvent();
  
  McTruthInfo& mc() {return mc_;}

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
  
  std::vector<Lepton*> e;
  std::vector<ZBoson*> Z;
  std::vector<HBoson*> H;
  */

  McTruthInfo mc_;
};

#endif
