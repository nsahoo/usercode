#ifndef HiggsAnalysis_DataFormats_McInfo_h
#define HiggsAnalysis_DataFormats_McInfo_h

#include "DataFormats/Math/interface/LorentzVector.h"
#include <vector>
#include <iostream>

using namespace std;



class McInfo 
{
 public:  
  McInfo():higgsMass_(-99),higgsP_(0,0,0,0)
    //:higgs_()
    //,z1_(),z2_()
    {}
  virtual ~McInfo();
  
  
  //void setHiggs(TLorentzVector& input){
  // cout << "=== in setHiggs: " << input.Perp() << endl;
  //higgs_ = input;}
  //void setZ1(TLorentzVector& input){z1_ = input;}
  //void setZ2(TLorentzVector& input){z2_ = input;}
  
  math::XYZTLorentzVector higgs_p4() const {return higgsP_;}
  double higgs_mass() const {return (higgsMass_);}
  //TLorentzVector& z1() {return z1_;}
  //TLorentzVector& z2() {return z2_;}
  
  //TLorentzVector higgs_;
  //z1_,z2_;
      
  //std::vector<Lepton*> e;
  //std::vector<ZBoson*> Z;
  //std::vector<HBoson*> H;
  
  double higgsMass_;
  math::XYZTLorentzVector higgsP_;
  
};
#endif
