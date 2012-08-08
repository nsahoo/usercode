#ifndef TreeModifiers_LeptSfProvider_h
#define TreeModifiers_LeptSfProvider_h

#include "TTree.h"
#include "TH2.h"
#include <string>
#include <iostream>

/** Class cache SF values from 2D histograms and returns them through convenient methods.
 */

class LeptSfProvider {
 public:
 LeptSfProvider() : is2011_(false)  {    }
  virtual ~LeptSfProvider(){};

  float getMuonIDSF(float pt, float eta);

  float getMuonIPSF(float pt, float eta);

  float getMuonIsoSF(float pt, float eta);

  float getMuonSF(float pt, float eta);

  float getElectronSF(float pt, float eta);
  
  float getPR(float pt, float eta, float id);

  float getSF(float pt, float eta, float id);

  void init(bool is2011){
    is2011_ = is2011;
    initMu(is2011);
    initEl(is2011);
  }

  void initMu(bool is2011); 

  void initEl(bool is2011); 
  

 private:  
  bool is2011_;
  Double_t* IP_eta; 
  Double_t* IP_pt_barrel; 
  Double_t* IP_pt_endcaps; 

  Double_t* ISO_eta; 
  Double_t* ISO_pt_barrel; 
  Double_t* ISO_pt_endcaps; 

  Double_t* ID_eta_2011A; 
  Double_t* ID_eta_2011B; 
  Double_t* ID_eta; 
  Double_t* ID_pt_barrel; 
  Double_t* ID_pt_endcaps; 


  TH2* electronscalefactors;
  std::vector<std::vector<float> > elesf;


};

#endif
