#ifndef TreeModifiers_XSecProvider_h
#define TreeModifiers_XSecProvider_h

#include "TTree.h"
#include "TH2.h"
#include <string>
#include <iostream>

/** Class cache XSection, BR, etc.. values and returns them through convenient methods.
 */

class XSecProvider {
 public:
 XSecProvider()  {    }
  virtual ~XSecProvider(){};

  float getXSec7TeVByID(int id);

  float getXSec8TeVByID(int id);

  float getEvt7TeV(int id){
    return evt_7TeV[id];
  }

  float getEvt8TeV(int id){
    return evt_8TeV[id];
  }

  float get7TeVggHXsec(float mass) {
    return xsec_7TeV_ggH[mass];
  }    

  float get8TeVggHXsec(float mass) {
    return xsec_8TeV_ggH[mass];
  }    

  float get7TeVVBFXsec(float mass) {
    return xsec_7TeV_VBF[mass];
  }    

  float get8TeVVBFXsec(float mass) {
    return xsec_8TeV_VBF[mass];
  }    

  float get7TeVggHXsecByChannel(float mass, int ch) {
    if (ch == 0 || ch == 1)return xsec_7TeV_ggH_SF[mass];
    else if (ch == 2) return xsec_7TeV_ggH_OF[mass];
    else return 0.;
  }

  float get8TeVggHXsecByChannel(float mass, int ch) {
    if (ch == 0 || ch == 1)return xsec_8TeV_ggH_SF[mass];
    else if (ch == 2) return xsec_8TeV_ggH_OF[mass];
    else return 0.;
  }

  float get7TeVVBFXsecByChannel(float mass, int ch) {
    if (ch == 0 || ch == 1)return xsec_7TeV_VBF_SF[mass];
    else if (ch == 2) return xsec_7TeV_VBF_OF[mass];
    else return 0.;
  }

  float get8TeVVBFXsecByChannel(float mass, int ch) {
    if (ch == 0 || ch == 1)return xsec_8TeV_VBF_SF[mass];
    else if (ch == 2) return xsec_8TeV_VBF_OF[mass];
    else return 0.;
  }

  float get7TeVggHQCDUp(float mass) {
    return qcdu_7TeV_ggH[mass];
  }    
  
  float get7TeVggHQCDDown(float mass) {
    return qcdd_7TeV_ggH[mass];
  }    

  float get7TeVVBFQCDUp(float mass) {
    return qcdu_7TeV_VBF[mass];
  }   
  
  float get7TeVVBFQCDDown(float mass) {
    return qcdd_7TeV_VBF[mass];
  }   

  float get8TeVggHQCDUp(float mass) {
    return qcdu_8TeV_ggH[mass];
  }    

  float get8TeVggHQCDDown(float mass) {
    return qcdd_8TeV_ggH[mass];
  }    

  float get8TeVVBFQCDUp(float mass) {
    return qcdu_8TeV_VBF[mass];
  }   

  float get8TeVVBFQCDDown(float mass) {
    return qcdd_8TeV_VBF[mass];
  }   

  float get7TeVggHPDFUp(float mass) {
    return qcdu_7TeV_ggH[mass];
  }   

  float get7TeVggHPDFDown(float mass) {
    return qcdd_7TeV_ggH[mass];
  }   

  float get7TeVVBFPDFUp(float mass) {
    return qcdu_7TeV_VBF[mass];
  }   

  float get7TeVVBFPDFDown(float mass) {
    return qcdd_7TeV_VBF[mass];
  }   

  float get8TeVggHPDFUp(float mass) {
    return qcdu_8TeV_ggH[mass];
  }

  float get8TeVggHPDFDown(float mass) {
    return qcdd_8TeV_ggH[mass];
  }

  float get8TeVVBFPDFUp(float mass) {
    return qcdu_8TeV_VBF[mass];
  }    

  float get8TeVVBFPDFDown(float mass) {
    return qcdd_8TeV_VBF[mass];
  }

  float get7TeVBkgXsec(int id) {
    return xsec_7TeV_bkg[id];
  }

  float get8TeVBkgXsec(int id) {
    return xsec_8TeV_bkg[id];
  }

  void initXsec();

  void initQCDScale();
    
  void initPDF();



 private:  
  std::map<int, float>   xsec_7TeV_bkg;
  std::map<int, float>   xsec_8TeV_bkg;

  std::map<float, float> xsec_7TeV_ggH;
  std::map<float, float> xsec_8TeV_ggH;
  std::map<float, float> xsec_7TeV_VBF;
  std::map<float, float> xsec_8TeV_VBF;

  std::map<float, float> xsec_7TeV_ggH_SF;
  std::map<float, float> xsec_8TeV_ggH_SF;
  std::map<float, float> xsec_7TeV_VBF_SF;
  std::map<float, float> xsec_8TeV_VBF_SF;
  
  std::map<float, float> xsec_7TeV_ggH_OF;
  std::map<float, float> xsec_8TeV_ggH_OF;
  std::map<float, float> xsec_7TeV_VBF_OF;
  std::map<float, float> xsec_8TeV_VBF_OF;
  
  std::map<float, float> qcdu_7TeV_ggH;
  std::map<float, float> qcdu_8TeV_ggH;
  std::map<float, float> qcdu_7TeV_VBF;
  std::map<float, float> qcdu_8TeV_VBF;

  std::map<float, float> qcdd_7TeV_ggH;
  std::map<float, float> qcdd_8TeV_ggH;
  std::map<float, float> qcdd_7TeV_VBF;
  std::map<float, float> qcdd_8TeV_VBF;

  std::map<float, float> pdfu_7TeV_ggH;
  std::map<float, float> pdfu_8TeV_ggH;
  std::map<float, float> pdfu_7TeV_VBF;
  std::map<float, float> pdfu_8TeV_VBF;

  std::map<float, float> pdfd_7TeV_ggH;
  std::map<float, float> pdfd_8TeV_ggH;
  std::map<float, float> pdfd_7TeV_VBF;
  std::map<float, float> pdfd_8TeV_VBF;

  std::map<int, float> evt_7TeV;
  std::map<int, float> evt_8TeV;

};

#endif
