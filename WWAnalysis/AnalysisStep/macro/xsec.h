#ifndef XSEC_H
#define XSEC_H

#include <map>
#include "WWAnalysis/TreeModifiers/interface/XSecProvider.h"
XSecProvider xsecProvider;

float evt_7TeV(int id){
  return xsecProvider.getEvt7TeV(id);
}

float evt_8TeV(int id){
  return xsecProvider.getEvt8TeV(id);
}

float get7TeVggHXsec(float mass) {
  return xsecProvider.get7TeVggHXsec(mass);
}    

float get8TeVggHXsec(float mass) {
  return xsecProvider.get8TeVggHXsec(mass);
}    

float get7TeVVBFXsec(float mass) {
  return xsecProvider.get7TeVVBFXsec(mass);
}    

float get8TeVVBFXsec(float mass) {
  return xsecProvider.get8TeVVBFXsec(mass);
}    

float get7TeVggHXsecByChannel(float mass, int ch) {
  return xsecProvider.get7TeVggHXsecByChannel(mass,ch);
}

float get8TeVggHXsecByChannel(float mass, int ch) {
  return xsecProvider.get8TeVggHXsecByChannel(mass,ch);
}

float get7TeVVBFXsecByChannel(float mass, int ch) {
  return xsecProvider.get7TeVVBFXsecByChannel(mass,ch);
}

float get8TeVVBFXsecByChannel(float mass, int ch) {
  return xsecProvider.get8TeVVBFXsecByChannel(mass,ch);
}

float get7TeVggHQCDUp(float mass) {
  return xsecProvider.get7TeVggHQCDUp(mass);
}    

float get7TeVggHQCDDown(float mass) {
  return xsecProvider.get7TeVggHQCDDown(mass);
}    

float get7TeVVBFQCDUp(float mass) {
  return xsecProvider.get7TeVVBFQCDUp(mass);
}   

float get7TeVVBFQCDDown(float mass) {
  return xsecProvider.get7TeVVBFQCDDown(mass);
}   

float get8TeVggHQCDUp(float mass) {
  return xsecProvider.get8TeVggHQCDUp(mass);
}    

float get8TeVggHQCDDown(float mass) {
  return xsecProvider.get8TeVggHQCDDown(mass);
}    

float get8TeVVBFQCDUp(float mass) {
  return xsecProvider.get8TeVVBFQCDUp(mass);
}   

float get8TeVVBFQCDDown(float mass) {
  return xsecProvider.get8TeVVBFQCDDown(mass);
}   

float get7TeVggHPDFUp(float mass) {
  return xsecProvider.get7TeVggHPDFUp(mass);
}   

float get7TeVggHPDFDown(float mass) {
  return xsecProvider.get7TeVggHPDFDown(mass);
}   

float get7TeVVBFPDFUp(float mass) {
  return xsecProvider.get7TeVVBFPDFUp(mass);
}   

float get7TeVVBFPDFDown(float mass) {
  return xsecProvider.get7TeVVBFPDFDown(mass);
}   

float get8TeVggHPDFUp(float mass) {
  return xsecProvider.get8TeVggHPDFUp(mass);
}

float get8TeVggHPDFDown(float mass) {
  return xsecProvider.get8TeVggHPDFDown(mass);
}

float get8TeVVBFPDFUp(float mass) {
  return xsecProvider.get8TeVVBFPDFUp(mass);
}    

float get8TeVVBFPDFDown(float mass) {
    return xsecProvider.get8TeVVBFPDFDown(mass);
}

float get7TeVBkgXsec(int id) {
  return xsecProvider.get7TeVBkgXsec(id);
}

float get8TeVBkgXsec(int id) {
  return xsecProvider.get8TeVBkgXsec(id);
}



#endif
