#ifndef SCALES_H
#define SCALES_H

#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <string>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include "xsec.h"
#include "pu.h"

#include "WWAnalysis/TreeModifiers/interface/LeptSfProvider.h"


/////////////////////// Lepton efficiency scale factors /////////////////////////////////////////////

bool is2011;
bool initialized;

LeptSfProvider sfProvider;

float getMuonIDSF(float pt, float eta) {
  return sfProvider.getMuonIDSF(pt,eta);
}

float getMuonIPSF(float pt, float eta) {
  return sfProvider.getMuonIDSF(pt,eta);
}

float getMuonIsoSF(float pt, float eta) {
  return  sfProvider.getMuonIsoSF(pt,eta);
}

float getMuonSF(float pt, float eta) {
  return sfProvider.getMuonSF(pt, eta);
}

float getElectronSF(float pt, float eta) {
  return sfProvider.getElectronSF(pt,eta);
}

float getPR(float pt, float eta, float id) {
  return  sfProvider.getSF(pt,eta,id);
}

float getSF(float pt, float eta, float id) {
  return  sfProvider.getSF(pt,eta,id);
}

///////////// Cross-sections and related uncertainties /////////////////////////////

float getBkgXsec(int id) {
    return is2011 ? get7TeVBkgXsec(id) : get8TeVBkgXsec(id);
}

float getXsecggH(float mass) {
    return is2011 ? get7TeVggHXsec(mass) : get8TeVggHXsec(mass);
}

float getXsecVBF(float mass) {
    return is2011 ? get7TeVVBFXsec(mass) : get8TeVVBFXsec(mass);
}

float getXsecWHi(float mass) {
    return is2011 ? get7TeVWHiXsec(mass) : get8TeVWHiXsec(mass);
}

float getXsecZHi(float mass) {
    return is2011 ? get7TeVZHiXsec(mass) : get8TeVZHiXsec(mass);
}

float getXsecttH(float mass) {
    return is2011 ? get7TeVttHXsec(mass) : get8TeVttHXsec(mass);
}

float getXsecggHByChannel(float mass, int ch) {
    return is2011 ? get7TeVggHXsecByChannel(mass, ch) : get8TeVggHXsecByChannel(mass, ch);
}

float getXsecVBFByChannel(float mass, int ch) {
    return is2011 ? get7TeVVBFXsecByChannel(mass, ch) : get8TeVVBFXsecByChannel(mass, ch);
}

float getXsecWHiByChannel(float mass, int ch) {
    return is2011 ? get7TeVWHiXsecByChannel(mass, ch) : get8TeVWHiXsecByChannel(mass, ch);
}

float getXsecZHiByChannel(float mass, int ch) {
    return is2011 ? get7TeVZHiXsecByChannel(mass, ch) : get8TeVZHiXsecByChannel(mass, ch);
}

float getXsecttHByChannel(float mass, int ch) {
    return is2011 ? get7TeVttHXsecByChannel(mass, ch) : get8TeVttHXsecByChannel(mass, ch);
}

float getggHQCDScaleUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVggHQCDUp(mass)   : get7TeVggHQCDUp(mass);
    else    return is2011 ? get7TeVggHQCDDown(mass) : get7TeVggHQCDDown(mass);

}

float getVBFQCDScaleUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVVBFQCDUp(mass)   : get7TeVVBFQCDUp(mass);
    else    return is2011 ? get7TeVVBFQCDDown(mass) : get7TeVVBFQCDDown(mass);

}

float getWHiQCDScaleUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVWHiQCDUp(mass)   : get7TeVWHiQCDUp(mass);
    else    return is2011 ? get7TeVWHiQCDDown(mass) : get7TeVWHiQCDDown(mass);

}

float getZHiQCDScaleUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVZHiQCDUp(mass)   : get7TeVZHiQCDUp(mass);
    else    return is2011 ? get7TeVZHiQCDDown(mass) : get7TeVZHiQCDDown(mass);

}

float getttHQCDScaleUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVttHQCDUp(mass)   : get7TeVttHQCDUp(mass);
    else    return is2011 ? get7TeVttHQCDDown(mass) : get7TeVttHQCDDown(mass);

}

float getggHPDFUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVggHPDFUp(mass)   : get7TeVggHPDFUp(mass);
    else    return is2011 ? get7TeVggHPDFDown(mass) : get7TeVggHPDFDown(mass);

}

float getVBFPDFUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVVBFPDFUp(mass)   : get7TeVVBFPDFUp(mass);
    else    return is2011 ? get7TeVVBFPDFDown(mass) : get7TeVVBFPDFDown(mass);

}

float getWHiPDFUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVWHiPDFUp(mass)   : get7TeVWHiPDFUp(mass);
    else    return is2011 ? get7TeVWHiPDFDown(mass) : get7TeVWHiPDFDown(mass);

}

float getZHiPDFUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVZHiPDFUp(mass)   : get7TeVZHiPDFUp(mass);
    else    return is2011 ? get7TeVZHiPDFDown(mass) : get7TeVZHiPDFDown(mass);

}

float getttHPDFUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVttHPDFUp(mass)   : get7TeVttHPDFUp(mass);
    else    return is2011 ? get7TeVttHPDFDown(mass) : get7TeVttHPDFDown(mass);

}

float getQQZZPDFUncertainty7TeV(float m4l) {
    return 1. + 0.0035*sqrt(m4l-30.);
}

float getGGZZPDFUncertainty7TeV(float m4l) {
    return 1. + 0.0066*sqrt(m4l-10.);
}

float getQQZZQCDScaleUncertainty7TeV(float m4l) {
    return 1. + 0.01*sqrt((m4l-20.)/13.);
}

float getGGZZQCDScaleUncertainty7TeV(float m4l) {
    return 1.04 + 0.1*sqrt((m4l+40.)/40.);
}



///////////// Pileup reweighting /////////////////////////////

float getPUWeight(float numsim, int mode=1) {
    return is2011 ? getPUWeight2011(numsim) : getPUWeight2012(numsim, mode);
}

///////////// Initialize all global variables /////////////////////////////

void init(bool is7TeV) {

    if (initialized) {
        if (is7TeV == is2011) return;            
        else {
            is2011 = is7TeV;
	    sfProvider.init(is2011);
        }
    }
    else {
        initialized = true;
        is2011 = is7TeV;
	sfProvider.init(is2011);
	xsecProvider.initXsec();
	xsecProvider.initQCDScale();
	xsecProvider.initPDF();
        initpuweights();
    }

}

#endif
