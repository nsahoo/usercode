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

/////////////////////// Lepton efficiency scale factors /////////////////////////////////////////////

bool is2011;
bool initialized;

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

void initmuonscalefactors() {

    std::string baseFolder(getenv("CMSSW_BASE"));
    std::string filepath = baseFolder + "/src/WWAnalysis/AnalysisStep/data/mu_sf.root";

    TFile file(filepath.c_str());

    if (is2011) {
        IP_eta        = ((TGraphAsymmErrors*)file.Get("IP_eta_2011"))->GetY();
        IP_pt_barrel  = ((TGraphAsymmErrors*)file.Get("IP_pt_barrel_2011"))->GetY();
        IP_pt_endcaps = ((TGraphAsymmErrors*)file.Get("IP_pt_endcaps_2011"))->GetY();
        
        ISO_eta       = ((TGraphAsymmErrors*)file.Get("ISO_eta_2011"))->GetY();
        ISO_pt_barrel = ((TGraphAsymmErrors*)file.Get("ISO_pt_barrel_2011"))->GetY();
        ISO_pt_endcaps= ((TGraphAsymmErrors*)file.Get("ISO_pt_endcaps_2011"))->GetY();
        
        ID_eta_2011A  = ((TGraphAsymmErrors*)file.Get("ID_eta_2011A"))->GetY();
        ID_eta_2011B  = ((TGraphAsymmErrors*)file.Get("ID_eta_2011B"))->GetY();
        ID_pt_barrel  = ((TGraphAsymmErrors*)file.Get("ID_pt_barrel_2011"))->GetY();
        ID_pt_endcaps = ((TGraphAsymmErrors*)file.Get("ID_pt_endcaps_2011"))->GetY();
    }

    else {
        IP_eta        = ((TGraphAsymmErrors*)file.Get("IP_eta_2012"))->GetY();
        IP_pt_barrel  = ((TGraphAsymmErrors*)file.Get("IP_pt_barrel_2012"))->GetY();
        IP_pt_endcaps = ((TGraphAsymmErrors*)file.Get("IP_pt_endcaps_2012"))->GetY();

        ISO_eta       = ((TGraphAsymmErrors*)file.Get("ISO_eta_2012"))->GetY();
        ISO_pt_barrel = ((TGraphAsymmErrors*)file.Get("ISO_pt_barrel_2012"))->GetY();
        ISO_pt_endcaps= ((TGraphAsymmErrors*)file.Get("ISO_pt_endcaps_2012"))->GetY();

        ID_eta        = ((TGraphAsymmErrors*)file.Get("ID_eta_2012"))->GetY();
        ID_pt_barrel  = ((TGraphAsymmErrors*)file.Get("ID_pt_barrel_2012"))->GetY();
        ID_pt_endcaps = ((TGraphAsymmErrors*)file.Get("ID_pt_endcaps_2012"))->GetY();
    }
        
}

void initelectronscalefactors() {
    std::string baseFolder(getenv("CMSSW_BASE"));
    std::string filepath = baseFolder + (is2011 ? "/src/WWAnalysis/AnalysisStep/data/hel_sf_2011.root" : "/src/WWAnalysis/AnalysisStep/data/hel_sf_2012.root");

    TFile file(filepath.c_str());

    electronscalefactors = (TH2*)file.Get("heff");

    if (elesf.size() == 0) {
        for (std::size_t i = 0; i < 7; i++) {
            elesf.push_back(std::vector<float>(5, 0.0));
        }
    }
    for (std::size_t i = 0; i < 7; i++) {
        for (std::size_t j = 0; j < 5; j++) {
            elesf[i][j] = electronscalefactors->GetBinContent(i+1, j+1);
        }
    }
}


float getMuonIDSF(float pt, float eta) {

    int bin = -1;

    if (pt > 20) {

        if (eta>-2.4 && eta<=-2.1) bin = 0;
        if (eta>-2.1 && eta<=-1.6) bin = 1;
        if (eta>-1.6 && eta<=-1.2) bin = 2;
        if (eta>-1.2 && eta<=-0.9) bin = 3;
        if (eta>-0.9 && eta<=-0.6) bin = 4;
        if (eta>-0.6 && eta<=-0.3) bin = 5;
        if (eta>-0.3 && eta<=-0.2) bin = 6;
        if (eta>-0.2 && eta<= 0.2) bin = 7;
        if (eta> 0.2 && eta<= 0.3) bin = 8;
        if (eta> 0.3 && eta<= 0.6) bin = 9;
        if (eta> 0.6 && eta<= 0.9) bin = 10;
        if (eta> 0.9 && eta<= 1.2) bin = 11;
        if (eta> 1.2 && eta<= 1.6) bin = 12;
        if (eta> 1.6 && eta<= 2.1) bin = 13;
        if (eta> 2.1 && eta<= 2.4) bin = 14;

        if (bin >= 0) {
            if (is2011) return 0.42*ID_eta_2011A[bin] + 0.58*ID_eta_2011B[bin];
            else return ID_eta[bin];
        }

    }

    else {
                
        if (fabs(eta) <= 1.2) {
            if (pt > 5.0 && pt <= 7.5) bin = 0;
            if (pt > 7.5 && pt <= 10.) bin = 1;
            if (pt > 10. && pt <= 15.) bin = 2;
            if (pt > 15. && pt <= 20.) bin = 3;

            if (bin >= 0) return ID_pt_barrel[bin];  
        }

        else {
            if (pt > 5.0 && pt <= 7.5) bin = 0;
            if (pt > 7.5 && pt <= 10.) bin = 1;
            if (pt > 10. && pt <= 15.) bin = 2;
            if (pt > 15. && pt <= 20.) bin = 3;

            if (bin >= 0) return ID_pt_endcaps[bin];  
        }

    }

    return 0.0;
}

float getMuonIPSF(float pt, float eta) {

    int bin = -1;

    if (pt > 20) {

        if (eta>-2.4 && eta<=-2.1) bin = 0;
        if (eta>-2.1 && eta<=-1.6) bin = 1;
        if (eta>-1.6 && eta<=-1.2) bin = 2;
        if (eta>-1.2 && eta<=-0.9) bin = 3;
        if (eta>-0.9 && eta<=-0.6) bin = 4;
        if (eta>-0.6 && eta<=-0.3) bin = 5;
        if (eta>-0.3 && eta<=-0.2) bin = 6;
        if (eta>-0.2 && eta<= 0.2) bin = 7;
        if (eta> 0.2 && eta<= 0.3) bin = 8;
        if (eta> 0.3 && eta<= 0.6) bin = 9;
        if (eta> 0.6 && eta<= 0.9) bin = 10;
        if (eta> 0.9 && eta<= 1.2) bin = 11;
        if (eta> 1.2 && eta<= 1.6) bin = 12;
        if (eta> 1.6 && eta<= 2.1) bin = 13;
        if (eta> 2.1 && eta<= 2.4) bin = 14;

        if (bin >= 0) return IP_eta[bin];

    }

    else {
                     
        if (fabs(eta) <= 1.2) {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return IP_pt_barrel[bin];  
        }

        else {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return IP_pt_endcaps[bin];
        }
        
    }

    return 0.0;

}

float getMuonIsoSF(float pt, float eta) {

    int bin = -1;

    if (pt > 20) {

        if (eta>-2.4 && eta<=-2.1) bin = 0;
        if (eta>-2.1 && eta<=-1.6) bin = 1;
        if (eta>-1.6 && eta<=-1.2) bin = 2;
        if (eta>-1.2 && eta<=-0.9) bin = 3;
        if (eta>-0.9 && eta<=-0.6) bin = 4;
        if (eta>-0.6 && eta<=-0.3) bin = 5;
        if (eta>-0.3 && eta<=-0.2) bin = 6;
        if (eta>-0.2 && eta<= 0.2) bin = 7;
        if (eta> 0.2 && eta<= 0.3) bin = 8;
        if (eta> 0.3 && eta<= 0.6) bin = 9;
        if (eta> 0.6 && eta<= 0.9) bin = 10;
        if (eta> 0.9 && eta<= 1.2) bin = 11;
        if (eta> 1.2 && eta<= 1.6) bin = 12;
        if (eta> 1.6 && eta<= 2.1) bin = 13;
        if (eta> 2.1 && eta<= 2.4) bin = 14;

        if (bin >= 0) return ISO_eta[bin];

    }

    else {

        if (fabs(eta) <= 1.2) {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return ISO_pt_barrel[bin];                               
        }

        else {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return ISO_pt_endcaps[bin];
        }

    }

    return 0.0;
}

float getMuonSF(float pt, float eta) {
    return getMuonIDSF(pt, eta) * getMuonIsoSF(pt, eta) * getMuonIPSF(pt, eta);
}

float getElectronSF(float pt, float eta) {
    int etabin = -1;
    int ptbin = -1;

    float abseta = fabs(eta);

    if (pt> 7 && pt<=10) ptbin = 0;
    if (pt>10 && pt<=15) ptbin = 1;
    if (pt>15 && pt<=20) ptbin = 2;
    if (pt>20 && pt<=30) ptbin = 3;
    if (pt>30 && pt<=40) ptbin = 4;
    if (pt>40 && pt<=50) ptbin = 5;
    if (pt>50 && pt<=80) ptbin = 6;
    if (pt>80)           ptbin = 6;

    if (abseta<=0.8)                     etabin = 0;
    if (abseta>0.8    && abseta<=1.4442) etabin = 1;
    if (abseta>1.4442 && abseta<=1.566)  etabin = 2;
    if (abseta>1.566  && abseta<=2.0)    etabin = 3;
    if (abseta>2.0    && abseta<=2.5)    etabin = 4;

    if (ptbin>=0 && etabin>=0) return elesf[ptbin][etabin];
    else return 0.0;
}

float getPR(float pt, float eta, float id) {

    if (id == 13 || id == -13) return getMuonIDSF(pt, eta) * getMuonIsoSF(pt, eta);
    else return getElectronSF(pt, eta);

}

float getSF(float pt, float eta, float id) {

    if (id == 13 || id == -13) return getMuonSF(pt, eta);

    else return getElectronSF(pt, eta);

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

float getggHQCDScaleUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVggHQCDUp(mass)   : get7TeVggHQCDUp(mass);
    else    return is2011 ? get7TeVggHQCDDown(mass) : get7TeVggHQCDDown(mass);

}

float getVBFQCDScaleUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVVBFQCDUp(mass)   : get7TeVVBFQCDUp(mass);
    else    return is2011 ? get7TeVVBFQCDDown(mass) : get7TeVVBFQCDDown(mass);

}


float getggHPDFUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVggHPDFUp(mass)   : get7TeVggHPDFUp(mass);
    else    return is2011 ? get7TeVggHPDFDown(mass) : get7TeVggHPDFDown(mass);

}

float getVBFPDFUncertainty(float mass, bool up) {

    if (up) return is2011 ? get7TeVVBFPDFUp(mass)   : get7TeVVBFPDFUp(mass);
    else    return is2011 ? get7TeVVBFPDFDown(mass) : get7TeVVBFPDFDown(mass);

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

float getPUWeight(int numsim) {
    return is2011 ? getPUWeight2011(numsim) : getPUWeight2012(numsim);
}

///////////// Initialize all global variables /////////////////////////////

void init(bool is7TeV) {

    if (initialized) {
        if (is7TeV == is2011) return;            
        else {
            is2011 = is7TeV;
            initmuonscalefactors();
            initelectronscalefactors();
        }
    }
    else {
        initialized = true;
        is2011 = is7TeV;
        initXsec();
        initQCDScale();
        initPDF();
        initpuweights();
        initmuonscalefactors();
        initelectronscalefactors();
    }

}

#endif
