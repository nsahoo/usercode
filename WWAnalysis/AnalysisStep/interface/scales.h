#ifndef SCALES_H
#define SCALES_H

#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <string>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>

std::vector<float> puweights_2011;
std::map<int, float> xsecweights;

Double_t* IP_eta_2011; 
Double_t* IP_pt_barrel_2011; 
Double_t* IP_pt_endcaps_2011; 

Double_t* ISO_eta_2011; 
Double_t* ISO_pt_barrel_2011; 
Double_t* ISO_pt_endcaps_2011; 

Double_t* ID_eta_2011A; 
Double_t* ID_eta_2011B; 
Double_t* ID_pt_barrel_2011; 
Double_t* ID_pt_endcaps_2011; 

TH2* electronscalefactors_2011;
std::vector<std::vector<float> > elesf_2011;

void initxsecweights() {
    xsecweights[100] = 1.74  / 524301.0;
    xsecweights[101] = 3.48  / 486993.0;
    xsecweights[102] = 15.34 / 499924.0;
    xsecweights[103] = 15.34 / 377430.0;
    xsecweights[104] = 15.34 / 499026.0;
    xsecweights[105] = 30.67 / 224340.0;
    xsecweights[106] = 30.67 / 493618.0;
    xsecweights[107] = 30.67 / 499927.0;
    xsecweights[111] = 3048000.0 / 36209629.0;
    xsecweights[112] = 17320.0 / 10189374.0;
    xsecweights[113] = 3190.0 / 259971.0;
    xsecweights[114] = 41920.0 / 3900171.0;
    xsecweights[115] = 7870.0 / 814390.0;
    xsecweights[116] = 1440.0 / 137980.0;
    xsecweights[117] = 22650.0 / 1944826.0;
    xsecweights[118] = 7870.0 / 323401.0;
    xsecweights[119] = 868.0 / 997759.0;
    xsecweights[120] = 4880.0 / 210667.0;

    xsecweights[121] = 0.1525  * 1000.0 * 0.3283 / 410728.;
    xsecweights[122] = 0.0659  * 1000.0 * 0.4131 / 397435.;
    xsecweights[123] = 0.1523  * 1000.0 * 0.3413 / 318885.;
    xsecweights[124] = 0.0664  * 1000.0 * 0.3764 / 478033.;
    xsecweights[125] = 0.1517  * 1000.0 * 0.3410 / 418014.;
    xsecweights[126] = 0.0661  * 1000.0 * 0.3765 / 479652.;
    
    xsecweights[200] = 18.13  * 1000.0 * 9.22e-5 / 296850.; 
    xsecweights[201] = 16.63  * 1000.0 * 1.67e-4 / 238879.;
    xsecweights[202] = 14.12  * 1000.0 * 4.14e-4 / 295356.;
    xsecweights[203] = 12.13  * 1000.0 * 7.11e-4 / 299973.;
    xsecweights[204] = 10.50  * 1000.0 * 8.50e-4 / 285262.;
    xsecweights[205] = 9.080  * 1000.0 * 4.25e-4 / 237830.;

    xsecweights[250] = 1.332  * 1000.0 * 9.22e-5 / 49939.;
    xsecweights[251] = 1.269  * 1000.0 * 1.67e-4 / 49935.;
    xsecweights[252] = 1.154  * 1000.0 * 4.14e-4 / 49948.;
    xsecweights[253] = 1.052  * 1000.0 * 7.11e-4 / 49950.;
    xsecweights[254] = 0.9617 * 1000.0 * 8.50e-4 / 45443.;
    xsecweights[255] = 0.8787 * 1000.0 * 4.25e-4 / 49955.;

}

void initpuweights() {
    puweights_2011.push_back(0.212929);
    puweights_2011.push_back(0.0208114);
    puweights_2011.push_back(0.0584048);
    puweights_2011.push_back(0.538898);
    puweights_2011.push_back(1.357);
    puweights_2011.push_back(1.49913);
    puweights_2011.push_back(1.42247);
    puweights_2011.push_back(1.35904);
    puweights_2011.push_back(1.29946);
    puweights_2011.push_back(1.27925);
    puweights_2011.push_back(1.37845);
    puweights_2011.push_back(1.71246);
    puweights_2011.push_back(1.5291);
    puweights_2011.push_back(1.35234);
    puweights_2011.push_back(1.22215);
    puweights_2011.push_back(1.0155);
    puweights_2011.push_back(1.01137);
    puweights_2011.push_back(0.395465);
    puweights_2011.push_back(0.230984);
    puweights_2011.push_back(0.109883);
    puweights_2011.push_back(0.0433739);
    puweights_2011.push_back(0.0111497);
    puweights_2011.push_back(0.00408801);
    puweights_2011.push_back(0.00115678);
    puweights_2011.push_back(0.000365505);
    puweights_2011.push_back(0.000112391);
    puweights_2011.push_back(3.83894e-05);
    puweights_2011.push_back(1.60651e-05);
    puweights_2011.push_back(4.81412e-06);
    puweights_2011.push_back(1.39717e-06);
    puweights_2011.push_back(1.92368e-06);
    puweights_2011.push_back(4.10748e-06);
    puweights_2011.push_back(2.33157e-05);
    puweights_2011.push_back(4.0181e-05);
    puweights_2011.push_back(4.87786e-05);
    puweights_2011.push_back(0.00194128);
    puweights_2011.push_back(8.97414e-05);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(0.000162709);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);
    puweights_2011.push_back(1);

}

void initmuonscalefactors() {

    std::string baseFolder(getenv("CMSSW_BASE"));
    std::string filepath = baseFolder + "/src/WWAnalysis/AnalysisStep/data/scale_factors.26-05-12.root";

    TFile file(filepath.c_str());

    IP_eta_2011        = ((TGraphAsymmErrors*)file.Get("IP_eta_2011"))->GetY();
    IP_pt_barrel_2011  = ((TGraphAsymmErrors*)file.Get("IP_pt_barrel_2011"))->GetY();
    IP_pt_endcaps_2011 = ((TGraphAsymmErrors*)file.Get("IP_pt_endcaps_2011"))->GetY();
    
    ISO_eta_2011       = ((TGraphAsymmErrors*)file.Get("ISO_eta_2011"))->GetY();
    ISO_pt_barrel_2011 = ((TGraphAsymmErrors*)file.Get("ISO_pt_barrel_2011"))->GetY();
    ISO_pt_endcaps_2011= ((TGraphAsymmErrors*)file.Get("ISO_pt_endcaps_2011"))->GetY();
    
    ID_eta_2011A       = ((TGraphAsymmErrors*)file.Get("ID_eta_2011A"))->GetY();
    ID_eta_2011B       = ((TGraphAsymmErrors*)file.Get("ID_eta_2011B"))->GetY();
    ID_pt_barrel_2011  = ((TGraphAsymmErrors*)file.Get("ID_pt_barrel_2011"))->GetY();
    ID_pt_endcaps_2011 = ((TGraphAsymmErrors*)file.Get("ID_pt_endcaps_2011"))->GetY();

        
}

void initelectronscalefactors() {
    std::string baseFolder(getenv("CMSSW_BASE"));
    std::string filepath = baseFolder + "/src/WWAnalysis/AnalysisStep/data/hel_sf_2011.root";

    TFile file(filepath.c_str());

    electronscalefactors_2011 = (TH2*)file.Get("heff");

    for (std::size_t i = 0; i < 7; i++) {
        elesf_2011.push_back(std::vector<float>(5, 0.0));
    }
    for (std::size_t i = 0; i < 7; i++) {
        for (std::size_t j = 0; j < 5; j++) {
            elesf_2011[i][j] = electronscalefactors_2011->GetBinContent(i+1, j+1);
        }
    }
}


void init() {

    initxsecweights();
    initpuweights();
    initmuonscalefactors();
    initelectronscalefactors();

}

float getPUWeight(int numsim) {
    if (numsim < 0 || (unsigned)numsim >= puweights_2011.size()) return 0;
    else return puweights_2011[(unsigned)numsim];
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

        if (bin >= 0) return 0.42*ID_eta_2011A[bin] + 0.58*ID_eta_2011B[bin];

    }

    else {
                
        if (fabs(eta) <= 1.2) {
            if (pt > 5.0 && pt <= 7.5) bin = 0;
            if (pt > 7.5 && pt <= 10.) bin = 1;
            if (pt > 10. && pt <= 15.) bin = 2;
            if (pt > 15. && pt <= 20.) bin = 3;

            if (bin >= 0) return ID_pt_barrel_2011[bin];  
        }

        else {
            if (pt > 5.0 && pt <= 7.5) bin = 0;
            if (pt > 7.5 && pt <= 10.) bin = 1;
            if (pt > 10. && pt <= 15.) bin = 2;
            if (pt > 15. && pt <= 20.) bin = 3;

            if (bin >= 0) return ID_pt_endcaps_2011[bin];  
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

        if (bin >= 0) return IP_eta_2011[bin];

    }

    else {
                     
        if (fabs(eta) <= 1.2) {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return IP_pt_barrel_2011[bin];  
        }

        else {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return IP_pt_endcaps_2011[bin];
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

        if (bin >= 0) return ISO_eta_2011[bin];

    }

    else {

        if (fabs(eta) <= 1.2) {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return ISO_pt_barrel_2011[bin];                               
        }

        else {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return ISO_pt_endcaps_2011[bin];
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

    if (ptbin>=0 && etabin>=0) return elesf_2011[ptbin][etabin];
    else return 0.0;
}

float getPR(float pt, float eta, float id) {

    return 1.0;

}

float getSF(float pt, float eta, float id) {

    if (id == 13 || id == -13) return getMuonSF(pt, eta);

    else return getElectronSF(pt, eta);

    //return 1.0;

}

#endif
