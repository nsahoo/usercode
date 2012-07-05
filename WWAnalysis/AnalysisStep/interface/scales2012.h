#ifndef SCALES2012_H
#define SCALES2012_H

#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include <string>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TH2F.h>

std::vector<float> puweights_2012;
TH1F* hist_puweights_2012;
std::map<int, float> xsecweights;

Double_t* IP_eta_2012; 
Double_t* IP_pt_barrel_2012; 
Double_t* IP_pt_endcaps_2012; 

Double_t* ISO_eta_2012; 
Double_t* ISO_pt_barrel_2012; 
Double_t* ISO_pt_endcaps_2012; 

Double_t* ID_eta_2012; 
Double_t* ID_pt_barrel_2012; 
Double_t* ID_pt_endcaps_2012; 

TH2* electronscalefactors_2012;
std::vector<std::vector<float> > elesf_2012;

void initxsecweights() {
    xsecweights[100] = 2.24  / 550862.0;
    xsecweights[101] = 4.47  / 400973.0;
    xsecweights[102] = 76.91 / 1489093.0;
    xsecweights[103] = 76.91 / 1476065.0;
    xsecweights[104] = 76.91 / 824466.0;
    xsecweights[105] = 176.7 / 1617951.0;
    xsecweights[106] = 176.7 / 823911.0;
    xsecweights[107] = 176.7 / 823922.0;
    xsecweights[111] = 3503700.0 / 17782545.0;
    xsecweights[112] = 23640.0 / 13948878.0;
    xsecweights[113] = 3890.0 / 497658.0;
    xsecweights[114] = 55530.0 / 493460.0;
    xsecweights[115] = 11770.0 / 23777.0;
    xsecweights[116] = 1760.0 / 1925072.0;
    xsecweights[117] = 30000.0 / 259961.0;
    xsecweights[118] = 117700.0 / 119974.0;
    xsecweights[119] = 1057.0 / 2170000.0;
    xsecweights[120] = 5995.0 / 1753235.0;

    xsecweights[1115] = 22.96  * 1000.0 * 9.22e-5 / 275378.; 
    xsecweights[1117] = 22.20  * 1000.0 * 1.18e-4 / 299982.;
    xsecweights[1119] = 21.48  * 1000.0 * 1.50e-4 / 274180.;
    xsecweights[1120] = 21.13  * 1000.0 * 1.67e-4 / 999921.;
    xsecweights[1121] = 20.80  * 1000.0 * 1.86e-4 / 297566.;
    xsecweights[1123] = 20.15  * 1000.0 * 2.29e-4 / 299977.;
    xsecweights[1124] = 19.83  * 1000.0 * 2.52e-4 / 299982.;
    xsecweights[1125] = 19.52  * 1000.0 * 2.76e-4 / 995117.;
    xsecweights[1126] = 19.22  * 1000.0 * 3.02e-4 / 299974.;
    xsecweights[1127] = 18.92  * 1000.0 * 3.28e-4 / 299971.;
    xsecweights[1145] = 14.59  * 1000.0 * 8.14e-4 / 299971.;
    xsecweights[1150] = 13.65  * 1000.0 * 8.50e-4 / 999314.;
    xsecweights[1180] = 8.874  * 1000.0 * 6.12e-4 / 297273.;
    xsecweights[1200] = 7.127  * 1000.0 * 2.59e-3 / 894023.;
    xsecweights[1220] = 6.038  * 1000.0 * 2.87e-3 / 299968.;

    xsecweights[2115] = 1.709  * 1000.0 * 9.22e-5 / 49937.;
    xsecweights[2117] = 1.678  * 1000.0 * 1.18e-4 / 49958.;
    xsecweights[2119] = 1.647  * 1000.0 * 1.50e-4 / 49940.;
    xsecweights[2120] = 1.632  * 1000.0 * 1.67e-4 / 49936.;
    xsecweights[2121] = 1.617  * 1000.0 * 1.86e-4 / 49941.;
    xsecweights[2123] = 1.588  * 1000.0 * 2.29e-4 / 49936.;
    xsecweights[2124] = 1.573  * 1000.0 * 2.52e-4 / 49941.;
    xsecweights[2125] = 1.559  * 1000.0 * 2.76e-4 / 49939.;
    xsecweights[2126] = 1.544  * 1000.0 * 3.02e-4 / 49945.;
    xsecweights[2127] = 1.531  * 1000.0 * 3.28e-4 / 49953.;
    xsecweights[2130] = 1.490  * 1000.0 * 4.14e-4 / 49946.;
    xsecweights[2150] = 1.251  * 1000.0 * 8.50e-4 / 49945.;
    xsecweights[2180] = 0.9820 * 1000.0 * 6.12e-4 / 49951.;
    xsecweights[2200] = 0.8441 * 1000.0 * 2.59e-3 / 49946.;
    xsecweights[2220] = 0.7229 * 1000.0 * 2.87e-3 / 49945.;

}

void initpuweights() {
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(0.222451);
    puweights_2012.push_back(0.0658851);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(0.150902);
    puweights_2012.push_back(0.202205);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1.50116);
    puweights_2012.push_back(2.79375);
    puweights_2012.push_back(0.198341);
    puweights_2012.push_back(0.246893);
    puweights_2012.push_back(0.28116);
    puweights_2012.push_back(0.449377);
    puweights_2012.push_back(0.553276);
    puweights_2012.push_back(1.48919);
    puweights_2012.push_back(2.15249);
    puweights_2012.push_back(3.62415);
    puweights_2012.push_back(4.33041);
    puweights_2012.push_back(3.57192);
    puweights_2012.push_back(4.99603);
    puweights_2012.push_back(7.79303);
    puweights_2012.push_back(8.04276);
    puweights_2012.push_back(8.05557);
    puweights_2012.push_back(12.9364);
    puweights_2012.push_back(9.9036);
    puweights_2012.push_back(14.6975);
    puweights_2012.push_back(13.3387);
    puweights_2012.push_back(10.9734);
    puweights_2012.push_back(12.6077);
    puweights_2012.push_back(11.5617);
    puweights_2012.push_back(10.8107);
    puweights_2012.push_back(14.5043);
    puweights_2012.push_back(17.8497);
    puweights_2012.push_back(11.8817);
    puweights_2012.push_back(9.6805);
    puweights_2012.push_back(12.2255);
    puweights_2012.push_back(10.1117);
    puweights_2012.push_back(10.2482);
    puweights_2012.push_back(11.5398);
    puweights_2012.push_back(9.35737);
    puweights_2012.push_back(9.90259);
    puweights_2012.push_back(9.19216);
    puweights_2012.push_back(7.57377);
    puweights_2012.push_back(7.94847);
    puweights_2012.push_back(7.15578);
    puweights_2012.push_back(5.63016);
    puweights_2012.push_back(5.35972);
    puweights_2012.push_back(5.05791);
    puweights_2012.push_back(3.35313);
    puweights_2012.push_back(3.60582);
    puweights_2012.push_back(3.35256);
    puweights_2012.push_back(2.49496);
    puweights_2012.push_back(2.28219);
    puweights_2012.push_back(2.21227);
    puweights_2012.push_back(1.76362);
    puweights_2012.push_back(1.68533);
    puweights_2012.push_back(1.62149);
    puweights_2012.push_back(1.34263);
    puweights_2012.push_back(1.30646);
    puweights_2012.push_back(1.21918);
    puweights_2012.push_back(1.10347);
    puweights_2012.push_back(1.08544);
    puweights_2012.push_back(1.0251);
    puweights_2012.push_back(0.907123);
    puweights_2012.push_back(0.905997);
    puweights_2012.push_back(0.869217);
    puweights_2012.push_back(0.816708);
    puweights_2012.push_back(0.76043);
    puweights_2012.push_back(0.714367);
    puweights_2012.push_back(0.679723);
    puweights_2012.push_back(0.665294);
    puweights_2012.push_back(0.609956);
    puweights_2012.push_back(0.586386);
    puweights_2012.push_back(0.548999);
    puweights_2012.push_back(0.521088);
    puweights_2012.push_back(0.4929);
    puweights_2012.push_back(0.453545);
    puweights_2012.push_back(0.44546);
    puweights_2012.push_back(0.406266);
    puweights_2012.push_back(0.378486);
    puweights_2012.push_back(0.347898);
    puweights_2012.push_back(0.337097);
    puweights_2012.push_back(0.313674);
    puweights_2012.push_back(0.291392);
    puweights_2012.push_back(0.283346);
    puweights_2012.push_back(0.25272);
    puweights_2012.push_back(0.244178);
    puweights_2012.push_back(0.228673);
    puweights_2012.push_back(0.211327);
    puweights_2012.push_back(0.19084);
    puweights_2012.push_back(0.179408);
    puweights_2012.push_back(0.169234);
    puweights_2012.push_back(0.157131);
    puweights_2012.push_back(0.143818);
    puweights_2012.push_back(0.140968);
    puweights_2012.push_back(0.124021);
    puweights_2012.push_back(0.118273);
    puweights_2012.push_back(0.109751);
    puweights_2012.push_back(0.0977754);
    puweights_2012.push_back(0.0967206);
    puweights_2012.push_back(0.0870401);
    puweights_2012.push_back(0.0826372);
    puweights_2012.push_back(0.0746777);
    puweights_2012.push_back(0.0698592);
    puweights_2012.push_back(0.0656062);
    puweights_2012.push_back(0.0601853);
    puweights_2012.push_back(0.057892);
    puweights_2012.push_back(0.0517871);
    puweights_2012.push_back(0.0512109);
    puweights_2012.push_back(0.0465423);
    puweights_2012.push_back(0.0403982);
    puweights_2012.push_back(0.0443631);
    puweights_2012.push_back(0.0399185);
    puweights_2012.push_back(0.0338933);
    puweights_2012.push_back(0.0354274);
    puweights_2012.push_back(0.0310775);
    puweights_2012.push_back(0.0261122);
    puweights_2012.push_back(0.0280878);
    puweights_2012.push_back(0.0264014);
    puweights_2012.push_back(0.021299);
    puweights_2012.push_back(0.0245197);
    puweights_2012.push_back(0.0221076);
    puweights_2012.push_back(0.0189236);
    puweights_2012.push_back(0.0202148);
    puweights_2012.push_back(0.0177248);
    puweights_2012.push_back(0.0163634);
    puweights_2012.push_back(0.0188307);
    puweights_2012.push_back(0.0144512);
    puweights_2012.push_back(0.0134599);
    puweights_2012.push_back(0.0143315);
    puweights_2012.push_back(0.0130668);
    puweights_2012.push_back(0.0108666);
    puweights_2012.push_back(0.0162516);
    puweights_2012.push_back(0.0126035);
    puweights_2012.push_back(0.0102154);
    puweights_2012.push_back(0.0154442);
    puweights_2012.push_back(0.00959973);
    puweights_2012.push_back(0.0106827);
    puweights_2012.push_back(0.0146624);
    puweights_2012.push_back(0.0155156);
    puweights_2012.push_back(0.00761674);
    puweights_2012.push_back(0.0187999);
    puweights_2012.push_back(0.0135013);
    puweights_2012.push_back(0.0160794);
    puweights_2012.push_back(0.0180586);
    puweights_2012.push_back(0.00905508);
    puweights_2012.push_back(0.00869858);
    puweights_2012.push_back(0.0193968);
    puweights_2012.push_back(0.0209201);
    puweights_2012.push_back(0.0084405);
    puweights_2012.push_back(0.0407657);
    puweights_2012.push_back(0.0109116);
    puweights_2012.push_back(0.0262218);
    puweights_2012.push_back(0.0104767);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(0.00658031);
    puweights_2012.push_back(0.0051814);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);
    puweights_2012.push_back(1);

    hist_puweights_2012 = new TH1F("hist_puweights_2012","",240,0.,80.);

    for(int k=0;k<240;k++){
        hist_puweights_2012->SetBinContent(k+1,puweights_2012[k]);
    }
}

void initmuonscalefactors() {

    std::string baseFolder(getenv("CMSSW_BASE"));
    std::string filepath = baseFolder + "/src/WWAnalysis/AnalysisStep/data/scale_factors.26-05-12.root";

    TFile file(filepath.c_str());

    IP_eta_2012        = ((TGraphAsymmErrors*)file.Get("IP_eta_2012"))->GetY();
    IP_pt_barrel_2012  = ((TGraphAsymmErrors*)file.Get("IP_pt_barrel_2012"))->GetY();
    IP_pt_endcaps_2012 = ((TGraphAsymmErrors*)file.Get("IP_pt_endcaps_2012"))->GetY();
    
    ISO_eta_2012       = ((TGraphAsymmErrors*)file.Get("ISO_eta_2012"))->GetY();
    ISO_pt_barrel_2012 = ((TGraphAsymmErrors*)file.Get("ISO_pt_barrel_2012"))->GetY();
    ISO_pt_endcaps_2012= ((TGraphAsymmErrors*)file.Get("ISO_pt_endcaps_2012"))->GetY();
    
    ID_eta_2012        = ((TGraphAsymmErrors*)file.Get("ID_eta_2012"))->GetY();
    ID_pt_barrel_2012  = ((TGraphAsymmErrors*)file.Get("ID_pt_barrel_2012"))->GetY();
    ID_pt_endcaps_2012 = ((TGraphAsymmErrors*)file.Get("ID_pt_endcaps_2012"))->GetY();

        
}

void initelectronscalefactors() {
    std::string baseFolder(getenv("CMSSW_BASE"));
    std::string filepath = baseFolder + "/src/WWAnalysis/AnalysisStep/data/hel_sf_2012.root";

    TFile file(filepath.c_str());

    electronscalefactors_2012 = (TH2*)file.Get("heff");

    for (std::size_t i = 0; i < 7; i++) {
        elesf_2012.push_back(std::vector<float>(5, 0.0));
    }
    for (std::size_t i = 0; i < 7; i++) {
        for (std::size_t j = 0; j < 5; j++) {
            elesf_2012[i][j] = electronscalefactors_2012->GetBinContent(i+1, j+1);
        }
    }
}


void init() {

    initxsecweights();
    initpuweights();
    initmuonscalefactors();
    initelectronscalefactors();

}

float getPUWeight(float numsim) {
    return hist_puweights_2012->GetBinContent(hist_puweights_2012->FindBin(numsim));
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

        if (bin >= 0) return ID_eta_2012[bin];

    }

    else {
                
        if (fabs(eta) <= 1.2) {
            if (pt > 5.0 && pt <= 7.5) bin = 0;
            if (pt > 7.5 && pt <= 10.) bin = 1;
            if (pt > 10. && pt <= 15.) bin = 2;
            if (pt > 15. && pt <= 20.) bin = 3;

            if (bin >= 0) return ID_pt_barrel_2012[bin];  
        }

        else {
            if (pt > 5.0 && pt <= 7.5) bin = 0;
            if (pt > 7.5 && pt <= 10.) bin = 1;
            if (pt > 10. && pt <= 15.) bin = 2;
            if (pt > 15. && pt <= 20.) bin = 3;

            if (bin >= 0) return ID_pt_endcaps_2012[bin];  
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

        if (bin >= 0) return IP_eta_2012[bin];

    }

    else {
                     
        if (fabs(eta) <= 1.2) {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return IP_pt_barrel_2012[bin];  
        }

        else {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return IP_pt_endcaps_2012[bin];
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

        if (bin >= 0) return ISO_eta_2012[bin];

    }

    else {

        if (fabs(eta) <= 1.2) {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return ISO_pt_barrel_2012[bin];                               
        }

        else {
            if (pt > 5.0 && pt <= 10.) bin = 0;
            if (pt > 10. && pt <= 15.) bin = 1;
            if (pt > 15. && pt <= 20.) bin = 2;

            if (bin >= 0) return ISO_pt_endcaps_2012[bin];
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

    if (ptbin>=0 && etabin>=0) return elesf_2012[ptbin][etabin];
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
