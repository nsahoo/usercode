#ifndef SCALES_H
#define SCALES_H

#include <vector>
#include <cmath>

float xsec_wjets       ;
float xsec_zjets       ;
float xsec_ttbar2l2nu2b;
float xsec_singletops  ;
float xsec_singletopt  ;
float xsec_singletopw  ;
float xsec_singletbars ;
float xsec_singletbart ;
float xsec_singletbarw ;
float xsec_wz3lnu      ;
float xsec_ww2l2nu     ;
float xsec_ggzz2l2l    ;
float xsec_ggzz4l      ;
float xsec_qqzz4e      ;
float xsec_qqzz4mu     ;
float xsec_qqzz4tau    ;
float xsec_qqzz2e2mu   ;
float xsec_qqzz2e2tau  ;
float xsec_qqzz2mu2tau ;
float xsec_gghzz120    ;
float xsec_vbfhzz120   ;
float xsec_vhzz120     ;

float evt_wjets        ; 
float evt_zjets        ; 
float evt_ttbar2l2nu2b ; 
float evt_singletops   ;
float evt_singletopt   ;
float evt_singletopw   ;
float evt_singletbars  ;
float evt_singletbart  ;
float evt_singletbarw  ;
float evt_wz3lnu       ; 
float evt_ww2l2nu      ;
float evt_ggzz2l2l     ;
float evt_ggzz4l       ;
float evt_qqzz4e       ;
float evt_qqzz4mu      ;
float evt_qqzz4tau     ;
float evt_qqzz2e2mu    ;
float evt_qqzz2e2tau   ;
float evt_qqzz2mu2tau  ;
float evt_ggzz4e       ;
float evt_ggzz4mu      ;
float evt_ggzz4tau     ;
float evt_gghzz120     ;
float evt_vbfhzz120    ;
float evt_vhzz120      ;

std::vector<float> puweights;

std::vector<float> fakemubar1;
std::vector<float> fakemuend1;
std::vector<float> fakemubar2;
std::vector<float> fakemuend2;

std::vector<float> fakeelbar1;
std::vector<float> fakeelend1;
std::vector<float> fakeelbar2;
std::vector<float> fakeelend2;

void initxsecscales() {
    xsec_wjets        = 31314.0 * 1000.0    ;
    xsec_zjets        = 3048.0  * 1000.0    ;
    xsec_ttbar2l2nu2b = 17.32   * 1000.0    ;
    xsec_singletops   = 3.19    * 1000.0    ;
    xsec_singletopt   = 41.92   * 1000.0    ;
    xsec_singletopw   = 7.87    * 1000.0    ;
    xsec_singletbars  = 1.44    * 1000.0    ;
    xsec_singletbart  = 22.65   * 1000.0    ;
    xsec_singletbarw  = 7.87    * 1000.0    ;
    xsec_wz3lnu       = 0.595   * 1000.0    ;
    xsec_ww2l2nu      = 4.88    * 1000.0    ;
    xsec_ggzz2l2l     = 3.48                ;
    xsec_ggzz4l       = 1.74                ;
    xsec_qqzz4e       = 12.9                ;
    xsec_qqzz4mu      = 12.9                ;
    xsec_qqzz4tau     = 12.9                ;
    xsec_qqzz2e2mu    = 24.42               ;
    xsec_qqzz2e2tau   = 24.42               ;
    xsec_qqzz2mu2tau  = 24.42               ;
    xsec_gghzz120     = 2.81                ;
    xsec_vbfhzz120    = 0.214               ;
    xsec_vhzz120      = 0.111               ;
    
    evt_wjets         = 81345381;
    evt_zjets         = 36209629;
    evt_ttbar2l2nu2b  = 10339374;
    evt_singletops    = 259971;
    evt_singletopt    = 3900171;
    evt_singletopw    = 814390;
    evt_singletbars   = 137980;
    evt_singletbart   = 1944826;
    evt_singletbarw   = 437798;
    evt_wz3lnu        = 1097759;
    evt_ww2l2nu       = 210667;
    evt_ggzz4l        = 524301;
    evt_ggzz2l2l      = 536993;
    evt_qqzz4e        = 499929;
    evt_qqzz4mu       = 499918;
    evt_qqzz4tau      = 473838;
    evt_qqzz2e2mu     = 499917;
    evt_qqzz2e2tau    = 484063;
    evt_qqzz2mu2tau   = 486266;
    evt_ggzz4e        = 174822; 
    evt_ggzz4mu       = 174602; 
    evt_ggzz4tau      = 174652; 
    evt_gghzz120      = 288879; 
    evt_vbfhzz120     = 49935 ; 
    evt_vhzz120       = 11480 ; 
}

void initpuweights() {
    puweights.push_back(0.232336075088);
    puweights.push_back(0.497862474791);
    puweights.push_back(0.747256129812);
    puweights.push_back(0.957615404034);
    puweights.push_back(1.12490413085);
    puweights.push_back(1.236084217);
    puweights.push_back(1.33301761958);
    puweights.push_back(1.41210447417);
    puweights.push_back(1.41306851511);
    puweights.push_back(1.40113107311);
    puweights.push_back(1.40143667483);
    puweights.push_back(1.36337918549);
    puweights.push_back(1.31087794774);
    puweights.push_back(1.23846099923);
    puweights.push_back(1.13781480338);
    puweights.push_back(1.03150871522);
    puweights.push_back(0.910180650156);
    puweights.push_back(0.778340121201);
    puweights.push_back(0.658396565454);
    puweights.push_back(0.546700073299);
    puweights.push_back(0.44247394088);
    puweights.push_back(0.351005629794);
    puweights.push_back(0.274366932731);
    puweights.push_back(0.211006902654);
    puweights.push_back(0.158984897593);
    puweights.push_back(0.118809864797);
    puweights.push_back(0.0874128086545);
    puweights.push_back(0.0636625788013);
    puweights.push_back(0.0460006001184);
    puweights.push_back(0.0326263540165);
    puweights.push_back(0.0228001990529);
    puweights.push_back(0.0160148490333);
    puweights.push_back(0.0111822575902);
    puweights.push_back(0.00764365686009);
    puweights.push_back(0.00523609984003);
}

void init() {

    initxsecscales();
    initpuweights();

}

float getPUWeight(int numsim) {
    if (numsim < 0 || (unsigned)numsim >= puweights.size()) return 0;
    else return puweights[(unsigned)numsim];
}

#endif
