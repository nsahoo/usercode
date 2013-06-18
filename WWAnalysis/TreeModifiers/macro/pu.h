#ifndef PU_H
#define PU_H

#include <vector>
#include <TH1F.h>


TH1F* hist_puweights_2011;
TH1F* hist_puweights_2012;
TH1F* hist_puweights_2012_52;

std::vector<float> puweights_2011;
std::vector<float> puweights_2012;
std::vector<float> puweights_2012_52;

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

    /*
    puweights_2012.push_back(0.409409);
    puweights_2012.push_back(0.527276);
    puweights_2012.push_back(0.39328);
    puweights_2012.push_back(0.507892);
    puweights_2012.push_back(0.48029);
    puweights_2012.push_back(0.787701);
    puweights_2012.push_back(0.632356);
    puweights_2012.push_back(0.618033);
    puweights_2012.push_back(0.806089);
    puweights_2012.push_back(1.14018);
    puweights_2012.push_back(1.5788);
    puweights_2012.push_back(1.93507);
    puweights_2012.push_back(1.957);
    puweights_2012.push_back(1.73004);
    puweights_2012.push_back(1.46737);
    puweights_2012.push_back(1.28278);
    puweights_2012.push_back(1.18189);
    puweights_2012.push_back(1.13388);
    puweights_2012.push_back(1.12578);
    puweights_2012.push_back(1.14415);
    puweights_2012.push_back(1.16048);
    puweights_2012.push_back(1.1618);
    puweights_2012.push_back(1.15318);
    puweights_2012.push_back(1.13405);
    puweights_2012.push_back(1.09239);
    puweights_2012.push_back(1.01915);
    puweights_2012.push_back(0.914837);
    puweights_2012.push_back(0.786744);
    puweights_2012.push_back(0.644879);
    puweights_2012.push_back(0.502039);
    puweights_2012.push_back(0.371688);
    puweights_2012.push_back(0.263586);
    puweights_2012.push_back(0.18067);
    puweights_2012.push_back(0.120472);
    puweights_2012.push_back(0.0780184);
    puweights_2012.push_back(0.0486113);
    puweights_2012.push_back(0.0289039);
    puweights_2012.push_back(0.0163367);
    puweights_2012.push_back(0.00879674);
    puweights_2012.push_back(0.00456046);
    puweights_2012.push_back(0.0023098);
    puweights_2012.push_back(0.00115977);
    puweights_2012.push_back(0.000583207);
    puweights_2012.push_back(0.000294815);
    puweights_2012.push_back(0.000149865);
    puweights_2012.push_back(7.62892e-05);
    puweights_2012.push_back(3.87537e-05);
    puweights_2012.push_back(1.96105e-05);
    puweights_2012.push_back(9.87744e-06);
    puweights_2012.push_back(4.95418e-06);
    puweights_2012.push_back(2.47913e-06);
    puweights_2012.push_back(1.23919e-06);
    puweights_2012.push_back(6.19751e-07);
    puweights_2012.push_back(3.10125e-07);
    puweights_2012.push_back(1.54934e-07);
    puweights_2012.push_back(7.71425e-08);
    puweights_2012.push_back(3.8182e-08);
    puweights_2012.push_back(1.87455e-08);
    puweights_2012.push_back(9.10765e-09);
    puweights_2012.push_back(9.19802e-09);
    */    

    // updated with 53X Jan22 rereco
    puweights_2012.push_back(0.242421);
    puweights_2012.push_back(0.314567);
    puweights_2012.push_back(0.328551);
    puweights_2012.push_back(0.341325);
    puweights_2012.push_back(0.311977);
    puweights_2012.push_back(0.560934);
    puweights_2012.push_back(0.443155);
    puweights_2012.push_back(0.444942);
    puweights_2012.push_back(0.625725);
    puweights_2012.push_back(0.940482);
    puweights_2012.push_back(1.33114);
    puweights_2012.push_back(1.67986);
    puweights_2012.push_back(1.7303);
    puweights_2012.push_back(1.54809);
    puweights_2012.push_back(1.32353);
    puweights_2012.push_back(1.15668);
    puweights_2012.push_back(1.07105);
    puweights_2012.push_back(1.04548);
    puweights_2012.push_back(1.06319);
    puweights_2012.push_back(1.10733);
    puweights_2012.push_back(1.14954);
    puweights_2012.push_back(1.17563);
    puweights_2012.push_back(1.188);
    puweights_2012.push_back(1.18718);
    puweights_2012.push_back(1.16671);
    puweights_2012.push_back(1.121);
    puweights_2012.push_back(1.04977);
    puweights_2012.push_back(0.956974);
    puweights_2012.push_back(0.84729);
    puweights_2012.push_back(0.727003);
    puweights_2012.push_back(0.603974);
    puweights_2012.push_back(0.485796);
    puweights_2012.push_back(0.377733);
    puweights_2012.push_back(0.283343);
    puweights_2012.push_back(0.204364);
    puweights_2012.push_back(0.14118);
    puweights_2012.push_back(0.0934506);
    puweights_2012.push_back(0.059445);
    puweights_2012.push_back(0.0365081);
    puweights_2012.push_back(0.0218306);
    puweights_2012.push_back(0.012844);
    puweights_2012.push_back(0.00753269);
    puweights_2012.push_back(0.00447223);
    puweights_2012.push_back(0.00273386);
    puweights_2012.push_back(0.00175157);
    puweights_2012.push_back(0.00118879);
    puweights_2012.push_back(0.000857334);
    puweights_2012.push_back(0.000653996);
    puweights_2012.push_back(0.000522478);
    puweights_2012.push_back(0.000432433);
    puweights_2012.push_back(0.000367567);
    puweights_2012.push_back(0.000318451);
    puweights_2012.push_back(0.000279865);
    puweights_2012.push_back(0.000248423);
    puweights_2012.push_back(0.000221711);
    puweights_2012.push_back(0.000198398);
    puweights_2012.push_back(0.000177509);
    puweights_2012.push_back(0.000158456);
    puweights_2012.push_back(0.000140801);
    puweights_2012.push_back(0.000261544);

    puweights_2012_52.push_back(0.0447136);     
    puweights_2012_52.push_back(0.11785);       
    puweights_2012_52.push_back(0.23825);
    puweights_2012_52.push_back(1.08447);
    puweights_2012_52.push_back(0.102575);
    puweights_2012_52.push_back(0.454605);
    puweights_2012_52.push_back(1.79761);
    puweights_2012_52.push_back(4.00271);
    puweights_2012_52.push_back(6.83281);
    puweights_2012_52.push_back(9.83701);
    puweights_2012_52.push_back(10.7966);
    puweights_2012_52.push_back(12.2356);
    puweights_2012_52.push_back(10.0247);
    puweights_2012_52.push_back(8.49395);
    puweights_2012_52.push_back(7.1125);
    puweights_2012_52.push_back(5.69527);
    puweights_2012_52.push_back(4.31256);
    puweights_2012_52.push_back(3.19305);
    puweights_2012_52.push_back(2.42035);
    puweights_2012_52.push_back(1.91666);
    puweights_2012_52.push_back(1.58485);
    puweights_2012_52.push_back(1.36297);
    puweights_2012_52.push_back(1.21166);
    puweights_2012_52.push_back(1.09466);
    puweights_2012_52.push_back(0.978941);
    puweights_2012_52.push_back(0.84653);
    puweights_2012_52.push_back(0.699235);
    puweights_2012_52.push_back(0.548996);
    puweights_2012_52.push_back(0.408673);
    puweights_2012_52.push_back(0.288194);
    puweights_2012_52.push_back(0.193367);
    puweights_2012_52.push_back(0.124653);
    puweights_2012_52.push_back(0.0781124);
    puweights_2012_52.push_back(0.0479268);
    puweights_2012_52.push_back(0.0287763);
    puweights_2012_52.push_back(0.0167744);
    puweights_2012_52.push_back(0.00941834);
    puweights_2012_52.push_back(0.00507877);
    puweights_2012_52.push_back(0.00264364);
    puweights_2012_52.push_back(0.00134612);
    puweights_2012_52.push_back(0.000682678);
    puweights_2012_52.push_back(0.000351412);
    puweights_2012_52.push_back(0.0001864);
    puweights_2012_52.push_back(0.00010259);
    puweights_2012_52.push_back(5.87818e-05);
    puweights_2012_52.push_back(3.5033e-05);
    puweights_2012_52.push_back(2.17116e-05);
    puweights_2012_52.push_back(1.39777e-05);
    puweights_2012_52.push_back(9.36123e-06);
    puweights_2012_52.push_back(6.53328e-06);
    puweights_2012_52.push_back(4.76598e-06);
    puweights_2012_52.push_back(3.64139e-06);
    puweights_2012_52.push_back(2.92018e-06);
    puweights_2012_52.push_back(2.4602e-06);
    puweights_2012_52.push_back(2.17291e-06);
    puweights_2012_52.push_back(2.01107e-06);
    puweights_2012_52.push_back(1.94392e-06);
    puweights_2012_52.push_back(1.9598e-06);
    puweights_2012_52.push_back(2.0583e-06);
    puweights_2012_52.push_back(2.24895e-06);


    hist_puweights_2012 = new TH1F("hist_puweights_2012","",60,0.,60.);

    for(int k=0;k<60;k++){
        hist_puweights_2012->SetBinContent(k+1,puweights_2012[k]);
    }

    hist_puweights_2012_52 = new TH1F("hist_puweights_2012_52","",60,0.,60.);

    for(int k=0;k<60;k++){
        hist_puweights_2012_52->SetBinContent(k+1,puweights_2012_52[k]);
    }

    hist_puweights_2011 = new TH1F("hist_puweights_2011","",50,0.,50.);

    for(int k=0;k<50;k++){
        hist_puweights_2011->SetBinContent(k+1,puweights_2011[k]);
    }
}

float getPUWeight2011(float numsim) {
    return hist_puweights_2011->GetBinContent(hist_puweights_2011->FindBin(numsim));
}

float getPUWeight2012(float numsim, int mode=1) {
    if (mode == 1) return hist_puweights_2012->GetBinContent(hist_puweights_2012->FindBin(numsim));
    else           return hist_puweights_2012_52->GetBinContent(hist_puweights_2012_52->FindBin(numsim));
}

#endif  
