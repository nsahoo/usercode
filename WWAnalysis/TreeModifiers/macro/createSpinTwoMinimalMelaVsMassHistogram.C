#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <utility>
#include <sstream>
#include <cmath>
#include <algorithm>

#include "YieldMaker.h"

using namespace RooFit;

float melacut           = 0.5;
float z1min             = 40.;
float z2min             = 12.;
float massLow           = 100.;
float massHigh          = 180.;
int   nBinsMass2D       = 40; //2GeV bins
int   nBinsMELA2D       = 30.; 
bool  do7TeV            = false; 
std::string treeFolder7 = "/data1/sixie/ntuples/step2/HZZ4L_42X_S1_V10_S2_V14/";
std::string treeFolder8 = "/data1/sixie/ntuples/step2/HZZ4L_53X_S1_V10_S2_V14/";

void fillMass(std::string rootfile, bool isSMHiggs, TH2F& h2D_sig_em, TH2F& h2D_sig_mm, TH2F& h2D_sig_ee) {
    std::string base_folder = do7TeV ? treeFolder7 : treeFolder8;
    ZZYieldMaker   ymaker;

    std::cout << "Filling " << rootfile << std::endl;

    ymaker.fill(base_folder+rootfile, 1.0, 0.0, isSMHiggs );

    ymaker.fillMelaSpinTwoMinimalVsMass2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_mm);
    ymaker.fillMelaSpinTwoMinimalVsMass2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_ee);
    ymaker.fillMelaSpinTwoMinimalVsMass2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_em);
    ymaker.fillMelaSpinTwoMinimalVsMass2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_em);
    
}

void fillBkg(TH2F& h2D_bkg_em, TH2F& h2D_bkg_mm, TH2F& h2D_bkg_ee) {
    ZZYieldMaker   ymaker_ggzz;
    ZZYieldMaker   ymaker_qqzz;

    if (do7TeV) {
        ymaker_qqzz.fill(treeFolder7+"/MC/7TeV/yesRegrYesCalib/hzzTree_id102.root" ,  getBkgXsec(102)/evt_7TeV(102), 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"/MC/7TeV/yesRegrYesCalib/hzzTree_id103.root" ,  getBkgXsec(103)/evt_7TeV(103), 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"/MC/7TeV/yesRegrYesCalib/hzzTree_id104.root" ,  getBkgXsec(104)/evt_7TeV(104), 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"/MC/7TeV/yesRegrYesCalib/hzzTree_id105.root" ,  getBkgXsec(105)/evt_7TeV(105), 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"/MC/7TeV/yesRegrYesCalib/hzzTree_id106.root" ,  getBkgXsec(106)/evt_7TeV(106), 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"/MC/7TeV/yesRegrYesCalib/hzzTree_id107.root" ,  getBkgXsec(107)/evt_7TeV(107), 0.0, false);
        ymaker_ggzz.fill(treeFolder7+"/MC/7TeV/yesRegrYesCalib/hzzTree_id101.root" ,  getBkgXsec(101)/evt_7TeV(101), 0.0, false);
        ymaker_ggzz.fill(treeFolder7+"/MC/7TeV/yesRegrYesCalib/hzzTree_id100.root" ,  getBkgXsec(100)/evt_7TeV(100), 0.0, false);
    }
    else {
        ymaker_qqzz.fill(treeFolder8+"/MC/8TeV/yesRegrYesCalib/hzzTree_id102.root" ,  getBkgXsec(102)/evt_8TeV(102), 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"/MC/8TeV/yesRegrYesCalib/hzzTree_id103.root" ,  getBkgXsec(103)/evt_8TeV(103), 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"/MC/8TeV/yesRegrYesCalib/hzzTree_id104.root" ,  getBkgXsec(104)/evt_8TeV(104), 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"/MC/8TeV/yesRegrYesCalib/hzzTree_id105.root" ,  getBkgXsec(105)/evt_8TeV(105), 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"/MC/8TeV/yesRegrYesCalib/hzzTree_id106.root" ,  getBkgXsec(106)/evt_8TeV(106), 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"/MC/8TeV/yesRegrYesCalib/hzzTree_id107.root" ,  getBkgXsec(107)/evt_8TeV(107), 0.0, false);
        ymaker_ggzz.fill(treeFolder8+"/MC/8TeV/yesRegrYesCalib/hzzTree_id101.root" ,  getBkgXsec(101)/evt_8TeV(101), 0.0, false);
        ymaker_ggzz.fill(treeFolder8+"/MC/8TeV/yesRegrYesCalib/hzzTree_id100.root" ,  getBkgXsec(100)/evt_8TeV(100), 0.0, false);
    }

    ymaker_ggzz.fillMelaSpinTwoMinimalVsMass2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_mm);
    ymaker_ggzz.fillMelaSpinTwoMinimalVsMass2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_ee);
    ymaker_ggzz.fillMelaSpinTwoMinimalVsMass2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
    ymaker_ggzz.fillMelaSpinTwoMinimalVsMass2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);

    ymaker_qqzz.fillMelaSpinTwoMinimalVsMass2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_mm);
    ymaker_qqzz.fillMelaSpinTwoMinimalVsMass2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_ee);
    ymaker_qqzz.fillMelaSpinTwoMinimalVsMass2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
    ymaker_qqzz.fillMelaSpinTwoMinimalVsMass2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
}

void fillZjets(TH2F& h2D_bkg_em, TH2F& h2D_bkg_mm, TH2F& h2D_bkg_ee, bool doSS) {
    std::string base_folder = do7TeV ? treeFolder7+"DATA/7TeV/yesRegrYesCalib/" : treeFolder8+"DATA/8TeV/yesRegrYesCalib/";

    FakeRateCalculator FR(base_folder+"hcp.root", do7TeV, 40, 120, 0.0, 0.0, true);

    ZXYieldMaker   ymaker_zxss;
    ymaker_zxss.fill(base_folder+"hcp.root", 1.0, FR, doSS);

    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_mm);
    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_ee);
    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);

    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_mm);
    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_ee);
    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
    ymaker_zxss.fillMelaSpinTwoMinimalVsMass2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
}

void smooth(TH2F* h, int binmin, int binmax, float arraysize) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_temp").c_str()));
    for (int i = binmin; i <= hist->GetNbinsX() && i <= binmax; ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float count = 0.;
            float val = 0;
            for (int m = i-arraysize; m <= i+arraysize; m++) {
                //for (int n = j-arraysize; n <= j+arraysize; n++) {
                for (int n = j-1; n <= j+1; n++) {
                    if (m >= binmin && m <= hist->GetNbinsX() && m <= binmax && n > 0 && n <= hist->GetNbinsY()) {
                        count += 1.0;
                        val += hist->GetBinContent(m, n);
                    }
                }
            }
            val /= count;
            if (val > 0.0) h->SetBinContent(i,j,val);
            else h->SetBinContent(i,j,0.0001);
        }
    }
}

void smoothSwiss(TH2F* h, float threshold) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_temp").c_str()));
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float val = hist->GetBinContent(i,j);
            float count = 0.;
            if (val < threshold) {
                if (i-1 != 0)                                                  val += hist->GetBinContent(i-1,j  );
                else                                                           count -= 1.0;;
                if (i+1 != hist->GetNbinsX()+1)                                val += hist->GetBinContent(i+1,j  );
                else                                                           count -= 1.0;
                if (j-1 != 0)                                                  val += hist->GetBinContent(i  ,j-1);
                else                                                           count -= 1.0;;
                if (j+1 != hist->GetNbinsY()+1)                                val += hist->GetBinContent(i  ,j+1);
                else                                                           count -= 1.0;
                val /= (9.0+count);
                h->SetBinContent(i,j,val);
            }
        }
    }
}

void smoothVertical(TH2F* h, float threshold) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_vert").c_str()));
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float val = hist->GetBinContent(i,j);
            float count = 0.;
            if (val < threshold) {
                if (j-1 != 0)                                                  val += hist->GetBinContent(i  ,j-1);
                else                                                           count -= 1.0;;
                if (j+1 != hist->GetNbinsY()+1)                                val += hist->GetBinContent(i  ,j+1);
                else                                                           count -= 1.0;
                val /= (3.0+count);
                h->SetBinContent(i,j,val);
            }
        }
    }
}

void smoothHorizontal(TH2F* h, float threshold) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_hori").c_str()));
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float val = hist->GetBinContent(i,j);
            float count = 0.;
            if (val < threshold) { 
                if (i-1 != 0)                                                  val += hist->GetBinContent(i-1,j  );
                else                                                           count -= 1.0;;
                if (i+1 != hist->GetNbinsX()+1)                                val += hist->GetBinContent(i+1,j  );
                else                                                           count -= 1.0;
                val /= (3.0+count);
                h->SetBinContent(i,j,val);
            }
        }
    }
}


void normalize(TH2F* hist) {
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        double histsum = 0;
        for(int j = 1; j <= hist->GetNbinsY(); ++j) histsum += hist->GetBinContent(i,j);
        if (histsum > 0) histsum = 1.0/histsum;
        for(int j = 1; j <= hist->GetNbinsY(); ++j) hist->SetBinContent(i,j, hist->GetBinContent(i,j) * histsum);
    }
}

void normalize(TH2F* hist1, TH2F* hist2) {
    for (int i = 1; i <= hist1->GetNbinsX(); ++i) {
        double hist1sum = 0;
        double hist2sum = 0;
        for(int j = 1; j <= hist1->GetNbinsY(); ++j) hist1sum += hist1->GetBinContent(i,j);
        for(int j = 1; j <= hist2->GetNbinsY(); ++j) hist2sum += hist2->GetBinContent(i,j);
        if (hist1sum > 0) hist1sum = hist2sum/hist1sum;
        for(int j = 1; j <= hist1->GetNbinsY(); ++j) hist1->SetBinContent(i,j, hist1->GetBinContent(i,j) * hist1sum);
    }
}

void createSpinTwoMinimalMelaVsMassHistogram() {
    TH2F h2D_sig_em("hist2D_sig_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_sig_mm("hist2D_sig_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_sig_ee("hist2D_sig_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);

    TH2F h2D_spinTwoMinimalsig_em("hist2D_spinTwoMinimalsig_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_spinTwoMinimalsig_mm("hist2D_spinTwoMinimalsig_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_spinTwoMinimalsig_ee("hist2D_spinTwoMinimalsig_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);

    TH2F h2D_bkg_em("hist2D_bkg_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_bkg_mm("hist2D_bkg_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_bkg_ee("hist2D_bkg_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);


    //***********************************************************************************************
    //Fill histogram for scalar higgs
    //***********************************************************************************************

//     do7TeV = true;
//     init(do7TeV);

//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1115.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1120.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1130.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1140.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1150.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1160.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1170.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1180.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1190.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1200.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1210.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1220.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1230.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1250.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1275.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1300.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1325.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1350.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1375.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1400.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1425.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1450.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1475.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
// //     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1500.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1525.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1550.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1575.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1600.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1650.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1700.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1750.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1800.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1900.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id1950.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
//     fillMass("/MC/7TeV/yesRegrYesCalib/hzzTree_id11000.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);

//     std::cout << "Filling background" << std::endl;
//     fillBkg (h2D_bkg_em, h2D_bkg_mm, h2D_bkg_ee);

    do7TeV = false;
    init(do7TeV);

    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1115.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1117.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1118.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1119.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1120.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1121.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1122.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1123.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1124.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
     fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1125.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1126.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1127.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1128.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1129.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1130.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1135.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1140.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1145.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1150.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1160.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1170.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1180.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1190.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1200.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1220.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1250.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1275.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1300.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1325.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1350.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1375.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1400.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1425.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1450.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1475.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1500.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1525.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1550.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1575.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1600.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1650.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1700.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1750.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1800.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1850.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1900.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id1950.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id11000.root", true, h2D_sig_em, h2D_sig_mm, h2D_sig_ee);

    //***********************************************************************************************
    //Fill histogram for pseudo-scalar higgs
    //***********************************************************************************************
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id8110.root", false, h2D_spinTwoMinimalsig_em, h2D_spinTwoMinimalsig_mm, h2D_spinTwoMinimalsig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id8120.root", false, h2D_spinTwoMinimalsig_em, h2D_spinTwoMinimalsig_mm, h2D_spinTwoMinimalsig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id8125.root", false, h2D_spinTwoMinimalsig_em, h2D_spinTwoMinimalsig_mm, h2D_spinTwoMinimalsig_ee);
    fillMass("/MC/8TeV/yesRegrYesCalib/hzzTree_id8130.root", false, h2D_spinTwoMinimalsig_em, h2D_spinTwoMinimalsig_mm, h2D_spinTwoMinimalsig_ee);


    //***********************************************************************************************
    //Fill histogram for Backgrounds
    //***********************************************************************************************
    std::cout << "Filling background" << std::endl;
    fillBkg (h2D_bkg_em, h2D_bkg_mm, h2D_bkg_ee);


    //***********************************************************************************************
    //Fill histograms for shape systematics
    //***********************************************************************************************
    std::cout << "Creating alternate shapes ...." << std::endl;
    TH2F h2D_bkg_em_up = *((TH2F*)h2D_bkg_em.Clone("hist2D_bkg_em_up"));
    TH2F h2D_bkg_mm_up = *((TH2F*)h2D_bkg_mm.Clone("hist2D_bkg_mm_up"));
    TH2F h2D_bkg_ee_up = *((TH2F*)h2D_bkg_ee.Clone("hist2D_bkg_ee_up"));

    TH2F h2D_bkg_em_dn = *((TH2F*)h2D_bkg_em.Clone("hist2D_bkg_em_dn"));
    TH2F h2D_bkg_mm_dn = *((TH2F*)h2D_bkg_mm.Clone("hist2D_bkg_mm_dn"));
    TH2F h2D_bkg_ee_dn = *((TH2F*)h2D_bkg_ee.Clone("hist2D_bkg_ee_dn"));

    nBinsMass2D = 6;
    Double_t systbins[] = {100., 120., 140., 160., 180., 300., 1000.};
    TH2F h2D_zz_em("hist2D_zz_em", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    TH2F h2D_zz_mm("hist2D_zz_mm", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    TH2F h2D_zz_ee("hist2D_zz_ee", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);

    TH2F h2D_zx_em("hist2D_zx_em", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    TH2F h2D_zx_mm("hist2D_zx_mm", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);
    TH2F h2D_zx_ee("hist2D_zx_ee", "", nBinsMass2D, systbins, nBinsMELA2D, 0, 1);

    h2D_zx_em.Sumw2();
    h2D_zx_mm.Sumw2();
    h2D_zx_ee.Sumw2();

//     do7TeV = true;
//     init(do7TeV);
//     fillBkg  (h2D_zz_em, h2D_zz_mm, h2D_zz_ee);
//     fillZjets(h2D_zx_em, h2D_zx_mm, h2D_zx_ee, true);

    do7TeV = false;
    init(do7TeV);
    fillBkg  (h2D_zz_em, h2D_zz_mm, h2D_zz_ee);
    fillZjets(h2D_zx_em, h2D_zx_mm, h2D_zx_ee, true);

    h2D_zz_em.Add(&h2D_zz_mm);
    h2D_zz_em.Add(&h2D_zz_ee);
    h2D_zx_em.Add(&h2D_zx_mm);
    h2D_zx_em.Add(&h2D_zx_ee);
    normalize(&h2D_zz_em, &h2D_zx_em);

    std::vector<float> fitp0;
    std::vector<float> fitp1;
    for (int i = 1; i <= nBinsMass2D; i++) {
        stringstream   histname;
        histname   << "hist";
        histname   << i;
        TH1F hist(histname.str().c_str()  , "", nBinsMELA2D, 0., 1.);
        TF1 f("f", "[0]+[1]*x", 0, 1000.);
        for (int j = 1; j <= nBinsMELA2D; j++) {
            hist.SetBinContent  (j, h2D_zx_em.GetBinContent(i,j)/h2D_zz_em.GetBinContent(i,j));
            hist.SetBinError    (j, h2D_zx_em.GetBinError  (i,j)/h2D_zz_em.GetBinContent(i,j));
        }
        hist.Fit(&f);
        fitp0.push_back(f.GetParameter(0));
        fitp1.push_back(f.GetParameter(1));
    }

    for (int i = 1; i <= h2D_bkg_em.GetNbinsX(); i++) {
        for (int j = 1; j <= h2D_bkg_em.GetNbinsY(); j++) {
            int k = nBinsMass2D;
            for (int m = 0; m < nBinsMass2D; m++) {
                if ((massLow + i*2) > systbins[m] && (massLow + i*2) <= systbins[m+1]) {
                    k = m;   
                }       
            }

            if (k < nBinsMass2D) {

                h2D_bkg_em_up.SetBinContent(i, j, h2D_bkg_em.GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
                h2D_bkg_ee_up.SetBinContent(i, j, h2D_bkg_ee.GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
                h2D_bkg_mm_up.SetBinContent(i, j, h2D_bkg_mm.GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));

                h2D_bkg_em_dn.SetBinContent(i, j, h2D_bkg_em.GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
                h2D_bkg_ee_dn.SetBinContent(i, j, h2D_bkg_ee.GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
                h2D_bkg_mm_dn.SetBinContent(i, j, h2D_bkg_mm.GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
            }
        }
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Smoothing ...." << std::endl;
    normalize(&h2D_sig_em);
    normalize(&h2D_sig_mm);
    normalize(&h2D_sig_ee);

    normalize(&h2D_spinTwoMinimalsig_em);
    normalize(&h2D_spinTwoMinimalsig_mm);
    normalize(&h2D_spinTwoMinimalsig_ee);

    normalize(&h2D_bkg_em);
    normalize(&h2D_bkg_mm);
    normalize(&h2D_bkg_ee);

    normalize(&h2D_bkg_em_up);
    normalize(&h2D_bkg_mm_up);
    normalize(&h2D_bkg_ee_up);

    normalize(&h2D_bkg_em_dn);
    normalize(&h2D_bkg_mm_dn);
    normalize(&h2D_bkg_ee_dn);


    int binmin1 = 1;
    int binmin2 = 41;
    int binmin3 = 61;

    int binmax1 = 40;
    int binmax2 = 60;
    int binmax3 = 450;

    smooth(&h2D_sig_em,    binmin1, binmax1, 1);
    smooth(&h2D_sig_mm,    binmin1, binmax1, 1);
    smooth(&h2D_sig_ee,    binmin1, binmax1, 1);
    
    smooth(&h2D_spinTwoMinimalsig_em,    binmin1, binmax1, 1);
    smooth(&h2D_spinTwoMinimalsig_mm,    binmin1, binmax1, 1);
    smooth(&h2D_spinTwoMinimalsig_ee,    binmin1, binmax1, 1);
    
    smooth(&h2D_bkg_em,    binmin1, binmax1, 1);
    smooth(&h2D_bkg_mm,    binmin1, binmax1, 1);
    smooth(&h2D_bkg_ee,    binmin1, binmax1, 1);

    smooth(&h2D_bkg_em_up, binmin1, binmax1, 1);
    smooth(&h2D_bkg_mm_up, binmin1, binmax1, 1);
    smooth(&h2D_bkg_ee_up, binmin1, binmax1, 1);

    smooth(&h2D_bkg_em_dn, binmin1, binmax1, 1);
    smooth(&h2D_bkg_mm_dn, binmin1, binmax1, 1);
    smooth(&h2D_bkg_ee_dn, binmin1, binmax1, 1);

    smooth(&h2D_sig_em,    binmin2, binmax2, 2);
    smooth(&h2D_sig_mm,    binmin2, binmax2, 2);
    smooth(&h2D_sig_ee,    binmin2, binmax2, 2);
    
    smooth(&h2D_spinTwoMinimalsig_em,    binmin2, binmax2, 2);
    smooth(&h2D_spinTwoMinimalsig_mm,    binmin2, binmax2, 2);
    smooth(&h2D_spinTwoMinimalsig_ee,    binmin2, binmax2, 2);
    
    smooth(&h2D_bkg_em,    binmin2, binmax2, 2);
    smooth(&h2D_bkg_mm,    binmin2, binmax2, 2);
    smooth(&h2D_bkg_ee,    binmin2, binmax2, 2);

    smooth(&h2D_bkg_em_up, binmin2, binmax2, 2);
    smooth(&h2D_bkg_mm_up, binmin2, binmax2, 2);
    smooth(&h2D_bkg_ee_up, binmin2, binmax2, 2);

    smooth(&h2D_bkg_em_dn, binmin2, binmax2, 2);
    smooth(&h2D_bkg_mm_dn, binmin2, binmax2, 2);
    smooth(&h2D_bkg_ee_dn, binmin2, binmax2, 2);

    smooth(&h2D_sig_em,    binmin3, binmax3, 3);
    smooth(&h2D_sig_mm,    binmin3, binmax3, 3);
    smooth(&h2D_sig_ee,    binmin3, binmax3, 3);
    
    smooth(&h2D_spinTwoMinimalsig_em,    binmin3, binmax3, 3);
    smooth(&h2D_spinTwoMinimalsig_mm,    binmin3, binmax3, 3);
    smooth(&h2D_spinTwoMinimalsig_ee,    binmin3, binmax3, 3);
    
    smooth(&h2D_bkg_em,    binmin3, binmax3, 3);
    smooth(&h2D_bkg_mm,    binmin3, binmax3, 3);
    smooth(&h2D_bkg_ee,    binmin3, binmax3, 3);

    smooth(&h2D_bkg_em_up, binmin3, binmax3, 3);
    smooth(&h2D_bkg_mm_up, binmin3, binmax3, 3);
    smooth(&h2D_bkg_ee_up, binmin3, binmax3, 3);

    smooth(&h2D_bkg_em_dn, binmin3, binmax3, 3);
    smooth(&h2D_bkg_mm_dn, binmin3, binmax3, 3);
    smooth(&h2D_bkg_ee_dn, binmin3, binmax3, 3);

//////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Writing to file ...." << std::endl;
    TFile histout("SpinTwoMinimalMelaVsMass2DShapes.root", "RECREATE");
    h2D_sig_em.Write();
    h2D_sig_mm.Write();
    h2D_sig_ee.Write();
    h2D_spinTwoMinimalsig_em.Write();
    h2D_spinTwoMinimalsig_mm.Write();
    h2D_spinTwoMinimalsig_ee.Write();
    h2D_bkg_em.Write();
    h2D_bkg_mm.Write();
    h2D_bkg_ee.Write();
    h2D_bkg_em_up.Write();
    h2D_bkg_mm_up.Write();
    h2D_bkg_ee_up.Write();
    h2D_bkg_em_dn.Write();
    h2D_bkg_mm_dn.Write();
    h2D_bkg_ee_dn.Write();
    histout.Close();
}
