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

float melacut           = -1.0;
float z1min             = 40.;
float z2min             = 12.;
float massLow           = 100.;
float massHigh          = 180.;
int   nBinsMass2D       = 40.;
int   nBinsMELA2D       = 30.; 
bool  do7TeV            = true; 
std::string treeFolder7 = "/home/avartak/CMS/Higgs/HZZ4L/CMSSW_4_2_8_patch7/src/WWAnalysis/AnalysisStep/trees/";
std::string treeFolder8 = "/home/avartak/CMS/Higgs/HZZ4L/CMSSW_5_2_4_patch4/src/WWAnalysis/AnalysisStep/trees/";

void fillMass(std::string ggh_rootfile, std::string vbf_rootfile, TH2F& h2D_sig_em, TH2F& h2D_sig_mm, TH2F& h2D_sig_ee) {
    std::string base_folder = do7TeV ? treeFolder7 : treeFolder8;
    ZZYieldMaker   ymaker_ghzz;
    ZZYieldMaker   ymaker_qhzz;

    std::cout << "Filling " << ggh_rootfile << " and " << vbf_rootfile << std::endl;

    ymaker_ghzz.fill(base_folder+ggh_rootfile, 1.0, 0.0, true );
    ymaker_qhzz.fill(base_folder+vbf_rootfile, 1.0, 0.0, true );

    ymaker_ghzz.get2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_mm);
    ymaker_ghzz.get2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_ee);
    ymaker_ghzz.get2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_em);
    ymaker_ghzz.get2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_em);
    
    ymaker_qhzz.get2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_mm);
    ymaker_qhzz.get2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_ee);
    ymaker_qhzz.get2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_em);
    ymaker_qhzz.get2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_sig_em);
}

void fillBkg(TH2F& h2D_bkg_em, TH2F& h2D_bkg_mm, TH2F& h2D_bkg_ee) {
    ZZYieldMaker   ymaker_ggzz;
    ZZYieldMaker   ymaker_qqzz;

    if (do7TeV) {
        ymaker_qqzz.fill(treeFolder7+"hzzTree_id121.root" , getBkgXsec(121)*5.05/evt_7TeV[121], 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"hzzTree_id122.root" , getBkgXsec(122)*5.05/evt_7TeV[122], 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"hzzTree_id123.root" , getBkgXsec(123)*5.05/evt_7TeV[123], 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"hzzTree_id124.root" , getBkgXsec(124)*5.05/evt_7TeV[124], 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"hzzTree_id125.root" , getBkgXsec(125)*5.05/evt_7TeV[125], 0.0, false);
        ymaker_qqzz.fill(treeFolder7+"hzzTree_id126.root" , getBkgXsec(126)*5.05/evt_7TeV[126], 0.0, false);
        ymaker_ggzz.fill(treeFolder7+"hzzTree_id101.root" , getBkgXsec(101)*5.05/evt_7TeV[101], 0.0, false);
        ymaker_ggzz.fill(treeFolder7+"hzzTree_id100.root" , getBkgXsec(100)*5.05/evt_7TeV[100], 0.0, false);
    }
    else {
        ymaker_qqzz.fill(treeFolder8+"hzzTree_id102.root" , getBkgXsec(102)*5.26/evt_8TeV[102], 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"hzzTree_id103.root" , getBkgXsec(103)*5.26/evt_8TeV[103], 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"hzzTree_id104.root" , getBkgXsec(104)*5.26/evt_8TeV[104], 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"hzzTree_id105.root" , getBkgXsec(105)*5.26/evt_8TeV[105], 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"hzzTree_id106.root" , getBkgXsec(106)*5.26/evt_8TeV[106], 0.0, false);
        ymaker_qqzz.fill(treeFolder8+"hzzTree_id107.root" , getBkgXsec(107)*5.26/evt_8TeV[107], 0.0, false);
        ymaker_ggzz.fill(treeFolder8+"hzzTree_id101.root" , getBkgXsec(101)*5.26/evt_8TeV[101], 0.0, false);
        ymaker_ggzz.fill(treeFolder8+"hzzTree_id100.root" , getBkgXsec(100)*5.26/evt_8TeV[100], 0.0, false);
    }

    ymaker_ggzz.get2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_mm);
    ymaker_ggzz.get2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_ee);
    ymaker_ggzz.get2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
    ymaker_ggzz.get2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);

    ymaker_qqzz.get2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_mm);
    ymaker_qqzz.get2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_ee);
    ymaker_qqzz.get2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
    ymaker_qqzz.get2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
}

void fillZjets(TH2F& h2D_bkg_em, TH2F& h2D_bkg_mm, TH2F& h2D_bkg_ee, bool doSS) {
    std::string base_folder = do7TeV ? treeFolder7 : treeFolder8;

    FakeRateCalculator FR(base_folder+"hzzTree.root", do7TeV, 40, 120, 0.0, 0.0, true);

    ZXYieldMaker   ymaker_zxss;
    ymaker_zxss.fill(base_folder+"hzzTree.root", 1.0, FR, doSS);

    ymaker_zxss.get2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_mm);
    ymaker_zxss.get2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_ee);
    ymaker_zxss.get2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
    ymaker_zxss.get2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);

    ymaker_zxss.get2DHist(0, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_mm);
    ymaker_zxss.get2DHist(1, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_ee);
    ymaker_zxss.get2DHist(2, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
    ymaker_zxss.get2DHist(3, z1min, z2min, massLow, massHigh, melacut, &h2D_bkg_em);
}

void smooth(TH2F* h, float threshold) {
    TH2F* hist = (TH2F*)(h->Clone((std::string(h->GetName())+"_temp").c_str()));
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        for(int j = 1; j <= hist->GetNbinsY(); ++j) {
            float count = 0.;
            float val = hist->GetBinContent(i,j);
            if (val<threshold) {
                if (i-1 != 0)                                                  val += hist->GetBinContent(i-1,j  );
                else                                                           count -= 1.;
                if (i+1 != hist->GetNbinsX()+1)                                val += hist->GetBinContent(i+1,j  );
                else                                                           count -= 1.;
                if (j-1 != 0)                                                  val += hist->GetBinContent(i  ,j-1);
                else                                                           count -= 1.;
                if (j+1 != hist->GetNbinsY()+1)                                val += hist->GetBinContent(i  ,j+1);
                else                                                           count -= 1.;
                if (i-1 != 0 && j-1 != 0)                                      val += hist->GetBinContent(i-1,j-1);
                else                                                           count -= 1.;
                if (i-1 != 0 && j+1 != hist->GetNbinsY()+1)                    val += hist->GetBinContent(i-1,j+1);
                else                                                           count -= 1.;
                if (i+1 != hist->GetNbinsX()+1 && j-1 != 0)                    val += hist->GetBinContent(i+1,j-1);
                else                                                           count -= 1.;
                if (i+1 != hist->GetNbinsX()+1 && j+1 != hist->GetNbinsY()+1)  val += hist->GetBinContent(i+1,j+1);
                else                                                           count -= 1.;
                val /= (9.0+count);
                h->SetBinContent(i,j,val);
            }
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
                val /= (5.0+count);
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

void create2DHistogram() {
    TH2F h2D_sig_em("hist2D_sig_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_sig_mm("hist2D_sig_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_sig_ee("hist2D_sig_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);

    TH2F h2D_bkg_em("hist2D_bkg_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_bkg_mm("hist2D_bkg_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_bkg_ee("hist2D_bkg_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);

    do7TeV = true;
    init(do7TeV);

    fillMass("hzzTree_id200.root", "hzzTree_id250.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id201.root", "hzzTree_id251.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id202.root", "hzzTree_id252.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id203.root", "hzzTree_id253.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id204.root", "hzzTree_id254.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id205.root", "hzzTree_id255.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id206.root", "hzzTree_id256.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id207.root", "hzzTree_id257.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);

    std::cout << "Filling background" << std::endl;
    fillBkg (h2D_bkg_em, h2D_bkg_mm, h2D_bkg_ee);

    do7TeV = false;
    init(do7TeV);

    fillMass("hzzTree_id1115.root", "hzzTree_id2115.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id1117.root", "hzzTree_id2117.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id1119.root", "hzzTree_id2119.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id1120.root", "hzzTree_id2120.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id1121.root", "hzzTree_id2121.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id1123.root", "hzzTree_id2123.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id1124.root", "hzzTree_id2124.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id1125.root", "hzzTree_id2125.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id1126.root", "hzzTree_id2126.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id1127.root", "hzzTree_id2127.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id1150.root", "hzzTree_id2150.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);
    fillMass("hzzTree_id1180.root", "hzzTree_id2180.root", h2D_sig_em, h2D_sig_mm, h2D_sig_ee);

    std::cout << "Filling background" << std::endl;
    fillBkg (h2D_bkg_em, h2D_bkg_mm, h2D_bkg_ee);

//////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Creating alternate shapes ...." << std::endl;
    TH2F h2D_bkg_em_up = *((TH2F*)h2D_bkg_em.Clone("hist2D_bkg_em_up"));
    TH2F h2D_bkg_mm_up = *((TH2F*)h2D_bkg_mm.Clone("hist2D_bkg_mm_up"));
    TH2F h2D_bkg_ee_up = *((TH2F*)h2D_bkg_ee.Clone("hist2D_bkg_ee_up"));

    TH2F h2D_bkg_em_dn = *((TH2F*)h2D_bkg_em.Clone("hist2D_bkg_em_dn"));
    TH2F h2D_bkg_mm_dn = *((TH2F*)h2D_bkg_mm.Clone("hist2D_bkg_mm_dn"));
    TH2F h2D_bkg_ee_dn = *((TH2F*)h2D_bkg_ee.Clone("hist2D_bkg_ee_dn"));

    nBinsMass2D = 4;
    TH2F h2D_zz_em("hist2D_zz_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_zz_mm("hist2D_zz_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_zz_ee("hist2D_zz_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);

    TH2F h2D_zx_em("hist2D_zx_em", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_zx_mm("hist2D_zx_mm", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);
    TH2F h2D_zx_ee("hist2D_zx_ee", "", nBinsMass2D, massLow, massHigh, nBinsMELA2D, 0, 1);

    h2D_zx_em.Sumw2();
    h2D_zx_mm.Sumw2();
    h2D_zx_ee.Sumw2();

    do7TeV = true;
    init(do7TeV);
    fillBkg  (h2D_zz_em, h2D_zz_mm, h2D_zz_ee);
    fillZjets(h2D_zx_em, h2D_zx_mm, h2D_zx_ee, false);

    do7TeV = false;
    init(do7TeV);
    fillBkg  (h2D_zz_em, h2D_zz_mm, h2D_zz_ee);
    fillZjets(h2D_zx_em, h2D_zx_mm, h2D_zx_ee, false);

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
            int k = (i*nBinsMass2D)/h2D_bkg_em.GetNbinsX();
            if (k == nBinsMass2D) k -= 1;
            h2D_bkg_em_up.SetBinContent(i, j, h2D_bkg_em.GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
            h2D_bkg_ee_up.SetBinContent(i, j, h2D_bkg_ee.GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
            h2D_bkg_mm_up.SetBinContent(i, j, h2D_bkg_mm.GetBinContent(i, j)*(fitp0[k]+fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));

            h2D_bkg_em_dn.SetBinContent(i, j, h2D_bkg_em.GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
            h2D_bkg_ee_dn.SetBinContent(i, j, h2D_bkg_ee.GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
            h2D_bkg_mm_dn.SetBinContent(i, j, h2D_bkg_mm.GetBinContent(i, j)*(2.0-fitp0[k]-fitp1[k]*((j-0.5)/h2D_bkg_em.GetNbinsY())));
        }
    }

//////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Smoothing ...." << std::endl;
    normalize(&h2D_sig_em);
    normalize(&h2D_sig_mm);
    normalize(&h2D_sig_ee);

    normalize(&h2D_bkg_em);
    normalize(&h2D_bkg_mm);
    normalize(&h2D_bkg_ee);

    normalize(&h2D_bkg_em_up);
    normalize(&h2D_bkg_mm_up);
    normalize(&h2D_bkg_ee_up);

    normalize(&h2D_bkg_em_dn);
    normalize(&h2D_bkg_mm_dn);
    normalize(&h2D_bkg_ee_dn);

    smooth(&h2D_sig_em, 1.0);
    smooth(&h2D_sig_mm, 1.0);
    smooth(&h2D_sig_ee, 1.0);
    
    smooth(&h2D_bkg_em, 1.0);
    smooth(&h2D_bkg_mm, 1.0);
    smooth(&h2D_bkg_ee, 1.0);

    smooth(&h2D_bkg_em_up, 1.0);
    smooth(&h2D_bkg_mm_up, 1.0);
    smooth(&h2D_bkg_ee_up, 1.0);

    smooth(&h2D_bkg_em_dn, 1.0);
    smooth(&h2D_bkg_mm_dn, 1.0);
    smooth(&h2D_bkg_ee_dn, 1.0);

    /*
    smoothVertical(&h2D_sig_em, 1.0);
    smoothVertical(&h2D_sig_mm, 1.0);
    smoothVertical(&h2D_sig_ee, 1.0);
    
    smoothVertical(&h2D_bkg_em, 1.0);
    smoothVertical(&h2D_bkg_mm, 1.0);
    smoothVertical(&h2D_bkg_ee, 1.0);

    smoothVertical(&h2D_bkg_em_up, 1.0);
    smoothVertical(&h2D_bkg_mm_up, 1.0);
    smoothVertical(&h2D_bkg_ee_up, 1.0);

    smoothVertical(&h2D_bkg_em_dn, 1.0);
    smoothVertical(&h2D_bkg_mm_dn, 1.0);
    smoothVertical(&h2D_bkg_ee_dn, 1.0);

    smoothHorizontal(&h2D_sig_em, 1.0);
    smoothHorizontal(&h2D_sig_mm, 1.0);
    smoothHorizontal(&h2D_sig_ee, 1.0);

    smoothHorizontal(&h2D_bkg_em, 1.0);
    smoothHorizontal(&h2D_bkg_mm, 1.0);
    smoothHorizontal(&h2D_bkg_ee, 1.0);

    smoothHorizontal(&h2D_bkg_em_up, 1.0);
    smoothHorizontal(&h2D_bkg_mm_up, 1.0);
    smoothHorizontal(&h2D_bkg_ee_up, 1.0);

    smoothHorizontal(&h2D_bkg_em_dn, 1.0);
    smoothHorizontal(&h2D_bkg_mm_dn, 1.0);
    smoothHorizontal(&h2D_bkg_ee_dn, 1.0);
    */

//////////////////////////////////////////////////////////////////////////////////////////////////////

    std::cout << "Writing to file ...." << std::endl;
    TFile histout("mela2DShapes.root", "RECREATE");
    h2D_sig_em.Write();
    h2D_sig_mm.Write();
    h2D_sig_ee.Write();
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
