#include <iostream>
#include <TH1F.h>
#include <THStack.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLatex.h>
#include <Math/Math.h>
#include <Math/QuantFuncMathCore.h>
#include "FakeRateCalculator.h"
#include "FitMaker.h"
#include "YieldMaker.h"

#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "HiggsAnalysis/CombinedLimit/interface/VerticalInterpHistPdf.h"
#include <RooDataHist.h>
#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h> 
#include <RooLandau.h> 
#include <RooBreitWigner.h>
#include <RooCBShape.h> 
#include <RooFFTConvPdf.h>
#include <RooFormulaVar.h>


TH1* getSigHist(const char* filename, int channel, const char* histname, int nbins, float xmin, float xmax, bool is7) {

    TFile* file = new TFile(filename);

    string wname = "w";
    string tevstr = is7 ? "_7TeV" : "_8TeV";
    string suffix = "";
    if (channel == 0) suffix += "_4mu";
    if (channel == 1) suffix += "_4e";
    if (channel == 2) suffix += "_2e2mu";
    string varname = "sig_ggh" + suffix + tevstr;

    RooFormulaVar var_mean   = *(RooFormulaVar*)(((RooWorkspace*)file->Get(wname.c_str()))->obj((varname+"_mean_CB").c_str()));
    RooFormulaVar var_sigma  = *(RooFormulaVar*)(((RooWorkspace*)file->Get(wname.c_str()))->obj((varname+"_sigma_CB").c_str()));
    RooFormulaVar var_alpha  = *(RooFormulaVar*)(((RooWorkspace*)file->Get(wname.c_str()))->obj((varname+"_alpha").c_str()));
    RooFormulaVar var_n      = *(RooFormulaVar*)(((RooWorkspace*)file->Get(wname.c_str()))->obj((varname+"_n").c_str()));
    RooFormulaVar ggh_norm   = *(RooFormulaVar*)(((RooWorkspace*)file->Get(wname.c_str()))->obj("sig_ggH_norm"));
    RooFormulaVar vbf_norm   = *(RooFormulaVar*)(((RooWorkspace*)file->Get(wname.c_str()))->obj("sig_VBF_norm"));

    RooRealVar mass("mass", "", xmin, xmin, xmax);


    RooRealVar mean_sig (("sig" + suffix + "_mean_sig" ).c_str(), "", var_mean.getVal()); 
    RooRealVar sigma_sig(("sig" + suffix + "_sigma_sig").c_str(), "", var_sigma.getVal()); 
    RooRealVar alpha    (("sig" + suffix + "_alpha"    ).c_str(), "", var_alpha.getVal()); 
    RooRealVar n        (("sig" + suffix + "_n"        ).c_str(), "", var_n.getVal()); 

    RooCBShape pdf("signalCB", "", mass, mean_sig,sigma_sig,alpha,n);


    TH1* hist = pdf.createHistogram(histname, mass, Binning(nbins, xmin, xmax));
    hist->Scale(1.0/hist->Integral());
    hist->Scale(ggh_norm.getVal()+vbf_norm.getVal());

    return hist;

}


TH1* getQQZZHist(const char* filename, int channel, const char* histname, int nbins, float xmin, float xmax, bool is7) {

    TFile* file = new TFile(filename);

    RooRealVar mass("mass", "", xmin, xmin, xmax);

    string suffix = "";
    string tevstr = is7 ? "_7TeV" : "_8TeV";
    if (channel == 0) suffix += "_4mu";
    if (channel == 1) suffix += "_4e";
    if (channel == 2) suffix += "_2e2mu";
    suffix += tevstr;    

    string wname = "w";
    string varname = "bkg_qqzz"+suffix;

    RooRealVar var_a0  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a0").c_str()));
    RooRealVar var_a1  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a1").c_str()));
    RooRealVar var_a2  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a2").c_str()));
    RooRealVar var_a3  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a3").c_str()));
    RooRealVar var_a4  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a4").c_str()));
    RooRealVar var_a5  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a5").c_str()));
    RooRealVar var_a6  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a6").c_str()));
    RooRealVar var_a7  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a7").c_str()));
    RooRealVar var_a8  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a8").c_str()));
    RooRealVar var_a9  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a9").c_str()));
    RooRealVar var_a10 = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a10").c_str()));
    RooRealVar var_a11 = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a11").c_str()));
    RooRealVar var_a12 = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a12").c_str()));
    RooRealVar var_a13 = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a13").c_str()));
    
    RooqqZZPdf_v2 pdf("pdf","",mass,
                       var_a0,
                       var_a1,
                       var_a2,
                       var_a3,
                       var_a4,
                       var_a5,
                       var_a6,
                       var_a7,
                       var_a8,
                       var_a9,   
                       var_a10,   
                       var_a11,   
                       var_a12,   
                       var_a13);   

    TH1* hist = pdf.createHistogram(histname, mass, Binning(nbins, xmin, xmax));

    return hist;
    

}

TH1* getGGZZHist(const char* filename, int channel, const char* histname, int nbins, float xmin, float xmax, bool is7) {

    TFile* file = new TFile(filename);

    RooRealVar mass("mass", "", xmin, xmin, xmax);

    string suffix = "";
    string tevstr = is7 ? "_7TeV" : "_8TeV";
    if (channel == 0) suffix += "_4mu";
    if (channel == 1) suffix += "_4e";
    if (channel == 2) suffix += "_2e2mu";
    suffix += tevstr;    

    string wname = "w";
    string varname = "bkg_ggzz"+suffix;

    RooRealVar var_a0  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a0").c_str()));
    RooRealVar var_a1  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a1").c_str()));
    RooRealVar var_a2  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a2").c_str()));
    RooRealVar var_a3  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a3").c_str()));
    RooRealVar var_a4  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a4").c_str()));
    RooRealVar var_a5  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a5").c_str()));
    RooRealVar var_a6  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a6").c_str()));
    RooRealVar var_a7  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a7").c_str()));
    RooRealVar var_a8  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a8").c_str()));
    RooRealVar var_a9  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_a9").c_str()));

    RooggZZPdf_v2 pdf("pdf","",mass,
                       var_a0,
                       var_a1,
                       var_a2,
                       var_a3,
                       var_a4,
                       var_a5,
                       var_a6,
                       var_a7,
                       var_a8,
                       var_a9);

    TH1* hist = pdf.createHistogram(histname, mass, Binning(nbins, xmin, xmax));

    return hist;


}

TH1* getZXSSHist(const char* filename, int channel, const char* histname, int nbins, float xmin, float xmax, bool is7) {

    TFile* file = new TFile(filename);

    RooRealVar mass("mass", "", xmin, xmin, xmax);

    string suffix = "";
    string tevstr = is7 ? "_7TeV" : "_8TeV";
    if (channel == 0) suffix += "_4mu";
    if (channel == 1) suffix += "_4e";
    if (channel == 2) suffix += "_2e2mu";
    suffix += tevstr;    

    string wname = "w";
    string varname = "bkg_zjets"+suffix;

    RooRealVar mean  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_mean_zx").c_str()));
    RooRealVar sigma = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_sigma_zx").c_str()));

    RooLandau pdf("pdf","",mass,mean,sigma);

    TH1* hist = pdf.createHistogram(histname, mass, Binning(nbins, xmin, xmax));

    return hist;

}

TH1* getDataHist(const char* filename, int channel, const char* histname, int nbins, float xmin, float xmax, DataYieldMaker& yd) {

        TH1F* hist = new TH1F("histname", "", nbins, xmin, xmax);
        hist->Sumw2();
        if (channel == 0) yd.get1DHist(0, 40., 12., xmin, xmax, -1., hist);
        if (channel == 1) yd.get1DHist(1, 40., 12., xmin, xmax, -1., hist);
        if (channel == 2) {
            yd.get1DHist(2, 40., 12., xmin, xmax, -1., hist); 
            yd.get1DHist(3, 40., 12., xmin, xmax, -1., hist); 
        }

        return hist;
}


void plotMass() {

    bool is7   = false;
    bool combine = false;
    //int nbins  = 27;
    //float xmin = 100.5;
    //float xmax = 181.5;
    int nbins  = 50;
    float xmin = 100.;
    float xmax = 600.;
    float lumi7 = 5.05;
    float lumi8 = 5.26;
    string base_folder7 = "/home/avartak/CMS/Higgs/HZZ4L/CMSSW_4_2_8_patch7/src/WWAnalysis/AnalysisStep/trees/";
    string base_folder8 = "/home/avartak/CMS/Higgs/HZZ4L/CMSSW_5_2_4_patch4/src/WWAnalysis/AnalysisStep/trees/";


    TH1* hist_hizz_4mu   = new TH1F("hist_hizz_4mu",   "", nbins, xmin, xmax);
    TH1* hist_qqzz_4mu   = new TH1F("hist_qqzz_4mu",   "", nbins, xmin, xmax);
    TH1* hist_ggzz_4mu   = new TH1F("hist_ggzz_4mu",   "", nbins, xmin, xmax);
    TH1* hist_zxss_4mu   = new TH1F("hist_zxss_4mu",   "", nbins, xmin, xmax);
    TH1* hist_data_4mu   = new TH1F("hist_data_4mu",   "", nbins, xmin, xmax);

    TH1* hist_hizz_4e    = new TH1F("hist_hizz_4e",    "", nbins, xmin, xmax);
    TH1* hist_qqzz_4e    = new TH1F("hist_qqzz_4e",    "", nbins, xmin, xmax);
    TH1* hist_ggzz_4e    = new TH1F("hist_ggzz_4e",    "", nbins, xmin, xmax);
    TH1* hist_zxss_4e    = new TH1F("hist_zxss_4e",    "", nbins, xmin, xmax);
    TH1* hist_data_4e    = new TH1F("hist_data_4e",    "", nbins, xmin, xmax);

    TH1* hist_hizz_2e2mu = new TH1F("hist_hizz_2e2mu", "", nbins, xmin, xmax);
    TH1* hist_qqzz_2e2mu = new TH1F("hist_qqzz_2e2mu", "", nbins, xmin, xmax);
    TH1* hist_ggzz_2e2mu = new TH1F("hist_ggzz_2e2mu", "", nbins, xmin, xmax);
    TH1* hist_zxss_2e2mu = new TH1F("hist_zxss_2e2mu", "", nbins, xmin, xmax);
    TH1* hist_data_2e2mu = new TH1F("hist_data_2e2mu", "", nbins, xmin, xmax);


    if (!is7 || combine) {
        DataYieldMaker ymaker_data;
        ZXYieldMaker   ymaker_zxss;
        ZZYieldMaker   ymaker_qqzz;
        ZZYieldMaker   ymaker_ggzz;

        FakeRateCalculator FR8(base_folder8+"hzzTree.root", false, 40, 120, 0.0, 0.0, true);
        ymaker_data.fill(base_folder8+"hzzTree.root");
        ymaker_zxss.fill(base_folder8+"hzzTree.root"       , 1.0, FR8, true);
        ymaker_qqzz.fill(base_folder8+"hzzTree_id102.root" , getBkgXsec(102)*lumi8/evt_8TeV[102], 0.0, false);
        ymaker_qqzz.fill(base_folder8+"hzzTree_id103.root" , getBkgXsec(103)*lumi8/evt_8TeV[103], 0.0, false);
        ymaker_qqzz.fill(base_folder8+"hzzTree_id104.root" , getBkgXsec(104)*lumi8/evt_8TeV[104], 0.0, false);
        ymaker_qqzz.fill(base_folder8+"hzzTree_id105.root" , getBkgXsec(105)*lumi8/evt_8TeV[105], 0.0, false);
        ymaker_qqzz.fill(base_folder8+"hzzTree_id106.root" , getBkgXsec(106)*lumi8/evt_8TeV[106], 0.0, false);
        ymaker_qqzz.fill(base_folder8+"hzzTree_id107.root" , getBkgXsec(107)*lumi8/evt_8TeV[107], 0.0, false);
        ymaker_ggzz.fill(base_folder8+"hzzTree_id101.root" , getBkgXsec(101)*lumi8/evt_8TeV[101], 0.0, false);
        ymaker_ggzz.fill(base_folder8+"hzzTree_id100.root" , getBkgXsec(100)*lumi8/evt_8TeV[100], 0.0, false);

        float yield_qqzz_4mu    = ymaker_qqzz.getYield(0, 40., 12., xmin, xmax, -1.);
        float yield_ggzz_4mu    = ymaker_ggzz.getYield(0, 40., 12., xmin, xmax, -1.);
        float yield_zxss_4mu    = ymaker_zxss.getYield(0, 40., 12., xmin, xmax, -1.);

        float yield_qqzz_4e     = ymaker_qqzz.getYield(1, 40., 12., xmin, xmax, -1.);
        float yield_ggzz_4e     = ymaker_ggzz.getYield(1, 40., 12., xmin, xmax, -1.);
        float yield_zxss_4e     = ymaker_zxss.getYield(1, 40., 12., xmin, xmax, -1.);

        float yield_qqzz_2e2mu  = ymaker_qqzz.getYield(2, 40., 12., xmin, xmax, -1.);
        float yield_ggzz_2e2mu  = ymaker_ggzz.getYield(2, 40., 12., xmin, xmax, -1.);
        float yield_zxss_2e2mu  = ymaker_zxss.getYield(2, 40., 12., xmin, xmax, -1.);

              yield_qqzz_2e2mu += ymaker_qqzz.getYield(3, 40., 12., xmin, xmax, -1.);
              yield_ggzz_2e2mu += ymaker_ggzz.getYield(3, 40., 12., xmin, xmax, -1.);
              yield_zxss_2e2mu += ymaker_zxss.getYield(3, 40., 12., xmin, xmax, -1.);

        hist_hizz_4mu   -> Add(getSigHist ("card_1D_m126_8TeV_4mu_workspace.root"  , 0, "hist_hizz_4mu_8TeV"  , nbins, xmin, xmax, false));
        hist_qqzz_4mu   -> Add(getQQZZHist("card_1D_m126_8TeV_4mu_workspace.root"  , 0, "hist_qqzz_4mu_8TeV"  , nbins, xmin, xmax, false), yield_qqzz_4mu);
        hist_ggzz_4mu   -> Add(getGGZZHist("card_1D_m126_8TeV_4mu_workspace.root"  , 0, "hist_ggzz_4mu_8TeV"  , nbins, xmin, xmax, false), yield_ggzz_4mu);
        hist_zxss_4mu   -> Add(getZXSSHist("card_1D_m126_8TeV_4mu_workspace.root"  , 0, "hist_zxss_4mu_8TeV"  , nbins, xmin, xmax, false), yield_zxss_4mu);
        hist_data_4mu   -> Add(getDataHist("card_1D_m126_8TeV_4mu_workspace.root"  , 0, "hist_data_4mu_8TeV"  , nbins, xmin, xmax, ymaker_data));
        
        hist_hizz_4e    -> Add(getSigHist ("card_1D_m126_8TeV_4e_workspace.root"   , 1, "hist_hizz_4e_8TeV"   , nbins, xmin, xmax, false));
        hist_qqzz_4e    -> Add(getQQZZHist("card_1D_m126_8TeV_4e_workspace.root"   , 1, "hist_qqzz_4e_8TeV"   , nbins, xmin, xmax, false), yield_qqzz_4e);
        hist_ggzz_4e    -> Add(getGGZZHist("card_1D_m126_8TeV_4e_workspace.root"   , 1, "hist_ggzz_4e_8TeV"   , nbins, xmin, xmax, false), yield_ggzz_4e);
        hist_zxss_4e    -> Add(getZXSSHist("card_1D_m126_8TeV_4e_workspace.root"   , 1, "hist_zxss_4e_8TeV"   , nbins, xmin, xmax, false), yield_zxss_4e);
        hist_data_4e    -> Add(getDataHist("card_1D_m126_8TeV_4e_workspace.root"   , 1, "hist_data_4e_8TeV"   , nbins, xmin, xmax, ymaker_data));
        
        hist_hizz_2e2mu -> Add(getSigHist ("card_1D_m126_8TeV_2e2mu_workspace.root", 2, "hist_hizz_2e2mu_8TeV", nbins, xmin, xmax, false));
        hist_qqzz_2e2mu -> Add(getQQZZHist("card_1D_m126_8TeV_2e2mu_workspace.root", 2, "hist_qqzz_2e2mu_8TeV", nbins, xmin, xmax, false), yield_qqzz_2e2mu);
        hist_ggzz_2e2mu -> Add(getGGZZHist("card_1D_m126_8TeV_2e2mu_workspace.root", 2, "hist_ggzz_2e2mu_8TeV", nbins, xmin, xmax, false), yield_ggzz_2e2mu);
        hist_zxss_2e2mu -> Add(getZXSSHist("card_1D_m126_8TeV_2e2mu_workspace.root", 2, "hist_zxss_2e2mu_8TeV", nbins, xmin, xmax, false), yield_zxss_2e2mu);
        hist_data_2e2mu -> Add(getDataHist("card_1D_m126_8TeV_2e2mu_workspace.root", 2, "hist_data_2e2mu_8TeV", nbins, xmin, xmax, ymaker_data));


    }
    if (is7 || combine) {
        DataYieldMaker ymaker_data;
        ZXYieldMaker   ymaker_zxss;
        ZZYieldMaker   ymaker_qqzz;
        ZZYieldMaker   ymaker_ggzz;

        FakeRateCalculator FR7(base_folder7+"hzzTree.root", true , 40, 120, 0.0, 0.0, true);
        ymaker_data.fill(base_folder7+"hzzTree.root");
        ymaker_zxss.fill(base_folder7+"hzzTree.root"       , 1.0, FR7, true);
        ymaker_qqzz.fill(base_folder7+"hzzTree_id121.root" , getBkgXsec(121)*lumi7/evt_7TeV[121], 0.0, false);
        ymaker_qqzz.fill(base_folder7+"hzzTree_id122.root" , getBkgXsec(122)*lumi7/evt_7TeV[122], 0.0, false);
        ymaker_qqzz.fill(base_folder7+"hzzTree_id123.root" , getBkgXsec(123)*lumi7/evt_7TeV[123], 0.0, false);
        ymaker_qqzz.fill(base_folder7+"hzzTree_id124.root" , getBkgXsec(124)*lumi7/evt_7TeV[124], 0.0, false);
        ymaker_qqzz.fill(base_folder7+"hzzTree_id125.root" , getBkgXsec(125)*lumi7/evt_7TeV[125], 0.0, false);
        ymaker_qqzz.fill(base_folder7+"hzzTree_id126.root" , getBkgXsec(126)*lumi7/evt_7TeV[126], 0.0, false);
        ymaker_ggzz.fill(base_folder7+"hzzTree_id101.root" , getBkgXsec(101)*lumi7/evt_7TeV[101], 0.0, false);
        ymaker_ggzz.fill(base_folder7+"hzzTree_id100.root" , getBkgXsec(100)*lumi7/evt_7TeV[100], 0.0, false);

        float yield_qqzz_4mu    = ymaker_qqzz.getYield(0, 40., 12., xmin, xmax, -1.);
        float yield_ggzz_4mu    = ymaker_ggzz.getYield(0, 40., 12., xmin, xmax, -1.);
        float yield_zxss_4mu    = ymaker_zxss.getYield(0, 40., 12., xmin, xmax, -1.);
        
        float yield_qqzz_4e     = ymaker_qqzz.getYield(1, 40., 12., xmin, xmax, -1.);
        float yield_ggzz_4e     = ymaker_ggzz.getYield(1, 40., 12., xmin, xmax, -1.);
        float yield_zxss_4e     = ymaker_zxss.getYield(1, 40., 12., xmin, xmax, -1.);
        
        float yield_qqzz_2e2mu  = ymaker_qqzz.getYield(2, 40., 12., xmin, xmax, -1.);
        float yield_ggzz_2e2mu  = ymaker_ggzz.getYield(2, 40., 12., xmin, xmax, -1.);
        float yield_zxss_2e2mu  = ymaker_zxss.getYield(2, 40., 12., xmin, xmax, -1.);
        
              yield_qqzz_2e2mu += ymaker_qqzz.getYield(3, 40., 12., xmin, xmax, -1.);
              yield_ggzz_2e2mu += ymaker_ggzz.getYield(3, 40., 12., xmin, xmax, -1.);
              yield_zxss_2e2mu += ymaker_zxss.getYield(3, 40., 12., xmin, xmax, -1.);


        hist_hizz_4mu   -> Add(getSigHist ("card_1D_m126_7TeV_4mu_workspace.root"  , 0, "hist_hizz_4mu_7TeV"  , nbins, xmin, xmax, true));
        hist_qqzz_4mu   -> Add(getQQZZHist("card_1D_m126_7TeV_4mu_workspace.root"  , 0, "hist_qqzz_4mu_7TeV"  , nbins, xmin, xmax, true), yield_qqzz_4mu);
        hist_ggzz_4mu   -> Add(getGGZZHist("card_1D_m126_7TeV_4mu_workspace.root"  , 0, "hist_ggzz_4mu_7TeV"  , nbins, xmin, xmax, true), yield_ggzz_4mu);
        hist_zxss_4mu   -> Add(getZXSSHist("card_1D_m126_7TeV_4mu_workspace.root"  , 0, "hist_zxss_4mu_7TeV"  , nbins, xmin, xmax, true), yield_zxss_4mu);
        hist_data_4mu   -> Add(getDataHist("card_1D_m126_7TeV_4mu_workspace.root"  , 0, "hist_data_4mu_7TeV"  , nbins, xmin, xmax, ymaker_data));
        
        hist_hizz_4e    -> Add(getSigHist ("card_1D_m126_7TeV_4e_workspace.root"   , 1, "hist_hizz_4e_7TeV"   , nbins, xmin, xmax, true));
        hist_qqzz_4e    -> Add(getQQZZHist("card_1D_m126_7TeV_4e_workspace.root"   , 1, "hist_qqzz_4e_7TeV"   , nbins, xmin, xmax, true), yield_qqzz_4e);
        hist_ggzz_4e    -> Add(getGGZZHist("card_1D_m126_7TeV_4e_workspace.root"   , 1, "hist_ggzz_4e_7TeV"   , nbins, xmin, xmax, true), yield_ggzz_4e);
        hist_zxss_4e    -> Add(getZXSSHist("card_1D_m126_7TeV_4e_workspace.root"   , 1, "hist_zxss_4e_7TeV"   , nbins, xmin, xmax, true), yield_zxss_4e);
        hist_data_4e    -> Add(getDataHist("card_1D_m126_7TeV_4e_workspace.root"   , 1, "hist_data_4e_7TeV"   , nbins, xmin, xmax, ymaker_data));
        
        hist_hizz_2e2mu -> Add(getSigHist ("card_1D_m126_7TeV_2e2mu_workspace.root", 2, "hist_hizz_2e2mu_7TeV", nbins, xmin, xmax, true));
        hist_qqzz_2e2mu -> Add(getQQZZHist("card_1D_m126_7TeV_2e2mu_workspace.root", 2, "hist_qqzz_2e2mu_7TeV", nbins, xmin, xmax, true), yield_qqzz_2e2mu);
        hist_ggzz_2e2mu -> Add(getGGZZHist("card_1D_m126_7TeV_2e2mu_workspace.root", 2, "hist_ggzz_2e2mu_7TeV", nbins, xmin, xmax, true), yield_ggzz_2e2mu);
        hist_zxss_2e2mu -> Add(getZXSSHist("card_1D_m126_7TeV_2e2mu_workspace.root", 2, "hist_zxss_2e2mu_7TeV", nbins, xmin, xmax, true), yield_zxss_2e2mu);
        hist_data_2e2mu -> Add(getDataHist("card_1D_m126_7TeV_2e2mu_workspace.root", 2, "hist_data_2e2mu_7TeV", nbins, xmin, xmax, ymaker_data));

    }

    hist_qqzz_4mu->Add(hist_qqzz_4e);
    hist_qqzz_4mu->Add(hist_qqzz_2e2mu);

    hist_ggzz_4mu->Add(hist_ggzz_4e);
    hist_ggzz_4mu->Add(hist_ggzz_2e2mu);

    hist_qqzz_4mu->Add(hist_ggzz_4mu);

    hist_hizz_4mu->Add(hist_hizz_4e);
    hist_hizz_4mu->Add(hist_hizz_2e2mu);

    hist_zxss_4mu->Add(hist_zxss_4e);
    hist_zxss_4mu->Add(hist_zxss_2e2mu);

    hist_data_4mu->Add(hist_data_4e);
    hist_data_4mu->Add(hist_data_2e2mu);

    hist_data_4mu->SetMarkerStyle(20);
    hist_data_4mu->SetMarkerSize(1.1);

    hist_qqzz_4mu->SetFillColor(kAzure-9);
    hist_zxss_4mu->SetFillColor(kGreen-5);
    hist_hizz_4mu->SetLineColor(kRed);

    hist_qqzz_4mu->SetLineWidth(2);
    hist_zxss_4mu->SetLineWidth(2);
    hist_hizz_4mu->SetLineWidth(2);
    hist_hizz_4mu->SetMarkerSize(0);

    THStack* hs = new THStack("hs", "");
    hs->Add(hist_zxss_4mu);
    hs->Add(hist_qqzz_4mu);
    hs->Add(hist_hizz_4mu);

    TLegend *leg = new TLegend(0.58,0.62,0.8,0.93);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hist_data_4mu, "Data");
    leg->AddEntry(hist_qqzz_4mu, "ZZ",  "F");
    leg->AddEntry(hist_zxss_4mu, "Z+X", "F");
    leg->AddEntry(hist_hizz_4mu, "m_{H} = 126 GeV");


    const Int_t n = hist_data_4mu->GetNbinsX();
    Double_t x[n],y[n],exl[n],eyl[n],exh[n],eyh[n];
    
    double q = (1-0.6827)/2.;
    
    for (Int_t i = 0; i < n; i++) {
        x[i] = hist_data_4mu->GetXaxis()->GetBinCenter(i+1);
        double N = hist_data_4mu->GetBinContent(i+1);
        
        exl[i] = 0; exh[i] = 0;
        
        if( N == 0 ){eyl[i] = 0; eyh[i] = 0; y[i] = -1;continue;}
        
        y[i] = N;
        eyl[i] = (N==0)?0:(N-ROOT::Math::chisquared_quantile_c(1-q,2*N)/2.);
        eyh[i] = ROOT::Math::chisquared_quantile_c(q,2*(N+1))/2.-N;
    }
    
    TGraphAsymmErrors *gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.1);
    gr->SetMarkerColor(kBlack);
    gr->SetLineWidth(1);


    hs->Draw();
    gr->Draw("SAME P");
    leg->Draw("SAME");

    TLatex* CP = new TLatex(99.,6.59, is7 ? "CMS Preliminary                          #sqrt{s} = 7 TeV, L = 5.05 fb^{-1}" : "CMS Preliminary                          #sqrt{s} = 8 TeV, L = 5.26 fb^{-1}");
    if (combine) CP->SetText(99., 6.59, "CMS Preliminary, #sqrt{s} = 7 TeV, L = 5.05 fb^{-1} & #sqrt{s} = 8 TeV, L = 5.26 fb^{-1}");
    CP->SetTextSize(0.035);
    CP->Draw("SAME");

    std::cout << hist_qqzz_4mu->Integral() << std::endl;
    std::cout << hist_zxss_4mu->Integral() << std::endl;
    std::cout << hist_hizz_4mu->Integral() << std::endl;
    std::cout << hist_data_4mu->Integral() << std::endl;
}

