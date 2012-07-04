#include <iostream>
#include <TH1F.h>
#include <THStack.h>
#include <TFile.h>
#include <TLegend.h>
#include <Math/Math.h>
#include <Math/QuantFuncMathCore.h>
#include "WWAnalysis/AnalysisStep/interface/scales.h"
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

TH1* getQQZZHist(const char* filename, int channel, const char* histname, int nbins, float xmin, float xmax) {

    TFile* file = new TFile(filename);

    RooRealVar mass("mass", "", xmin, xmin, xmax);

    string suffix = "";
    if (channel == 0) suffix += "_4mu";
    if (channel == 1) suffix += "_4e";
    if (channel == 2) suffix += "_2e2mu";

    string wname = "w"+suffix;
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

    TH1* hist = pdf.createHistogram(histname, mass, Binning(50, 100., 600.));

    return hist;
    

}

TH1* getGGZZHist(const char* filename, int channel, const char* histname, int nbins, float xmin, float xmax) {

    TFile* file = new TFile(filename);

    RooRealVar mass("mass", "", xmin, xmin, xmax);

    string suffix = "";
    if (channel == 0) suffix += "_4mu";
    if (channel == 1) suffix += "_4e";
    if (channel == 2) suffix += "_2e2mu";

    string wname = "w"+suffix;
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

    TH1* hist = pdf.createHistogram(histname, mass, Binning(50, 100., 600.));

    return hist;


}

TH1* getZXSSHist(const char* filename, int channel, const char* histname, int nbins, float xmin, float xmax) {

    TFile* file = new TFile(filename);

    RooRealVar mass("mass", "", xmin, xmin, xmax);

    string suffix = "";
    if (channel == 0) suffix += "_4mu";
    if (channel == 1) suffix += "_4e";
    if (channel == 2) suffix += "_2e2mu";

    string wname = "w"+suffix;
    string varname = "bkg_zjets"+suffix;

    RooRealVar mean  = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_mean_zx").c_str()));
    RooRealVar sigma = *(((RooWorkspace*)file->Get(wname.c_str()))->var((varname+"_sigma_zx").c_str()));

    RooLandau pdf("pdf","",mass,mean,sigma);

    TH1* hist = pdf.createHistogram(histname, mass, Binning(50, 100., 600.));

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

    init();

    string base_folder = "/home/avartak/CMS/Higgs/HZZ4L/CMSSW_5_2_4_patch4/src/WWAnalysis/AnalysisStep/trees/7TeV/";

    DataYieldMaker ymaker_data;
    ZXYieldMaker   ymaker_zxss;
    ZZYieldMaker   ymaker_qqzz;
    ZZYieldMaker   ymaker_ggzz;
    //ZZYieldMaker   ymaker_ghzz;
    //ZZYieldMaker   ymaker_qhzz;

    FakeRateCalculator FR(base_folder+"hzzTree.root", true, 40, 120, 0.0, 0.0, true);

    ymaker_data.fill(base_folder+"hzzTree.root");
    ymaker_zxss.fill(base_folder+"hzzTree.root"       , 1.0, FR, true);
    ymaker_qqzz.fill(base_folder+"hzzTree_id121.root" , xsecweights[121]*5.05, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id122.root" , xsecweights[122]*5.05, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id123.root" , xsecweights[123]*5.05, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id124.root" , xsecweights[124]*5.05, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id125.root" , xsecweights[125]*5.05, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id126.root" , xsecweights[126]*5.05, 0.0, false);
    ymaker_ggzz.fill(base_folder+"hzzTree_id101.root" , xsecweights[101]*5.05, 0.0, false);
    ymaker_ggzz.fill(base_folder+"hzzTree_id100.root" , xsecweights[100]*5.05, 0.0, false);
    //ymaker_ghzz.fill(base_folder+ggh_rootfile         , xsecweights[hinfo.gghid]*hinfo.lumi, 0.0, true );
    //ymaker_qhzz.fill(base_folder+vbf_rootfile         , xsecweights[hinfo.vbfid]*hinfo.lumi, 0.0, true );

    TH1* hist_qqzz_4mu   = getQQZZHist("hzz_m120_4mu_workspace.root",   0, "hist_qqzz_4mu",   50, 100., 600.);
    TH1* hist_ggzz_4mu   = getGGZZHist("hzz_m120_4mu_workspace.root",   0, "hist_ggzz_4mu",   50, 100., 600.);
    TH1* hist_zxss_4mu   = getZXSSHist("hzz_m120_4mu_workspace.root",   0, "hist_zxss_4mu",   50, 100., 600.);
    TH1* hist_data_4mu   = getDataHist("hzz_m120_4mu_workspace.root",   0, "hist_data_4mu",   50, 100., 600., ymaker_data);

    TH1* hist_qqzz_4e    = getQQZZHist("hzz_m120_4e_workspace.root",    1, "hist_qqzz_4e",    50, 100., 600.);
    TH1* hist_ggzz_4e    = getGGZZHist("hzz_m120_4e_workspace.root",    1, "hist_ggzz_4e",    50, 100., 600.);
    TH1* hist_zxss_4e    = getZXSSHist("hzz_m120_4e_workspace.root",    1, "hist_zxss_4e",    50, 100., 600.);
    TH1* hist_data_4e    = getDataHist("hzz_m120_4e_workspace.root",    1, "hist_data_4e",    50, 100., 600., ymaker_data);

    TH1* hist_qqzz_2e2mu = getQQZZHist("hzz_m120_2e2mu_workspace.root", 2, "hist_qqzz_2e2mu", 50, 100., 600.);
    TH1* hist_ggzz_2e2mu = getGGZZHist("hzz_m120_2e2mu_workspace.root", 2, "hist_ggzz_2e2mu", 50, 100., 600.);
    TH1* hist_zxss_2e2mu = getZXSSHist("hzz_m120_2e2mu_workspace.root", 2, "hist_zxss_2e2mu", 50, 100., 600.);
    TH1* hist_data_2e2mu = getDataHist("hzz_m120_2e2mu_workspace.root", 2, "hist_data_2e2mu", 50, 100., 600., ymaker_data);

    float yield_qqzz_4mu    = ymaker_qqzz.getYield(0, 40., 12., 100., 600., -1.);
    float yield_ggzz_4mu    = ymaker_ggzz.getYield(0, 40., 12., 100., 600., -1.);
    float yield_zxss_4mu    = ymaker_zxss.getYield(0, 40., 12., 100., 600., -1.);

    float yield_qqzz_4e     = ymaker_qqzz.getYield(1, 40., 12., 100., 600., -1.);
    float yield_ggzz_4e     = ymaker_ggzz.getYield(1, 40., 12., 100., 600., -1.);
    float yield_zxss_4e     = ymaker_zxss.getYield(1, 40., 12., 100., 600., -1.);

    float yield_qqzz_2e2mu  = ymaker_qqzz.getYield(2, 40., 12., 100., 600., -1.);
    float yield_ggzz_2e2mu  = ymaker_ggzz.getYield(2, 40., 12., 100., 600., -1.);
    float yield_zxss_2e2mu  = ymaker_zxss.getYield(2, 40., 12., 100., 600., -1.);

          yield_qqzz_2e2mu += ymaker_qqzz.getYield(3, 40., 12., 100., 600., -1.);
          yield_ggzz_2e2mu += ymaker_ggzz.getYield(3, 40., 12., 100., 600., -1.);
          yield_zxss_2e2mu += ymaker_zxss.getYield(3, 40., 12., 100., 600., -1.);

    hist_qqzz_4mu->Scale(yield_qqzz_4mu);
    hist_ggzz_4mu->Scale(yield_ggzz_4mu);
    hist_zxss_4mu->Scale(yield_zxss_4mu);

    hist_qqzz_4e->Scale(yield_qqzz_4e);
    hist_ggzz_4e->Scale(yield_ggzz_4e);
    hist_zxss_4e->Scale(yield_zxss_4e);

    hist_qqzz_2e2mu->Scale(yield_qqzz_2e2mu);
    hist_ggzz_2e2mu->Scale(yield_ggzz_2e2mu);
    hist_zxss_2e2mu->Scale(yield_zxss_2e2mu);

    hist_qqzz_4mu->Add(hist_qqzz_4e);
    hist_qqzz_4mu->Add(hist_qqzz_2e2mu);

    hist_ggzz_4mu->Add(hist_ggzz_4e);
    hist_ggzz_4mu->Add(hist_ggzz_2e2mu);

    hist_qqzz_4mu->Add(hist_ggzz_4mu);

    hist_zxss_4mu->Add(hist_zxss_4e);
    hist_zxss_4mu->Add(hist_zxss_2e2mu);

    hist_data_4mu->Add(hist_data_4e);
    hist_data_4mu->Add(hist_data_2e2mu);

    hist_data_4mu->SetMarkerStyle(20);
    hist_data_4mu->SetMarkerSize(1.1);

    hist_qqzz_4mu->SetFillColor(kAzure-9);
    hist_zxss_4mu->SetFillColor(kGreen-5);

    hist_qqzz_4mu->SetLineWidth(2);
    hist_zxss_4mu->SetLineWidth(2);

    THStack* hs = new THStack("hs", "");
    hs->Add(hist_zxss_4mu);
    hs->Add(hist_qqzz_4mu);

    TLegend *leg = new TLegend(0.58,0.62,0.8,0.93);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(hist_data_4mu, "Data");
    leg->AddEntry(hist_qqzz_4mu, "ZZ",  "F");
    leg->AddEntry(hist_zxss_4mu, "Z+X", "F");


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

    std::cout << hist_qqzz_4mu->Integral() << std::endl;
    std::cout << hist_zxss_4mu->Integral() << std::endl;
    std::cout << hist_data_4mu->Integral() << std::endl;
}

