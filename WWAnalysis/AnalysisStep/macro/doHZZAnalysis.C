#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TBranch.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <utility>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h> 
#include <RooFormulaVar.h> 
#include <RooWorkspace.h> 
#include <RooLandau.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h> 
#include <RooFFTConvPdf.h>
#include <RooProdPdf.h>

#include "FakeRateCalculator.h"
#include "CardTemplate.h"
#include "FitMaker.h"
#include "YieldMaker.h"
#include "findAndReplace.h"

using namespace RooFit;

struct SignalCardInfo {

    float mean_CB;
    float sigma_CB;
    float alpha;
    float n;
    float mean_BW;
    float gamma_BW;
    float higgsmass;
    float yield;

    SignalCardInfo() {
        mean_CB   =  0.0;
        sigma_CB  =  0.0;
        alpha     =  0.0;
        n         =  0.0;
        mean_BW   =  0.0;
        gamma_BW  =  0.0;
        higgsmass =  0.0;
        yield     =  0.0;
    }

    SignalCardInfo(SignalFitMaker sfm) {
        mean_CB   =  sfm.getVarMeanCB()  ;
        sigma_CB  =  sfm.getVarSigmaCB() ;
        alpha     =  sfm.getVarAlphaCB() ;
        n         =  sfm.getVarNCB()     ;
        mean_BW   =  sfm.getVarMeanBW()  ;
        gamma_BW  =  sfm.getVarGammaBW() ;
        higgsmass =  sfm.getVarMass()    ;
        yield     =  0.0;    
    }

    SignalCardInfo(SignalFitMaker sfm, float y) {
        mean_CB   =  sfm.getVarMeanCB()  ;
        sigma_CB  =  sfm.getVarSigmaCB() ;
        alpha     =  sfm.getVarAlphaCB() ;
        n         =  sfm.getVarNCB()     ;
        mean_BW   =  sfm.getVarMeanBW()  ;
        gamma_BW  =  sfm.getVarGammaBW() ;
        higgsmass =  sfm.getVarMass()    ;
        yield     =  y                   ;
    }    

};

struct QQZZCardInfo {

    float a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13;

    QQZZCardInfo() {
        a0  = 0.0;
        a1  = 0.0;
        a2  = 0.0;
        a3  = 0.0;
        a4  = 0.0;
        a5  = 0.0;
        a6  = 0.0;
        a7  = 0.0;
        a8  = 0.0;
        a9  = 0.0;
        a10 = 0.0;
        a11 = 0.0;
        a12 = 0.0;
        a13 = 0.0;
    }


    QQZZCardInfo(QQZZFitMaker fm) {
        a0  = fm.getVarA0();
        a1  = fm.getVarA1();
        a2  = fm.getVarA2();
        a3  = fm.getVarA3();
        a4  = fm.getVarA4();
        a5  = fm.getVarA5();
        a6  = fm.getVarA6();
        a7  = fm.getVarA7();
        a8  = fm.getVarA8();
        a9  = fm.getVarA9();
        a10 = fm.getVarA10();
        a11 = fm.getVarA11();
        a12 = fm.getVarA12();
        a13 = fm.getVarA13();
    }

};

struct GGZZCardInfo {

    float a0, a1, a2, a3, a4, a5, a6, a7, a8, a9;

    GGZZCardInfo() {
        a0  = 0.0;
        a1  = 0.0;
        a2  = 0.0;
        a3  = 0.0;
        a4  = 0.0;
        a5  = 0.0;
        a6  = 0.0;
        a7  = 0.0;
        a8  = 0.0;
        a9  = 0.0;
    }


    GGZZCardInfo(GGZZFitMaker fm) {
        a0  = fm.getVarA0();
        a1  = fm.getVarA1();
        a2  = fm.getVarA2();
        a3  = fm.getVarA3();
        a4  = fm.getVarA4();
        a5  = fm.getVarA5();
        a6  = fm.getVarA6();
        a7  = fm.getVarA7();
        a8  = fm.getVarA8();
        a9  = fm.getVarA9();
    }

};

struct ZXCardInfo {

    float mean, sigma;

    ZXCardInfo() {
        mean   = 0.0;
        sigma  = 0.0;
    }

    ZXCardInfo(ZXFitMaker fm) {
        mean   = fm.getVarMean()  ;
        sigma  = fm.getVarSigma() ;
    }

};

struct HiggsMassPointInfo {

    float lumi;
    float melacut;
    float z1min;
    float z2min;
    float massLowBkgFit;
    float massHighBkgFit;
    bool  do1D;
    bool  doSS;
    bool  do7TeV;
    std::string treeFolder;
    std::string melafilename;
    
    std::vector<SignalCardInfo> signalpoints_2e2mu;
    std::vector<SignalCardInfo> signalpoints_4e;
    std::vector<SignalCardInfo> signalpoints_4mu;
    
    QQZZCardInfo qqzzbkg_2e2mu;
    QQZZCardInfo qqzzbkg_4e;
    QQZZCardInfo qqzzbkg_4mu;
    
    GGZZCardInfo ggzzbkg_2e2mu;
    GGZZCardInfo ggzzbkg_4e;
    GGZZCardInfo ggzzbkg_4mu;
    
    ZXCardInfo zxbkg_2e2mu;
    ZXCardInfo zxbkg_4e;
    ZXCardInfo zxbkg_4mu;


    float CBMeanP0_2e2mu;
    float CBMeanP1_2e2mu;
    float CBSigmaP0_2e2mu;
    float CBSigmaP1_2e2mu;
    float CBAlphaP0_2e2mu;
    float CBNP0_2e2mu;
    float EfficiencyP0_2e2mu;
    float EfficiencyP1_2e2mu;
    float EfficiencyP2_2e2mu;
    float EfficiencyP3_2e2mu;

    float CBMeanP0_4e;
    float CBMeanP1_4e;
    float CBSigmaP0_4e;
    float CBSigmaP1_4e;
    float CBAlphaP0_4e;
    float CBNP0_4e;
    float EfficiencyP0_4e;
    float EfficiencyP1_4e;
    float EfficiencyP2_4e;
    float EfficiencyP3_4e;

    float CBMeanP0_4mu;
    float CBMeanP1_4mu;
    float CBSigmaP0_4mu;
    float CBSigmaP1_4mu;
    float CBAlphaP0_4mu;
    float CBNP0_4mu;
    float EfficiencyP0_4mu;
    float EfficiencyP1_4mu;
    float EfficiencyP2_4mu;
    float EfficiencyP3_4mu;

    FakeRateCalculator FR;

    DataYieldMaker ymaker_data;
    ZXYieldMaker   ymaker_zxss;
    ZZYieldMaker   ymaker_qqzz;
    ZZYieldMaker   ymaker_ggzz;

    float getSignalCBMean(float m, int ch) {
        if      (ch == 0) return CBMeanP0_4mu   + m*CBMeanP1_4mu;
        else if (ch == 1) return CBMeanP0_4e    + m*CBMeanP1_4e;
        else if (ch == 2) return CBMeanP0_2e2mu + m*CBMeanP1_2e2mu;
        else return 0.0;
    }

    float getSignalCBSigma(float m, int ch) {
        if      (ch == 0) return CBSigmaP0_4mu   + m*CBSigmaP1_4mu;
        else if (ch == 1) return CBSigmaP0_4e    + m*CBSigmaP1_4e;
        else if (ch == 2) return CBSigmaP0_2e2mu + m*CBSigmaP1_2e2mu;
        else return 0.0;
    }

    float getSignalCBAlpha(float m, int ch) {
        if      (ch == 0) return CBAlphaP0_4mu   ;
        else if (ch == 1) return CBAlphaP0_4e    ;
        else if (ch == 2) return CBAlphaP0_2e2mu ;
        else return 0.0;
    }

    float getSignalCBN(float m, int ch) {
        if      (ch == 0) return CBNP0_4mu   ;
        else if (ch == 1) return CBNP0_4e    ;
        else if (ch == 2) return CBNP0_2e2mu ;
        else return 0.0;
    }

    std::string getSignalCBMeanString(float m, int ch) {
        stringstream ss;
        if      (ch == 0) ss << CBMeanP0_4mu   << " + @0*" << CBMeanP1_4mu    << " + @0*@1";
        else if (ch == 1) ss << CBMeanP0_4e    << " + @0*" << CBMeanP1_4e     << " + @0*@1";
        else if (ch == 2) ss << CBMeanP0_2e2mu << " + @0*" << CBMeanP1_2e2mu  << " + @0*@1";
        return ss.str();
    }

    std::string getSignalCBSigmaString(float m, int ch) {
        stringstream ss;
        if      (ch == 0) ss << "(" << CBSigmaP0_4mu   << " + @0*" << CBSigmaP1_4mu    << ") * (1+@1)";
        else if (ch == 1) ss << "(" << CBSigmaP0_4e    << " + @0*" << CBSigmaP1_4e     << ") * (1+@1)";
        else if (ch == 2) ss << "(" << CBSigmaP0_2e2mu << " + @0*" << CBSigmaP1_2e2mu  << ") * (1+@1)";
        return ss.str();
    }

    float getEfficiency(float m, int ch) {
        if      (ch == 0) return EfficiencyP0_4mu   + EfficiencyP1_4mu*m   + EfficiencyP2_4mu*m*m   + EfficiencyP3_4mu*m*m*m   ;
        else if (ch == 1) return EfficiencyP0_4e    + EfficiencyP1_4e*m    + EfficiencyP2_4e*m*m    + EfficiencyP3_4e*m*m*m    ;
        else if (ch == 2) return EfficiencyP0_2e2mu + EfficiencyP1_2e2mu*m + EfficiencyP2_2e2mu*m*m + EfficiencyP3_2e2mu*m*m*m ;
        else return 0.0;
    }

    void printFit(TGraph gr, TF1 f, std::string filename) {
        TCanvas c("c", "", 400, 400);
        gr.Draw("LAP");
        f.Draw("SAME");
        c.Print(filename.c_str());
    }

    void fitCBMean() {
        if (signalpoints_4mu.size() == 0) {
            CBMeanP0_4mu = 0.0;
            CBMeanP1_4mu = 0.0;
        }
        else {
            const int n = signalpoints_4mu.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_4mu[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_4mu[i].mean_CB;
            TGraph gr(n, x, y);
            TF1 f("f", "[0]+[1]*x", 0, 1000.);
            gr.Fit(&f);
            CBMeanP0_4mu = f.GetParameter(0);
            CBMeanP1_4mu = f.GetParameter(1);
            printFit(gr, f, "CBMean_4mu.pdf");
        }

        if (signalpoints_4e.size() == 0) {
            CBMeanP0_4e = 0.0;
            CBMeanP1_4e = 0.0;
        }
        else {
            const int n = signalpoints_4e.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_4e[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_4e[i].mean_CB;
            TGraph gr(n, x, y);
            TF1 f("f", "[0]+[1]*x", 0, 1000.);
            gr.Fit(&f);
            CBMeanP0_4e = f.GetParameter(0);
            CBMeanP1_4e = f.GetParameter(1);
            printFit(gr, f, "CBMean_4e.pdf");
        }

        if (signalpoints_2e2mu.size() == 0) {
            CBMeanP0_2e2mu = 0.0;
            CBMeanP1_2e2mu = 0.0;
        }
        else {
            const int n = signalpoints_2e2mu.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_2e2mu[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_2e2mu[i].mean_CB;
            TGraph gr(n, x, y);
            TF1 f("f", "[0]+[1]*x", 0, 1000.);
            gr.Fit(&f);
            CBMeanP0_2e2mu = f.GetParameter(0);
            CBMeanP1_2e2mu = f.GetParameter(1);
            printFit(gr, f, "CBMean_2e2mu.pdf");
        }

    }


    void fitCBSigma() {
        if (signalpoints_4mu.size() == 0) {
            CBSigmaP0_4mu = 0.0;
            CBSigmaP1_4mu = 0.0;
        }
        else {
            const int n = signalpoints_4mu.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_4mu[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_4mu[i].sigma_CB;
            TGraph gr(n, x, y);
            TF1 f("f", "[0]+[1]*x", 0, 1000.);
            gr.Fit(&f);
            CBSigmaP0_4mu = f.GetParameter(0);
            CBSigmaP1_4mu = f.GetParameter(1);
            printFit(gr, f, "CBSigma_4mu.pdf");
        }

        if (signalpoints_4e.size() == 0) {
            CBSigmaP0_4e = 0.0;
            CBSigmaP1_4e = 0.0;
        }
        else {
            const int n = signalpoints_4e.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_4e[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_4e[i].sigma_CB;
            TGraph gr(n, x, y);
            TF1 f("f", "[0]+[1]*x", 0, 1000.);
            gr.Fit(&f);
            CBSigmaP0_4e = f.GetParameter(0);
            CBSigmaP1_4e = f.GetParameter(1);
            printFit(gr, f, "CBSigma_4e.pdf");
        }

        if (signalpoints_2e2mu.size() == 0) {
            CBSigmaP0_2e2mu = 0.0;
            CBSigmaP1_2e2mu = 0.0;
        }
        else {
            const int n = signalpoints_2e2mu.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_2e2mu[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_2e2mu[i].sigma_CB;
            TGraph gr(n, x, y);
            TF1 f("f", "[0]+[1]*x", 0, 1000.);
            gr.Fit(&f);
            CBSigmaP0_2e2mu = f.GetParameter(0);
            CBSigmaP1_2e2mu = f.GetParameter(1);
            printFit(gr, f, "CBSigma_2e2mu.pdf");
        }

    }

    void fitCBAlpha() {
        if (signalpoints_4mu.size() == 0) {
            CBAlphaP0_4mu = 0.0;
        }
        else {
            const int n = signalpoints_4mu.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_4mu[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_4mu[i].alpha;
            TGraph gr(n, x, y);
            TF1 f("f", "[0]", 0, 1000.);
            gr.Fit(&f);
            CBAlphaP0_4mu = f.GetParameter(0);
            printFit(gr, f, "CBAlpha_4mu.pdf");
        }

        if (signalpoints_4e.size() == 0) {
            CBAlphaP0_4e = 0.0;
        }
        else {
            const int n = signalpoints_4e.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_4e[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_4e[i].alpha;
            TGraph gr(n, x, y);
            TF1 f("f", "[0]", 0, 1000.);
            gr.Fit(&f);
            CBAlphaP0_4e = f.GetParameter(0);
            printFit(gr, f, "CBAlpha_4e.pdf");
        }

        if (signalpoints_2e2mu.size() == 0) {
            CBAlphaP0_2e2mu = 0.0;
        }
        else {
            const int n = signalpoints_2e2mu.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_2e2mu[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_2e2mu[i].alpha;
            TGraph gr(n, x, y);
            TF1 f("f", "[0]", 0, 1000.);
            gr.Fit(&f);
            CBAlphaP0_2e2mu = f.GetParameter(0);
            printFit(gr, f, "CBAlpha_2e2mu.pdf");
        }
    }

    void fitCBN() {
        if (signalpoints_4mu.size() == 0) {
            CBNP0_4mu = 0.0;
        }
        else {
            const int n = signalpoints_4mu.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_4mu[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_4mu[i].n;
            TGraph gr(n, x, y);
            TF1 f("f", "[0]", 0, 1000.);
            gr.Fit(&f);
            CBNP0_4mu = f.GetParameter(0);
            printFit(gr, f, "CBAlpha_4mu.pdf");
        }

        if (signalpoints_4e.size() == 0) {
            CBNP0_4e = 0.0;
        }
        else {
            const int n = signalpoints_4e.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_4e[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_4e[i].n;
            TGraph gr(n, x, y);
            TF1 f("f", "[0]", 0, 1000.);
            gr.Fit(&f);
            CBNP0_4e = f.GetParameter(0);
            printFit(gr, f, "CBAlpha_4e.pdf");
        }

        if (signalpoints_2e2mu.size() == 0) {
            CBNP0_2e2mu = 0.0;
        }
        else {
            const int n = signalpoints_2e2mu.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_2e2mu[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_2e2mu[i].n;
            TGraph gr(n, x, y);
            TF1 f("f", "[0]", 0, 1000.);
            gr.Fit(&f);
            CBNP0_2e2mu = f.GetParameter(0);
            printFit(gr, f, "CBAlpha_2e2mu.pdf");
        }
    }


    void fitEfficiency() {
        if (signalpoints_4mu.size() == 0) {
            EfficiencyP0_4mu = 0.0;
            EfficiencyP1_4mu = 0.0;
            EfficiencyP2_4mu = 0.0;
            EfficiencyP3_4mu = 0.0;
        }
        else {
            const int n = signalpoints_4mu.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_4mu[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_4mu[i].yield / (do7TeV ? (get7TeVggHXsec(x[i]) + get7TeVVBFXsec(x[i])) : (get8TeVggHXsec(x[i]) + get8TeVVBFXsec(x[i])));
            TGraph gr(n, x, y);
            TF1 f("f", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0, 1000.);
            gr.Fit(&f);
            EfficiencyP0_4mu = f.GetParameter(0);
            EfficiencyP1_4mu = f.GetParameter(1);
            EfficiencyP2_4mu = f.GetParameter(2);
            EfficiencyP3_4mu = f.GetParameter(3);
            gr.GetYaxis()->SetRangeUser(0.1,1.0);
            printFit(gr, f, "Efficiency_4mu.pdf");
        }

        if (signalpoints_4e.size() == 0) {
            EfficiencyP0_4e = 0.0;
            EfficiencyP1_4e = 0.0;
            EfficiencyP2_4e = 0.0;
            EfficiencyP3_4e = 0.0;
        }
        else {
            const int n = signalpoints_4e.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_4e[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_4e[i].yield / (do7TeV ? (get7TeVggHXsec(x[i]) + get7TeVVBFXsec(x[i])) : (get8TeVggHXsec(x[i]) + get8TeVVBFXsec(x[i])));
            TGraph gr(n, x, y);
            TF1 f("f", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0, 1000.);
            gr.Fit(&f); 
            EfficiencyP0_4e = f.GetParameter(0);
            EfficiencyP1_4e = f.GetParameter(1);
            EfficiencyP2_4e = f.GetParameter(2);
            EfficiencyP3_4e = f.GetParameter(3);
            gr.GetYaxis()->SetRangeUser(0.1,1.0);
            printFit(gr, f, "Efficiency_4e.pdf");
        }

        if (signalpoints_2e2mu.size() == 0) {
            EfficiencyP0_2e2mu = 0.0;
            EfficiencyP1_2e2mu = 0.0;
            EfficiencyP2_2e2mu = 0.0;
            EfficiencyP3_2e2mu = 0.0;
        }
        else {
            const int n = signalpoints_2e2mu.size();
            float x[n];
            float y[n];

            for (int i = 0; i < n; i++) x[i] = signalpoints_2e2mu[i].higgsmass;
            for (int i = 0; i < n; i++) y[i] = signalpoints_2e2mu[i].yield / (do7TeV ? (get7TeVggHXsec(x[i]) + get7TeVVBFXsec(x[i])) : (get8TeVggHXsec(x[i]) + get8TeVVBFXsec(x[i])));
            TGraph gr(n, x, y);
            TF1 f("f", "[0]+[1]*x+[2]*x*x+[3]*x*x*x", 0, 1000.);
            gr.Fit(&f); 
            EfficiencyP0_2e2mu = f.GetParameter(0);
            EfficiencyP1_2e2mu = f.GetParameter(1);
            EfficiencyP2_2e2mu = f.GetParameter(2);
            EfficiencyP3_2e2mu = f.GetParameter(3);
            gr.GetYaxis()->SetRangeUser(0.1,1.0);
            printFit(gr, f, "Efficiency_2e2mu.pdf");
        }


    }

    void analyze(float mass, float massLowSigFit, float massHighSigFit, int igghid, int ivbfid, bool fitSig, bool fitBkg) {

        stringstream mass_str_ss;
        stringstream gghid_ss;
        stringstream vbfid_ss;
        
        mass_str_ss << mass;
        gghid_ss << igghid;
        vbfid_ss << ivbfid;
        
        std::string mass_str = mass_str_ss.str();
        std::string gghid = gghid_ss.str();
        std::string vbfid = vbfid_ss.str();
        
        std::cout << "Analyzing " << mass_str << " GeV mass point ... " << std::endl;
        
        std::string card_name   = std::string("hzz_m")+mass_str;
        std::string workspace_2e2mu = card_name+"_2e2mu_workspace.root";
        std::string workspace_4e    = card_name+"_4e_workspace.root";
        std::string workspace_4mu   = card_name+"_4mu_workspace.root";
        
        std::string ggh_rootfile = "hzzTree_id";
        ggh_rootfile += gghid;
        ggh_rootfile += ".root";
        
        std::string vbf_rootfile = "hzzTree_id";
        vbf_rootfile += vbfid;
        vbf_rootfile += ".root";
        
        ZZYieldMaker   ymaker_ghzz;
        
        ymaker_data.fill(treeFolder+"hzzTree.root");
        ymaker_zxss.fill(treeFolder+"hzzTree.root"       , 1.0, FR, doSS);
        if (do7TeV) {
            ymaker_ghzz.fill(treeFolder+ggh_rootfile         , getXsecggH(mass)*lumi/evt_7TeV[igghid], 0.0, true );
            ymaker_ghzz.fill(treeFolder+vbf_rootfile         , getXsecVBF(mass)*lumi/evt_7TeV[ivbfid], 0.0, true );
        }
        else {
            ymaker_ghzz.fill(treeFolder+ggh_rootfile         , getXsecggH(mass)*lumi/evt_8TeV[igghid], 0.0, true );
            ymaker_ghzz.fill(treeFolder+vbf_rootfile         , getXsecVBF(mass)*lumi/evt_8TeV[ivbfid], 0.0, true );
        }
        
        QQZZFitMaker fitmaker_2e2mu    ("bkg_qqzz_2e2mu" , massLowBkgFit, massHighBkgFit);
        QQZZFitMaker fitmaker_4e       ("bkg_qqzz_4e"    , massLowBkgFit, massHighBkgFit);
        QQZZFitMaker fitmaker_4mu      ("bkg_qqzz_4mu"   , massLowBkgFit, massHighBkgFit);
        
        GGZZFitMaker fitmaker_gz2e2mu  ("bkg_ggzz_2e2mu" , massLowBkgFit, massHighBkgFit);
        GGZZFitMaker fitmaker_gz4e     ("bkg_ggzz_4e"    , massLowBkgFit, massHighBkgFit);
        GGZZFitMaker fitmaker_gz4mu    ("bkg_ggzz_4mu"   , massLowBkgFit, massHighBkgFit);
        
        ZXFitMaker   fitmaker_zx2e2mu  ("bkg_zjets_2e2mu", massLowBkgFit, massHighBkgFit);
        ZXFitMaker   fitmaker_zx4e     ("bkg_zjets_4e"   , massLowBkgFit, massHighBkgFit);
        ZXFitMaker   fitmaker_zx4mu    ("bkg_zjets_4mu"  , massLowBkgFit, massHighBkgFit);
        
        SignalFitMaker  fitmaker_g2e2mu("sig_ggH_2e2mu"  , mass   , massLowSigFit,  massHighSigFit);
        SignalFitMaker  fitmaker_g4e   ("sig_ggH_4e"     , mass   , massLowSigFit,  massHighSigFit);
        SignalFitMaker  fitmaker_g4mu  ("sig_ggH_4mu"    , mass   , massLowSigFit,  massHighSigFit);
        
        fitmaker_gz4mu.add  (ymaker_ggzz.getFitDataSet(0, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        fitmaker_gz4e.add   (ymaker_ggzz.getFitDataSet(1, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        fitmaker_gz2e2mu.add(ymaker_ggzz.getFitDataSet(2, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        fitmaker_gz2e2mu.add(ymaker_ggzz.getFitDataSet(3, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        
        fitmaker_4mu.add    (ymaker_qqzz.getFitDataSet(0, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        fitmaker_4e.add     (ymaker_qqzz.getFitDataSet(1, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        fitmaker_2e2mu.add  (ymaker_qqzz.getFitDataSet(2, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        fitmaker_2e2mu.add  (ymaker_qqzz.getFitDataSet(3, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        
        fitmaker_zx4mu.add  (ymaker_zxss.getFitDataSet(0, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        fitmaker_zx4e.add   (ymaker_zxss.getFitDataSet(1, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        fitmaker_zx2e2mu.add(ymaker_zxss.getFitDataSet(2, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        fitmaker_zx2e2mu.add(ymaker_zxss.getFitDataSet(3, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        
        fitmaker_g4mu.add   (ymaker_ghzz.getFitDataSet(0, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        fitmaker_g4e.add    (ymaker_ghzz.getFitDataSet(1, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        fitmaker_g2e2mu.add (ymaker_ghzz.getFitDataSet(2, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        fitmaker_g2e2mu.add (ymaker_ghzz.getFitDataSet(3, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut));
        
        if (fitSig) {
            fitmaker_g2e2mu.fit();
            fitmaker_g4e.fit();
            fitmaker_g4mu.fit();
            
            fitmaker_g2e2mu.print("higgs_2e2mu_"+mass_str+".pdf");
            fitmaker_g4e.print   ("higgs_4e_"   +mass_str+".pdf");
            fitmaker_g4mu.print  ("higgs_4mu_"  +mass_str+".pdf");

            float yield_ggh_mm  = ymaker_ghzz.getYield(0, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut);
            float yield_ggh_ee  = ymaker_ghzz.getYield(1, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut);
            float yield_ggh_em  = ymaker_ghzz.getYield(2, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut);
                  yield_ggh_em += ymaker_ghzz.getYield(3, z1min, z2min, massLowBkgFit, massHighBkgFit, melacut);

            SignalCardInfo sci_2e2mu(fitmaker_g2e2mu, yield_ggh_em);
            SignalCardInfo sci_4e   (fitmaker_g4e   , yield_ggh_ee);
            SignalCardInfo sci_4mu  (fitmaker_g4mu  , yield_ggh_mm);
        
            signalpoints_2e2mu.push_back(sci_2e2mu); 
            signalpoints_4e.push_back(sci_4e); 
            signalpoints_4mu.push_back(sci_4mu); 
        }
        
        if (fitBkg) {
            fitmaker_zx2e2mu.fit();
            fitmaker_zx4e.fit();
            fitmaker_zx4mu.fit();
            
            fitmaker_zx2e2mu.print("zjets_2e2mu_"+mass_str+".pdf", 15);
            fitmaker_zx4e.print   ("zjets_4e_"   +mass_str+".pdf", 15);
            fitmaker_zx4mu.print  ("zjets_4mu_"  +mass_str+".pdf", 15);
            
            ZXCardInfo zci_2e2mu(fitmaker_zx2e2mu);
            ZXCardInfo zci_4e   (fitmaker_zx4e   );
            ZXCardInfo zci_4mu  (fitmaker_zx4mu  );

            zxbkg_2e2mu = zci_2e2mu;
            zxbkg_4e    = zci_4e;
            zxbkg_4mu   = zci_4mu;

            fitmaker_2e2mu.fit();
            fitmaker_4e.fit();
            fitmaker_4mu.fit();
            
            fitmaker_2e2mu.print("qqZZ_2e2mu_"+mass_str+".pdf");
            fitmaker_4e.print   ("qqZZ_4e_"   +mass_str+".pdf");
            fitmaker_4mu.print  ("qqZZ_4mu_"  +mass_str+".pdf");
            
            QQZZCardInfo qzci_2e2mu(fitmaker_2e2mu);
            QQZZCardInfo qzci_4e   (fitmaker_4e   );
            QQZZCardInfo qzci_4mu  (fitmaker_4mu  );

            qqzzbkg_2e2mu = qzci_2e2mu;
            qqzzbkg_4e    = qzci_4e;
            qqzzbkg_4mu   = qzci_4mu;

            fitmaker_gz2e2mu.fit();
            fitmaker_gz4e.fit();
            fitmaker_gz4mu.fit();
            
            fitmaker_gz2e2mu.print("ggZZ_2e2mu_"+mass_str+".pdf");
            fitmaker_gz4e.print   ("ggZZ_4e_"   +mass_str+".pdf");
            fitmaker_gz4mu.print  ("ggZZ_4mu_"  +mass_str+".pdf");

            GGZZCardInfo gzci_2e2mu(fitmaker_gz2e2mu);
            GGZZCardInfo gzci_4e   (fitmaker_gz4e   );
            GGZZCardInfo gzci_4mu  (fitmaker_gz4mu  );

            ggzzbkg_2e2mu = gzci_2e2mu;
            ggzzbkg_4e    = gzci_4e;
            ggzzbkg_4mu   = gzci_4mu;

        }    
       
    } 


    void createCard(float mass, float massLow, float massHigh) {

        float gghxsec = getXsecggH(mass); 
        float vbfxsec = getXsecVBF(mass); 
        
        stringstream mass_str_ss;
        mass_str_ss << mass;
        std::string mass_str = mass_str_ss.str();
        
        std::cout << "Creating datacard for " << mass_str << " GeV mass point ... " << std::endl;
        
        std::string card_name   = do1D ? (std::string("card_1D_m")+mass_str+(do7TeV ? "_7TeV" : "_8TeV")) : (std::string("card_2D_m")+mass_str+(do7TeV ? "_7TeV" : "_8TeV"));
        std::string workspace_2e2mu = card_name+"_2e2mu_workspace.root";
        std::string workspace_4e    = card_name+"_4e_workspace.root";
        std::string workspace_4mu   = card_name+"_4mu_workspace.root";
        
        float yield_dt_mm  = ymaker_data.getYield(0, z1min, z2min, massLow, massHigh, melacut);
        float yield_dt_ee  = ymaker_data.getYield(1, z1min, z2min, massLow, massHigh, melacut);
        float yield_dt_em  = ymaker_data.getYield(2, z1min, z2min, massLow, massHigh, melacut);
              yield_dt_em += ymaker_data.getYield(3, z1min, z2min, massLow, massHigh, melacut);
         
        float yield_qq_mm  = ymaker_qqzz.getYield(0, z1min, z2min, massLow, massHigh, melacut);
        float yield_qq_ee  = ymaker_qqzz.getYield(1, z1min, z2min, massLow, massHigh, melacut);
        float yield_qq_em  = ymaker_qqzz.getYield(2, z1min, z2min, massLow, massHigh, melacut);
              yield_qq_em += ymaker_qqzz.getYield(3, z1min, z2min, massLow, massHigh, melacut);
         
        float yield_gg_mm  = ymaker_ggzz.getYield(0, z1min, z2min, massLow, massHigh, melacut);
        float yield_gg_ee  = ymaker_ggzz.getYield(1, z1min, z2min, massLow, massHigh, melacut);
        float yield_gg_em  = ymaker_ggzz.getYield(2, z1min, z2min, massLow, massHigh, melacut);
              yield_gg_em += ymaker_ggzz.getYield(3, z1min, z2min, massLow, massHigh, melacut);
         
        float yield_zj_mm  = ymaker_zxss.getYield(0, z1min, z2min, massLow, massHigh, melacut);
        float yield_zj_ee  = ymaker_zxss.getYield(1, z1min, z2min, massLow, massHigh, melacut);
        float yield_zj_em  = ymaker_zxss.getYield(2, z1min, z2min, massLow, massHigh, melacut);
              yield_zj_em += ymaker_zxss.getYield(3, z1min, z2min, massLow, massHigh, melacut);
         
        std::string card_4mu   = createCardTemplate(do7TeV, 0, do1D, workspace_4mu.c_str());
        std::string card_4e    = createCardTemplate(do7TeV, 1, do1D, workspace_4e.c_str());
        std::string card_2e2mu = createCardTemplate(do7TeV, 2, do1D, workspace_2e2mu.c_str());
        
        card_2e2mu = findAndReplace(card_2e2mu, "SIG_GGH_NAME"  ,  "sig_ggH_2e2mu");
        card_4e    = findAndReplace(card_4e,    "SIG_GGH_NAME"  ,  "sig_ggH_4e");
        card_4mu   = findAndReplace(card_4mu  , "SIG_GGH_NAME"  ,  "sig_ggH_4mu");
         
        card_2e2mu = findAndReplace(card_2e2mu, "SIG_VBF_NAME"  ,  "sig_VBF_2e2mu");
        card_4e    = findAndReplace(card_4e,    "SIG_VBF_NAME"  ,  "sig_VBF_4e");
        card_4mu   = findAndReplace(card_4mu,   "SIG_VBF_NAME"  ,  "sig_VBF_4mu");
         
        card_2e2mu = findAndReplace(card_2e2mu, "BKG_QQZZ_NAME" ,  "bkg_qqzz_2e2mu");
        card_4e    = findAndReplace(card_4e,    "BKG_QQZZ_NAME" ,  "bkg_qqzz_4e");
        card_4mu   = findAndReplace(card_4mu,   "BKG_QQZZ_NAME" ,  "bkg_qqzz_4mu");
         
        card_2e2mu = findAndReplace(card_2e2mu, "BKG_GGZZ_NAME" ,  "bkg_ggzz_2e2mu");
        card_4e    = findAndReplace(card_4e,    "BKG_GGZZ_NAME" ,  "bkg_ggzz_4e");
        card_4mu   = findAndReplace(card_4mu,   "BKG_GGZZ_NAME" ,  "bkg_ggzz_4mu");
         
        card_2e2mu = findAndReplace(card_2e2mu, "BKG_ZJETS_NAME",  "bkg_zjets_2e2mu");
        card_4e    = findAndReplace(card_4e,    "BKG_ZJETS_NAME",  "bkg_zjets_4e");
        card_4mu   = findAndReplace(card_4mu,   "BKG_ZJETS_NAME",  "bkg_zjets_4mu");
        
        card_2e2mu = findAndReplace(card_2e2mu, "GGZZ_PDF"      ,  getGGZZPDFUncertainty7TeV(mass));
        card_4e    = findAndReplace(card_4e,    "GGZZ_PDF"      ,  getGGZZPDFUncertainty7TeV(mass));
        card_4mu   = findAndReplace(card_4mu  , "GGZZ_PDF"      ,  getGGZZPDFUncertainty7TeV(mass));
         
        card_2e2mu = findAndReplace(card_2e2mu, "QQZZ_PDF"      ,  getQQZZPDFUncertainty7TeV(mass));
        card_4e    = findAndReplace(card_4e,    "QQZZ_PDF"      ,  getQQZZPDFUncertainty7TeV(mass));
        card_4mu   = findAndReplace(card_4mu  , "QQZZ_PDF"      ,  getQQZZPDFUncertainty7TeV(mass));
         
        card_2e2mu = findAndReplace(card_2e2mu, "GGZZ_QCD"      ,  getGGZZQCDScaleUncertainty7TeV(mass));
        card_4e    = findAndReplace(card_4e,    "GGZZ_QCD"      ,  getGGZZQCDScaleUncertainty7TeV(mass));
        card_4mu   = findAndReplace(card_4mu  , "GGZZ_QCD"      ,  getGGZZQCDScaleUncertainty7TeV(mass));
         
        card_2e2mu = findAndReplace(card_2e2mu, "QQZZ_QCD"      ,  getQQZZQCDScaleUncertainty7TeV(mass));
        card_4e    = findAndReplace(card_4e,    "QQZZ_QCD"      ,  getQQZZQCDScaleUncertainty7TeV(mass));
        card_4mu   = findAndReplace(card_4mu  , "QQZZ_QCD"      ,  getQQZZQCDScaleUncertainty7TeV(mass));
         
        card_2e2mu = findAndReplace(card_2e2mu, "GGH_PDF"       ,  getggHPDFUncertainty(mass, true), getggHPDFUncertainty(mass, false));
        card_4e    = findAndReplace(card_4e,    "GGH_PDF"       ,  getggHPDFUncertainty(mass, true), getggHPDFUncertainty(mass, false));
        card_4mu   = findAndReplace(card_4mu  , "GGH_PDF"       ,  getggHPDFUncertainty(mass, true), getggHPDFUncertainty(mass, false));
         
        card_2e2mu = findAndReplace(card_2e2mu, "VBF_PDF"       ,  getVBFPDFUncertainty(mass, true), getVBFPDFUncertainty(mass, false));
        card_4e    = findAndReplace(card_4e,    "VBF_PDF"       ,  getVBFPDFUncertainty(mass, true), getVBFPDFUncertainty(mass, false));
        card_4mu   = findAndReplace(card_4mu  , "VBF_PDF"       ,  getVBFPDFUncertainty(mass, true), getVBFPDFUncertainty(mass, false));
         
        card_2e2mu = findAndReplace(card_2e2mu, "GGH_QCD"       ,  getggHQCDScaleUncertainty(mass, true), getggHQCDScaleUncertainty(mass, false));
        card_4e    = findAndReplace(card_4e,    "GGH_QCD"       ,  getggHQCDScaleUncertainty(mass, true), getggHQCDScaleUncertainty(mass, false));
        card_4mu   = findAndReplace(card_4mu  , "GGH_QCD"       ,  getggHQCDScaleUncertainty(mass, true), getggHQCDScaleUncertainty(mass, false));
         
        card_2e2mu = findAndReplace(card_2e2mu, "VBF_QCD"       ,  getVBFQCDScaleUncertainty(mass, true), getVBFQCDScaleUncertainty(mass, false));
        card_4e    = findAndReplace(card_4e,    "VBF_QCD"       ,  getVBFQCDScaleUncertainty(mass, true), getVBFQCDScaleUncertainty(mass, false));
        card_4mu   = findAndReplace(card_4mu  , "VBF_QCD"       ,  getVBFQCDScaleUncertainty(mass, true), getVBFQCDScaleUncertainty(mass, false));
         
        card_2e2mu = findAndReplace(card_2e2mu, "SIG_GGH_YIELD" ,  gghxsec*getEfficiency(mass, 2));
        card_4e    = findAndReplace(card_4e,    "SIG_GGH_YIELD" ,  gghxsec*getEfficiency(mass, 1));
        card_4mu   = findAndReplace(card_4mu  , "SIG_GGH_YIELD" ,  gghxsec*getEfficiency(mass, 0));
         
        card_2e2mu = findAndReplace(card_2e2mu, "SIG_VBF_YIELD" ,  vbfxsec*getEfficiency(mass, 2));
        card_4e    = findAndReplace(card_4e,    "SIG_VBF_YIELD" ,  vbfxsec*getEfficiency(mass, 1));
        card_4mu   = findAndReplace(card_4mu,   "SIG_VBF_YIELD" ,  vbfxsec*getEfficiency(mass, 0));
         
        card_2e2mu = findAndReplace(card_2e2mu, "BKG_QQZZ_YIELD",  yield_qq_em);
        card_4e    = findAndReplace(card_4e,    "BKG_QQZZ_YIELD",  yield_qq_ee);
        card_4mu   = findAndReplace(card_4mu,   "BKG_QQZZ_YIELD",  yield_qq_mm);
         
        card_2e2mu = findAndReplace(card_2e2mu, "BKG_GGZZ_YIELD",  yield_gg_em);
        card_4e    = findAndReplace(card_4e,    "BKG_GGZZ_YIELD",  yield_gg_ee);
        card_4mu   = findAndReplace(card_4mu,   "BKG_GGZZ_YIELD",  yield_gg_mm);
         
        card_2e2mu = findAndReplace(card_2e2mu, "BKG_ZJETS_YIELD", yield_zj_em);
        card_4e    = findAndReplace(card_4e,    "BKG_ZJETS_YIELD", yield_zj_ee);
        card_4mu   = findAndReplace(card_4mu,   "BKG_ZJETS_YIELD", yield_zj_mm);
         
        card_2e2mu = findAndReplace(card_2e2mu, "MEAN_0"         , 0);
        card_4e    = findAndReplace(card_4e,    "MEAN_0"         , 0);
        card_4mu   = findAndReplace(card_4mu,   "MEAN_0"         , 0);
         
        card_2e2mu = findAndReplace(card_2e2mu, "MEAN_1"         , 0.005);
        card_4e    = findAndReplace(card_4e,    "MEAN_1"         , 0.004);
        card_4mu   = findAndReplace(card_4mu,   "MEAN_1"         , 0.005);
         
        card_2e2mu = findAndReplace(card_2e2mu, "SIGMA_0"        , 0);
        card_4e    = findAndReplace(card_4e,    "SIGMA_0"        , 0);
        card_4mu   = findAndReplace(card_4mu,   "SIGMA_0"        , 0);
         
        card_2e2mu = findAndReplace(card_2e2mu, "SIGMA_1"        , 0.3);
        card_4e    = findAndReplace(card_4e,    "SIGMA_1"        , 0.3);
        card_4mu   = findAndReplace(card_4mu,   "SIGMA_1"        , 0.3);
         
        card_2e2mu = findAndReplace(card_2e2mu, "BIN",   "a3");
        card_4e    = findAndReplace(card_4e,    "BIN",   "a2");
        card_4mu   = findAndReplace(card_4mu,   "BIN",   "a1");
        
        card_2e2mu = findAndReplace(card_2e2mu, "OBS",      yield_dt_em);
        card_4e    = findAndReplace(card_4e,    "OBS",      yield_dt_ee);
        card_4mu   = findAndReplace(card_4mu,   "OBS",      yield_dt_mm);
        
        ofstream file_2e2mu;
        file_2e2mu.open ((card_name +"_2e2mu.txt").c_str());
        file_2e2mu << card_2e2mu;
        file_2e2mu.close();
        
        ofstream file_4e;
        file_4e.open ((card_name +"_4e.txt").c_str());
        file_4e << card_4e;
        file_4e.close();
        
        ofstream file_4mu;
        file_4mu.open ((card_name +"_4mu.txt").c_str());
        file_4mu << card_4mu;
        file_4mu.close();
        
        
        RooWorkspace w_2e2mu("w_2e2mu", "");
        RooWorkspace w_4e   ("w_4e"   , "");
        RooWorkspace w_4mu  ("w_4mu"  , "");
        
        RooRealVar CMS_zz4l_melaLD ("CMS_zz4l_melaLD" , "MELA" ,   0,             1,            "");
        RooRealVar CMS_zz4l_mass_1D("CMS_zz4l_mass_1D", "M(4l)", massLow, massHigh, "GeV/c^{2}");
        CMS_zz4l_mass_1D.setBins(100000, "fft");
        
        if (do1D) {
            RooArgSet argset_obs_2e2mu(CMS_zz4l_mass_1D, "argset_obs");
            RooArgSet argset_obs_4mu  (CMS_zz4l_mass_1D, "argset_obs");
            RooArgSet argset_obs_4e   (CMS_zz4l_mass_1D, "argset_obs");
            RooDataSet data_obs_2e2mu("data_obs", "data_obs", argset_obs_2e2mu);
            RooDataSet data_obs_4mu  ("data_obs", "data_obs", argset_obs_4mu);
            RooDataSet data_obs_4e   ("data_obs", "data_obs", argset_obs_4e);
        
            ymaker_data.getDataSet1D(0, z1min, z2min, massLow, massHigh, melacut, data_obs_4mu  , CMS_zz4l_mass_1D);
            ymaker_data.getDataSet1D(1, z1min, z2min, massLow, massHigh, melacut, data_obs_4e   , CMS_zz4l_mass_1D);
            ymaker_data.getDataSet1D(2, z1min, z2min, massLow, massHigh, melacut, data_obs_2e2mu, CMS_zz4l_mass_1D);
            ymaker_data.getDataSet1D(3, z1min, z2min, massLow, massHigh, melacut, data_obs_2e2mu, CMS_zz4l_mass_1D);
            
            w_2e2mu.import(data_obs_2e2mu);
            w_4e.import   (data_obs_4e);
            w_4mu.import  (data_obs_4mu);
        }
        
        else {
            RooArgSet argset_obs_2e2mu(CMS_zz4l_mass_1D, CMS_zz4l_melaLD, "argset_obs");
            RooArgSet argset_obs_4mu  (CMS_zz4l_mass_1D, CMS_zz4l_melaLD, "argset_obs");
            RooArgSet argset_obs_4e   (CMS_zz4l_mass_1D, CMS_zz4l_melaLD, "argset_obs");
            RooDataSet data_obs_2e2mu("data_obs", "data_obs", argset_obs_2e2mu);
            RooDataSet data_obs_4mu  ("data_obs", "data_obs", argset_obs_4mu);
            RooDataSet data_obs_4e   ("data_obs", "data_obs", argset_obs_4e);
            
            ymaker_data.getDataSet2D(0, z1min, z2min, massLow, massHigh, melacut, data_obs_4mu  , CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
            ymaker_data.getDataSet2D(1, z1min, z2min, massLow, massHigh, melacut, data_obs_4e   , CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
            ymaker_data.getDataSet2D(2, z1min, z2min, massLow, massHigh, melacut, data_obs_2e2mu, CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
            ymaker_data.getDataSet2D(3, z1min, z2min, massLow, massHigh, melacut, data_obs_2e2mu, CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
            
            w_2e2mu.import(data_obs_2e2mu);
            w_4e.import   (data_obs_4e);
            w_4mu.import  (data_obs_4mu);
        }
        
        
        ///////////////////// Define parameters //////////////////////////////////
        std::string tevstr = do7TeV ? "_7TeV" : "_8TeV"; 
        RooRealVar qqzz_2e2mu_a0 ("bkg_qqzz_2e2mu_a0" , "", qqzzbkg_2e2mu.a0 );
        RooRealVar qqzz_2e2mu_a1 ("bkg_qqzz_2e2mu_a1" , "", qqzzbkg_2e2mu.a1 );
        RooRealVar qqzz_2e2mu_a2 ("bkg_qqzz_2e2mu_a2" , "", qqzzbkg_2e2mu.a2 );
        RooRealVar qqzz_2e2mu_a3 ("bkg_qqzz_2e2mu_a3" , "", qqzzbkg_2e2mu.a3 );
        RooRealVar qqzz_2e2mu_a4 ("bkg_qqzz_2e2mu_a4" , "", qqzzbkg_2e2mu.a4 );
        RooRealVar qqzz_2e2mu_a5 ("bkg_qqzz_2e2mu_a5" , "", qqzzbkg_2e2mu.a5 );
        RooRealVar qqzz_2e2mu_a6 ("bkg_qqzz_2e2mu_a6" , "", qqzzbkg_2e2mu.a6 );
        RooRealVar qqzz_2e2mu_a7 ("bkg_qqzz_2e2mu_a7" , "", qqzzbkg_2e2mu.a7 );
        RooRealVar qqzz_2e2mu_a8 ("bkg_qqzz_2e2mu_a8" , "", qqzzbkg_2e2mu.a8 );
        RooRealVar qqzz_2e2mu_a9 ("bkg_qqzz_2e2mu_a9" , "", qqzzbkg_2e2mu.a9 );
        RooRealVar qqzz_2e2mu_a10("bkg_qqzz_2e2mu_a10", "", qqzzbkg_2e2mu.a10);
        RooRealVar qqzz_2e2mu_a11("bkg_qqzz_2e2mu_a11", "", qqzzbkg_2e2mu.a11);
        RooRealVar qqzz_2e2mu_a12("bkg_qqzz_2e2mu_a12", "", qqzzbkg_2e2mu.a12);
        RooRealVar qqzz_2e2mu_a13("bkg_qqzz_2e2mu_a13", "", qqzzbkg_2e2mu.a13);
        
        RooRealVar ggzz_2e2mu_a0 ("bkg_ggzz_2e2mu_a0" , "", ggzzbkg_2e2mu.a0 );
        RooRealVar ggzz_2e2mu_a1 ("bkg_ggzz_2e2mu_a1" , "", ggzzbkg_2e2mu.a1 );
        RooRealVar ggzz_2e2mu_a2 ("bkg_ggzz_2e2mu_a2" , "", ggzzbkg_2e2mu.a2 );
        RooRealVar ggzz_2e2mu_a3 ("bkg_ggzz_2e2mu_a3" , "", ggzzbkg_2e2mu.a3 );
        RooRealVar ggzz_2e2mu_a4 ("bkg_ggzz_2e2mu_a4" , "", ggzzbkg_2e2mu.a4 );
        RooRealVar ggzz_2e2mu_a5 ("bkg_ggzz_2e2mu_a5" , "", ggzzbkg_2e2mu.a5 );
        RooRealVar ggzz_2e2mu_a6 ("bkg_ggzz_2e2mu_a6" , "", ggzzbkg_2e2mu.a6 );
        RooRealVar ggzz_2e2mu_a7 ("bkg_ggzz_2e2mu_a7" , "", ggzzbkg_2e2mu.a7 );
        RooRealVar ggzz_2e2mu_a8 ("bkg_ggzz_2e2mu_a8" , "", ggzzbkg_2e2mu.a8 );
        RooRealVar ggzz_2e2mu_a9 ("bkg_ggzz_2e2mu_a9" , "", ggzzbkg_2e2mu.a9 );
        
        RooRealVar zx_2e2mu_mean ("bkg_zjets_2e2mu_mean_zx" , "", zxbkg_2e2mu.mean);
        RooRealVar zx_2e2mu_sigma("bkg_zjets_2e2mu_sigma_zx", "", zxbkg_2e2mu.sigma);
        
        RooRealVar sig_2e2mu_mean_err    (("sig_2e2mu_mean_err"  +tevstr).c_str()  , "", 0., -10., 10.);
        RooRealVar sig_2e2mu_sigma_err   (("sig_2e2mu_sigma_err" +tevstr).c_str()  , "", 0., -10., 10.);

        RooRealVar ggh_2e2mu_alpha       ("sig_ggh_2e2mu_alpha"    , "", getSignalCBAlpha(mass, 2));
        RooRealVar ggh_2e2mu_n           ("sig_ggh_2e2mu_n"        , "", getSignalCBN    (mass, 2));
        RooRealVar ggh_2e2mu_mean_BW     ("sig_ggh_2e2mu_mean_BW"  , "", mass);
        RooRealVar ggh_2e2mu_gamma_BW    ("sig_ggh_2e2mu_gamma_BW" , "", 1.0);
        RooFormulaVar ggh_2e2mu_mean_CB  ("sig_ggh_2e2mu_mean_CB" ,  getSignalCBMeanString(mass, 2).c_str()  , RooArgList(ggh_2e2mu_mean_BW, sig_2e2mu_mean_err));
        RooFormulaVar ggh_2e2mu_sigma_CB ("sig_ggh_2e2mu_sigma_CB",  getSignalCBSigmaString(mass, 2).c_str() , RooArgList(ggh_2e2mu_mean_BW, sig_2e2mu_sigma_err));
        
        RooRealVar vbf_2e2mu_alpha       ("sig_vbf_2e2mu_alpha"    , "", getSignalCBAlpha(mass, 2));
        RooRealVar vbf_2e2mu_n           ("sig_vbf_2e2mu_n"        , "", getSignalCBN    (mass, 2));
        RooRealVar vbf_2e2mu_mean_BW     ("sig_vbf_2e2mu_mean_BW"  , "", mass);
        RooRealVar vbf_2e2mu_gamma_BW    ("sig_vbf_2e2mu_gamma_BW" , "", 1.0);
        RooFormulaVar vbf_2e2mu_mean_CB  ("sig_vbf_2e2mu_mean_CB"  , getSignalCBMeanString(mass, 2).c_str()  , RooArgList(vbf_2e2mu_mean_BW, sig_2e2mu_mean_err));
        RooFormulaVar vbf_2e2mu_sigma_CB ("sig_vbf_2e2mu_sigma_CB" , getSignalCBSigmaString(mass, 2).c_str() , RooArgList(vbf_2e2mu_mean_BW, sig_2e2mu_sigma_err));
        
        RooRealVar qqzz_4e_a0 ("bkg_qqzz_4e_a0" , "", qqzzbkg_4e.a0 );
        RooRealVar qqzz_4e_a1 ("bkg_qqzz_4e_a1" , "", qqzzbkg_4e.a1 );
        RooRealVar qqzz_4e_a2 ("bkg_qqzz_4e_a2" , "", qqzzbkg_4e.a2 );
        RooRealVar qqzz_4e_a3 ("bkg_qqzz_4e_a3" , "", qqzzbkg_4e.a3 );
        RooRealVar qqzz_4e_a4 ("bkg_qqzz_4e_a4" , "", qqzzbkg_4e.a4 );
        RooRealVar qqzz_4e_a5 ("bkg_qqzz_4e_a5" , "", qqzzbkg_4e.a5 );
        RooRealVar qqzz_4e_a6 ("bkg_qqzz_4e_a6" , "", qqzzbkg_4e.a6 );
        RooRealVar qqzz_4e_a7 ("bkg_qqzz_4e_a7" , "", qqzzbkg_4e.a7 );
        RooRealVar qqzz_4e_a8 ("bkg_qqzz_4e_a8" , "", qqzzbkg_4e.a8 );
        RooRealVar qqzz_4e_a9 ("bkg_qqzz_4e_a9" , "", qqzzbkg_4e.a9 );
        RooRealVar qqzz_4e_a10("bkg_qqzz_4e_a10", "", qqzzbkg_4e.a10);
        RooRealVar qqzz_4e_a11("bkg_qqzz_4e_a11", "", qqzzbkg_4e.a11);
        RooRealVar qqzz_4e_a12("bkg_qqzz_4e_a12", "", qqzzbkg_4e.a12);
        RooRealVar qqzz_4e_a13("bkg_qqzz_4e_a13", "", qqzzbkg_4e.a13);
        
        RooRealVar ggzz_4e_a0 ("bkg_ggzz_4e_a0" , "", ggzzbkg_4e.a0 );
        RooRealVar ggzz_4e_a1 ("bkg_ggzz_4e_a1" , "", ggzzbkg_4e.a1 );
        RooRealVar ggzz_4e_a2 ("bkg_ggzz_4e_a2" , "", ggzzbkg_4e.a2 );
        RooRealVar ggzz_4e_a3 ("bkg_ggzz_4e_a3" , "", ggzzbkg_4e.a3 );
        RooRealVar ggzz_4e_a4 ("bkg_ggzz_4e_a4" , "", ggzzbkg_4e.a4 );
        RooRealVar ggzz_4e_a5 ("bkg_ggzz_4e_a5" , "", ggzzbkg_4e.a5 );
        RooRealVar ggzz_4e_a6 ("bkg_ggzz_4e_a6" , "", ggzzbkg_4e.a6 );
        RooRealVar ggzz_4e_a7 ("bkg_ggzz_4e_a7" , "", ggzzbkg_4e.a7 );
        RooRealVar ggzz_4e_a8 ("bkg_ggzz_4e_a8" , "", ggzzbkg_4e.a8 );
        RooRealVar ggzz_4e_a9 ("bkg_ggzz_4e_a9" , "", ggzzbkg_4e.a9 );
        
        RooRealVar zx_4e_mean ("bkg_zjets_4e_mean_zx" , "", zxbkg_4e.mean);
        RooRealVar zx_4e_sigma("bkg_zjets_4e_sigma_zx", "", zxbkg_4e.sigma);
        
        RooRealVar sig_4e_mean_err    (("sig_4e_mean_err"  +tevstr).c_str()  , "", 0., -10., 10.);
        RooRealVar sig_4e_sigma_err   (("sig_4e_sigma_err" +tevstr).c_str()  , "", 0., -10., 10.);

        RooRealVar ggh_4e_alpha       ("sig_ggh_4e_alpha"    , "", getSignalCBAlpha(mass, 1));
        RooRealVar ggh_4e_n           ("sig_ggh_4e_n"        , "", getSignalCBN    (mass, 1));
        RooRealVar ggh_4e_mean_BW     ("sig_ggh_4e_mean_BW"  , "", mass);
        RooRealVar ggh_4e_gamma_BW    ("sig_ggh_4e_gamma_BW" , "", 1.0);
        RooFormulaVar ggh_4e_mean_CB  ("sig_ggh_4e_mean_CB" ,  getSignalCBMeanString(mass, 1).c_str()  , RooArgList(ggh_4e_mean_BW, sig_4e_mean_err));
        RooFormulaVar ggh_4e_sigma_CB ("sig_ggh_4e_sigma_CB",  getSignalCBSigmaString(mass, 1).c_str() , RooArgList(ggh_4e_mean_BW, sig_4e_sigma_err));
        
        RooRealVar vbf_4e_alpha       ("sig_vbf_4e_alpha"    , "", getSignalCBAlpha(mass, 1));
        RooRealVar vbf_4e_n           ("sig_vbf_4e_n"        , "", getSignalCBN    (mass, 1));
        RooRealVar vbf_4e_mean_BW     ("sig_vbf_4e_mean_BW"  , "", mass);
        RooRealVar vbf_4e_gamma_BW    ("sig_vbf_4e_gamma_BW" , "", 1.0);
        RooFormulaVar vbf_4e_mean_CB  ("sig_vbf_4e_mean_CB"  , getSignalCBMeanString(mass, 1).c_str()  , RooArgList(vbf_4e_mean_BW, sig_4e_mean_err));
        RooFormulaVar vbf_4e_sigma_CB ("sig_vbf_4e_sigma_CB" , getSignalCBSigmaString(mass, 1).c_str() , RooArgList(vbf_4e_mean_BW, sig_4e_sigma_err));
        
        RooRealVar qqzz_4mu_a0 ("bkg_qqzz_4mu_a0" , "", qqzzbkg_4mu.a0 ); 
        RooRealVar qqzz_4mu_a1 ("bkg_qqzz_4mu_a1" , "", qqzzbkg_4mu.a1 ); 
        RooRealVar qqzz_4mu_a2 ("bkg_qqzz_4mu_a2" , "", qqzzbkg_4mu.a2 ); 
        RooRealVar qqzz_4mu_a3 ("bkg_qqzz_4mu_a3" , "", qqzzbkg_4mu.a3 ); 
        RooRealVar qqzz_4mu_a4 ("bkg_qqzz_4mu_a4" , "", qqzzbkg_4mu.a4 ); 
        RooRealVar qqzz_4mu_a5 ("bkg_qqzz_4mu_a5" , "", qqzzbkg_4mu.a5 ); 
        RooRealVar qqzz_4mu_a6 ("bkg_qqzz_4mu_a6" , "", qqzzbkg_4mu.a6 ); 
        RooRealVar qqzz_4mu_a7 ("bkg_qqzz_4mu_a7" , "", qqzzbkg_4mu.a7 ); 
        RooRealVar qqzz_4mu_a8 ("bkg_qqzz_4mu_a8" , "", qqzzbkg_4mu.a8 ); 
        RooRealVar qqzz_4mu_a9 ("bkg_qqzz_4mu_a9" , "", qqzzbkg_4mu.a9 ); 
        RooRealVar qqzz_4mu_a10("bkg_qqzz_4mu_a10", "", qqzzbkg_4mu.a10);
        RooRealVar qqzz_4mu_a11("bkg_qqzz_4mu_a11", "", qqzzbkg_4mu.a11);
        RooRealVar qqzz_4mu_a12("bkg_qqzz_4mu_a12", "", qqzzbkg_4mu.a12);
        RooRealVar qqzz_4mu_a13("bkg_qqzz_4mu_a13", "", qqzzbkg_4mu.a13);
        
        RooRealVar ggzz_4mu_a0 ("bkg_ggzz_4mu_a0" , "", ggzzbkg_4mu.a0 ); 
        RooRealVar ggzz_4mu_a1 ("bkg_ggzz_4mu_a1" , "", ggzzbkg_4mu.a1 ); 
        RooRealVar ggzz_4mu_a2 ("bkg_ggzz_4mu_a2" , "", ggzzbkg_4mu.a2 ); 
        RooRealVar ggzz_4mu_a3 ("bkg_ggzz_4mu_a3" , "", ggzzbkg_4mu.a3 ); 
        RooRealVar ggzz_4mu_a4 ("bkg_ggzz_4mu_a4" , "", ggzzbkg_4mu.a4 ); 
        RooRealVar ggzz_4mu_a5 ("bkg_ggzz_4mu_a5" , "", ggzzbkg_4mu.a5 ); 
        RooRealVar ggzz_4mu_a6 ("bkg_ggzz_4mu_a6" , "", ggzzbkg_4mu.a6 ); 
        RooRealVar ggzz_4mu_a7 ("bkg_ggzz_4mu_a7" , "", ggzzbkg_4mu.a7 ); 
        RooRealVar ggzz_4mu_a8 ("bkg_ggzz_4mu_a8" , "", ggzzbkg_4mu.a8 ); 
        RooRealVar ggzz_4mu_a9 ("bkg_ggzz_4mu_a9" , "", ggzzbkg_4mu.a9 ); 
        
        RooRealVar zx_4mu_mean ("bkg_zjets_4mu_mean_zx" , "", zxbkg_4mu.mean);
        RooRealVar zx_4mu_sigma("bkg_zjets_4mu_sigma_zx", "", zxbkg_4mu.sigma);
        
        RooRealVar sig_4mu_mean_err    (("sig_4mu_mean_err"  +tevstr).c_str()  , "", 0., -10., 10.);
        RooRealVar sig_4mu_sigma_err   (("sig_4mu_sigma_err" +tevstr).c_str()  , "", 0., -10., 10.);

        RooRealVar ggh_4mu_alpha       ("sig_ggh_4mu_alpha"    , "", getSignalCBAlpha(mass, 0));
        RooRealVar ggh_4mu_n           ("sig_ggh_4mu_n"        , "", getSignalCBN    (mass, 0));
        RooRealVar ggh_4mu_mean_BW     ("sig_ggh_4mu_mean_BW"  , "", mass);
        RooRealVar ggh_4mu_gamma_BW    ("sig_ggh_4mu_gamma_BW" , "", 1.0);
        RooFormulaVar ggh_4mu_mean_CB  ("sig_ggh_4mu_mean_CB"  , getSignalCBMeanString(mass, 0).c_str()  , RooArgList(ggh_4mu_mean_BW, sig_4mu_mean_err));
        RooFormulaVar ggh_4mu_sigma_CB ("sig_ggh_4mu_sigma_CB" , getSignalCBSigmaString(mass, 0).c_str() , RooArgList(ggh_4mu_mean_BW, sig_4mu_sigma_err));
        
        RooRealVar vbf_4mu_alpha       ("sig_vbf_4mu_alpha"    , "", getSignalCBAlpha(mass, 0));
        RooRealVar vbf_4mu_n           ("sig_vbf_4mu_n"        , "", getSignalCBN    (mass, 0));
        RooRealVar vbf_4mu_mean_BW     ("sig_vbf_4mu_mean_BW"  , "", mass);
        RooRealVar vbf_4mu_gamma_BW    ("sig_vbf_4mu_gamma_BW" , "", 1.0);
        RooFormulaVar vbf_4mu_mean_CB  ("sig_vbf_4mu_mean_CB"  , getSignalCBMeanString(mass, 0).c_str()  , RooArgList(vbf_4mu_mean_BW, sig_4mu_mean_err));
        RooFormulaVar vbf_4mu_sigma_CB ("sig_vbf_4mu_sigma_CB" , getSignalCBSigmaString(mass, 0).c_str() , RooArgList(vbf_4mu_mean_BW, sig_4mu_sigma_err));
        
        /////////// Set parameters to constant //////////////////
        qqzz_2e2mu_a0      .setConstant(kTRUE);
        qqzz_2e2mu_a1      .setConstant(kTRUE);
        qqzz_2e2mu_a2      .setConstant(kTRUE);
        qqzz_2e2mu_a3      .setConstant(kTRUE);
        qqzz_2e2mu_a4      .setConstant(kTRUE);
        qqzz_2e2mu_a5      .setConstant(kTRUE);
        qqzz_2e2mu_a6      .setConstant(kTRUE);
        qqzz_2e2mu_a7      .setConstant(kTRUE);
        qqzz_2e2mu_a8      .setConstant(kTRUE);
        qqzz_2e2mu_a9      .setConstant(kTRUE);
        qqzz_2e2mu_a10     .setConstant(kTRUE);
        qqzz_2e2mu_a11     .setConstant(kTRUE);
        qqzz_2e2mu_a12     .setConstant(kTRUE);
        qqzz_2e2mu_a13     .setConstant(kTRUE);
        
        ggzz_2e2mu_a0      .setConstant(kTRUE);
        ggzz_2e2mu_a1      .setConstant(kTRUE);
        ggzz_2e2mu_a2      .setConstant(kTRUE);
        ggzz_2e2mu_a3      .setConstant(kTRUE);
        ggzz_2e2mu_a4      .setConstant(kTRUE);
        ggzz_2e2mu_a5      .setConstant(kTRUE);
        ggzz_2e2mu_a6      .setConstant(kTRUE);
        ggzz_2e2mu_a7      .setConstant(kTRUE);
        ggzz_2e2mu_a8      .setConstant(kTRUE);
        ggzz_2e2mu_a9      .setConstant(kTRUE);
        
        zx_2e2mu_mean      .setConstant(kTRUE);
        zx_2e2mu_sigma     .setConstant(kTRUE);
       
        ggh_2e2mu_alpha    .setConstant(kTRUE);
        ggh_2e2mu_n        .setConstant(kTRUE);
        ggh_2e2mu_mean_BW  .setConstant(kTRUE);
        ggh_2e2mu_gamma_BW .setConstant(kTRUE);
        
        vbf_2e2mu_alpha    .setConstant(kTRUE);
        vbf_2e2mu_n        .setConstant(kTRUE);
        vbf_2e2mu_mean_BW  .setConstant(kTRUE);
        vbf_2e2mu_gamma_BW .setConstant(kTRUE);
        
        qqzz_4e_a0      .setConstant(kTRUE);
        qqzz_4e_a1      .setConstant(kTRUE);
        qqzz_4e_a2      .setConstant(kTRUE);
        qqzz_4e_a3      .setConstant(kTRUE);
        qqzz_4e_a4      .setConstant(kTRUE);
        qqzz_4e_a5      .setConstant(kTRUE);
        qqzz_4e_a6      .setConstant(kTRUE);
        qqzz_4e_a7      .setConstant(kTRUE);
        qqzz_4e_a8      .setConstant(kTRUE);
        qqzz_4e_a9      .setConstant(kTRUE);
        qqzz_4e_a10     .setConstant(kTRUE);
        qqzz_4e_a11     .setConstant(kTRUE);
        qqzz_4e_a12     .setConstant(kTRUE);
        qqzz_4e_a13     .setConstant(kTRUE);
        
        ggzz_4e_a0      .setConstant(kTRUE);
        ggzz_4e_a1      .setConstant(kTRUE);
        ggzz_4e_a2      .setConstant(kTRUE);
        ggzz_4e_a3      .setConstant(kTRUE);
        ggzz_4e_a4      .setConstant(kTRUE);
        ggzz_4e_a5      .setConstant(kTRUE);
        ggzz_4e_a6      .setConstant(kTRUE);
        ggzz_4e_a7      .setConstant(kTRUE);
        ggzz_4e_a8      .setConstant(kTRUE);
        ggzz_4e_a9      .setConstant(kTRUE);
        
        zx_4e_mean      .setConstant(kTRUE);
        zx_4e_sigma     .setConstant(kTRUE);
        
        ggh_4e_alpha    .setConstant(kTRUE);
        ggh_4e_n        .setConstant(kTRUE);
        ggh_4e_mean_BW  .setConstant(kTRUE);
        ggh_4e_gamma_BW .setConstant(kTRUE);
 
        vbf_4e_alpha    .setConstant(kTRUE);
        vbf_4e_n        .setConstant(kTRUE);
        vbf_4e_mean_BW  .setConstant(kTRUE);
        vbf_4e_gamma_BW .setConstant(kTRUE);
        
        qqzz_4mu_a0      .setConstant(kTRUE);
        qqzz_4mu_a1      .setConstant(kTRUE);
        qqzz_4mu_a2      .setConstant(kTRUE);
        qqzz_4mu_a3      .setConstant(kTRUE);
        qqzz_4mu_a4      .setConstant(kTRUE);
        qqzz_4mu_a5      .setConstant(kTRUE);
        qqzz_4mu_a6      .setConstant(kTRUE);
        qqzz_4mu_a7      .setConstant(kTRUE);
        qqzz_4mu_a8      .setConstant(kTRUE);
        qqzz_4mu_a9      .setConstant(kTRUE);
        qqzz_4mu_a10     .setConstant(kTRUE);
        qqzz_4mu_a11     .setConstant(kTRUE);
        qqzz_4mu_a12     .setConstant(kTRUE);
        qqzz_4mu_a13     .setConstant(kTRUE);
        
        ggzz_4mu_a0      .setConstant(kTRUE);
        ggzz_4mu_a1      .setConstant(kTRUE);
        ggzz_4mu_a2      .setConstant(kTRUE);
        ggzz_4mu_a3      .setConstant(kTRUE);
        ggzz_4mu_a4      .setConstant(kTRUE);
        ggzz_4mu_a5      .setConstant(kTRUE);
        ggzz_4mu_a6      .setConstant(kTRUE);
        ggzz_4mu_a7      .setConstant(kTRUE);
        ggzz_4mu_a8      .setConstant(kTRUE);
        ggzz_4mu_a9      .setConstant(kTRUE);
        
        zx_4mu_mean      .setConstant(kTRUE);
        zx_4mu_sigma     .setConstant(kTRUE);
        
        ggh_4mu_alpha    .setConstant(kTRUE);
        ggh_4mu_n        .setConstant(kTRUE);
        ggh_4mu_mean_BW  .setConstant(kTRUE);
        ggh_4mu_gamma_BW .setConstant(kTRUE);
        
        vbf_4mu_alpha    .setConstant(kTRUE);
        vbf_4mu_n        .setConstant(kTRUE);
        vbf_4mu_mean_BW  .setConstant(kTRUE);
        vbf_4mu_gamma_BW .setConstant(kTRUE);
        
        ////////////////// Define the PDFs /////////////////////////////////
        
        const char* bkg_qqzz_2e2mu_pdf_name  = do1D ? "bkg_qqzz_2e2mu"  : "bkg_qqzz_2e2mu_1D";
        const char* bkg_qqzz_4e_pdf_name     = do1D ? "bkg_qqzz_4e"     : "bkg_qqzz_4e_1D"   ;
        const char* bkg_qqzz_4mu_pdf_name    = do1D ? "bkg_qqzz_4mu"    : "bkg_qqzz_4mu_1D"  ;
        
        const char* bkg_ggzz_2e2mu_pdf_name  = do1D ? "bkg_ggzz_2e2mu"  : "bkg_ggzz_2e2mu_1D";
        const char* bkg_ggzz_4e_pdf_name     = do1D ? "bkg_ggzz_4e"     : "bkg_ggzz_4e_1D"   ;
        const char* bkg_ggzz_4mu_pdf_name    = do1D ? "bkg_ggzz_4mu"    : "bkg_ggzz_4mu_1D"  ;
        
        const char* bkg_zjets_2e2mu_pdf_name = do1D ? "bkg_zjets_2e2mu" : "bkg_zjets_2e2mu_1D";
        const char* bkg_zjets_4e_pdf_name    = do1D ? "bkg_zjets_4e"    : "bkg_zjets_4e_1D"   ;
        const char* bkg_zjets_4mu_pdf_name   = do1D ? "bkg_zjets_4mu"   : "bkg_zjets_4mu_1D"  ;
        
        const char* sig_ggH_2e2mu_pdf_name   = do1D ? "sig_ggH_2e2mu"   : "sig_ggH_2e2mu_1D" ;
        const char* sig_ggH_4e_pdf_name      = do1D ? "sig_ggH_4e"      : "sig_ggH_4e_1D"    ;
        const char* sig_ggH_4mu_pdf_name     = do1D ? "sig_ggH_4mu"     : "sig_ggH_4mu_1D"   ;
        
        const char* sig_VBF_2e2mu_pdf_name   = do1D ? "sig_VBF_2e2mu"   : "sig_VBF_2e2mu_1D" ;
        const char* sig_VBF_4e_pdf_name      = do1D ? "sig_VBF_4e"      : "sig_VBF_4e_1D"    ;
        const char* sig_VBF_4mu_pdf_name     = do1D ? "sig_VBF_4mu"     : "sig_VBF_4mu_1D"   ;
        
        RooqqZZPdf_v2 bkg_qqzz_2e2mu_pdf(bkg_qqzz_2e2mu_pdf_name,"",CMS_zz4l_mass_1D,
                           qqzz_2e2mu_a0,
                           qqzz_2e2mu_a1,
                           qqzz_2e2mu_a2,
                           qqzz_2e2mu_a3,
                           qqzz_2e2mu_a4,
                           qqzz_2e2mu_a5,
                           qqzz_2e2mu_a6,
                           qqzz_2e2mu_a7,
                           qqzz_2e2mu_a8,
                           qqzz_2e2mu_a9,
                           qqzz_2e2mu_a10,
                           qqzz_2e2mu_a11,
                           qqzz_2e2mu_a12,
                           qqzz_2e2mu_a13);
        
        RooqqZZPdf_v2 bkg_qqzz_4e_pdf(bkg_qqzz_4e_pdf_name,"",CMS_zz4l_mass_1D,
                           qqzz_4e_a0,
                           qqzz_4e_a1,
                           qqzz_4e_a2,
                           qqzz_4e_a3,
                           qqzz_4e_a4,
                           qqzz_4e_a5,
                           qqzz_4e_a6,
                           qqzz_4e_a7,
                           qqzz_4e_a8,
                           qqzz_4e_a9,
                           qqzz_4e_a10,
                           qqzz_4e_a11,
                           qqzz_4e_a12,
                           qqzz_4e_a13);
        
        RooqqZZPdf_v2 bkg_qqzz_4mu_pdf(bkg_qqzz_4mu_pdf_name,"",CMS_zz4l_mass_1D,
                           qqzz_4mu_a0,
                           qqzz_4mu_a1,
                           qqzz_4mu_a2,
                           qqzz_4mu_a3,
                           qqzz_4mu_a4,
                           qqzz_4mu_a5,
                           qqzz_4mu_a6,
                           qqzz_4mu_a7,
                           qqzz_4mu_a8,
                           qqzz_4mu_a9,
                           qqzz_4mu_a10,
                           qqzz_4mu_a11,
                           qqzz_4mu_a12,
                           qqzz_4mu_a13);
        
        RooggZZPdf_v2 bkg_ggzz_2e2mu_pdf(bkg_ggzz_2e2mu_pdf_name,"",CMS_zz4l_mass_1D,
                           ggzz_2e2mu_a0,
                           ggzz_2e2mu_a1,
                           ggzz_2e2mu_a2,
                           ggzz_2e2mu_a3,
                           ggzz_2e2mu_a4,
                           ggzz_2e2mu_a5,
                           ggzz_2e2mu_a6,
                           ggzz_2e2mu_a7,
                           ggzz_2e2mu_a8,
                           ggzz_2e2mu_a9);
        
        RooggZZPdf_v2 bkg_ggzz_4e_pdf(bkg_ggzz_4e_pdf_name,"",CMS_zz4l_mass_1D,
                           ggzz_4e_a0,
                           ggzz_4e_a1,
                           ggzz_4e_a2,
                           ggzz_4e_a3,
                           ggzz_4e_a4,
                           ggzz_4e_a5,
                           ggzz_4e_a6,
                           ggzz_4e_a7,
                           ggzz_4e_a8,
                           ggzz_4e_a9);
        
        RooggZZPdf_v2 bkg_ggzz_4mu_pdf(bkg_ggzz_4mu_pdf_name,"",CMS_zz4l_mass_1D,
                           ggzz_4mu_a0,
                           ggzz_4mu_a1,
                           ggzz_4mu_a2,
                           ggzz_4mu_a3,
                           ggzz_4mu_a4,
                           ggzz_4mu_a5,
                           ggzz_4mu_a6,
                           ggzz_4mu_a7,
                           ggzz_4mu_a8,
                           ggzz_4mu_a9);
        
        RooLandau bkg_zjets_2e2mu_pdf(bkg_zjets_2e2mu_pdf_name, "",CMS_zz4l_mass_1D,zx_2e2mu_mean,zx_2e2mu_sigma);
        RooLandau bkg_zjets_4e_pdf   (bkg_zjets_4e_pdf_name   , "",CMS_zz4l_mass_1D,zx_4e_mean   ,zx_4e_sigma   );
        RooLandau bkg_zjets_4mu_pdf  (bkg_zjets_4mu_pdf_name  , "",CMS_zz4l_mass_1D,zx_4mu_mean  ,zx_4mu_sigma  );
        
        RooCBShape      signalCB_ggH_2e2mu   ("signalCB_ggH_2e2mu", "", CMS_zz4l_mass_1D, ggh_2e2mu_mean_CB,ggh_2e2mu_sigma_CB,ggh_2e2mu_alpha,ggh_2e2mu_n);
        RooRelBWUFParam signalBW_ggH_2e2mu   ("signalBW_ggH_2e2mu", "", CMS_zz4l_mass_1D, ggh_2e2mu_mean_BW,ggh_2e2mu_gamma_BW);
        RooFFTConvPdf   sig_ggH_2e2mu_pdf(sig_ggH_2e2mu_pdf_name, "", CMS_zz4l_mass_1D, signalBW_ggH_2e2mu,signalCB_ggH_2e2mu,2);
        sig_ggH_2e2mu_pdf.setBufferFraction(0.2);
        
        RooCBShape      signalCB_ggH_4e   ("signalCB_ggH_4e", "", CMS_zz4l_mass_1D, ggh_4e_mean_CB,ggh_4e_sigma_CB,ggh_4e_alpha,ggh_4e_n);
        RooRelBWUFParam signalBW_ggH_4e   ("signalBW_ggH_4e", "", CMS_zz4l_mass_1D, ggh_4e_mean_BW,ggh_4e_gamma_BW);
        RooFFTConvPdf   sig_ggH_4e_pdf(sig_ggH_4e_pdf_name, "", CMS_zz4l_mass_1D, signalBW_ggH_4e,signalCB_ggH_4e,2);
        sig_ggH_4e_pdf.setBufferFraction(0.2);
        
        RooCBShape      signalCB_ggH_4mu   ("signalCB_ggH_4mu", "", CMS_zz4l_mass_1D, ggh_4mu_mean_CB,ggh_4mu_sigma_CB,ggh_4mu_alpha,ggh_4mu_n);
        RooRelBWUFParam signalBW_ggH_4mu   ("signalBW_ggH_4mu", "", CMS_zz4l_mass_1D, ggh_4mu_mean_BW,ggh_4mu_gamma_BW);
        RooFFTConvPdf   sig_ggH_4mu_pdf(sig_ggH_4mu_pdf_name, "", CMS_zz4l_mass_1D, signalBW_ggH_4mu,signalCB_ggH_4mu,2);
        sig_ggH_4mu_pdf.setBufferFraction(0.2);
        
        RooCBShape      signalCB_VBF_2e2mu   ("signalCB_VBF_2e2mu", "", CMS_zz4l_mass_1D, vbf_2e2mu_mean_CB,vbf_2e2mu_sigma_CB,vbf_2e2mu_alpha,vbf_2e2mu_n);
        RooRelBWUFParam signalBW_VBF_2e2mu   ("signalBW_VBF_2e2mu", "", CMS_zz4l_mass_1D, vbf_2e2mu_mean_BW,vbf_2e2mu_gamma_BW);
        RooFFTConvPdf   sig_VBF_2e2mu_pdf(sig_VBF_2e2mu_pdf_name, "", CMS_zz4l_mass_1D, signalBW_VBF_2e2mu,signalCB_VBF_2e2mu,2);
        sig_VBF_2e2mu_pdf.setBufferFraction(0.2);
        
        RooCBShape      signalCB_VBF_4e   ("signalCB_VBF_4e", "", CMS_zz4l_mass_1D, vbf_4e_mean_CB,vbf_4e_sigma_CB,vbf_4e_alpha,vbf_4e_n);
        RooRelBWUFParam signalBW_VBF_4e   ("signalBW_VBF_4e", "", CMS_zz4l_mass_1D, vbf_4e_mean_BW,vbf_4e_gamma_BW);
        RooFFTConvPdf   sig_VBF_4e_pdf(sig_VBF_4e_pdf_name, "", CMS_zz4l_mass_1D, signalBW_VBF_4e,signalCB_VBF_4e,2);
        sig_VBF_4e_pdf.setBufferFraction(0.2);
        
        RooCBShape      signalCB_VBF_4mu   ("signalCB_VBF_4mu", "", CMS_zz4l_mass_1D, vbf_4mu_mean_CB,vbf_4mu_sigma_CB,vbf_4mu_alpha,vbf_4mu_n);
        RooRelBWUFParam signalBW_VBF_4mu   ("signalBW_VBF_4mu", "", CMS_zz4l_mass_1D, vbf_4mu_mean_BW,vbf_4mu_gamma_BW);
        RooFFTConvPdf   sig_VBF_4mu_pdf(sig_VBF_4mu_pdf_name, "", CMS_zz4l_mass_1D, signalBW_VBF_4mu,signalCB_VBF_4mu,2);
        sig_VBF_4mu_pdf.setBufferFraction(0.2);
        
        if (do1D) { 
            w_2e2mu.import(bkg_qqzz_2e2mu_pdf);
            w_2e2mu.import(bkg_ggzz_2e2mu_pdf);
            w_2e2mu.import(bkg_zjets_2e2mu_pdf);
            w_2e2mu.import(sig_ggH_2e2mu_pdf);
            w_2e2mu.import(sig_VBF_2e2mu_pdf);
            
            w_4mu.import(bkg_qqzz_4mu_pdf);
            w_4mu.import(bkg_ggzz_4mu_pdf);
            w_4mu.import(bkg_zjets_4mu_pdf);
            w_4mu.import(sig_ggH_4mu_pdf);
            w_4mu.import(sig_VBF_4mu_pdf);
            
            w_4e.import(bkg_qqzz_4e_pdf);
            w_4e.import(bkg_ggzz_4e_pdf);
            w_4e.import(bkg_zjets_4e_pdf);
            w_4e.import(sig_ggH_4e_pdf);
            w_4e.import(sig_VBF_4e_pdf);
        }

        else {
            TFile melafile(melafilename.c_str());
            RooArgList v2dList(CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
            RooArgSet  v2dSet (CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
            
            TH2F* melashape_zz_2e2mu = (TH2F*)(melafile.Get("hist2D_bkg_em"));
            TH2F* melashape_zz_4e    = (TH2F*)(melafile.Get("hist2D_bkg_ee"));
            TH2F* melashape_zz_4mu   = (TH2F*)(melafile.Get("hist2D_bkg_mm"));
            
            TH2F* melashape_si_2e2mu = (TH2F*)(melafile.Get("hist2D_sig_em"));
            TH2F* melashape_si_4e    = (TH2F*)(melafile.Get("hist2D_sig_ee"));
            TH2F* melashape_si_4mu   = (TH2F*)(melafile.Get("hist2D_sig_mm"));
            
            TH2F* melashape_zu_2e2mu = (TH2F*)(melafile.Get("hist2D_bkg_em_up"));
            TH2F* melashape_zu_4e    = (TH2F*)(melafile.Get("hist2D_bkg_ee_up"));
            TH2F* melashape_zu_4mu   = (TH2F*)(melafile.Get("hist2D_bkg_mm_up"));
            
            TH2F* melashape_zd_2e2mu = (TH2F*)(melafile.Get("hist2D_bkg_em_dn"));
            TH2F* melashape_zd_4e    = (TH2F*)(melafile.Get("hist2D_bkg_ee_dn"));
            TH2F* melashape_zd_4mu   = (TH2F*)(melafile.Get("hist2D_bkg_mm_dn"));
            
            RooDataHist rhist_qqzz_2e2mu ("rhist_qqzz_2e2mu" , "", v2dList, melashape_zz_2e2mu);
            RooDataHist rhist_qqzz_4e    ("rhist_qqzz_4e"    , "", v2dList, melashape_zz_4e   );
            RooDataHist rhist_qqzz_4mu   ("rhist_qqzz_4mu"   , "", v2dList, melashape_zz_4mu  );
            
            RooDataHist rhist_ggzz_2e2mu ("rhist_ggzz_2e2mu" , "", v2dList, melashape_zz_2e2mu);
            RooDataHist rhist_ggzz_4e    ("rhist_ggzz_4e"    , "", v2dList, melashape_zz_4e   );
            RooDataHist rhist_ggzz_4mu   ("rhist_ggzz_4mu"   , "", v2dList, melashape_zz_4mu  );
            
            RooDataHist rhist_zjets_2e2mu("rhist_zjets_2e2mu", "", v2dList, melashape_zz_2e2mu);
            RooDataHist rhist_zjets_4e   ("rhist_zjets_4e"   , "", v2dList, melashape_zz_4e   );
            RooDataHist rhist_zjets_4mu  ("rhist_zjets_4mu"  , "", v2dList, melashape_zz_4mu  );
            
            RooDataHist rhist_ggH_2e2mu  ("rhist_ggH_2e2mu"  , "", v2dList, melashape_si_2e2mu);
            RooDataHist rhist_ggH_4e     ("rhist_ggH_4e"     , "", v2dList, melashape_si_4e   );
            RooDataHist rhist_ggH_4mu    ("rhist_ggH_4mu"    , "", v2dList, melashape_si_4mu  );
            
            RooDataHist rhist_VBF_2e2mu  ("rhist_VBF_2e2mu"  , "", v2dList, melashape_si_2e2mu);
            RooDataHist rhist_VBF_4e     ("rhist_VBF_4e"     , "", v2dList, melashape_si_4e   );
            RooDataHist rhist_VBF_4mu    ("rhist_VBF_4mu"    , "", v2dList, melashape_si_4mu  );
            
            RooDataHist rhist_zjup_2e2mu ("rhist_zjup_2e2mu" , "", v2dList, melashape_zu_2e2mu);
            RooDataHist rhist_zjup_4e    ("rhist_zjup_4e"    , "", v2dList, melashape_zu_4e   );
            RooDataHist rhist_zjup_4mu   ("rhist_zjup_4mu"   , "", v2dList, melashape_zu_4mu  );
            
            RooDataHist rhist_zjdn_2e2mu ("rhist_zjdn_2e2mu" , "", v2dList, melashape_zd_2e2mu);
            RooDataHist rhist_zjdn_4e    ("rhist_zjdn_4e"    , "", v2dList, melashape_zd_4e   );
            RooDataHist rhist_zjdn_4mu   ("rhist_zjdn_4mu"   , "", v2dList, melashape_zd_4mu  );
            
            RooHistPdf rpdf_qqzz_2e2mu ("bkg_qqzz_2e2mu_mela2D_pdf" , "", v2dSet , rhist_qqzz_2e2mu);
            RooHistPdf rpdf_qqzz_4e    ("bkg_qqzz_4e_mela2D_pdf"    , "", v2dSet , rhist_qqzz_4e   );
            RooHistPdf rpdf_qqzz_4mu   ("bkg_qqzz_4mu_mela2D_pdf"   , "", v2dSet , rhist_qqzz_4mu  );
            
            RooHistPdf rpdf_ggzz_2e2mu ("bkg_ggzz_2e2mu_mela2D_pdf" , "", v2dSet , rhist_ggzz_2e2mu);
            RooHistPdf rpdf_ggzz_4e    ("bkg_ggzz_4e_mela2D_pdf"    , "", v2dSet , rhist_ggzz_4e   );
            RooHistPdf rpdf_ggzz_4mu   ("bkg_ggzz_4mu_mela2D_pdf"   , "", v2dSet , rhist_ggzz_4mu  );
            
            RooHistPdf rpdf_zjets_2e2mu("bkg_zjets_2e2mu_mela2D_pdf", "", v2dSet , rhist_zjets_2e2mu);
            RooHistPdf rpdf_zjets_4e   ("bkg_zjets_4e_mela2D_pdf"   , "", v2dSet , rhist_zjets_4e   );
            RooHistPdf rpdf_zjets_4mu  ("bkg_zjets_4mu_mela2D_pdf"  , "", v2dSet , rhist_zjets_4mu  );
            
            RooHistPdf rpdf_ggH_2e2mu  ("sig_ggH_2e2mu_mela2D_pdf"  , "", v2dSet , rhist_ggH_2e2mu);
            RooHistPdf rpdf_ggH_4e     ("sig_ggH_4e_mela2D_pdf"     , "", v2dSet , rhist_ggH_4e   );
            RooHistPdf rpdf_ggH_4mu    ("sig_ggH_4mu_mela2D_pdf"    , "", v2dSet , rhist_ggH_4mu  );
            
            RooHistPdf rpdf_VBF_2e2mu  ("sig_VBF_2e2mu_mela2D_pdf"  , "", v2dSet , rhist_VBF_2e2mu);
            RooHistPdf rpdf_VBF_4e     ("sig_VBF_4e_mela2D_pdf"     , "", v2dSet , rhist_VBF_4e   );
            RooHistPdf rpdf_VBF_4mu    ("sig_VBF_4mu_mela2D_pdf"    , "", v2dSet , rhist_VBF_4mu  );
            
            RooHistPdf rpdf_zjup_2e2mu ("bkg_zjup_2e2mu_mela2D_pdf" , "", v2dSet , rhist_zjup_2e2mu);
            RooHistPdf rpdf_zjup_4e    ("bkg_zjup_4e_mela2D_pdf"    , "", v2dSet , rhist_zjup_4e   );
            RooHistPdf rpdf_zjup_4mu   ("bkg_zjup_4mu_mela2D_pdf"   , "", v2dSet , rhist_zjup_4mu  );
            
            RooHistPdf rpdf_zjdn_2e2mu ("bkg_zjdn_2e2mu_mela2D_pdf" , "", v2dSet , rhist_zjdn_2e2mu);
            RooHistPdf rpdf_zjdn_4e    ("bkg_zjdn_4e_mela2D_pdf"    , "", v2dSet , rhist_zjdn_4e   );
            RooHistPdf rpdf_zjdn_4mu   ("bkg_zjdn_4mu_mela2D_pdf"   , "", v2dSet , rhist_zjdn_4mu  );
           
            RooArgList list_zjets_2e2mu;
            RooArgList list_zjets_4e;
            RooArgList list_zjets_4mu;
            list_zjets_2e2mu.add(rpdf_zjets_2e2mu);
            list_zjets_2e2mu.add(rpdf_zjup_2e2mu);
            list_zjets_2e2mu.add(rpdf_zjdn_2e2mu);
            list_zjets_4e.add(rpdf_zjets_4e);
            list_zjets_4e.add(rpdf_zjup_4e);
            list_zjets_4e.add(rpdf_zjdn_4e);
            list_zjets_4mu.add(rpdf_zjets_4mu);
            list_zjets_4mu.add(rpdf_zjup_4mu);
            list_zjets_4mu.add(rpdf_zjdn_4mu);
            RooRealVar CMS_zz4l_bkgMELA("CMS_zz4l_bkgMELA" ,"" ,0,-10,10); 
 
            FastVerticalInterpHistPdf2D plpdf_qqzz_2e2mu ("bkg_qqzz_2e2mu_FVIHP" ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_qqzz_2e2mu) ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_qqzz_4e    ("bkg_qqzz_4e_FVIHP"    ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_qqzz_4e)    ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_qqzz_4mu   ("bkg_qqzz_4mu_FVIHP"   ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_qqzz_4mu)   ,RooArgList()                ,1.0,1);
            
            FastVerticalInterpHistPdf2D plpdf_ggzz_2e2mu ("bkg_ggzz_2e2mu_FVIHP" ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ggzz_2e2mu) ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_ggzz_4e    ("bkg_ggzz_4e_FVIHP"    ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ggzz_4e)    ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_ggzz_4mu   ("bkg_ggzz_4mu_FVIHP"   ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ggzz_4mu)   ,RooArgList()                ,1.0,1);
            
            FastVerticalInterpHistPdf2D plpdf_ggH_2e2mu  ("sig_ggH_2e2mu_FVIHP"  ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ggH_2e2mu)  ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_ggH_4e     ("sig_ggH_4e_FVIHP"     ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ggH_4e)     ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_ggH_4mu    ("sig_ggH_4mu_FVIHP"    ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ggH_4mu)    ,RooArgList()                ,1.0,1);
            
            FastVerticalInterpHistPdf2D plpdf_VBF_2e2mu  ("sig_VBF_2e2mu_FVIHP"  ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_VBF_2e2mu)  ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_VBF_4e     ("sig_VBF_4e_FVIHP"     ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_VBF_4e)     ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_VBF_4mu    ("sig_VBF_4mu_FVIHP"    ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_VBF_4mu)    ,RooArgList()                ,1.0,1);
            
            FastVerticalInterpHistPdf2D plpdf_zjets_2e2mu("bkg_zjets_2e2mu_FVIHP","",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,list_zjets_2e2mu            ,RooArgList(CMS_zz4l_bkgMELA),1.0,1);
            FastVerticalInterpHistPdf2D plpdf_zjets_4e   ("bkg_zjets_4e_FVIHP"   ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,list_zjets_4e               ,RooArgList(CMS_zz4l_bkgMELA),1.0,1);
            FastVerticalInterpHistPdf2D plpdf_zjets_4mu  ("bkg_zjets_4mu_FVIHP"  ,"",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,list_zjets_4mu              ,RooArgList(CMS_zz4l_bkgMELA),1.0,1);
            

            RooProdPdf bkg_qqzz_2e2mu_pdf_2D ("bkg_qqzz_2e2mu" , "", bkg_qqzz_2e2mu_pdf ,Conditional(plpdf_qqzz_2e2mu , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf bkg_qqzz_4e_pdf_2D    ("bkg_qqzz_4e"    , "", bkg_qqzz_4e_pdf    ,Conditional(plpdf_qqzz_4e    , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf bkg_qqzz_4mu_pdf_2D   ("bkg_qqzz_4mu"   , "", bkg_qqzz_4mu_pdf   ,Conditional(plpdf_qqzz_4mu   , RooArgSet(CMS_zz4l_melaLD)));
            
            RooProdPdf bkg_ggzz_2e2mu_pdf_2D ("bkg_ggzz_2e2mu" , "", bkg_ggzz_2e2mu_pdf ,Conditional(plpdf_ggzz_2e2mu , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf bkg_ggzz_4e_pdf_2D    ("bkg_ggzz_4e"    , "", bkg_ggzz_4e_pdf    ,Conditional(plpdf_ggzz_4e    , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf bkg_ggzz_4mu_pdf_2D   ("bkg_ggzz_4mu"   , "", bkg_ggzz_4mu_pdf   ,Conditional(plpdf_ggzz_4mu   , RooArgSet(CMS_zz4l_melaLD)));
            
            RooProdPdf bkg_zjets_2e2mu_pdf_2D("bkg_zjets_2e2mu", "", bkg_zjets_2e2mu_pdf,Conditional(plpdf_zjets_2e2mu, RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf bkg_zjets_4e_pdf_2D   ("bkg_zjets_4e"   , "", bkg_zjets_4e_pdf   ,Conditional(plpdf_zjets_4e   , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf bkg_zjets_4mu_pdf_2D  ("bkg_zjets_4mu"  , "", bkg_zjets_4mu_pdf  ,Conditional(plpdf_zjets_4mu  , RooArgSet(CMS_zz4l_melaLD)));
            
            RooProdPdf sig_ggH_2e2mu_pdf_2D  ("sig_ggH_2e2mu"  , "", sig_ggH_2e2mu_pdf  ,Conditional(plpdf_ggH_2e2mu  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ggH_4e_pdf_2D     ("sig_ggH_4e"     , "", sig_ggH_4e_pdf     ,Conditional(plpdf_ggH_4e     , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ggH_4mu_pdf_2D    ("sig_ggH_4mu"    , "", sig_ggH_4mu_pdf    ,Conditional(plpdf_ggH_4mu    , RooArgSet(CMS_zz4l_melaLD)));
            
            RooProdPdf sig_VBF_2e2mu_pdf_2D  ("sig_VBF_2e2mu"  , "", sig_VBF_2e2mu_pdf  ,Conditional(plpdf_VBF_2e2mu  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_VBF_4e_pdf_2D     ("sig_VBF_4e"     , "", sig_VBF_4e_pdf     ,Conditional(plpdf_VBF_4e     , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_VBF_4mu_pdf_2D    ("sig_VBF_4mu"    , "", sig_VBF_4mu_pdf    ,Conditional(plpdf_VBF_4mu    , RooArgSet(CMS_zz4l_melaLD)));
            
            w_2e2mu.import(bkg_qqzz_2e2mu_pdf_2D); 
            w_2e2mu.import(bkg_ggzz_2e2mu_pdf_2D); 
            w_2e2mu.import(bkg_zjets_2e2mu_pdf_2D); 
            w_2e2mu.import(sig_ggH_2e2mu_pdf_2D); 
            w_2e2mu.import(sig_VBF_2e2mu_pdf_2D);
            
            w_4e.import(bkg_qqzz_4e_pdf_2D); 
            w_4e.import(bkg_ggzz_4e_pdf_2D); 
            w_4e.import(bkg_zjets_4e_pdf_2D); 
            w_4e.import(sig_ggH_4e_pdf_2D); 
            w_4e.import(sig_VBF_4e_pdf_2D);
            
            w_4mu.import(bkg_qqzz_4mu_pdf_2D); 
            w_4mu.import(bkg_ggzz_4mu_pdf_2D); 
            w_4mu.import(bkg_zjets_4mu_pdf_2D); 
            w_4mu.import(sig_ggH_4mu_pdf_2D); 
            w_4mu.import(sig_VBF_4mu_pdf_2D);
 
        }
        
        w_2e2mu.writeToFile(workspace_2e2mu.c_str());
        w_4e.writeToFile(workspace_4e.c_str());
        w_4mu.writeToFile(workspace_4mu.c_str());

    }


    void writeParameters(std::string filename) {
        std::cout << "Creating " << filename << " ...." << std::endl;
        
        stringstream ss (stringstream::in | stringstream::out);
        ss << "#ifndef WRITEPARAMETERS_H\n"; 
        ss << "#define WRITEPARAMETERS_H\n"; 
        
        if (do7TeV) {
            ss << "float CBMeanP0_2e2mu_7TeV = "     << CBMeanP0_2e2mu << "\n"; 
            ss << "float CBMeanP1_2e2mu_7TeV = "     << CBMeanP1_2e2mu << "\n";
            ss << "float CBSigmaP0_2e2mu_7TeV = "    << CBSigmaP0_2e2mu << "\n";
            ss << "float CBSigmaP1_2e2mu_7TeV = "    << CBSigmaP1_2e2mu << "\n";
            ss << "float CBAlphaP0_2e2mu_7TeV = "    << CBAlphaP0_2e2mu << "\n";
            ss << "float CBNP0_2e2mu_7TeV = "        << CBNP0_2e2mu << "\n";
            ss << "float EfficiencyP0_2e2mu_7TeV = " << EfficiencyP0_2e2mu << "\n";
            ss << "float EfficiencyP1_2e2mu_7TeV = " << EfficiencyP1_2e2mu << "\n";
            ss << "float EfficiencyP2_2e2mu_7TeV = " << EfficiencyP2_2e2mu << "\n";
            ss << "float EfficiencyP3_2e2mu_7TeV = " << EfficiencyP3_2e2mu << "\n";
            ss << "\n";
            
            ss << "float CBMeanP0_4e_7TeV = "        << CBMeanP0_4e << "\n";
            ss << "float CBMeanP1_4e_7TeV = "        << CBMeanP1_4e << "\n";
            ss << "float CBSigmaP0_4e_7TeV = "       << CBSigmaP0_4e << "\n";
            ss << "float CBSigmaP1_4e_7TeV = "       << CBSigmaP1_4e << "\n";
            ss << "float CBAlphaP0_4e_7TeV = "       << CBAlphaP0_4e << "\n";
            ss << "float CBNP0_4e_7TeV = "           << CBNP0_4e << "\n";
            ss << "float EfficiencyP0_4e_7TeV = "    << EfficiencyP0_4e << "\n";
            ss << "float EfficiencyP1_4e_7TeV = "    << EfficiencyP1_4e << "\n";
            ss << "float EfficiencyP2_4e_7TeV = "    << EfficiencyP2_4e << "\n";
            ss << "float EfficiencyP3_4e_7TeV = "    << EfficiencyP3_4e << "\n";
            ss << "\n";
            
            ss << "float CBMeanP0_4mu_7TeV = "       << CBMeanP0_4mu << "\n";
            ss << "float CBMeanP1_4mu_7TeV = "       << CBMeanP1_4mu << "\n";
            ss << "float CBSigmaP0_4mu_7TeV = "      << CBSigmaP0_4mu << "\n";
            ss << "float CBSigmaP1_4mu_7TeV = "      << CBSigmaP1_4mu << "\n";
            ss << "float CBAlphaP0_4mu_7TeV = "      << CBAlphaP0_4mu << "\n";
            ss << "float CBNP0_4mu_7TeV = "          << CBNP0_4mu << "\n";
            ss << "float EfficiencyP0_4mu_7TeV = "   << EfficiencyP0_4mu << "\n";
            ss << "float EfficiencyP1_4mu_7TeV = "   << EfficiencyP1_4mu << "\n";
            ss << "float EfficiencyP2_4mu_7TeV = "   << EfficiencyP2_4mu << "\n";
            ss << "float EfficiencyP3_4mu_7TeV = "   << EfficiencyP3_4mu << "\n";
            ss << "\n";
        }
        
        else {
            ss << "float CBMeanP0_2e2mu_8TeV = "     << CBMeanP0_2e2mu << "\n";
            ss << "float CBMeanP1_2e2mu_8TeV = "     << CBMeanP1_2e2mu << "\n";
            ss << "float CBSigmaP0_2e2mu_8TeV = "    << CBSigmaP0_2e2mu << "\n";
            ss << "float CBSigmaP1_2e2mu_8TeV = "    << CBSigmaP1_2e2mu << "\n";
            ss << "float CBAlphaP0_2e2mu_8TeV = "    << CBAlphaP0_2e2mu << "\n";
            ss << "float CBNP0_2e2mu_8TeV = "        << CBNP0_2e2mu << "\n";
            ss << "float EfficiencyP0_2e2mu_8TeV = " << EfficiencyP0_2e2mu << "\n";
            ss << "float EfficiencyP1_2e2mu_8TeV = " << EfficiencyP1_2e2mu << "\n";
            ss << "float EfficiencyP2_2e2mu_8TeV = " << EfficiencyP2_2e2mu << "\n";
            ss << "float EfficiencyP3_2e2mu_8TeV = " << EfficiencyP3_2e2mu << "\n";
            ss << "\n";
        
            ss << "float CBMeanP0_4e_8TeV = "        << CBMeanP0_4e << "\n";
            ss << "float CBMeanP1_4e_8TeV = "        << CBMeanP1_4e << "\n";
            ss << "float CBSigmaP0_4e_8TeV = "       << CBSigmaP0_4e << "\n";
            ss << "float CBSigmaP1_4e_8TeV = "       << CBSigmaP1_4e << "\n";
            ss << "float CBAlphaP0_4e_8TeV = "       << CBAlphaP0_4e << "\n";
            ss << "float CBNP0_4e_8TeV = "           << CBNP0_4e << "\n";
            ss << "float EfficiencyP0_4e_8TeV = "    << EfficiencyP0_4e << "\n";
            ss << "float EfficiencyP1_4e_8TeV = "    << EfficiencyP1_4e << "\n";
            ss << "float EfficiencyP2_4e_8TeV = "    << EfficiencyP2_4e << "\n";
            ss << "float EfficiencyP3_4e_8TeV = "    << EfficiencyP3_4e << "\n";
            ss << "\n";
        
            ss << "float CBMeanP0_4mu_8TeV = "       << CBMeanP0_4mu << "\n";
            ss << "float CBMeanP1_4mu_8TeV = "       << CBMeanP1_4mu << "\n";
            ss << "float CBSigmaP0_4mu_8TeV = "      << CBSigmaP0_4mu << "\n";
            ss << "float CBSigmaP1_4mu_8TeV = "      << CBSigmaP1_4mu << "\n";
            ss << "float CBAlphaP0_4mu_8TeV = "      << CBAlphaP0_4mu << "\n";
            ss << "float CBNP0_4mu_8TeV = "          << CBNP0_4mu << "\n";
            ss << "float EfficiencyP0_4mu_8TeV = "   << EfficiencyP0_4mu << "\n";
            ss << "float EfficiencyP1_4mu_8TeV = "   << EfficiencyP1_4mu << "\n";
            ss << "float EfficiencyP2_4mu_8TeV = "   << EfficiencyP2_4mu << "\n";
            ss << "float EfficiencyP3_4mu_8TeV = "   << EfficiencyP3_4mu << "\n";
        }
        
        ss << "#endif\n"; 
        
        ofstream file;
        file.open (filename.c_str());
        file << ss.str();
        file.close();

    }



};


void doHZZAnalysis() {

    HiggsMassPointInfo hmpi7;
    hmpi7.lumi = 5.05;
    hmpi7.z1min = 40.;
    hmpi7.z2min = 12.;
    hmpi7.massLowBkgFit = 100.;
    hmpi7.massHighBkgFit = 600.;
    hmpi7.melacut = -1.0;
    hmpi7.do1D = true;
    hmpi7.doSS = true;
    hmpi7.do7TeV = true;
    hmpi7.treeFolder = "/home/avartak/CMS/Higgs/HZZ4L/CMSSW_4_2_8_patch7/src/WWAnalysis/AnalysisStep/trees/";
    hmpi7.melafilename = "mela2DShapes.root";

    init(hmpi7.do7TeV);

    FakeRateCalculator FR_7TeV(hmpi7.treeFolder+"hzzTree.root", hmpi7.do7TeV, 40, 120, 0.0, 0.0, true);
    hmpi7.FR = FR_7TeV;
    
    hmpi7.ymaker_data.fill(hmpi7.treeFolder+"hzzTree.root");
    hmpi7.ymaker_zxss.fill(hmpi7.treeFolder+"hzzTree.root"       , 1.0, hmpi7.FR, hmpi7.doSS);
    hmpi7.ymaker_qqzz.fill(hmpi7.treeFolder+"hzzTree_id121.root" , getBkgXsec(121)*hmpi7.lumi/evt_7TeV[121], 0.0, false);
    hmpi7.ymaker_qqzz.fill(hmpi7.treeFolder+"hzzTree_id122.root" , getBkgXsec(122)*hmpi7.lumi/evt_7TeV[122], 0.0, false);
    hmpi7.ymaker_qqzz.fill(hmpi7.treeFolder+"hzzTree_id123.root" , getBkgXsec(123)*hmpi7.lumi/evt_7TeV[123], 0.0, false);
    hmpi7.ymaker_qqzz.fill(hmpi7.treeFolder+"hzzTree_id124.root" , getBkgXsec(124)*hmpi7.lumi/evt_7TeV[124], 0.0, false);
    hmpi7.ymaker_qqzz.fill(hmpi7.treeFolder+"hzzTree_id125.root" , getBkgXsec(125)*hmpi7.lumi/evt_7TeV[125], 0.0, false);
    hmpi7.ymaker_qqzz.fill(hmpi7.treeFolder+"hzzTree_id126.root" , getBkgXsec(126)*hmpi7.lumi/evt_7TeV[126], 0.0, false);
    hmpi7.ymaker_ggzz.fill(hmpi7.treeFolder+"hzzTree_id101.root" , getBkgXsec(101)*hmpi7.lumi/evt_7TeV[101], 0.0, false);
    hmpi7.ymaker_ggzz.fill(hmpi7.treeFolder+"hzzTree_id100.root" , getBkgXsec(100)*hmpi7.lumi/evt_7TeV[100], 0.0, false);

    hmpi7.analyze(115, 100, 130, 200, 250, false, false);
    hmpi7.analyze(115, 100, 130, 200, 250, true , false);
    hmpi7.analyze(120, 100, 135, 201, 251, true , true );
    hmpi7.analyze(130, 110, 145, 202, 252, true , false);
    hmpi7.analyze(140, 120, 155, 203, 253, true , false);
    hmpi7.analyze(150, 130, 165, 204, 254, true , false);
    hmpi7.analyze(160, 140, 175, 205, 255, true , false);
    hmpi7.analyze(170, 150, 185, 206, 256, true , false);
    hmpi7.analyze(180, 160, 195, 207, 257, true , false);
    hmpi7.analyze(190, 170, 205, 208, 258, true , false);
    hmpi7.analyze(200, 180, 215, 209, 259, true , false);
    hmpi7.analyze(210, 190, 225, 210, 260, true , false);
    hmpi7.analyze(220, 200, 235, 211, 261, true , false);

    hmpi7.fitCBMean();
    hmpi7.fitCBSigma();
    hmpi7.fitCBAlpha();
    hmpi7.fitCBN();
    hmpi7.fitEfficiency();

    for (float i = 114.; i <= 160.; i += 1.) hmpi7.createCard(i, std::max<float>(i-20., 100.), i+15.);
    hmpi7.createCard(162., 142., 177.);
    hmpi7.createCard(164., 144., 179.);

    hmpi7.do1D = false;
    for (float i = 114.; i <= 160.; i += 1.) hmpi7.createCard(i, std::max<float>(i-20., 100.), i+15.);
    hmpi7.createCard(162., 142., 177.);
    hmpi7.createCard(164., 144., 179.);
    

    HiggsMassPointInfo hmpi8;
    hmpi8.lumi = 5.26; 
    hmpi8.z1min = 40.;
    hmpi8.z2min = 12.;
    hmpi8.massLowBkgFit = 100.;
    hmpi8.massHighBkgFit = 600.;
    hmpi8.melacut = -1.0;
    hmpi8.do1D = true;
    hmpi8.doSS = true;    
    hmpi8.do7TeV = false;
    hmpi8.treeFolder = "/home/avartak/CMS/Higgs/HZZ4L/CMSSW_5_2_4_patch4/src/WWAnalysis/AnalysisStep/trees/";
    hmpi8.melafilename = "mela2DShapes.root";

    init(hmpi8.do7TeV);


    FakeRateCalculator FR_8TeV(hmpi8.treeFolder+"hzzTree.root", hmpi8.do7TeV, 40, 120, 0.0, 0.0, true);
    hmpi8.FR = FR_8TeV;
    
    DataYieldMaker ymaker_data_8TeV;
    ZXYieldMaker   ymaker_zxss_8TeV;
    ZZYieldMaker   ymaker_qqzz_8TeV;
    ZZYieldMaker   ymaker_ggzz_8TeV;
    
    hmpi8.ymaker_data.fill(hmpi8.treeFolder+"hzzTree.root");
    hmpi8.ymaker_zxss.fill(hmpi8.treeFolder+"hzzTree.root"       , 1.0, hmpi8.FR, hmpi8.doSS);
    hmpi8.ymaker_qqzz.fill(hmpi8.treeFolder+"hzzTree_id102.root" , getBkgXsec(102)*hmpi8.lumi/evt_8TeV[102], 0.0, false);
    hmpi8.ymaker_qqzz.fill(hmpi8.treeFolder+"hzzTree_id103.root" , getBkgXsec(103)*hmpi8.lumi/evt_8TeV[103], 0.0, false);
    hmpi8.ymaker_qqzz.fill(hmpi8.treeFolder+"hzzTree_id104.root" , getBkgXsec(104)*hmpi8.lumi/evt_8TeV[104], 0.0, false);
    hmpi8.ymaker_qqzz.fill(hmpi8.treeFolder+"hzzTree_id105.root" , getBkgXsec(105)*hmpi8.lumi/evt_8TeV[105], 0.0, false);
    hmpi8.ymaker_qqzz.fill(hmpi8.treeFolder+"hzzTree_id106.root" , getBkgXsec(106)*hmpi8.lumi/evt_8TeV[106], 0.0, false);
    hmpi8.ymaker_qqzz.fill(hmpi8.treeFolder+"hzzTree_id107.root" , getBkgXsec(107)*hmpi8.lumi/evt_8TeV[107], 0.0, false);
    hmpi8.ymaker_ggzz.fill(hmpi8.treeFolder+"hzzTree_id101.root" , getBkgXsec(101)*hmpi8.lumi/evt_8TeV[101], 0.0, false);
    hmpi8.ymaker_ggzz.fill(hmpi8.treeFolder+"hzzTree_id100.root" , getBkgXsec(100)*hmpi8.lumi/evt_8TeV[100], 0.0, false);

    hmpi8.analyze(115, 100, 130, 1115, 2115, true, false);
    hmpi8.analyze(117, 100, 132, 1117, 2117, true, false);
    hmpi8.analyze(119, 100, 134, 1119, 2119, true, false);
    hmpi8.analyze(120, 100, 135, 1120, 2120, true, true );
    hmpi8.analyze(121, 101, 136, 1121, 2121, true, false);
    hmpi8.analyze(123, 103, 138, 1123, 2123, true, false);
    hmpi8.analyze(124, 104, 139, 1124, 2124, true, false);
    hmpi8.analyze(125, 105, 140, 1125, 2125, true, false);
    hmpi8.analyze(126, 106, 141, 1126, 2126, true, false);
    hmpi8.analyze(127, 107, 142, 1127, 2127, true, false);
    hmpi8.analyze(150, 130, 165, 1150, 2150, true, false);
    hmpi8.analyze(180, 160, 195, 1180, 2180, true, false);
    hmpi8.analyze(200, 180, 215, 1200, 2200, true, false);
    hmpi8.analyze(220, 200, 235, 1220, 2220, true, false);


    hmpi8.fitCBMean();
    hmpi8.fitCBSigma();
    hmpi8.fitCBAlpha();
    hmpi8.fitCBN();
    hmpi8.fitEfficiency();

    for (float i = 114.; i <= 160.; i += 1.) hmpi8.createCard(i, std::max<float>(i-20., 100.), i+15.);
    hmpi8.createCard(162., 142., 177.);
    hmpi8.createCard(164., 144., 179.);

    hmpi8.do1D = false;
    for (float i = 114.; i <= 160.; i += 1.) hmpi8.createCard(i, std::max<float>(i-20., 100.), i+15.);
    hmpi8.createCard(162., 142., 177.);
    hmpi8.createCard(164., 144., 179.);
}


