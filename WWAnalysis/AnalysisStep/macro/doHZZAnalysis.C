#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TGraphErrors.h>
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
#include <RooHistFunc.h>

#include "FakeRateCalculator.h"
#include "CardTemplate.h"
#include "FitMaker.h"
#include "YieldMaker.h"
#include "findAndReplace.h"

using namespace RooFit;

struct HiggsMassPointInfo {

    float lumi;
    float base_lumi;
    float melacut;
    float z1min;
    float z2min;
    float massLowBkgFit;
    float massHighBkgFit;
    bool  do1D;
    bool  do7TeV;
    bool  doFFT;
    std::string treeFolder;
    std::string melafilename;
    
    std::map<std::string, std::vector<float> > interpolationmap;

    DataYieldMaker ymaker_data;

    std::string getSignalCBMeanString(float m, int ch) {
        stringstream fss;
        fss << "( ";  

        if (!doFFT) fss << "@0 + ";

        if (do7TeV) {
            if (ch == 0) fss << "0.2500 - 0.002641*@0";
            if (ch == 1) fss << "0.2340 - 0.002254*@0";
            if (ch == 2) fss << "0.1671 - 0.001355*@0";
        }
        else {
            if (ch == 0) fss << "0.03633 - 0.000925*@0";
            if (ch == 1) fss << "0.3118  - 0.001246*@0";
            if (ch == 2) fss << "0.03999 + 0.000621*@0";
        }
        fss << " ) + @0*@1";
        return fss.str();
    }
    
    std::string getSignalCBSigmaString(float m, int ch) {
        stringstream fss;
        fss << "( ";  

        if (do7TeV) {
            if (ch == 0) fss << "-0.1151 + 0.01108*@0";
            if (ch == 1) fss << "0.9918  + 0.01221*@0";
            if (ch == 2) fss << "0.5971  + 0.01041*@0";
        }
        else {
            if (ch == 0) fss << "-0.09006+ 0.01053*@0";
            if (ch == 1) fss << "1.374   + 0.01468*@0";
            if (ch == 2) fss << "0.8079  + 0.01151*@0";
        }
        fss << " ) * (1+@1)";
        return fss.str();
    }
   
    std::string getSignalCBAlphaString(float m, int ch) {
        stringstream fss;
        
        if (do7TeV) {
            if (ch == 0) fss << "1.419";
            if (ch == 1) fss << "0.8779";
            if (ch == 2) fss << "0.9944";
        }
        else {
            if (ch == 0) fss << "1.44";
            if (ch == 1) fss << "1.012";
            if (ch == 2) fss << "1.069";
        }
        return fss.str();
    }

    std::string getSignalCBNString(float m, int ch) {
        stringstream fss;

        if (do7TeV) {
            if (ch == 0) fss << "2.58  - 0.004141*@0";
            if (ch == 1) fss << "6.87  - 0.01909*@0";
            if (ch == 2) fss << "4.364 - 0.007203*@0";
        }
        else {
            if (ch == 0) fss << "2.226 - 0.002341*@0";
            if (ch == 1) fss << "9.419 - 0.02975*@0";
            if (ch == 2) fss << "4.721 - 0.009716*@0";
        }
        return fss.str();
    }

    std::string getYieldEfficiencyString(float m, int ch) {
        stringstream fss;

        if (do7TeV) {
            if      (ch == 0) return "-0.3257 + 0.005203*@0 + (2.046e-5)*@0*@0 - (1.099e-7)*@0*@0*@0";
            else if (ch == 1) return "-0.4718 + 0.004796*@0 + (2.008e-5)*@0*@0 - (1.047e-7)*@0*@0*@0";
            else              return "-0.3608 + 0.004852*@0 + (1.930e-5)*@0*@0 - (1.044e-7)*@0*@0*@0";
        }
        else {
            if      (ch == 0) return "-0.4888 + 0.005325*@0 + (2.980e-5)*@0*@0 - (1.419e-7)*@0*@0*@0";
            else if (ch == 1) return "-0.5458 + 0.004162*@0 + (2.651e-5)*@0*@0 - (1.179e-7)*@0*@0*@0";
            else              return "-0.5142 + 0.004792*@0 + (2.806e-5)*@0*@0 - (1.323e-7)*@0*@0*@0";
        }

        return fss.str();
    }


    void createCard(float mass, float massLow, float massHigh, int ch) {

        std::string chstr;
        if (ch == 0) chstr = "4mu";
        if (ch == 1) chstr = "4e";
        if (ch == 2) chstr = "2e2mu";

        stringstream mass_str_ss;
        mass_str_ss << mass;
        std::string mass_str = mass_str_ss.str();
        
        std::cout << "Creating datacard for " << mass_str << " GeV mass point and channel " << chstr << " ... " << std::endl;
        
        std::string card_name   = do1D ? (std::string("card_1D_m")+mass_str+"_"+(do7TeV ? "7TeV_" : "8TeV_")) : (std::string("card_2D_m")+mass_str+"_"+(do7TeV ? "7TeV_" : "8TeV_"));
        card_name += chstr;
        std::string workspace = card_name+"_workspace.root";
        
        float yield_dt  = ymaker_data.getYield(ch, z1min, z2min, massLow, massHigh, melacut);
        if (ch == 2) yield_dt += ymaker_data.getYield(3, z1min, z2min, massLow, massHigh, melacut);
         
        RooWorkspace w("w", "");
        
        RooRealVar CMS_zz4l_melaLD ("CMS_zz4l_melaLD" , "MELA" ,   0,             1,            "");
        RooRealVar CMS_zz4l_mass_1D("CMS_zz4l_mass_1D", "M(4l)", massLow, massHigh, "GeV/c^{2}");
        if (doFFT) CMS_zz4l_mass_1D.setBins(100000, "fft");
        
        if (do1D) {
            RooArgSet argset_obs(CMS_zz4l_mass_1D, "argset_obs");
            RooDataSet data_obs("data_obs", "data_obs", argset_obs);
        
            ymaker_data.getDataSet1D(ch, z1min, z2min, massLow, massHigh, melacut, data_obs, CMS_zz4l_mass_1D);
            if (ch == 2) ymaker_data.getDataSet1D(3 , z1min, z2min, massLow, massHigh, melacut, data_obs, CMS_zz4l_mass_1D);
            
            w.import(data_obs);
        }
        
        else {
            RooArgSet argset_obs(CMS_zz4l_mass_1D, CMS_zz4l_melaLD, "argset_obs");
            RooDataSet data_obs("data_obs", "data_obs", argset_obs);
            
            ymaker_data.getDataSet2D(ch, z1min, z2min, massLow, massHigh, melacut, data_obs, CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
            if (ch == 2) ymaker_data.getDataSet2D(3, z1min, z2min, massLow, massHigh, melacut, data_obs, CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
            
            w.import(data_obs);
        }
        
        
        ///////////////////// Define parameters //////////////////////////////////

        float qqa0  = 0.;
        float qqa1  = 0.;
        float qqa2  = 0.;
        float qqa3  = 0.;
        float qqa4  = 0.;
        float qqa5  = 0.;
        float qqa6  = 0.;
        float qqa7  = 0.;
        float qqa8  = 0.;
        float qqa9  = 0.;
        float qqa10 = 0.;
        float qqa11 = 0.;
        float qqa12 = 0.;
        float qqa13 = 0.;

        float gga0  = 0.;
        float gga1  = 0.;
        float gga2  = 0.;
        float gga3  = 0.;
        float gga4  = 0.;
        float gga5  = 0.;
        float gga6  = 0.;
        float gga7  = 0.;
        float gga8  = 0.;
        float gga9  = 0.;

        float zxme  = 0.;
        float zxsi  = 0.;

        if (do7TeV) {

            if (ch == 0) {
                qqa0  = 107.581;
                qqa1  = 8.49052;
                qqa2  = 111.682;
                qqa3  = 0.0450007;
                qqa4  = 185.303;
                qqa5  = 8.64641;
                qqa6  = 36.2561;
                qqa7  = 0.0956393;
                qqa8  = 57.3221;
                qqa9  = 0.0490451;
                qqa10 = 100.677;
                qqa11 = -8.81871;
                qqa12 = 1770.75;
                qqa13 = 0.0693848;

                gga0  = 158.179;
                gga1  = 82.4557;
                gga2  = 113.057;
                gga3  = 0.0496534;
                gga4  = 183.759;
                gga5  = 9.2115;
                gga6  = 43.0377;
                gga7  = 0.456734;
                gga8  = 44.0894;
                gga9  = -0.186183;

                zxme  = 142.11;
                zxsi  = 26.5834;

            }

            else if (ch == 1) {
                qqa1  = 17.7472;
                qqa2  = 115.855;
                qqa3  = 0.0403639;
                qqa4  = 184.6;
                qqa5  = 11.4167;
                qqa6  = 35.0533;
                qqa7  = 0.109481;
                qqa8  = 59.8558;
                qqa9  = 0.0587938;
                qqa10 = 98.0922;
                qqa11 = -8.0287;
                qqa12 = 2358.66;
                qqa13 = 0.0666297;

                gga0  = 165.586;
                gga1  = 66.0532;
                gga2  = 107.971;
                gga3  = 0.0507375;
                gga4  = 184.216;
                gga5  = 11.3268;
                gga6  = 42.9324;
                gga7  = 0.45907;
                gga8  = 44.2638;
                gga9  = -0.187391;

                zxme  = 157.293;
                zxsi  = 22.9188;

            }

            else {
                qqa0  = 111.288;
                qqa1  = 13.9389;
                qqa2  = 114.437;
                qqa3  = 0.0410222;
                qqa4  = 185.424;
                qqa5  = 10.302;
                qqa6  = 33.866;
                qqa7  = 0.10935;
                qqa8  = 61.3538;
                qqa9  = 0.0536563;
                qqa10 = 96.897;
                qqa11 = -6.84111;
                qqa12 = 4310.44;
                qqa13 = 0.0754481;

                gga0  = 162.617;
                gga1  = 53.0264;
                gga2  = 107.535;
                gga3  = 0.0529968;
                gga4  = 184.481;
                gga5  = 10.0783;
                gga6  = 42.8566;
                gga7  = 0.437623;
                gga8  = 44.636;
                gga9  = -0.197371;

                zxme  = 153.855;
                zxsi  = 22.46;

            }

        }

        else {
            if (ch == 0) {
                qqa0  = 108.004;
                qqa1  = 8.95585;
                qqa2  = 116.539;
                qqa3  = 0.0446413;
                qqa4  = 185.882;
                qqa5  = 9.87494;
                qqa6  = 35.1431;
                qqa7  = 0.095469;
                qqa8  = 58.2072;
                qqa9  = 0.0481988;
                qqa10 = 100.493;
                qqa11 = -7.34781;
                qqa12 = 1847.28;
                qqa13 = 0.0768976;

                gga0  = 148.424;
                gga1  = 91.4164;
                gga2  = 119.695;
                gga3  = 0.0500115;
                gga4  = 184.823;
                gga5  = 8.09216;
                gga6  = 43.819;
                gga7  = 0.437734;
                gga8  = 42.8844;
                gga9  = -0.184409;

                zxme  = 135.702;
                zxsi  = 27.1273;

            }

            else if (ch == 1) {
                qqa0  = 113.928;
                qqa1  = 22.0867;
                qqa2  = 124.569;
                qqa3  = 0.0371985;
                qqa4  = 185.846;
                qqa5  = 12.9931;
                qqa6  = 34.8594;
                qqa7  = 0.124639;
                qqa8  = 59.9443;
                qqa9  = 0.0641501;
                qqa10 = 99.3032;
                qqa11 = -7.55575;
                qqa12 = 2092.07;
                qqa13 = 0.0671605;

                gga0  = 164.221;
                gga1  = 84.3258;
                gga2  = 129.772;
                gga3  = 0.0458068;
                gga4  = 185.791;
                gga5  = 11.7362;
                gga6  = 43.0303;
                gga7  = 0.493722;
                gga8  = 44.3683;
                gga9  = -0.169254;

                zxme  = 154.622;
                zxsi  = 24.2448;

            }
        
            else {
                qqa0  = 111.724;
                qqa1  = 12.8939;
                qqa2  = 119.094;
                qqa3  = 0.0404216;
                qqa4  = 185.375;
                qqa5  = 9.11556;
                qqa6  = 35.4589;
                qqa7  = 0.108189;
                qqa8  = 63.0395;
                qqa9  = 0.059775;
                qqa10 = 96.0745;
                qqa11 = -8.06897;
                qqa12 = 2131.03;
                qqa13 = 0.0811445;

                gga0  = 162.743;
                gga1  = 54.6051;
                gga2  = 124.825;
                gga3  = 0.0482855;
                gga4  = 187.618;
                gga5  = 11.3918;
                gga6  = 42.8017;
                gga7  = 0.470558;
                gga8  = 45.3276;
                gga9  = -0.18372;

                zxme  = 153.393;
                zxsi  = 23.4266;

            }
        }





        std::string tevstr = do7TeV ? "_7TeV" : "_8TeV"; 
        stringstream lumiss;
        lumiss << lumi;
        std::string lumistr = lumiss.str(); 
        RooRealVar qqzz_a0 (("bkg_qqzz_"+chstr+tevstr+"_a0" ).c_str(), "", qqa0 );
        RooRealVar qqzz_a1 (("bkg_qqzz_"+chstr+tevstr+"_a1" ).c_str(), "", qqa1 );
        RooRealVar qqzz_a2 (("bkg_qqzz_"+chstr+tevstr+"_a2" ).c_str(), "", qqa2 );
        RooRealVar qqzz_a3 (("bkg_qqzz_"+chstr+tevstr+"_a3" ).c_str(), "", qqa3 );
        RooRealVar qqzz_a4 (("bkg_qqzz_"+chstr+tevstr+"_a4" ).c_str(), "", qqa4 );
        RooRealVar qqzz_a5 (("bkg_qqzz_"+chstr+tevstr+"_a5" ).c_str(), "", qqa5 );
        RooRealVar qqzz_a6 (("bkg_qqzz_"+chstr+tevstr+"_a6" ).c_str(), "", qqa6 );
        RooRealVar qqzz_a7 (("bkg_qqzz_"+chstr+tevstr+"_a7" ).c_str(), "", qqa7 );
        RooRealVar qqzz_a8 (("bkg_qqzz_"+chstr+tevstr+"_a8" ).c_str(), "", qqa8 );
        RooRealVar qqzz_a9 (("bkg_qqzz_"+chstr+tevstr+"_a9" ).c_str(), "", qqa9 );
        RooRealVar qqzz_a10(("bkg_qqzz_"+chstr+tevstr+"_a10").c_str(), "", qqa10);
        RooRealVar qqzz_a11(("bkg_qqzz_"+chstr+tevstr+"_a11").c_str(), "", qqa11);
        RooRealVar qqzz_a12(("bkg_qqzz_"+chstr+tevstr+"_a12").c_str(), "", qqa12);
        RooRealVar qqzz_a13(("bkg_qqzz_"+chstr+tevstr+"_a13").c_str(), "", qqa13);
        
        RooRealVar ggzz_a0 (("bkg_ggzz_"+chstr+tevstr+"_a0" ).c_str(), "", gga0 );
        RooRealVar ggzz_a1 (("bkg_ggzz_"+chstr+tevstr+"_a1" ).c_str(), "", gga1 );
        RooRealVar ggzz_a2 (("bkg_ggzz_"+chstr+tevstr+"_a2" ).c_str(), "", gga2 );
        RooRealVar ggzz_a3 (("bkg_ggzz_"+chstr+tevstr+"_a3" ).c_str(), "", gga3 );
        RooRealVar ggzz_a4 (("bkg_ggzz_"+chstr+tevstr+"_a4" ).c_str(), "", gga4 );
        RooRealVar ggzz_a5 (("bkg_ggzz_"+chstr+tevstr+"_a5" ).c_str(), "", gga5 );
        RooRealVar ggzz_a6 (("bkg_ggzz_"+chstr+tevstr+"_a6" ).c_str(), "", gga6 );
        RooRealVar ggzz_a7 (("bkg_ggzz_"+chstr+tevstr+"_a7" ).c_str(), "", gga7 );
        RooRealVar ggzz_a8 (("bkg_ggzz_"+chstr+tevstr+"_a8" ).c_str(), "", gga8 );
        RooRealVar ggzz_a9 (("bkg_ggzz_"+chstr+tevstr+"_a9" ).c_str(), "", gga9 );
        
        RooRealVar zx_mean (("bkg_zjets_"+chstr+tevstr+"_mean_zx" ).c_str(), "", zxme);
        RooRealVar zx_sigma(("bkg_zjets_"+chstr+tevstr+"_sigma_zx").c_str(), "", zxsi);
        
        RooRealVar sig_mean_err    (("sig_"+chstr+"_mean_err"  +tevstr).c_str()  , "", 0., -10., 10.);
        RooRealVar sig_sigma_err   (("sig_"+chstr+"_sigma_err" +tevstr).c_str()  , "", 0., -10., 10.);

        RooRealVar masshiggs       ("MH", "", mass);

        std::string cs_scale_str = "";
        if (do7TeV) cs_scale_str += "0.5+0.5*TMath::Erf((@0 - 80.85)/50.42)";
        else cs_scale_str += "1";    
        RooFormulaVar cs_scale_z2_ggh  (("cs_scale_z2_ggH"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar cs_scale_z2_vbf  (("cs_scale_z2_VBF"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar cs_scale_z2_whi  (("cs_scale_z2_WHi"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar cs_scale_z2_zhi  (("cs_scale_z2_ZHi"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar cs_scale_z2_tth  (("cs_scale_z2_ttH"+tevstr).c_str()       , cs_scale_str.c_str()                      , RooArgList(masshiggs));
        RooFormulaVar yield_var_ggh    (("yield_eff_ggH_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch).c_str(), RooArgList(masshiggs));
        RooFormulaVar yield_var_vbf    (("yield_eff_VBF_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch).c_str(), RooArgList(masshiggs));
        RooFormulaVar yield_var_whi    (("yield_eff_WHi_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch).c_str(), RooArgList(masshiggs));
        RooFormulaVar yield_var_zhi    (("yield_eff_ZHi_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch).c_str(), RooArgList(masshiggs));
        RooFormulaVar yield_var_tth    (("yield_eff_ttH_"+chstr+tevstr).c_str()  , getYieldEfficiencyString(mass, ch).c_str(), RooArgList(masshiggs));

        TH1F* histxsecbrggh;
        TH1F* histxsecbrvbf;
        TH1F* histxsecbrwhi;
        TH1F* histxsecbrzhi;
        TH1F* histxsecbrtth;


        if (mass >= 100 && mass < 160) {
            Float_t gghxsecbry[60];
            Float_t vbfxsecbry[60];
            Float_t whixsecbry[60];
            Float_t zhixsecbry[60];
            Float_t tthxsecbry[60];
            
            
            for (float i = 100.; i < 160.; i = i+1.) {
                gghxsecbry[(int)(i-100.)] = getXsecggHByChannel(i, ch); 
                vbfxsecbry[(int)(i-100.)] = getXsecVBFByChannel(i, ch); 
                whixsecbry[(int)(i-100.)] = getXsecWHiByChannel(i, ch); 
                zhixsecbry[(int)(i-100.)] = getXsecZHiByChannel(i, ch); 
                tthxsecbry[(int)(i-100.)] = getXsecttHByChannel(i, ch); 
            }
            
            histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 60, 100., 160.);
            histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 60, 100., 160.);
            histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 60, 100., 160.);
            histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 60, 100., 160.);
            histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 60, 100., 160.);
            
            for (int i = 1; i <= 60; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
            for (int i = 1; i <= 60; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
            for (int i = 1; i <= 60; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
            for (int i = 1; i <= 60; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
            for (int i = 1; i <= 60; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
        }

        else if (mass >= 160 && mass < 290) {
            Float_t gghxsecbry[65];
            Float_t vbfxsecbry[65];
            Float_t whixsecbry[65];
            Float_t zhixsecbry[65];
            Float_t tthxsecbry[65];


            for (float i = 160.; i < 290.; i = i+2.) {
                gghxsecbry[((int)(i-160.))/2] = getXsecggHByChannel(i, ch);
                vbfxsecbry[((int)(i-160.))/2] = getXsecVBFByChannel(i, ch);
                whixsecbry[((int)(i-160.))/2] = getXsecWHiByChannel(i, ch);
                zhixsecbry[((int)(i-160.))/2] = getXsecZHiByChannel(i, ch);
                tthxsecbry[((int)(i-160.))/2] = getXsecttHByChannel(i, ch);
            }

            histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 65, 160., 290.);
            histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 65, 160., 290.);
            histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 65, 160., 290.);
            histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 65, 160., 290.);
            histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 65, 160., 290.);

            for (int i = 1; i <= 65; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
            for (int i = 1; i <= 65; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
            for (int i = 1; i <= 65; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
            for (int i = 1; i <= 65; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
            for (int i = 1; i <= 65; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
        }

        else if (mass >= 290 && mass < 350) {
            Float_t gghxsecbry[12];
            Float_t vbfxsecbry[12];
            Float_t whixsecbry[12];
            Float_t zhixsecbry[12];
            Float_t tthxsecbry[12];


            for (float i = 290.; i < 350.; i = i+5.) {
                gghxsecbry[((int)(i-290.))/5] = getXsecggHByChannel(i, ch);
                vbfxsecbry[((int)(i-290.))/5] = getXsecVBFByChannel(i, ch);
                whixsecbry[((int)(i-290.))/5] = (i <= 300) ? getXsecWHiByChannel(i, ch) : 0.0;
                zhixsecbry[((int)(i-290.))/5] = (i <= 300) ? getXsecZHiByChannel(i, ch) : 0.0;
                tthxsecbry[((int)(i-290.))/5] = (i <= 300) ? getXsecttHByChannel(i, ch) : 0.0;
            }

            histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 12, 290., 350);
            histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 12, 290., 350);
            histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 12, 290., 350);
            histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 12, 290., 350);
            histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 12, 290., 350);

            for (int i = 1; i <= 12; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
            for (int i = 1; i <= 12; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
            for (int i = 1; i <= 12; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
            for (int i = 1; i <= 12; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
            for (int i = 1; i <= 12; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
        }

        else if (mass >= 350 && mass < 400) {
            Float_t gghxsecbry[5];
            Float_t vbfxsecbry[5];
            Float_t whixsecbry[5];
            Float_t zhixsecbry[5];
            Float_t tthxsecbry[5];


            for (float i = 350.; i < 400.; i = i+10.) {
                gghxsecbry[((int)(i-350.))/10] = getXsecggHByChannel(i, ch);
                vbfxsecbry[((int)(i-350.))/10] = getXsecVBFByChannel(i, ch);
                whixsecbry[((int)(i-350.))/10] = 0.0;
                zhixsecbry[((int)(i-350.))/10] = 0.0;
                tthxsecbry[((int)(i-350.))/10] = 0.0;
            }

            histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 5, 350., 400.);
            histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 5, 350., 400.);
            histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 5, 350., 400.);
            histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 5, 350., 400.);
            histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 5, 350., 400.);

            for (int i = 1; i <= 5; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
            for (int i = 1; i <= 5; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
            for (int i = 1; i <= 5; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
            for (int i = 1; i <= 5; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
            for (int i = 1; i <= 5; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
        }

        else  {
            Float_t gghxsecbry[30];
            Float_t vbfxsecbry[30];
            Float_t whixsecbry[30];
            Float_t zhixsecbry[30];
            Float_t tthxsecbry[30];


            for (float i = 400.; i < 1000.; i = i+20.) {
                gghxsecbry[((int)(i-400.))/20] = getXsecggHByChannel(i, ch);
                vbfxsecbry[((int)(i-400.))/20] = getXsecVBFByChannel(i, ch);
                whixsecbry[((int)(i-400.))/20] = 0.;
                zhixsecbry[((int)(i-400.))/20] = 0.;
                tthxsecbry[((int)(i-400.))/20] = 0.;
            }

            histxsecbrggh = new TH1F(("histxsecbrggh_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);
            histxsecbrvbf = new TH1F(("histxsecbrvbf_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);
            histxsecbrwhi = new TH1F(("histxsecbrwhi_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);
            histxsecbrzhi = new TH1F(("histxsecbrzhi_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);
            histxsecbrtth = new TH1F(("histxsecbrtth_"+chstr+tevstr).c_str(), "", 30, 400., 1000.);

            for (int i = 1; i <= 30; i++) histxsecbrggh->SetBinContent(i, gghxsecbry[i-1]);
            for (int i = 1; i <= 30; i++) histxsecbrvbf->SetBinContent(i, vbfxsecbry[i-1]);
            for (int i = 1; i <= 30; i++) histxsecbrwhi->SetBinContent(i, whixsecbry[i-1]);
            for (int i = 1; i <= 30; i++) histxsecbrzhi->SetBinContent(i, zhixsecbry[i-1]);
            for (int i = 1; i <= 30; i++) histxsecbrtth->SetBinContent(i, tthxsecbry[i-1]);
        }


        RooDataHist dhistxsecbrggh(("rdhistxsecbr_ggH_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), histxsecbrggh);
        RooDataHist dhistxsecbrvbf(("rdhistxsecbr_VBF_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), histxsecbrvbf);
        RooDataHist dhistxsecbrwhi(("rdhistxsecbr_WHi_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), histxsecbrwhi);
        RooDataHist dhistxsecbrzhi(("rdhistxsecbr_ZHi_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), histxsecbrzhi);
        RooDataHist dhistxsecbrtth(("rdhistxsecbr_ttH_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), histxsecbrtth);

        RooHistFunc ggh_xsecbr(("xsecbr_ggH_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrggh, 1);
        RooHistFunc vbf_xsecbr(("xsecbr_VBF_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrvbf, 1);
        RooHistFunc whi_xsecbr(("xsecbr_WHi_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrwhi, 1);
        RooHistFunc zhi_xsecbr(("xsecbr_ZHi_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrzhi, 1);
        RooHistFunc tth_xsecbr(("xsecbr_ttH_"+chstr+tevstr).c_str(), "", RooArgList(masshiggs), dhistxsecbrtth, 1);


        RooRealVar ggh_gamma_BW    (("sig_ggh_"+chstr+tevstr+"_gamma_BW" ).c_str(), "", 1.0);
        RooFormulaVar ggh_mean_CB  (("sig_ggh_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString (mass, ch).c_str() , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar ggh_sigma_CB (("sig_ggh_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar ggh_alpha    (("sig_ggh_"+chstr+tevstr+"_alpha"    ).c_str(), getSignalCBAlphaString(mass, ch).c_str() , RooArgList(masshiggs));
        RooFormulaVar ggh_n        (("sig_ggh_"+chstr+tevstr+"_n"        ).c_str(), getSignalCBNString(mass, ch).c_str()     , RooArgList(masshiggs));
        RooFormulaVar ggh_norm     ("ggH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_ggh, ggh_xsecbr, yield_var_ggh));
        
        RooRealVar vbf_gamma_BW    (("sig_vbf_"+chstr+tevstr+"_gamma_BW" ).c_str(), "", 1.0);
        RooFormulaVar vbf_mean_CB  (("sig_vbf_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString(mass, ch).c_str()  , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar vbf_sigma_CB (("sig_vbf_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar vbf_alpha    (("sig_vbf_"+chstr+tevstr+"_alpha"    ).c_str(), getSignalCBAlphaString(mass, ch).c_str() , RooArgList(masshiggs));
        RooFormulaVar vbf_n        (("sig_vbf_"+chstr+tevstr+"_n"        ).c_str(), getSignalCBNString(mass, ch).c_str()     , RooArgList(masshiggs));
        RooFormulaVar vbf_norm     ("qqH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_vbf, vbf_xsecbr, yield_var_vbf));
      
        RooRealVar whi_gamma_BW    (("sig_whi_"+chstr+tevstr+"_gamma_BW" ).c_str(), "", 1.0);
        RooFormulaVar whi_mean_CB  (("sig_whi_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString (mass, ch).c_str() , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar whi_sigma_CB (("sig_whi_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar whi_alpha    (("sig_whi_"+chstr+tevstr+"_alpha"    ).c_str(), getSignalCBAlphaString(mass, ch).c_str() , RooArgList(masshiggs));
        RooFormulaVar whi_n        (("sig_whi_"+chstr+tevstr+"_n"        ).c_str(), getSignalCBNString(mass, ch).c_str()     , RooArgList(masshiggs));
        RooFormulaVar whi_norm     ("WH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_whi, whi_xsecbr, yield_var_whi));
 
        RooRealVar zhi_gamma_BW    (("sig_zhi_"+chstr+tevstr+"_gamma_BW" ).c_str(), "", 1.0);
        RooFormulaVar zhi_mean_CB  (("sig_zhi_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString(mass, ch).c_str()  , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar zhi_sigma_CB (("sig_zhi_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar zhi_alpha    (("sig_zhi_"+chstr+tevstr+"_alpha"    ).c_str(), getSignalCBAlphaString(mass, ch).c_str() , RooArgList(masshiggs));
        RooFormulaVar zhi_n        (("sig_zhi_"+chstr+tevstr+"_n"        ).c_str(), getSignalCBNString(mass, ch).c_str()     , RooArgList(masshiggs));
        RooFormulaVar zhi_norm     ("ZH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_zhi, zhi_xsecbr, yield_var_zhi));
 
        RooRealVar tth_gamma_BW    (("sig_tth_"+chstr+tevstr+"_gamma_BW" ).c_str(), "", 1.0);
        RooFormulaVar tth_mean_CB  (("sig_tth_"+chstr+tevstr+"_mean_CB"  ).c_str(), getSignalCBMeanString(mass, ch).c_str()  , RooArgList(masshiggs, sig_mean_err));
        RooFormulaVar tth_sigma_CB (("sig_tth_"+chstr+tevstr+"_sigma_CB" ).c_str(), getSignalCBSigmaString(mass, ch).c_str() , RooArgList(masshiggs, sig_sigma_err));
        RooFormulaVar tth_alpha    (("sig_tth_"+chstr+tevstr+"_alpha"    ).c_str(), getSignalCBAlphaString(mass, ch).c_str() , RooArgList(masshiggs));
        RooFormulaVar tth_n        (("sig_tth_"+chstr+tevstr+"_n"        ).c_str(), getSignalCBNString(mass, ch).c_str()     , RooArgList(masshiggs));
        RooFormulaVar tth_norm     ("ttH_norm"                                , ("@0*@1*@2*"+lumistr).c_str()            , RooArgList(cs_scale_z2_tth, tth_xsecbr, yield_var_tth));

        /////////// Set parameters to constant //////////////////
        qqzz_a0      .setConstant(kTRUE);
        qqzz_a1      .setConstant(kTRUE);
        qqzz_a2      .setConstant(kTRUE);
        qqzz_a3      .setConstant(kTRUE);
        qqzz_a4      .setConstant(kTRUE);
        qqzz_a5      .setConstant(kTRUE);
        qqzz_a6      .setConstant(kTRUE);
        qqzz_a7      .setConstant(kTRUE);
        qqzz_a8      .setConstant(kTRUE);
        qqzz_a9      .setConstant(kTRUE);
        qqzz_a10     .setConstant(kTRUE);
        qqzz_a11     .setConstant(kTRUE);
        qqzz_a12     .setConstant(kTRUE);
        qqzz_a13     .setConstant(kTRUE);
        
        ggzz_a0      .setConstant(kTRUE);
        ggzz_a1      .setConstant(kTRUE);
        ggzz_a2      .setConstant(kTRUE);
        ggzz_a3      .setConstant(kTRUE);
        ggzz_a4      .setConstant(kTRUE);
        ggzz_a5      .setConstant(kTRUE);
        ggzz_a6      .setConstant(kTRUE);
        ggzz_a7      .setConstant(kTRUE);
        ggzz_a8      .setConstant(kTRUE);
        ggzz_a9      .setConstant(kTRUE);
        
        zx_mean      .setConstant(kTRUE);
        zx_sigma     .setConstant(kTRUE);
       
        masshiggs    .setConstant(kTRUE);
        ggh_gamma_BW .setConstant(kTRUE);
        vbf_gamma_BW .setConstant(kTRUE);
        whi_gamma_BW .setConstant(kTRUE);
        zhi_gamma_BW .setConstant(kTRUE);
        tth_gamma_BW .setConstant(kTRUE);
        sig_mean_err .setConstant(kTRUE);
        sig_sigma_err.setConstant(kTRUE);
        
        
        ////////////////// Define the PDFs /////////////////////////////////
        
        const char* bkg_qqzz_pdf_name  = do1D ? "bkg_qqzz"  : "bkg_qqzz_1D" ;
        const char* bkg_ggzz_pdf_name  = do1D ? "bkg_ggzz"  : "bkg_ggzz_1D" ;
        const char* bkg_zjets_pdf_name = do1D ? "bkg_zjets" : "bkg_zjets_1D";
        const char* sig_ggH_pdf_name   = do1D ? "ggH"       : "ggH_1D"  ;
        const char* sig_VBF_pdf_name   = do1D ? "qqH"       : "qqH_1D"  ;
        const char* sig_WHi_pdf_name   = do1D ? "WH"        : "WH_1D"  ;
        const char* sig_ZHi_pdf_name   = do1D ? "ZH"        : "ZH_1D"  ;
        const char* sig_ttH_pdf_name   = do1D ? "ttH"       : "ttH_1D"  ;
        
        RooqqZZPdf_v2 bkg_qqzz_pdf(bkg_qqzz_pdf_name,"",CMS_zz4l_mass_1D,
                           qqzz_a0,
                           qqzz_a1,
                           qqzz_a2,
                           qqzz_a3,
                           qqzz_a4,
                           qqzz_a5,
                           qqzz_a6,
                           qqzz_a7,
                           qqzz_a8,
                           qqzz_a9,
                           qqzz_a10,
                           qqzz_a11,
                           qqzz_a12,
                           qqzz_a13);
        
        RooggZZPdf_v2 bkg_ggzz_pdf(bkg_ggzz_pdf_name,"",CMS_zz4l_mass_1D,
                           ggzz_a0,
                           ggzz_a1,
                           ggzz_a2,
                           ggzz_a3,
                           ggzz_a4,
                           ggzz_a5,
                           ggzz_a6,
                           ggzz_a7,
                           ggzz_a8,
                           ggzz_a9);
        
        RooLandau bkg_zjets_pdf(bkg_zjets_pdf_name, "",CMS_zz4l_mass_1D,zx_mean,zx_sigma);
        
        RooCBShape      signalCB_ggH   (doFFT ? "signalCB_ggH" : sig_ggH_pdf_name, "", CMS_zz4l_mass_1D, ggh_mean_CB,ggh_sigma_CB,ggh_alpha,ggh_n);
        RooRelBWUFParam signalBW_ggH   ("signalBW_ggH", "", CMS_zz4l_mass_1D, masshiggs,ggh_gamma_BW);
        RooFFTConvPdf   sig_ggH_pdf(sig_ggH_pdf_name, "", CMS_zz4l_mass_1D, signalBW_ggH,signalCB_ggH,2);
        sig_ggH_pdf.setBufferFraction(0.2);
        
        RooCBShape      signalCB_VBF   (doFFT ? "signalCB_VBF" : sig_VBF_pdf_name, "", CMS_zz4l_mass_1D, vbf_mean_CB,vbf_sigma_CB,vbf_alpha,vbf_n);
        RooRelBWUFParam signalBW_VBF   ("signalBW_VBF_", "", CMS_zz4l_mass_1D, masshiggs,vbf_gamma_BW);
        RooFFTConvPdf   sig_VBF_pdf(sig_VBF_pdf_name, "", CMS_zz4l_mass_1D, signalBW_VBF,signalCB_VBF,2);
        sig_VBF_pdf.setBufferFraction(0.2);
       
        RooCBShape      signalCB_WHi   (doFFT ? "signalCB_WHi" : sig_WHi_pdf_name, "", CMS_zz4l_mass_1D, whi_mean_CB,whi_sigma_CB,whi_alpha,whi_n);
        RooRelBWUFParam signalBW_WHi   ("signalBW_WHi", "", CMS_zz4l_mass_1D, masshiggs,whi_gamma_BW);
        RooFFTConvPdf   sig_WHi_pdf(sig_WHi_pdf_name, "", CMS_zz4l_mass_1D, signalBW_WHi,signalCB_WHi,2);
        sig_WHi_pdf.setBufferFraction(0.2);

        RooCBShape      signalCB_ZHi   (doFFT ? "signalCB_ZHi" : sig_ZHi_pdf_name, "", CMS_zz4l_mass_1D, zhi_mean_CB,zhi_sigma_CB,zhi_alpha,zhi_n);
        RooRelBWUFParam signalBW_ZHi   ("signalBW_ZHi_", "", CMS_zz4l_mass_1D, masshiggs,zhi_gamma_BW);
        RooFFTConvPdf   sig_ZHi_pdf(sig_ZHi_pdf_name, "", CMS_zz4l_mass_1D, signalBW_ZHi,signalCB_ZHi,2);
        sig_ZHi_pdf.setBufferFraction(0.2);
 
        RooCBShape      signalCB_ttH   (doFFT ? "signalCB_ttH" : sig_ttH_pdf_name, "", CMS_zz4l_mass_1D, tth_mean_CB,tth_sigma_CB,tth_alpha,tth_n);
        RooRelBWUFParam signalBW_ttH   ("signalBW_ttH_", "", CMS_zz4l_mass_1D, masshiggs,tth_gamma_BW);
        RooFFTConvPdf   sig_ttH_pdf(sig_ttH_pdf_name, "", CMS_zz4l_mass_1D, signalBW_ttH,signalCB_ttH,2);
        sig_ttH_pdf.setBufferFraction(0.2);


        //////////////// Yields and printing cards /////////////////////////

        RooRealVar CMS_zz4l_mass_norm("CMS_zz4l_mass_norm", "M(4l)", massLowBkgFit, massHighBkgFit, "GeV/c^{2}");
        CMS_zz4l_mass_norm.setRange("card", massLow      , massHigh);
        CMS_zz4l_mass_norm.setRange("full", massLowBkgFit, massHighBkgFit);

        RooqqZZPdf_v2 bkg_qqzz_norm((std::string(bkg_qqzz_pdf_name)+"_normalization").c_str(),"",CMS_zz4l_mass_norm,
                           qqzz_a0,
                           qqzz_a1,
                           qqzz_a2,
                           qqzz_a3,
                           qqzz_a4,
                           qqzz_a5,
                           qqzz_a6,
                           qqzz_a7,
                           qqzz_a8,
                           qqzz_a9,
                           qqzz_a10,
                           qqzz_a11,
                           qqzz_a12,
                           qqzz_a13);

        RooggZZPdf_v2 bkg_ggzz_norm((std::string(bkg_ggzz_pdf_name)+"_normalization").c_str(),"",CMS_zz4l_mass_norm,
                           ggzz_a0,
                           ggzz_a1,
                           ggzz_a2,
                           ggzz_a3,
                           ggzz_a4,
                           ggzz_a5,
                           ggzz_a6,
                           ggzz_a7,
                           ggzz_a8,
                           ggzz_a9);

        RooLandau bkg_zjet_norm((std::string(bkg_zjets_pdf_name)+"_normalization").c_str(), "",CMS_zz4l_mass_norm,zx_mean,zx_sigma);
        
        double qqzz_fullyield = 0.0;
        double ggzz_fullyield = 0.0;
        double zjet_fullyield = 0.0;

        if (ch == 0) {
            qqzz_fullyield = do7TeV ? 21.2 : 25.0;
            ggzz_fullyield = do7TeV ? 1.17 : 0.63;
            zjet_fullyield = do7TeV ? 0.72 : 1.30;
        }

        else if (ch == 1) {
            qqzz_fullyield = do7TeV ? 13.7 : 14.4;
            ggzz_fullyield = do7TeV ? 0.83 : 0.41;
            zjet_fullyield = do7TeV ? 2.20 : 2.70;
        }

        else {
            qqzz_fullyield = do7TeV ? 33.0 : 37.2;
            ggzz_fullyield = do7TeV ? 1.60 : 1.00;
            zjet_fullyield = do7TeV ? 2.90 : 3.80;
        }

        float yield_qq = qqzz_fullyield * ((bkg_qqzz_norm.createIntegral(RooArgSet(CMS_zz4l_mass_norm), Range("card")))->getVal()) / ((bkg_qqzz_norm.createIntegral(RooArgSet(CMS_zz4l_mass_norm), Range("full")))->getVal());
        float yield_gg = ggzz_fullyield * ((bkg_ggzz_norm.createIntegral(RooArgSet(CMS_zz4l_mass_norm), Range("card")))->getVal()) / ((bkg_ggzz_norm.createIntegral(RooArgSet(CMS_zz4l_mass_norm), Range("full")))->getVal());
        float yield_zj = zjet_fullyield * ((bkg_zjet_norm.createIntegral(RooArgSet(CMS_zz4l_mass_norm), Range("card")))->getVal()) / ((bkg_zjet_norm.createIntegral(RooArgSet(CMS_zz4l_mass_norm), Range("full")))->getVal());


        std::string card   = createCardTemplate(do7TeV, ch, do1D, workspace.c_str());

        std::string binname;
        if (ch == 0) binname = "a1";
        if (ch == 1) binname = "a2";
        if (ch == 2) binname = "a3";

        card = findAndReplace(card, "GGZZ_PDF"       , getGGZZPDFUncertainty7TeV(mass));
        card = findAndReplace(card, "QQZZ_PDF"       , getQQZZPDFUncertainty7TeV(mass));
        card = findAndReplace(card, "GGZZ_QCD"       , getGGZZQCDScaleUncertainty7TeV(mass));
        card = findAndReplace(card, "QQZZ_QCD"       , getQQZZQCDScaleUncertainty7TeV(mass));
        card = findAndReplace(card, "GGH_PDF"        , getggHPDFUncertainty(mass, true), getggHPDFUncertainty(mass, false));
        card = findAndReplace(card, "VBF_PDF"        , getVBFPDFUncertainty(mass, true), getVBFPDFUncertainty(mass, false));
        card = findAndReplace(card, "WHI_PDF"        , getWHiPDFUncertainty(mass, true), getWHiPDFUncertainty(mass, false));
        card = findAndReplace(card, "ZHI_PDF"        , getZHiPDFUncertainty(mass, true), getZHiPDFUncertainty(mass, false));
        card = findAndReplace(card, "TTH_PDF"        , getttHPDFUncertainty(mass, true), getttHPDFUncertainty(mass, false));
        card = findAndReplace(card, "GGH_QCD"        , getggHQCDScaleUncertainty(mass, true), getggHQCDScaleUncertainty(mass, false));
        card = findAndReplace(card, "VBF_QCD"        , getVBFQCDScaleUncertainty(mass, true), getVBFQCDScaleUncertainty(mass, false));
        card = findAndReplace(card, "WHI_QCD"        , getWHiQCDScaleUncertainty(mass, true), getWHiQCDScaleUncertainty(mass, false));
        card = findAndReplace(card, "ZHI_QCD"        , getZHiQCDScaleUncertainty(mass, true), getZHiQCDScaleUncertainty(mass, false));
        card = findAndReplace(card, "TTH_QCD"        , getttHQCDScaleUncertainty(mass, true), getttHQCDScaleUncertainty(mass, false));
        card = findAndReplace(card, "HIGH_MH"        , mass < 200 ? 1 : 1+(mass/1000.));
        card = findAndReplace(card, "SIG_GGH_YIELD"  , 1);
        card = findAndReplace(card, "SIG_VBF_YIELD"  , 1);
        card = findAndReplace(card, "SIG_WHI_YIELD"  , 1);
        card = findAndReplace(card, "SIG_ZHI_YIELD"  , 1);
        card = findAndReplace(card, "SIG_TTH_YIELD"  , 1);
        card = findAndReplace(card, "BKG_QQZZ_YIELD" , yield_qq);
        card = findAndReplace(card, "BKG_GGZZ_YIELD" , yield_gg);
        card = findAndReplace(card, "BKG_ZJETS_YIELD", yield_zj);
        card = findAndReplace(card, "BIN"            , binname);
        card = findAndReplace(card, "OBS"            , yield_dt);

        ofstream file;
        file.open ((card_name +".txt").c_str());
        file << card;
        file.close();



 
        w.import(ggh_norm);
        w.import(vbf_norm);
        w.import(whi_norm);
        w.import(zhi_norm);
        w.import(tth_norm);
        if (do1D) { 
            w.import(bkg_qqzz_pdf);
            w.import(bkg_ggzz_pdf);
            w.import(bkg_zjets_pdf);
            if (doFFT) {
                w.import(sig_ggH_pdf);
                w.import(sig_VBF_pdf);
                w.import(sig_WHi_pdf);
                w.import(sig_ZHi_pdf);
                w.import(sig_ttH_pdf);
            }        
            else {
                w.import(signalCB_ggH);
                w.import(signalCB_VBF);
                w.import(signalCB_WHi);
                w.import(signalCB_ZHi);
                w.import(signalCB_ttH);
            }        
        }

        else {
            TFile melafile(melafilename.c_str());
            RooArgList v2dList(CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
            RooArgSet  v2dSet (CMS_zz4l_mass_1D, CMS_zz4l_melaLD);
            
            std::string chstrsmall;
            if (ch == 0) chstrsmall = "mm";
            if (ch == 1) chstrsmall = "ee";
            if (ch == 2) chstrsmall = "em";

            TH2F* melashape_zz = (TH2F*)(melafile.Get(("hist2D_bkg_"+chstrsmall).c_str()));
            TH2F* melashape_si = (TH2F*)(melafile.Get(("hist2D_sig_"+chstrsmall).c_str()));
            TH2F* melashape_zu = (TH2F*)(melafile.Get(("hist2D_bkg_"+chstrsmall+"_up").c_str()));
            TH2F* melashape_zd = (TH2F*)(melafile.Get(("hist2D_bkg_"+chstrsmall+"_dn").c_str()));
            
            RooDataHist rhist_qqzz (("rhist_qqzz_" +chstr+tevstr).c_str(), "", v2dList, melashape_zz);
            RooDataHist rhist_ggzz (("rhist_ggzz_" +chstr+tevstr).c_str(), "", v2dList, melashape_zz);
            RooDataHist rhist_zjets(("rhist_zjets_"+chstr+tevstr).c_str(), "", v2dList, melashape_zz);
            RooDataHist rhist_ggH  (("rhist_ggH_"  +chstr+tevstr).c_str(), "", v2dList, melashape_si);
            RooDataHist rhist_VBF  (("rhist_VBF_"  +chstr+tevstr).c_str(), "", v2dList, melashape_si);
            RooDataHist rhist_WHi  (("rhist_WHi_"  +chstr+tevstr).c_str(), "", v2dList, melashape_si);
            RooDataHist rhist_ZHi  (("rhist_ZHi_"  +chstr+tevstr).c_str(), "", v2dList, melashape_si);
            RooDataHist rhist_ttH  (("rhist_ttH_"  +chstr+tevstr).c_str(), "", v2dList, melashape_si);
            RooDataHist rhist_zjup (("rhist_zjup_" +chstr+tevstr).c_str(), "", v2dList, melashape_zu);
            RooDataHist rhist_zjdn (("rhist_zjdn_" +chstr+tevstr).c_str(), "", v2dList, melashape_zd);
            RooHistPdf rpdf_qqzz (("bkg_qqzz_mela2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_qqzz);
            RooHistPdf rpdf_ggzz (("bkg_ggzz_mela2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_ggzz);
            RooHistPdf rpdf_zjets(("bkg_zjets_mela2D_pdf_"+chstr+tevstr).c_str(), "", v2dSet , rhist_zjets);
            RooHistPdf rpdf_ggH  (("sig_ggH_mela2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet , rhist_ggH);
            RooHistPdf rpdf_VBF  (("sig_VBF_mela2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet , rhist_VBF);
            RooHistPdf rpdf_WHi  (("sig_WHi_mela2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet , rhist_WHi);
            RooHistPdf rpdf_ZHi  (("sig_ZHi_mela2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet , rhist_ZHi);
            RooHistPdf rpdf_ttH  (("sig_ttH_mela2D_pdf_"  +chstr+tevstr).c_str(), "", v2dSet , rhist_ttH);
            RooHistPdf rpdf_zjup (("bkg_zjup_mela2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_zjup);
            RooHistPdf rpdf_zjdn (("bkg_zjdn_mela2D_pdf_" +chstr+tevstr).c_str(), "", v2dSet , rhist_zjdn);

            RooArgList list_zjets;
            list_zjets.add(rpdf_zjets);
            list_zjets.add(rpdf_zjup);
            list_zjets.add(rpdf_zjdn);

            RooRealVar CMS_zz4l_bkgMELA("CMS_zz4l_bkgMELA" ,"" ,0,-10,10); 
 
            FastVerticalInterpHistPdf2D plpdf_qqzz (("bkg_qqzz_FVIHP_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_qqzz) ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_ggzz (("bkg_ggzz_FVIHP_" +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ggzz) ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_ggH  (("sig_ggH_FVIHP_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ggH)  ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_VBF  (("sig_VBF_FVIHP_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_VBF)  ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_WHi  (("sig_WHi_FVIHP_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_WHi)  ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_ZHi  (("sig_ZHi_FVIHP_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ZHi)  ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_ttH  (("sig_ttH_FVIHP_"  +chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,RooArgList(rpdf_ttH)  ,RooArgList()                ,1.0,1);
            FastVerticalInterpHistPdf2D plpdf_zjets(("bkg_zjets_FVIHP_"+chstr+tevstr).c_str(), "",CMS_zz4l_mass_1D,CMS_zz4l_melaLD,true,list_zjets            ,RooArgList(CMS_zz4l_bkgMELA),1.0,1);

            RooProdPdf bkg_qqzz_pdf_2D ("bkg_qqzz" , "", bkg_qqzz_pdf  ,Conditional(plpdf_qqzz , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf bkg_ggzz_pdf_2D ("bkg_ggzz" , "", bkg_ggzz_pdf  ,Conditional(plpdf_ggzz , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf bkg_zjets_pdf_2D("bkg_zjets", "", bkg_zjets_pdf ,Conditional(plpdf_zjets, RooArgSet(CMS_zz4l_melaLD))); 
           
            if (doFFT) { 
            RooProdPdf sig_ggH_pdf_2D  ("ggH", "", sig_ggH_pdf   ,Conditional(plpdf_ggH  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_VBF_pdf_2D  ("qqH", "", sig_VBF_pdf   ,Conditional(plpdf_VBF  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_WHi_pdf_2D  ("WH",  "", sig_WHi_pdf   ,Conditional(plpdf_WHi  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ZHi_pdf_2D  ("ZH",  "", sig_ZHi_pdf   ,Conditional(plpdf_ZHi  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ttH_pdf_2D  ("ttH", "", sig_ttH_pdf   ,Conditional(plpdf_ttH  , RooArgSet(CMS_zz4l_melaLD))); 

            w.import(bkg_qqzz_pdf_2D); 
            w.import(bkg_ggzz_pdf_2D); 
            w.import(bkg_zjets_pdf_2D); 
            w.import(sig_ggH_pdf_2D); 
            w.import(sig_VBF_pdf_2D);
            w.import(sig_WHi_pdf_2D); 
            w.import(sig_ZHi_pdf_2D);
            w.import(sig_ttH_pdf_2D);
            }
            else {
            RooProdPdf sig_ggH_pdf_2D  ("ggH", "", signalCB_ggH  ,Conditional(plpdf_ggH  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_VBF_pdf_2D  ("qqH", "", signalCB_VBF  ,Conditional(plpdf_VBF  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_WHi_pdf_2D  ("WH",  "", signalCB_WHi  ,Conditional(plpdf_WHi  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ZHi_pdf_2D  ("ZH",  "", signalCB_ZHi  ,Conditional(plpdf_ZHi  , RooArgSet(CMS_zz4l_melaLD))); 
            RooProdPdf sig_ttH_pdf_2D  ("ttH", "", signalCB_ttH  ,Conditional(plpdf_ttH  , RooArgSet(CMS_zz4l_melaLD))); 
            w.import(bkg_qqzz_pdf_2D); 
            w.import(bkg_ggzz_pdf_2D); 
            w.import(bkg_zjets_pdf_2D); 
            w.import(sig_ggH_pdf_2D); 
            w.import(sig_VBF_pdf_2D);
            w.import(sig_WHi_pdf_2D); 
            w.import(sig_ZHi_pdf_2D);
            w.import(sig_ttH_pdf_2D);
            }

        }
        
        w.writeToFile(workspace.c_str());

        delete histxsecbrggh;
        delete histxsecbrvbf;
        delete histxsecbrwhi;
        delete histxsecbrzhi;
        delete histxsecbrtth;

    }

    float getMassCut(float cardmass, bool low) {
        std::map<float, float> higgswidth;
        higgswidth[90.0] = 2.22E-03;
        higgswidth[95.0] = 2.35E-03;
        higgswidth[100.0] = 2.49E-03;
        higgswidth[105.0] = 2.65E-03;
        higgswidth[110.0] = 2.85E-03;
        higgswidth[110.5] = 2.87E-03;
        higgswidth[111.0] = 2.90E-03;
        higgswidth[111.5] = 2.92E-03;
        higgswidth[112.0] = 2.95E-03;
        higgswidth[112.5] = 2.98E-03;
        higgswidth[113.0] = 3.00E-03;
        higgswidth[113.5] = 3.03E-03;
        higgswidth[114.0] = 3.06E-03;
        higgswidth[114.5] = 3.09E-03;
        higgswidth[115.0] = 3.12E-03;
        higgswidth[115.5] = 3.16E-03;
        higgswidth[116.0] = 3.19E-03;
        higgswidth[116.5] = 3.23E-03;
        higgswidth[117.0] = 3.26E-03;
        higgswidth[117.5] = 3.30E-03;
        higgswidth[118.0] = 3.34E-03;
        higgswidth[118.5] = 3.37E-03;
        higgswidth[119.0] = 3.41E-03;
        higgswidth[119.5] = 3.46E-03;
        higgswidth[120.0] = 3.50E-03;
        higgswidth[120.5] = 3.55E-03;
        higgswidth[121.0] = 3.60E-03;
        higgswidth[121.5] = 3.65E-03;
        higgswidth[122.0] = 3.70E-03;
        higgswidth[122.5] = 3.76E-03;
        higgswidth[123.0] = 3.82E-03;
        higgswidth[123.5] = 3.88E-03;
        higgswidth[124.0] = 3.94E-03;
        higgswidth[124.5] = 4.00E-03;
        higgswidth[125.0] = 4.07E-03;
        higgswidth[125.5] = 4.14E-03;
        higgswidth[126.0] = 4.21E-03;
        higgswidth[126.5] = 4.29E-03;
        higgswidth[127.0] = 4.37E-03;
        higgswidth[127.5] = 4.45E-03;
        higgswidth[128.0] = 4.53E-03;
        higgswidth[128.5] = 4.62E-03;
        higgswidth[129.0] = 4.71E-03;
        higgswidth[129.5] = 4.81E-03;
        higgswidth[130.0] = 4.91E-03;
        higgswidth[130.5] = 5.01E-03;
        higgswidth[131.0] = 5.12E-03;
        higgswidth[131.5] = 5.23E-03;
        higgswidth[132.0] = 5.35E-03;
        higgswidth[132.5] = 5.47E-03;
        higgswidth[133.0] = 5.60E-03;
        higgswidth[133.5] = 5.74E-03;
        higgswidth[134.0] = 5.88E-03;
        higgswidth[134.5] = 6.03E-03;
        higgswidth[135.0] = 6.18E-03;
        higgswidth[135.5] = 6.34E-03;
        higgswidth[136.0] = 6.51E-03;
        higgswidth[136.5] = 6.69E-03;
        higgswidth[137.0] = 6.87E-03;
        higgswidth[137.5] = 7.06E-03;
        higgswidth[138.0] = 7.27E-03;
        higgswidth[138.5] = 7.48E-03;
        higgswidth[139.0] = 7.70E-03;
        higgswidth[139.5] = 7.93E-03;
        higgswidth[140.0] = 8.18E-03;
        higgswidth[141.0] = 8.70E-03;
        higgswidth[142.0] = 9.28E-03;
        higgswidth[143.0] = 9.93E-03;
        higgswidth[144.0] = 1.07E-02;
        higgswidth[145.0] = 1.15E-02;
        higgswidth[146.0] = 1.23E-02;
        higgswidth[147.0] = 1.34E-02;
        higgswidth[148.0] = 1.45E-02;
        higgswidth[149.0] = 1.58E-02;
        higgswidth[150.0] = 1.73E-02;
        higgswidth[151.0] = 1.91E-02;
        higgswidth[152.0] = 2.11E-02;
        higgswidth[153.0] = 2.36E-02;
        higgswidth[154.0] = 2.66E-02;
        higgswidth[155.0] = 3.03E-02;
        higgswidth[156.0] = 3.51E-02;
        higgswidth[157.0] = 4.14E-02;
        higgswidth[158.0] = 5.02E-02;
        higgswidth[159.0] = 6.32E-02;
        higgswidth[160.0] = 8.31E-02;
        higgswidth[162.0] = 1.47E-01;
        higgswidth[164.0] = 2.15E-01;
        higgswidth[165.0] = 2.46E-01;
        higgswidth[166.0] = 2.76E-01;
        higgswidth[168.0] = 3.30E-01;
        higgswidth[170.0] = 3.80E-01;
        higgswidth[172.0] = 4.29E-01;
        higgswidth[174.0] = 4.77E-01;
        higgswidth[175.0] = 5.01E-01;
        higgswidth[176.0] = 5.25E-01;
        higgswidth[178.0] = 5.75E-01;
        higgswidth[180.0] = 6.31E-01;
        higgswidth[182.0] = 7.00E-01;
        higgswidth[184.0] = 7.88E-01;
        higgswidth[185.0] = 8.32E-01;
        higgswidth[186.0] = 8.76E-01;
        higgswidth[188.0] = 9.60E-01;
        higgswidth[190.0] = 1.04E+00;
        higgswidth[192.0] = 1.12E+00;
        higgswidth[194.0] = 1.20E+00;
        higgswidth[195.0] = 1.24E+00;
        higgswidth[196.0] = 1.28E+00;
        higgswidth[198.0] = 1.35E+00;
        higgswidth[200.0] = 1.43E+00;
        higgswidth[202.0] = 1.51E+00;
        higgswidth[204.0] = 1.59E+00;
        higgswidth[206.0] = 1.68E+00;
        higgswidth[208.0] = 1.76E+00;
        higgswidth[210.0] = 1.85E+00;
        higgswidth[212.0] = 1.93E+00;
        higgswidth[214.0] = 2.02E+00;
        higgswidth[216.0] = 2.12E+00;
        higgswidth[218.0] = 2.21E+00;
        higgswidth[220.0] = 2.31E+00;
        higgswidth[222.0] = 2.40E+00;
        higgswidth[224.0] = 2.50E+00;
        higgswidth[226.0] = 2.61E+00;
        higgswidth[228.0] = 2.71E+00;
        higgswidth[230.0] = 2.82E+00;
        higgswidth[232.0] = 2.93E+00;
        higgswidth[234.0] = 3.04E+00;
        higgswidth[236.0] = 3.16E+00;
        higgswidth[238.0] = 3.27E+00;
        higgswidth[240.0] = 3.40E+00;
        higgswidth[242.0] = 3.52E+00;
        higgswidth[244.0] = 3.64E+00;
        higgswidth[246.0] = 3.77E+00;
        higgswidth[248.0] = 3.91E+00;
        higgswidth[250.0] = 4.04E+00;
        higgswidth[252.0] = 4.18E+00;
        higgswidth[254.0] = 4.32E+00;
        higgswidth[256.0] = 4.46E+00;
        higgswidth[258.0] = 4.61E+00;
        higgswidth[260.0] = 4.76E+00;
        higgswidth[262.0] = 4.91E+00;
        higgswidth[264.0] = 5.07E+00;
        higgswidth[266.0] = 5.23E+00;
        higgswidth[268.0] = 5.39E+00;
        higgswidth[270.0] = 5.55E+00;
        higgswidth[272.0] = 5.72E+00;
        higgswidth[274.0] = 5.89E+00;
        higgswidth[276.0] = 6.07E+00;
        higgswidth[278.0] = 6.25E+00;
        higgswidth[280.0] = 6.43E+00;
        higgswidth[282.0] = 6.61E+00;
        higgswidth[284.0] = 6.80E+00;
        higgswidth[286.0] = 6.99E+00;
        higgswidth[288.0] = 7.19E+00;
        higgswidth[290.0] = 7.39E+00;
        higgswidth[295.0] = 7.90E+00;
        higgswidth[300.0] = 8.43E+00;
        higgswidth[305.0] = 8.99E+00;
        higgswidth[310.0] = 9.57E+00;
        higgswidth[315.0] = 1.02E+01;
        higgswidth[320.0] = 1.08E+01;
        higgswidth[325.0] = 1.14E+01;
        higgswidth[330.0] = 1.21E+01;
        higgswidth[335.0] = 1.28E+01;
        higgswidth[340.0] = 1.35E+01;
        higgswidth[345.0] = 1.42E+01;
        higgswidth[350.0] = 1.52E+01;
        higgswidth[360.0] = 1.76E+01;
        higgswidth[370.0] = 2.02E+01;
        higgswidth[380.0] = 2.31E+01;
        higgswidth[390.0] = 2.61E+01;
        higgswidth[400.0] = 2.92E+01;
        higgswidth[410.0] = 3.25E+01;
        higgswidth[420.0] = 3.59E+01;
        higgswidth[430.0] = 3.94E+01;
        higgswidth[440.0] = 4.30E+01;
        higgswidth[450.0] = 4.68E+01;
        higgswidth[460.0] = 5.08E+01;
        higgswidth[470.0] = 5.49E+01;
        higgswidth[480.0] = 5.91E+01;
        higgswidth[490.0] = 6.35E+01;
        higgswidth[500.0] = 6.80E+01;
        higgswidth[510.0] = 7.27E+01;
        higgswidth[520.0] = 7.75E+01;
        higgswidth[530.0] = 8.25E+01;
        higgswidth[540.0] = 8.77E+01;
        higgswidth[550.0] = 9.30E+01;
        higgswidth[560.0] = 9.86E+01;
        higgswidth[570.0] = 1.04E+02;
        higgswidth[580.0] = 1.10E+02;
        higgswidth[590.0] = 1.16E+02;
        higgswidth[600.0] = 1.23E+02;
        higgswidth[610.0] = 1.29E+02;
        higgswidth[620.0] = 1.36E+02;
        higgswidth[630.0] = 1.43E+02;
        higgswidth[640.0] = 1.50E+02;
        higgswidth[650.0] = 1.58E+02;
        higgswidth[660.0] = 1.65E+02;
        higgswidth[670.0] = 1.73E+02;
        higgswidth[680.0] = 1.82E+02;
        higgswidth[690.0] = 1.90E+02;
        higgswidth[700.0] = 1.99E+02;
        higgswidth[710.0] = 2.08E+02;
        higgswidth[720.0] = 2.17E+02;
        higgswidth[730.0] = 2.27E+02;
        higgswidth[740.0] = 2.37E+02;
        higgswidth[750.0] = 2.47E+02;
        higgswidth[760.0] = 2.58E+02;
        higgswidth[770.0] = 2.69E+02;
        higgswidth[780.0] = 2.80E+02;
        higgswidth[790.0] = 2.92E+02;
        higgswidth[800.0] = 3.04E+02;
        higgswidth[810.0] = 3.17E+02;
        higgswidth[820.0] = 3.30E+02;
        higgswidth[830.0] = 3.43E+02;
        higgswidth[840.0] = 3.57E+02;
        higgswidth[850.0] = 3.71E+02;
        higgswidth[860.0] = 3.86E+02;
        higgswidth[870.0] = 4.01E+02;
        higgswidth[880.0] = 4.16E+02;
        higgswidth[890.0] = 4.32E+02;
        higgswidth[900.0] = 4.49E+02;
        higgswidth[910.0] = 4.66E+02;
        higgswidth[920.0] = 4.84E+02;
        higgswidth[930.0] = 5.02E+02;
        higgswidth[940.0] = 5.21E+02;
        higgswidth[950.0] = 5.40E+02;
        higgswidth[960.0] = 5.60E+02;
        higgswidth[970.0] = 5.81E+02;
        higgswidth[980.0] = 6.02E+02;
        higgswidth[990.0] = 6.24E+02;
        higgswidth[1000.0] = 6.47E+02;


        double windowVal = max(higgswidth[cardmass], float(1.));
        double lowside = (cardmass >= 275) ? 180. : 100.;
        if (low) return std::max((cardmass - 20.*windowVal), lowside);
        else return std::min((cardmass + 15.*windowVal), 800.);
    }


    void makeCards() {
        for (float i = 114.; i <= 160.; i += 1.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
        }
        for (float i = 162.; i <= 290.; i += 2.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
        }
        for (float i = 295.; i <= 350.; i += 4.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
        }
        for (float i = 360.; i <= 600.; i += 5.) {
            createCard(i, getMassCut(i, true), getMassCut(i, false), 0);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 1);
            createCard(i, getMassCut(i, true), getMassCut(i, false), 2);
        }
    }
};

void doHZZAnalysis() {

    HiggsMassPointInfo hmpi7;
    hmpi7.lumi = 5.05;
    hmpi7.base_lumi = 5.05;
    hmpi7.z1min = 40.;
    hmpi7.z2min = 12.;
    hmpi7.massLowBkgFit = 100.;
    hmpi7.massHighBkgFit = 800.;
    hmpi7.melacut = -1.0;
    hmpi7.do1D = true;
    hmpi7.do7TeV = true;
    hmpi7.doFFT = false;
    hmpi7.treeFolder = "/home/avartak/CMS/Higgs/HZZ4L/CMSSW_4_2_8_patch7/src/WWAnalysis/AnalysisStep/trees/";
    hmpi7.melafilename = "mela2DShapes.root";

    init(hmpi7.do7TeV);

    hmpi7.ymaker_data.fill(hmpi7.treeFolder+"hzzTree.root");

    hmpi7.makeCards();
    hmpi7.do1D = false;
    hmpi7.makeCards();

    HiggsMassPointInfo hmpi8;
    hmpi8.lumi = 5.26; 
    hmpi8.base_lumi = 5.26; 
    hmpi8.z1min = 40.;
    hmpi8.z2min = 12.;
    hmpi8.massLowBkgFit = 100.;
    hmpi8.massHighBkgFit = 800.;
    hmpi8.melacut = -1.0;
    hmpi8.do1D = true;
    hmpi8.do7TeV = false;
    hmpi8.doFFT = false;
    hmpi8.treeFolder = "/home/avartak/CMS/Higgs/HZZ4L/CMSSW_5_2_4_patch4/src/WWAnalysis/AnalysisStep/trees/";
    hmpi8.melafilename = "mela2DShapes.root";

    init(hmpi8.do7TeV);


    hmpi8.ymaker_data.fill(hmpi8.treeFolder+"hzzTree.root");

    hmpi8.makeCards();
    hmpi8.do1D = false;
    hmpi8.makeCards();

}


