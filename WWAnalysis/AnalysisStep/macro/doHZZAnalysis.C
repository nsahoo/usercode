#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TBranch.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <utility>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "WWAnalysis/AnalysisStep/interface/scales.h"
#include "WWAnalysis/AnalysisStep/interface/findAndReplace.h"
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
#include <RooWorkspace.h> 
#include <RooLandau.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h> 
#include <RooFFTConvPdf.h>

#include "WWAnalysis/AnalysisStep/macro/FakeRateCalculator.h"
#include "WWAnalysis/AnalysisStep/macro/CardTemplate.h"
#include "WWAnalysis/AnalysisStep/macro/FitMaker.h"
#include "WWAnalysis/AnalysisStep/macro/YieldMaker.h"

using namespace RooFit;

struct HiggsMassPointInfo {

    static float melacut;
    static float z1min;
    static float z2min;
    static bool doSS;
    static bool doShapeAnalysis;
    static bool do1D;

    int mass;
    int massLow;
    int massHigh;
    int gghid;
    int vbfid;

    HiggsMassPointInfo(int m, int mLow, int mHigh, int gid, int vid):
        mass(m),
        massLow(mLow),
        massHigh(mHigh),
        gghid(gid),
        vbfid(vid)
    {}
};

void analysisEngine(HiggsMassPointInfo hinfo) {

    stringstream mass_str_ss;
    stringstream gghid_ss;
    stringstream vbfid_ss;

    mass_str_ss << hinfo.mass;
    gghid_ss << hinfo.gghid;
    vbfid_ss << hinfo.vbfid;

    std::string mass_str = mass_str_ss.str();
    std::string gghid = gghid_ss.str();
    std::string vbfid = vbfid_ss.str();

    std::cout << "Analyzing " << mass_str << " GeV mass point ... " << std::endl;

    std::string base_folder = "/home/avartak/CMS/Higgs/CMSSW_5_2_3/src/WWAnalysis/AnalysisStep/trees/fsr/";
    std::string card_name   = std::string("hzz_m")+mass_str;
    std::string workspace_2e2mu = card_name+"_2e2mu_workspace.root";
    std::string workspace_4e    = card_name+"_4e_workspace.root";
    std::string workspace_4mu   = card_name+"_4mu_workspace.root";

    init();

    ofstream result_out;
    result_out.open ((std::string("m")+mass_str +"_yields.txt").c_str());

    FakeRateCalculator FR(base_folder+"hzzTree.root", true, 40, 120, 0.0, 0.0, true);
    FR.printMuonFakeRate(result_out);
    FR.printElectronFakeRate(result_out);

    std::string ggh_rootfile = "hzzTree_id";
    ggh_rootfile += gghid;
    ggh_rootfile += ".root";

    std::string vbf_rootfile = "hzzTree_id";
    vbf_rootfile += vbfid;
    vbf_rootfile += ".root";

    ZXYieldMaker   ymaker_zxss;
    ZZYieldMaker   ymaker_qqzz;
    ZZYieldMaker   ymaker_ggzz;
    ZZYieldMaker   ymaker_ghzz;
    ZZYieldMaker   ymaker_qhzz;

    TH2F h2D_zxss_em("h2D_zxss_em", "", 40, 100, 180, 30, 0, 1);
    TH2F h2D_zxss_mm("h2D_zxss_mm", "", 40, 100, 180, 30, 0, 1);
    TH2F h2D_zxss_ee("h2D_zxss_ee", "", 40, 100, 180, 30, 0, 1);

    TH2F h2D_qqzz_em("h2D_qqzz_em", "", 40, 100, 180, 30, 0, 1);
    TH2F h2D_qqzz_mm("h2D_qqzz_mm", "", 40, 100, 180, 30, 0, 1);
    TH2F h2D_qqzz_ee("h2D_qqzz_ee", "", 40, 100, 180, 30, 0, 1);

    TH2F h2D_ggzz_em("h2D_ggzz_em", "", 40, 100, 180, 30, 0, 1);
    TH2F h2D_ggzz_mm("h2D_ggzz_mm", "", 40, 100, 180, 30, 0, 1);
    TH2F h2D_ggzz_ee("h2D_ggzz_ee", "", 40, 100, 180, 30, 0, 1);

    TH2F h2D_ghzz_em("h2D_ghzz_em", "", 40, 100, 180, 30, 0, 1);
    TH2F h2D_ghzz_mm("h2D_ghzz_mm", "", 40, 100, 180, 30, 0, 1);
    TH2F h2D_ghzz_ee("h2D_ghzz_ee", "", 40, 100, 180, 30, 0, 1);

    TH2F h2D_qhzz_em("h2D_qhzz_em", "", 40, 100, 180, 30, 0, 1);
    TH2F h2D_qhzz_mm("h2D_qhzz_mm", "", 40, 100, 180, 30, 0, 1);
    TH2F h2D_qhzz_ee("h2D_qhzz_ee", "", 40, 100, 180, 30, 0, 1);

    ymaker_zxss.fill(base_folder+"hzzTree.root"       , 1.0, FR, hinfo.doSS);
    ymaker_qqzz.fill(base_folder+"hzzTree_id121.root" , xsecweights[121]        *5.05, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id122.root" , xsecweights[122]        *5.05, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id123.root" , xsecweights[123]        *5.05, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id124.root" , xsecweights[124]        *5.05, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id125.root" , xsecweights[125]        *5.05, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id126.root" , xsecweights[126]        *5.05, 0.0, false);
    ymaker_ggzz.fill(base_folder+"hzzTree_id101.root" , xsecweights[101]        *5.05, 0.0, false);
    ymaker_ggzz.fill(base_folder+"hzzTree_id100.root" , xsecweights[100]        *5.05, 0.0, false);
    ymaker_ghzz.fill(base_folder+ggh_rootfile         , xsecweights[hinfo.gghid]*5.05, 0.0, true );
    ymaker_qhzz.fill(base_folder+vbf_rootfile         , xsecweights[hinfo.vbfid]*5.05, 0.0, true );

    float yield_ggh_mm  = ymaker_ghzz.getYield(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_ggh_ee  = ymaker_ghzz.getYield(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_ggh_em  = ymaker_ghzz.getYield(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_ggh_em += ymaker_ghzz.getYield(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_vbh_mm  = ymaker_qhzz.getYield(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_vbh_ee  = ymaker_qhzz.getYield(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_vbh_em  = ymaker_qhzz.getYield(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_vbh_em += ymaker_qhzz.getYield(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_qq_mm   = ymaker_qqzz.getYield(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_qq_ee   = ymaker_qqzz.getYield(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_qq_em   = ymaker_qqzz.getYield(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_qq_em  += ymaker_qqzz.getYield(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_gg_mm   = ymaker_ggzz.getYield(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_gg_ee   = ymaker_ggzz.getYield(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_gg_em   = ymaker_ggzz.getYield(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_gg_em  += ymaker_ggzz.getYield(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_zj_mm   = ymaker_zxss.getYield(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_zj_ee   = ymaker_zxss.getYield(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_zj_em   = ymaker_zxss.getYield(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_zj_em  += ymaker_zxss.getYield(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_zj_mm_count  = ymaker_zxss.getCount(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_zj_ee_count  = ymaker_zxss.getCount(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_zj_em_count  = ymaker_zxss.getCount(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_zj_em_count += ymaker_zxss.getCount(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_zj_mm_error  = ymaker_zxss.getYieldError(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_zj_ee_error  = ymaker_zxss.getYieldError(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_zj_em_error  = ymaker_zxss.getYieldError(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_zj_em_error += ymaker_zxss.getYieldError(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);

    std::string card_4mu   = createCardTemplate(0, hinfo.doShapeAnalysis, workspace_4mu.c_str());
    std::string card_4e    = createCardTemplate(1, hinfo.doShapeAnalysis, workspace_4e.c_str());
    std::string card_2e2mu = createCardTemplate(2, hinfo.doShapeAnalysis, workspace_2e2mu.c_str());

    card_2e2mu = findAndReplace(card_2e2mu, "SIG_GGH",  yield_ggh_em);
    card_4e    = findAndReplace(card_4e,    "SIG_GGH",  yield_ggh_ee);
    card_4mu   = findAndReplace(card_4mu  , "SIG_GGH",  yield_ggh_mm);
     
    card_2e2mu = findAndReplace(card_2e2mu, "SIG_VBF",  yield_vbh_em);
    card_4e    = findAndReplace(card_4e,    "SIG_VBF",  yield_vbh_ee);
    card_4mu   = findAndReplace(card_4mu,   "SIG_VBF",  yield_vbh_mm);
     
    card_2e2mu = findAndReplace(card_2e2mu, "BKG_QQZZ", yield_qq_em);
    card_4e    = findAndReplace(card_4e,    "BKG_QQZZ", yield_qq_ee);
    card_4mu   = findAndReplace(card_4mu,   "BKG_QQZZ", yield_qq_mm);
     
    card_2e2mu = findAndReplace(card_2e2mu, "BKG_GGZZ", yield_gg_em);
    card_4e    = findAndReplace(card_4e,    "BKG_GGZZ", yield_gg_ee);
    card_4mu   = findAndReplace(card_4mu,   "BKG_GGZZ", yield_gg_mm);
     
    card_2e2mu = findAndReplace(card_2e2mu, "BKG_ZJETS",yield_zj_em);
    card_4e    = findAndReplace(card_4e,    "BKG_ZJETS",yield_zj_ee);
    card_4mu   = findAndReplace(card_4mu,   "BKG_ZJETS",yield_zj_mm);
     
    card_2e2mu = findAndReplace(card_2e2mu, "ZJEVT",    yield_zj_em_count);
    card_4e    = findAndReplace(card_4e,    "ZJEVT",    yield_zj_ee_count);
    card_4mu   = findAndReplace(card_4mu,   "ZJEVT",    yield_zj_mm_count);
     
    card_2e2mu = findAndReplace(card_2e2mu, "ZJALPHA",  yield_zj_em / yield_zj_em_count);
    card_4e    = findAndReplace(card_4e,    "ZJALPHA",  yield_zj_ee / yield_zj_ee_count);
    card_4mu   = findAndReplace(card_4mu,   "ZJALPHA",  yield_zj_mm / yield_zj_mm_count);
     
    card_2e2mu = findAndReplace(card_2e2mu, "ZJSYST",   1.0+(yield_zj_em_error / yield_zj_em));
    card_4e    = findAndReplace(card_4e,    "ZJSYST",   1.0+(yield_zj_ee_error / yield_zj_ee));
    card_4mu   = findAndReplace(card_4mu,   "ZJSYST",   1.0+(yield_zj_mm_error / yield_zj_mm));

    card_2e2mu = findAndReplace(card_2e2mu, "BIN",   "a3");
    card_4e    = findAndReplace(card_4e,    "BIN",   "a2");
    card_4mu   = findAndReplace(card_4mu,   "BIN",   "a1");

    if (hinfo.doShapeAnalysis) {
        card_2e2mu = findAndReplace(card_2e2mu, "OBS",      2.0);
        card_4e    = findAndReplace(card_4e,    "OBS",      2.0);
        card_4mu   = findAndReplace(card_4mu,   "OBS",      2.0);
    }     
    else {
        card_2e2mu = findAndReplace(card_2e2mu, "OBS",      yield_zj_em + yield_gg_em + yield_qq_em);
        card_4e    = findAndReplace(card_4e,    "OBS",      yield_zj_ee + yield_gg_ee + yield_qq_ee);
        card_4mu   = findAndReplace(card_4mu,   "OBS",      yield_zj_mm + yield_gg_mm + yield_qq_mm);
    }
  
    result_out << "---- ZX Yields ----" << std::endl;
    result_out << "4e     : " << yield_zj_ee << " +/- " << yield_zj_ee_error << " (" << yield_zj_ee_count << " events in CR)" << std::endl;
    result_out << "4mu    : " << yield_zj_mm << " +/- " << yield_zj_mm_error << " (" << yield_zj_mm_count << " events in CR)" << std::endl;
    result_out << "2e2mu  : " << yield_zj_em << " +/- " << yield_zj_em_error << " (" << yield_zj_em_count << " events in CR)" << std::endl;

    result_out << std::endl;

    result_out << "---- ZZ Yields ----" << std::endl;

    result_out << "qq->ZZ->4e     : " << yield_qq_ee << std::endl;
    result_out << "qq->ZZ->4mu    : " << yield_qq_mm << std::endl;
    result_out << "qq->ZZ->2e2mu  : " << yield_qq_em << std::endl;

    result_out << std::endl;

    result_out << "gg->ZZ->4e     : " << yield_gg_ee << std::endl;
    result_out << "gg->ZZ->4mu    : " << yield_gg_mm << std::endl;
    result_out << "gg->ZZ->2e2mu  : " << yield_gg_em << std::endl;

    result_out << std::endl;

    result_out << "ZZ->4e         : " << yield_qq_ee+yield_gg_ee<< std::endl;
    result_out << "ZZ->4mu        : " << yield_qq_mm+yield_gg_mm<< std::endl;
    result_out << "ZZ->2e2mu      : " << yield_qq_em+yield_gg_em<< std::endl;

    result_out << std::endl;
    result_out << std::endl;

    result_out << "---- Higgs Yields ----" << std::endl;

    result_out << "ggH->ZZ->4e    : " << yield_ggh_ee << std::endl;
    result_out << "ggH->ZZ->4mu   : " << yield_ggh_mm << std::endl;
    result_out << "ggH->ZZ->2e2mu : " << yield_ggh_em << std::endl;

    result_out << std::endl;

    result_out << "qqH->ZZ->4e    : " << yield_vbh_ee << std::endl;
    result_out << "qqH->ZZ->4mu   : " << yield_vbh_mm << std::endl;
    result_out << "qqH->ZZ->2e2mu : " << yield_vbh_em << std::endl;

    result_out << std::endl;

    result_out << "H->ZZ->4e      : " << yield_ggh_ee + yield_vbh_ee << std::endl;
    result_out << "H->ZZ->4mu     : " << yield_ggh_mm + yield_vbh_mm << std::endl;
    result_out << "H->ZZ->2e2mu   : " << yield_ggh_em + yield_vbh_em << std::endl;

    result_out << std::endl;

    result_out << std::endl;
    result_out.close();

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


 
    QQZZFitMaker fitmaker_2e2mu    ("bkg_qqzz",  100, 95., 605.);
    QQZZFitMaker fitmaker_4e       ("bkg_qqzz",  100, 95., 605.);
    QQZZFitMaker fitmaker_4mu      ("bkg_qqzz",  100, 95., 605.);
    
    GGZZFitMaker fitmaker_gz2e2mu  ("bkg_ggzz",  100, 95., 605.);
    GGZZFitMaker fitmaker_gz4e     ("bkg_ggzz",  100, 95., 605.);
    GGZZFitMaker fitmaker_gz4mu    ("bkg_ggzz",  100, 95., 605.);
    
    ZXFitMaker   fitmaker_zx2e2mu  ("bkg_zjets", 100, 95., 605.);
    ZXFitMaker   fitmaker_zx4e     ("bkg_zjets", 100, 95., 605.);
    ZXFitMaker   fitmaker_zx4mu    ("bkg_zjets", 100, 95., 605.);
    
    SignalFitMaker  fitmaker_g2e2mu("sig_ggH", hinfo.mass, hinfo.mass, hinfo.massLow, hinfo.massHigh);
    SignalFitMaker  fitmaker_g4e   ("sig_ggH", hinfo.mass, hinfo.mass, hinfo.massLow, hinfo.massHigh);
    SignalFitMaker  fitmaker_g4mu  ("sig_ggH", hinfo.mass, hinfo.mass, hinfo.massLow, hinfo.massHigh);
    
    SignalFitMaker  fitmaker_q2e2mu("sig_VBF", hinfo.mass, hinfo.mass, hinfo.massLow, hinfo.massHigh);
    SignalFitMaker  fitmaker_q4e   ("sig_VBF", hinfo.mass, hinfo.mass, hinfo.massLow, hinfo.massHigh);
    SignalFitMaker  fitmaker_q4mu  ("sig_VBF", hinfo.mass, hinfo.mass, hinfo.massLow, hinfo.massHigh);
        
    if (hinfo.doShapeAnalysis) {  
        fitmaker_gz4mu.add  (ymaker_ggzz.getFitDataSet(0, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_gz4e.add   (ymaker_ggzz.getFitDataSet(1, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_gz2e2mu.add(ymaker_ggzz.getFitDataSet(2, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_gz2e2mu.add(ymaker_ggzz.getFitDataSet(3, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        
        fitmaker_4mu.add    (ymaker_qqzz.getFitDataSet(0, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_4e.add     (ymaker_qqzz.getFitDataSet(1, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_2e2mu.add  (ymaker_qqzz.getFitDataSet(2, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_2e2mu.add  (ymaker_qqzz.getFitDataSet(3, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        
        fitmaker_zx4mu.add  (ymaker_zxss.getFitDataSet(0, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_zx4e.add   (ymaker_zxss.getFitDataSet(1, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_zx2e2mu.add(ymaker_zxss.getFitDataSet(2, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_zx2e2mu.add(ymaker_zxss.getFitDataSet(3, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        
        fitmaker_g4mu.add   (ymaker_ghzz.getFitDataSet(0, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_g4e.add    (ymaker_ghzz.getFitDataSet(1, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_g2e2mu.add (ymaker_ghzz.getFitDataSet(2, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_g2e2mu.add (ymaker_ghzz.getFitDataSet(3, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        
        fitmaker_q4mu.add   (ymaker_qhzz.getFitDataSet(0, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_q4e.add    (ymaker_qhzz.getFitDataSet(1, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_q2e2mu.add (ymaker_qhzz.getFitDataSet(2, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));
        fitmaker_q2e2mu.add (ymaker_qhzz.getFitDataSet(3, hinfo.z1min, hinfo.z2min, 95., 605., hinfo.melacut));

        fitmaker_g2e2mu.fit();
        fitmaker_g4e.fit();
        fitmaker_g4mu.fit();
        
        fitmaker_g2e2mu.print("ggH_2e2mu_"+mass_str+".pdf");
        fitmaker_g4e.print   ("ggH_4e_"   +mass_str+".pdf");
        fitmaker_g4mu.print  ("ggH_4mu_"  +mass_str+".pdf");
        
        fitmaker_q2e2mu.fit();
        fitmaker_q4e.fit();
        fitmaker_q4mu.fit();
        
        fitmaker_q2e2mu.print("VBF_2e2mu_"+mass_str+".pdf");
        fitmaker_q4e.print   ("VBF_4e_"   +mass_str+".pdf");
        fitmaker_q4mu.print  ("VBF_4mu_"  +mass_str+".pdf");
        
        fitmaker_zx2e2mu.fit();
        fitmaker_zx4e.fit();
        fitmaker_zx4mu.fit();
        
        fitmaker_zx2e2mu.print("zjets_2e2mu_"+mass_str+".pdf");
        fitmaker_zx4e.print   ("zjets_4e_"   +mass_str+".pdf");
        fitmaker_zx4mu.print  ("zjets_4mu_"  +mass_str+".pdf");
        
        fitmaker_2e2mu.fit();
        fitmaker_4e.fit();
        fitmaker_4mu.fit();
        
        fitmaker_2e2mu.print("qqZZ_2e2mu_"+mass_str+".pdf");
        fitmaker_4e.print   ("qqZZ_4e_"   +mass_str+".pdf");
        fitmaker_4mu.print  ("qqZZ_4mu_"  +mass_str+".pdf");
        
        fitmaker_gz2e2mu.fit();
        fitmaker_gz4e.fit();
        fitmaker_gz4mu.fit();
        
        fitmaker_gz2e2mu.print("ggZZ_2e2mu_"+mass_str+".pdf");
        fitmaker_gz4e.print   ("ggZZ_4e_"   +mass_str+".pdf");
        fitmaker_gz4mu.print  ("ggZZ_4mu_"  +mass_str+".pdf");
        
        RooWorkspace w_2e2mu("w_2e2mu", "");
        RooWorkspace w_4e   ("w_4e"   , "");
        RooWorkspace w_4mu  ("w_4mu"  , "");

        if (hinfo.do1D) { 
            RooRealVar CMS_zz4l_mass_1D("CMS_zz4l_mass_1D", "M(4l)", hinfo.massLow, hinfo.massHigh, "GeV/c^{2}");
            fitmaker_g2e2mu.makeWorkspace1D(w_2e2mu , CMS_zz4l_mass_1D);
            fitmaker_g4e.makeWorkspace1D   (w_4e    , CMS_zz4l_mass_1D);
            fitmaker_g4mu.makeWorkspace1D  (w_4mu   , CMS_zz4l_mass_1D);
            
            fitmaker_q2e2mu.makeWorkspace1D(w_2e2mu , CMS_zz4l_mass_1D);
            fitmaker_q4e.makeWorkspace1D   (w_4e    , CMS_zz4l_mass_1D);
            fitmaker_q4mu.makeWorkspace1D  (w_4mu   , CMS_zz4l_mass_1D);

            fitmaker_zx2e2mu.makeWorkspace1D(w_2e2mu , CMS_zz4l_mass_1D);
            fitmaker_zx4e.makeWorkspace1D   (w_4e    , CMS_zz4l_mass_1D);
            fitmaker_zx4mu.makeWorkspace1D  (w_4mu   , CMS_zz4l_mass_1D);
            
            fitmaker_2e2mu.makeWorkspace1D  (w_2e2mu , CMS_zz4l_mass_1D);
            fitmaker_4e.makeWorkspace1D     (w_4e    , CMS_zz4l_mass_1D);
            fitmaker_4mu.makeWorkspace1D    (w_4mu   , CMS_zz4l_mass_1D);
            
            fitmaker_gz2e2mu.makeWorkspace1D(w_2e2mu , CMS_zz4l_mass_1D);
            fitmaker_gz4e.makeWorkspace1D   (w_4e    , CMS_zz4l_mass_1D);
            fitmaker_gz4mu.makeWorkspace1D  (w_4mu   , CMS_zz4l_mass_1D);

            RooRealVar weight("weight", "weight", 0.,  0.,  10.);
            RooArgSet argset_obs(CMS_zz4l_mass_1D, weight, "argset_obs");
            RooDataSet data_obs("data_obs", "data_obs", argset_obs, RooFit::WeightVar("weight"));
            
            argset_obs.setRealValue("CMS_zz4l_mass_1D", hinfo.massLow+1.0);
            argset_obs.setRealValue("weight", 1.0);
            data_obs.add(argset_obs, 1.0);
            
            argset_obs.setRealValue("mass", hinfo.massHigh-1.0);
            argset_obs.setRealValue("weight", 1.0);
            data_obs.add(argset_obs, 1.0);
            
            w_2e2mu.import(data_obs);
            w_4e.import   (data_obs);
            w_4mu.import  (data_obs);

            TH1* shapehist_zx2e2mu = fitmaker_zx2e2mu.getShapeHistogram("shapehist_zx2e2mu", 31, 99, 161);
            TH1* shapehist_zx4e    = fitmaker_zx4e.getShapeHistogram   ("shapehist_zx4e"   , 31, 99, 161);
            TH1* shapehist_zx4mu   = fitmaker_zx4mu.getShapeHistogram  ("shapehist_zx4mu"  , 31, 99, 161);
            TH1* shapehist_gz2e2mu = fitmaker_gz2e2mu.getShapeHistogram("shapehist_gz2e2mu", 31, 99, 161);
            TH1* shapehist_gz4e    = fitmaker_gz4e.getShapeHistogram   ("shapehist_gz4e"   , 31, 99, 161);
            TH1* shapehist_gz4mu   = fitmaker_gz4mu.getShapeHistogram  ("shapehist_gz4mu"  , 31, 99, 161);
            TH1* shapehist_qz2e2mu = fitmaker_2e2mu.getShapeHistogram  ("shapehist_qz2e2mu", 31, 99, 161);
            TH1* shapehist_qz4e    = fitmaker_4e.getShapeHistogram     ("shapehist_qz4e"   , 31, 99, 161);
            TH1* shapehist_qz4mu   = fitmaker_4mu.getShapeHistogram    ("shapehist_qz4mu"  , 31, 99, 161);
            TH1* shapehist_gh2e2mu = fitmaker_g2e2mu.getShapeHistogram ("shapehist_gh2e2mu", 31, 99, 161);
            TH1* shapehist_gh4e    = fitmaker_g4e.getShapeHistogram    ("shapehist_gh4e"   , 31, 99, 161);
            TH1* shapehist_gh4mu   = fitmaker_g4mu.getShapeHistogram   ("shapehist_gh4mu"  , 31, 99, 161);
            TH1* shapehist_qh2e2mu = fitmaker_q2e2mu.getShapeHistogram ("shapehist_qh2e2mu", 31, 99, 161);
            TH1* shapehist_qh4e    = fitmaker_q4e.getShapeHistogram    ("shapehist_qh4e"   , 31, 99, 161);
            TH1* shapehist_qh4mu   = fitmaker_q4mu.getShapeHistogram   ("shapehist_qh4mu"  , 31, 99, 161);

            TFile* shapefile = new TFile(("shapehist_"+mass_str+".root").c_str(), "RECREATE");

            shapehist_zx2e2mu ->Scale(ymaker_zxss.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_zx4e    ->Scale(ymaker_zxss.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_zx4mu   ->Scale(ymaker_zxss.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_gz2e2mu ->Scale(ymaker_ggzz.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_gz4e    ->Scale(ymaker_ggzz.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_gz4mu   ->Scale(ymaker_ggzz.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_qz2e2mu ->Scale(ymaker_qqzz.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_qz4e    ->Scale(ymaker_qqzz.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_qz4mu   ->Scale(ymaker_qqzz.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_gh2e2mu ->Scale(ymaker_ghzz.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_gh4e    ->Scale(ymaker_ghzz.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_gh4mu   ->Scale(ymaker_ghzz.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_qh2e2mu ->Scale(ymaker_qhzz.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_qh4e    ->Scale(ymaker_qhzz.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));
            shapehist_qh4mu   ->Scale(ymaker_qhzz.getYield(0, hinfo.z1min, hinfo.z2min, 99, 161., hinfo.melacut));


            shapehist_gz2e2mu ->Write();
            shapehist_gz4e    ->Write();
            shapehist_gz4mu   ->Write();
            shapehist_qz2e2mu ->Write();
            shapehist_qz4e    ->Write();
            shapehist_qz4mu   ->Write();
            shapehist_zx2e2mu ->Write();
            shapehist_zx4e    ->Write();
            shapehist_zx4mu   ->Write();
            shapehist_gh2e2mu ->Write();
            shapehist_gh4e    ->Write();
            shapehist_gh4mu   ->Write();
            shapehist_qh2e2mu ->Write();
            shapehist_qh4e    ->Write();
            shapehist_qh4mu   ->Write();
                            
            shapefile->Close();
            delete shapefile;
        }      
        else {
            RooRealVar CMS_zz4l_mass_2D("CMS_zz4l_mass_2D", "M(4l)", 100,           180,            "GeV/c^{2}");
            RooRealVar CMS_zz4l_melaLD ("CMS_zz4l_melaLD" , "MELA" ,   0,             1,            "");
            CMS_zz4l_melaLD.setBins(30);            

            ymaker_ghzz.get2DHist(0, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_ghzz_mm);
            ymaker_ghzz.get2DHist(1, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_ghzz_ee);
            ymaker_ghzz.get2DHist(2, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_ghzz_em);
            ymaker_ghzz.get2DHist(3, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_ghzz_em);
            
            ymaker_qhzz.get2DHist(0, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_qhzz_mm);
            ymaker_qhzz.get2DHist(1, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_qhzz_ee);
            ymaker_qhzz.get2DHist(2, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_qhzz_em);
            ymaker_qhzz.get2DHist(3, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_qhzz_em);
            
            ymaker_qqzz.get2DHist(0, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_qqzz_mm);
            ymaker_qqzz.get2DHist(1, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_qqzz_ee);
            ymaker_qqzz.get2DHist(2, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_qqzz_em);
            ymaker_qqzz.get2DHist(3, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_qqzz_em);
            
            ymaker_ggzz.get2DHist(0, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_ggzz_mm);
            ymaker_ggzz.get2DHist(1, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_ggzz_ee);
            ymaker_ggzz.get2DHist(2, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_ggzz_em);
            ymaker_ggzz.get2DHist(3, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_ggzz_em);
            
            ymaker_zxss.get2DHist(0, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_zxss_mm);
            ymaker_zxss.get2DHist(1, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_zxss_ee);
            ymaker_zxss.get2DHist(2, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_zxss_em);
            ymaker_zxss.get2DHist(3, hinfo.z1min, hinfo.z2min, 100, 180, 0.0, &h2D_zxss_em);

            fitmaker_g2e2mu.makeWorkspace2D(w_2e2mu , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_ghzz_em);
            fitmaker_g4e.makeWorkspace2D   (w_4e    , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_ghzz_ee);
            fitmaker_g4mu.makeWorkspace2D  (w_4mu   , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_ghzz_mm);
            
            fitmaker_q2e2mu.makeWorkspace2D(w_2e2mu , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_qhzz_em);
            fitmaker_q4e.makeWorkspace2D   (w_4e    , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_qhzz_ee);
            fitmaker_q4mu.makeWorkspace2D  (w_4mu   , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_qhzz_mm);

            fitmaker_zx2e2mu.makeWorkspace2D(w_2e2mu , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_zxss_em);
            fitmaker_zx4e.makeWorkspace2D   (w_4e    , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_zxss_em);
            fitmaker_zx4mu.makeWorkspace2D  (w_4mu   , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_zxss_em);
            
            fitmaker_2e2mu.makeWorkspace2D  (w_2e2mu , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_qqzz_em);
            fitmaker_4e.makeWorkspace2D     (w_4e    , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_qqzz_em);
            fitmaker_4mu.makeWorkspace2D    (w_4mu   , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_qqzz_em);
            
            fitmaker_gz2e2mu.makeWorkspace2D(w_2e2mu , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_ggzz_em);
            fitmaker_gz4e.makeWorkspace2D   (w_4e    , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_ggzz_em);
            fitmaker_gz4mu.makeWorkspace2D  (w_4mu   , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_ggzz_em);
        
            RooRealVar weight("weight", "weight", 0.,  0.,  10.);
            RooArgSet argset_obs(CMS_zz4l_mass_2D, CMS_zz4l_melaLD, weight, "argset_obs");
            RooDataSet data_obs("data_obs", "data_obs", argset_obs, RooFit::WeightVar("weight"));
        
            argset_obs.setRealValue("CMS_zz4l_mass_1D", 120);
            argset_obs.setRealValue("CMS_zz4l_mass_2D", 0.1);
            argset_obs.setRealValue("weight", 1.0);
            data_obs.add(argset_obs, 1.0);
        
            argset_obs.setRealValue("CMS_zz4l_mass_1D", 160);
            argset_obs.setRealValue("CMS_zz4l_mass_2D", 0.9);
            argset_obs.setRealValue("weight", 1.0);
            data_obs.add(argset_obs, 1.0);
        
            w_2e2mu.import(data_obs);
            w_4e.import   (data_obs);
            w_4mu.import  (data_obs);
        }
        
        w_2e2mu.writeToFile(workspace_2e2mu.c_str());
        w_4e.writeToFile(workspace_4e.c_str());
        w_4mu.writeToFile(workspace_4mu.c_str());
    }

}


float HiggsMassPointInfo::z1min = 40.;
float HiggsMassPointInfo::z2min = 12.;
float HiggsMassPointInfo::melacut = -1.0;
bool  HiggsMassPointInfo::doShapeAnalysis = true;
bool  HiggsMassPointInfo::do1D = true;
bool  HiggsMassPointInfo::doSS = true;

void doHZZAnalysis() {

    HiggsMassPointInfo h115(115, 100, 130, 200, 250);
    HiggsMassPointInfo h120(120, 100, 135, 201, 251);
    HiggsMassPointInfo h130(130, 110, 145, 202, 252);
    HiggsMassPointInfo h140(140, 120, 155, 203, 253);
    HiggsMassPointInfo h150(150, 130, 165, 204, 254);
    HiggsMassPointInfo h160(160, 140, 175, 205, 255);

    //HiggsMassPointInfo h115(115, 100, 800, 200, 250);
    //HiggsMassPointInfo h120(120, 100, 800, 201, 251);
    //HiggsMassPointInfo h130(130, 100, 800, 202, 252);
    //HiggsMassPointInfo h140(140, 100, 800, 203, 253);
    //HiggsMassPointInfo h150(150, 100, 800, 204, 254);
    //HiggsMassPointInfo h160(160, 100, 800, 205, 255);

    
    bool hinfodoshape = HiggsMassPointInfo::doShapeAnalysis;
    HiggsMassPointInfo::doShapeAnalysis = false;
    analysisEngine(h115);
    HiggsMassPointInfo::doShapeAnalysis = hinfodoshape;

    if (HiggsMassPointInfo::doShapeAnalysis) {
        //analysisEngine(h115);
        analysisEngine(h120);
        //analysisEngine(h130);
        //analysisEngine(h140);
        //analysisEngine(h150);
        //analysisEngine(h160);
    }

    else {
        //analysisEngine(h115);
        analysisEngine(h120);
        analysisEngine(h130);
        //analysisEngine(h140);
        //analysisEngine(h150);
        //analysisEngine(h160);
    }
}
