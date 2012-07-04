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

    static float lumi;
    static float melacut;
    static float z1min;
    static float z2min;
    static float massLow;
    static float massHigh;
    static float massLowBkgFit;
    static float massHighBkgFit;
    static int   nBinsMass2D;
    static int   nBinsMELA2D;
    static bool  doSS;
    static bool  doShapeAnalysis;
    static bool  do1D;
    static std::string treeFolder;

    int mass;
    int massLowSigFit;
    int massHighSigFit;
    int gghid;
    int vbfid;

    HiggsMassPointInfo(int m, int mLow, int mHigh, int gid, int vid):
        mass(m),
        massLowSigFit(mLow),
        massHighSigFit(mHigh),
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

    std::string base_folder = hinfo.treeFolder;
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

    DataYieldMaker ymaker_data;
    ZXYieldMaker   ymaker_zxss;
    ZZYieldMaker   ymaker_qqzz;
    ZZYieldMaker   ymaker_ggzz;
    ZZYieldMaker   ymaker_ghzz;
    ZZYieldMaker   ymaker_qhzz;

    TH2F h2D_bkg_em("h2D_bkg_em", "", hinfo.nBinsMass2D, hinfo.massLow, hinfo.massHigh, hinfo.nBinsMELA2D, 0, 1);
    TH2F h2D_bkg_mm("h2D_bkg_mm", "", hinfo.nBinsMass2D, hinfo.massLow, hinfo.massHigh, hinfo.nBinsMELA2D, 0, 1);
    TH2F h2D_bkg_ee("h2D_bkg_ee", "", hinfo.nBinsMass2D, hinfo.massLow, hinfo.massHigh, hinfo.nBinsMELA2D, 0, 1);

    TH2F h2D_sig_em("h2D_sig_em", "", hinfo.nBinsMass2D, hinfo.massLow, hinfo.massHigh, hinfo.nBinsMELA2D, 0, 1);
    TH2F h2D_sig_mm("h2D_sig_mm", "", hinfo.nBinsMass2D, hinfo.massLow, hinfo.massHigh, hinfo.nBinsMELA2D, 0, 1);
    TH2F h2D_sig_ee("h2D_sig_ee", "", hinfo.nBinsMass2D, hinfo.massLow, hinfo.massHigh, hinfo.nBinsMELA2D, 0, 1);

    ymaker_data.fill(base_folder+"hzzTree.root");
    ymaker_zxss.fill(base_folder+"hzzTree.root"       , 1.0, FR, hinfo.doSS);
    ymaker_qqzz.fill(base_folder+"hzzTree_id121.root" , xsecweights[121]        *hinfo.lumi, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id122.root" , xsecweights[122]        *hinfo.lumi, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id123.root" , xsecweights[123]        *hinfo.lumi, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id124.root" , xsecweights[124]        *hinfo.lumi, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id125.root" , xsecweights[125]        *hinfo.lumi, 0.0, false);
    ymaker_qqzz.fill(base_folder+"hzzTree_id126.root" , xsecweights[126]        *hinfo.lumi, 0.0, false);
    ymaker_ggzz.fill(base_folder+"hzzTree_id101.root" , xsecweights[101]        *hinfo.lumi, 0.0, false);
    ymaker_ggzz.fill(base_folder+"hzzTree_id100.root" , xsecweights[100]        *hinfo.lumi, 0.0, false);
    ymaker_ghzz.fill(base_folder+ggh_rootfile         , xsecweights[hinfo.gghid]*hinfo.lumi, 0.0, true );
    ymaker_qhzz.fill(base_folder+vbf_rootfile         , xsecweights[hinfo.vbfid]*hinfo.lumi, 0.0, true );

    float yield_dt_mm        = ymaker_data.getYield     (0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_dt_ee        = ymaker_data.getYield     (1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_dt_em        = ymaker_data.getYield     (2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_dt_em       += ymaker_data.getYield     (3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_ggh_mm       = ymaker_ghzz.getYield     (0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_ggh_ee       = ymaker_ghzz.getYield     (1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_ggh_em       = ymaker_ghzz.getYield     (2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_ggh_em      += ymaker_ghzz.getYield     (3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_vbh_mm       = ymaker_qhzz.getYield     (0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_vbh_ee       = ymaker_qhzz.getYield     (1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_vbh_em       = ymaker_qhzz.getYield     (2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_vbh_em      += ymaker_qhzz.getYield     (3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_qq_mm        = ymaker_qqzz.getYield     (0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_qq_ee        = ymaker_qqzz.getYield     (1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_qq_em        = ymaker_qqzz.getYield     (2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_qq_em       += ymaker_qqzz.getYield     (3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_gg_mm        = ymaker_ggzz.getYield     (0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_gg_ee        = ymaker_ggzz.getYield     (1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_gg_em        = ymaker_ggzz.getYield     (2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_gg_em       += ymaker_ggzz.getYield     (3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_zj_mm        = ymaker_zxss.getYield     (0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_zj_ee        = ymaker_zxss.getYield     (1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_zj_em        = ymaker_zxss.getYield     (2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_zj_em       += ymaker_zxss.getYield     (3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_zj_mm_count  = ymaker_zxss.getCount     (0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_zj_ee_count  = ymaker_zxss.getCount     (1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_zj_em_count  = ymaker_zxss.getCount     (2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_zj_em_count += ymaker_zxss.getCount     (3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
     
    float yield_zj_mm_error  = ymaker_zxss.getYieldError(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_zj_ee_error  = ymaker_zxss.getYieldError(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
    float yield_zj_em_error  = ymaker_zxss.getYieldError(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);
          yield_zj_em_error += ymaker_zxss.getYieldError(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut);

    std::string card_4mu   = createCardTemplate(0, hinfo.doShapeAnalysis, workspace_4mu.c_str());
    std::string card_4e    = createCardTemplate(1, hinfo.doShapeAnalysis, workspace_4e.c_str());
    std::string card_2e2mu = createCardTemplate(2, hinfo.doShapeAnalysis, workspace_2e2mu.c_str());

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

    card_2e2mu = findAndReplace(card_2e2mu, "GGZZ_PDF"      ,  getGGZZPDFUncertainty7TeV(hinfo.mass));
    card_4e    = findAndReplace(card_4e,    "GGZZ_PDF"      ,  getGGZZPDFUncertainty7TeV(hinfo.mass));
    card_4mu   = findAndReplace(card_4mu  , "GGZZ_PDF"      ,  getGGZZPDFUncertainty7TeV(hinfo.mass));
     
    card_2e2mu = findAndReplace(card_2e2mu, "QQZZ_PDF"      ,  getQQZZPDFUncertainty7TeV(hinfo.mass));
    card_4e    = findAndReplace(card_4e,    "QQZZ_PDF"      ,  getQQZZPDFUncertainty7TeV(hinfo.mass));
    card_4mu   = findAndReplace(card_4mu  , "QQZZ_PDF"      ,  getQQZZPDFUncertainty7TeV(hinfo.mass));
     
    card_2e2mu = findAndReplace(card_2e2mu, "GGZZ_QCD"      ,  getGGZZQCDScaleUncertainty7TeV(hinfo.mass));
    card_4e    = findAndReplace(card_4e,    "GGZZ_QCD"      ,  getGGZZQCDScaleUncertainty7TeV(hinfo.mass));
    card_4mu   = findAndReplace(card_4mu  , "GGZZ_QCD"      ,  getGGZZQCDScaleUncertainty7TeV(hinfo.mass));
     
    card_2e2mu = findAndReplace(card_2e2mu, "QQZZ_QCD"      ,  getQQZZQCDScaleUncertainty7TeV(hinfo.mass));
    card_4e    = findAndReplace(card_4e,    "QQZZ_QCD"      ,  getQQZZQCDScaleUncertainty7TeV(hinfo.mass));
    card_4mu   = findAndReplace(card_4mu  , "QQZZ_QCD"      ,  getQQZZQCDScaleUncertainty7TeV(hinfo.mass));
     
    card_2e2mu = findAndReplace(card_2e2mu, "GGH_PDF"       ,  getGGHPDFUncertaintyUp7TeV(hinfo.mass), getGGHPDFUncertaintyDown7TeV(hinfo.mass));
    card_4e    = findAndReplace(card_4e,    "GGH_PDF"       ,  getGGHPDFUncertaintyUp7TeV(hinfo.mass), getGGHPDFUncertaintyDown7TeV(hinfo.mass));
    card_4mu   = findAndReplace(card_4mu  , "GGH_PDF"       ,  getGGHPDFUncertaintyUp7TeV(hinfo.mass), getGGHPDFUncertaintyDown7TeV(hinfo.mass));
     
    card_2e2mu = findAndReplace(card_2e2mu, "VBF_PDF"       ,  getVBFPDFUncertaintyUp7TeV(hinfo.mass), getVBFPDFUncertaintyDown7TeV(hinfo.mass));
    card_4e    = findAndReplace(card_4e,    "VBF_PDF"       ,  getVBFPDFUncertaintyUp7TeV(hinfo.mass), getVBFPDFUncertaintyDown7TeV(hinfo.mass));
    card_4mu   = findAndReplace(card_4mu  , "VBF_PDF"       ,  getVBFPDFUncertaintyUp7TeV(hinfo.mass), getVBFPDFUncertaintyDown7TeV(hinfo.mass));
     
    card_2e2mu = findAndReplace(card_2e2mu, "GGH_QCD"       ,  getGGHQCDScaleUncertaintyUp7TeV(hinfo.mass), getGGHQCDScaleUncertaintyDown7TeV(hinfo.mass));
    card_4e    = findAndReplace(card_4e,    "GGH_QCD"       ,  getGGHQCDScaleUncertaintyUp7TeV(hinfo.mass), getGGHQCDScaleUncertaintyDown7TeV(hinfo.mass));
    card_4mu   = findAndReplace(card_4mu  , "GGH_QCD"       ,  getGGHQCDScaleUncertaintyUp7TeV(hinfo.mass), getGGHQCDScaleUncertaintyDown7TeV(hinfo.mass));
     
    card_2e2mu = findAndReplace(card_2e2mu, "VBF_QCD"       ,  getVBFQCDScaleUncertaintyUp7TeV(hinfo.mass), getVBFQCDScaleUncertaintyDown7TeV(hinfo.mass));
    card_4e    = findAndReplace(card_4e,    "VBF_QCD"       ,  getVBFQCDScaleUncertaintyUp7TeV(hinfo.mass), getVBFQCDScaleUncertaintyDown7TeV(hinfo.mass));
    card_4mu   = findAndReplace(card_4mu  , "VBF_QCD"       ,  getVBFQCDScaleUncertaintyUp7TeV(hinfo.mass), getVBFQCDScaleUncertaintyDown7TeV(hinfo.mass));
     
    card_2e2mu = findAndReplace(card_2e2mu, "SIG_GGH_YIELD" ,  yield_ggh_em);
    card_4e    = findAndReplace(card_4e,    "SIG_GGH_YIELD" ,  yield_ggh_ee);
    card_4mu   = findAndReplace(card_4mu  , "SIG_GGH_YIELD" ,  yield_ggh_mm);
     
    card_2e2mu = findAndReplace(card_2e2mu, "SIG_VBF_YIELD" ,  yield_vbh_em);
    card_4e    = findAndReplace(card_4e,    "SIG_VBF_YIELD" ,  yield_vbh_ee);
    card_4mu   = findAndReplace(card_4mu,   "SIG_VBF_YIELD" ,  yield_vbh_mm);
     
    card_2e2mu = findAndReplace(card_2e2mu, "BKG_QQZZ_YIELD",  yield_qq_em);
    card_4e    = findAndReplace(card_4e,    "BKG_QQZZ_YIELD",  yield_qq_ee);
    card_4mu   = findAndReplace(card_4mu,   "BKG_QQZZ_YIELD",  yield_qq_mm);
     
    card_2e2mu = findAndReplace(card_2e2mu, "BKG_GGZZ_YIELD",  yield_gg_em);
    card_4e    = findAndReplace(card_4e,    "BKG_GGZZ_YIELD",  yield_gg_ee);
    card_4mu   = findAndReplace(card_4mu,   "BKG_GGZZ_YIELD",  yield_gg_mm);
     
    card_2e2mu = findAndReplace(card_2e2mu, "BKG_ZJETS_YIELD", yield_zj_em);
    card_4e    = findAndReplace(card_4e,    "BKG_ZJETS_YIELD", yield_zj_ee);
    card_4mu   = findAndReplace(card_4mu,   "BKG_ZJETS_YIELD", yield_zj_mm);
     
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

    card_2e2mu = findAndReplace(card_2e2mu, "OBS",      yield_dt_em);
    card_4e    = findAndReplace(card_4e,    "OBS",      yield_dt_ee);
    card_4mu   = findAndReplace(card_4mu,   "OBS",      yield_dt_mm);
  
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

 
    QQZZFitMaker fitmaker_2e2mu    ("bkg_qqzz_2e2mu" , hinfo.massLowBkgFit, hinfo.massHighBkgFit);
    QQZZFitMaker fitmaker_4e       ("bkg_qqzz_4e"    , hinfo.massLowBkgFit, hinfo.massHighBkgFit);
    QQZZFitMaker fitmaker_4mu      ("bkg_qqzz_4mu"   , hinfo.massLowBkgFit, hinfo.massHighBkgFit);
    
    GGZZFitMaker fitmaker_gz2e2mu  ("bkg_ggzz_2e2mu" , hinfo.massLowBkgFit, hinfo.massHighBkgFit);
    GGZZFitMaker fitmaker_gz4e     ("bkg_ggzz_4e"    , hinfo.massLowBkgFit, hinfo.massHighBkgFit);
    GGZZFitMaker fitmaker_gz4mu    ("bkg_ggzz_4mu"   , hinfo.massLowBkgFit, hinfo.massHighBkgFit);
    
    ZXFitMaker   fitmaker_zx2e2mu  ("bkg_zjets_2e2mu", hinfo.massLowBkgFit, hinfo.massHighBkgFit);
    ZXFitMaker   fitmaker_zx4e     ("bkg_zjets_4e"   , hinfo.massLowBkgFit, hinfo.massHighBkgFit);
    ZXFitMaker   fitmaker_zx4mu    ("bkg_zjets_4mu"  , hinfo.massLowBkgFit, hinfo.massHighBkgFit);
    
    SignalFitMaker  fitmaker_g2e2mu("sig_ggH_2e2mu"  , hinfo.mass   , hinfo.massLowSigFit,  hinfo.massHighSigFit);
    SignalFitMaker  fitmaker_g4e   ("sig_ggH_4e"     , hinfo.mass   , hinfo.massLowSigFit,  hinfo.massHighSigFit);
    SignalFitMaker  fitmaker_g4mu  ("sig_ggH_4mu"    , hinfo.mass   , hinfo.massLowSigFit,  hinfo.massHighSigFit);
    
    SignalFitMaker  fitmaker_q2e2mu("sig_VBF_2e2mu"  , hinfo.mass   , hinfo.massLowSigFit,  hinfo.massHighSigFit);
    SignalFitMaker  fitmaker_q4e   ("sig_VBF_4e"     , hinfo.mass   , hinfo.massLowSigFit,  hinfo.massHighSigFit);
    SignalFitMaker  fitmaker_q4mu  ("sig_VBF_4mu"    , hinfo.mass   , hinfo.massLowSigFit,  hinfo.massHighSigFit);
        
    if (hinfo.doShapeAnalysis) {  
        fitmaker_gz4mu.add  (ymaker_ggzz.getFitDataSet(0, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_gz4e.add   (ymaker_ggzz.getFitDataSet(1, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_gz2e2mu.add(ymaker_ggzz.getFitDataSet(2, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_gz2e2mu.add(ymaker_ggzz.getFitDataSet(3, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        
        fitmaker_4mu.add    (ymaker_qqzz.getFitDataSet(0, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_4e.add     (ymaker_qqzz.getFitDataSet(1, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_2e2mu.add  (ymaker_qqzz.getFitDataSet(2, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_2e2mu.add  (ymaker_qqzz.getFitDataSet(3, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        
        fitmaker_zx4mu.add  (ymaker_zxss.getFitDataSet(0, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_zx4e.add   (ymaker_zxss.getFitDataSet(1, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_zx2e2mu.add(ymaker_zxss.getFitDataSet(2, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_zx2e2mu.add(ymaker_zxss.getFitDataSet(3, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        
        fitmaker_g4mu.add   (ymaker_ghzz.getFitDataSet(0, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_g4e.add    (ymaker_ghzz.getFitDataSet(1, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_g2e2mu.add (ymaker_ghzz.getFitDataSet(2, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_g2e2mu.add (ymaker_ghzz.getFitDataSet(3, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        
        fitmaker_g4mu.add   (ymaker_qhzz.getFitDataSet(0, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_g4e.add    (ymaker_qhzz.getFitDataSet(1, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_g2e2mu.add (ymaker_qhzz.getFitDataSet(2, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_g2e2mu.add (ymaker_qhzz.getFitDataSet(3, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        
        fitmaker_q4mu.add   (ymaker_ghzz.getFitDataSet(0, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_q4e.add    (ymaker_ghzz.getFitDataSet(1, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_q2e2mu.add (ymaker_ghzz.getFitDataSet(2, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_q2e2mu.add (ymaker_ghzz.getFitDataSet(3, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));

        fitmaker_q4mu.add   (ymaker_qhzz.getFitDataSet(0, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_q4e.add    (ymaker_qhzz.getFitDataSet(1, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_q2e2mu.add (ymaker_qhzz.getFitDataSet(2, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));
        fitmaker_q2e2mu.add (ymaker_qhzz.getFitDataSet(3, hinfo.z1min, hinfo.z2min, hinfo.massLowBkgFit, hinfo.massHighBkgFit, hinfo.melacut));

        fitmaker_g2e2mu.fit();
        fitmaker_g4e.fit();
        fitmaker_g4mu.fit();
        
        fitmaker_q2e2mu.fit();
        fitmaker_q4e.fit();
        fitmaker_q4mu.fit();
        
        fitmaker_g2e2mu.print("higgs_2e2mu_"+mass_str+".pdf");
        fitmaker_g4e.print   ("higgs_4e_"   +mass_str+".pdf");
        fitmaker_g4mu.print  ("higgs_4mu_"  +mass_str+".pdf");
        
        fitmaker_zx2e2mu.fit();
        fitmaker_zx4e.fit();
        fitmaker_zx4mu.fit();
        
        fitmaker_zx2e2mu.print("zjets_2e2mu_"+mass_str+".pdf", 15);
        fitmaker_zx4e.print   ("zjets_4e_"   +mass_str+".pdf", 15);
        fitmaker_zx4mu.print  ("zjets_4mu_"  +mass_str+".pdf", 15);
        
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
            fitmaker_g2e2mu.makeWorkspace1D (w_2e2mu , CMS_zz4l_mass_1D);
            fitmaker_g4e.makeWorkspace1D    (w_4e    , CMS_zz4l_mass_1D);
            fitmaker_g4mu.makeWorkspace1D   (w_4mu   , CMS_zz4l_mass_1D);
            
            fitmaker_q2e2mu.makeWorkspace1D (w_2e2mu , CMS_zz4l_mass_1D);
            fitmaker_q4e.makeWorkspace1D    (w_4e    , CMS_zz4l_mass_1D);
            fitmaker_q4mu.makeWorkspace1D   (w_4mu   , CMS_zz4l_mass_1D);

            fitmaker_zx2e2mu.makeWorkspace1D(w_2e2mu , CMS_zz4l_mass_1D);
            fitmaker_zx4e.makeWorkspace1D   (w_4e    , CMS_zz4l_mass_1D);
            fitmaker_zx4mu.makeWorkspace1D  (w_4mu   , CMS_zz4l_mass_1D);
            
            fitmaker_2e2mu.makeWorkspace1D  (w_2e2mu , CMS_zz4l_mass_1D);
            fitmaker_4e.makeWorkspace1D     (w_4e    , CMS_zz4l_mass_1D);
            fitmaker_4mu.makeWorkspace1D    (w_4mu   , CMS_zz4l_mass_1D);
            
            fitmaker_gz2e2mu.makeWorkspace1D(w_2e2mu , CMS_zz4l_mass_1D);
            fitmaker_gz4e.makeWorkspace1D   (w_4e    , CMS_zz4l_mass_1D);
            fitmaker_gz4mu.makeWorkspace1D  (w_4mu   , CMS_zz4l_mass_1D);

            RooArgSet argset_obs(CMS_zz4l_mass_1D, "argset_obs");
            RooDataSet data_obs("data_obs", "data_obs", argset_obs);
           
            ymaker_data.getDataSet1D(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, data_obs, CMS_zz4l_mass_1D);
            ymaker_data.getDataSet1D(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, data_obs, CMS_zz4l_mass_1D);
            ymaker_data.getDataSet1D(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, data_obs, CMS_zz4l_mass_1D);
            ymaker_data.getDataSet1D(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, data_obs, CMS_zz4l_mass_1D);
            
            w_2e2mu.import(data_obs);
            w_4e.import   (data_obs);
            w_4mu.import  (data_obs);

        }      
        else {
            RooRealVar CMS_zz4l_mass_2D("CMS_zz4l_mass_2D", "M(4l)", hinfo.massLow, hinfo.massHigh, "GeV/c^{2}");
            RooRealVar CMS_zz4l_melaLD ("CMS_zz4l_melaLD" , "MELA" ,   0,             1,            "");
            CMS_zz4l_melaLD.setBins(30);            

            ymaker_ghzz.get2DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_sig_mm);
            ymaker_ghzz.get2DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_sig_ee);
            ymaker_ghzz.get2DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_sig_em);
            ymaker_ghzz.get2DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_sig_em);
            
            ymaker_qhzz.get2DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_sig_mm);
            ymaker_qhzz.get2DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_sig_ee);
            ymaker_qhzz.get2DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_sig_em);
            ymaker_qhzz.get2DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_sig_em);
            
            ymaker_qqzz.get2DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_bkg_mm);
            ymaker_qqzz.get2DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_bkg_ee);
            ymaker_qqzz.get2DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_bkg_em);
            ymaker_qqzz.get2DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_bkg_em);
            
            ymaker_ggzz.get2DHist(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_bkg_mm);
            ymaker_ggzz.get2DHist(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_bkg_ee);
            ymaker_ggzz.get2DHist(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_bkg_em);
            ymaker_ggzz.get2DHist(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, &h2D_bkg_em);
            
            fitmaker_g2e2mu.makeWorkspace2D(w_2e2mu  , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_sig_em);
            fitmaker_g4e.makeWorkspace2D   (w_4e     , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_sig_ee);
            fitmaker_g4mu.makeWorkspace2D  (w_4mu    , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_sig_mm);
            
            fitmaker_q2e2mu.makeWorkspace2D(w_2e2mu  , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_sig_em);
            fitmaker_q4e.makeWorkspace2D   (w_4e     , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_sig_ee);
            fitmaker_q4mu.makeWorkspace2D  (w_4mu    , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_sig_mm);

            fitmaker_zx2e2mu.makeWorkspace2D(w_2e2mu , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_bkg_em);
            fitmaker_zx4e.makeWorkspace2D   (w_4e    , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_bkg_em);
            fitmaker_zx4mu.makeWorkspace2D  (w_4mu   , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_bkg_em);
            
            fitmaker_2e2mu.makeWorkspace2D  (w_2e2mu , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_bkg_em);
            fitmaker_4e.makeWorkspace2D     (w_4e    , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_bkg_em);
            fitmaker_4mu.makeWorkspace2D    (w_4mu   , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_bkg_em);
            
            fitmaker_gz2e2mu.makeWorkspace2D(w_2e2mu , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_bkg_em);
            fitmaker_gz4e.makeWorkspace2D   (w_4e    , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_bkg_em);
            fitmaker_gz4mu.makeWorkspace2D  (w_4mu   , CMS_zz4l_mass_2D, CMS_zz4l_melaLD, &h2D_bkg_em);
        
            RooArgSet argset_obs(CMS_zz4l_mass_2D, CMS_zz4l_melaLD, "argset_obs");
            RooDataSet data_obs("data_obs", "data_obs", argset_obs);
        
            ymaker_data.getDataSet2D(0, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, data_obs, CMS_zz4l_mass_2D, CMS_zz4l_melaLD);              
            ymaker_data.getDataSet2D(1, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, data_obs, CMS_zz4l_mass_2D, CMS_zz4l_melaLD);
            ymaker_data.getDataSet2D(2, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, data_obs, CMS_zz4l_mass_2D, CMS_zz4l_melaLD);
            ymaker_data.getDataSet2D(3, hinfo.z1min, hinfo.z2min, hinfo.massLow, hinfo.massHigh, hinfo.melacut, data_obs, CMS_zz4l_mass_2D, CMS_zz4l_melaLD);
 
            w_2e2mu.import(data_obs);
            w_4e.import   (data_obs);
            w_4mu.import  (data_obs);
        }

        w_2e2mu.writeToFile(workspace_2e2mu.c_str());
        w_4e.writeToFile(workspace_4e.c_str());
        w_4mu.writeToFile(workspace_4mu.c_str());
    }

}


float       HiggsMassPointInfo::lumi = 5.05;
float       HiggsMassPointInfo::z1min = 40.;
float       HiggsMassPointInfo::z2min = 12.;
float       HiggsMassPointInfo::massLow = 100.;
float       HiggsMassPointInfo::massHigh = 300.;
float       HiggsMassPointInfo::massLowBkgFit = 100.;
float       HiggsMassPointInfo::massHighBkgFit = 600.;
float       HiggsMassPointInfo::melacut = -1.0;
int         HiggsMassPointInfo::nBinsMass2D = 100;
int         HiggsMassPointInfo::nBinsMELA2D = 30;
bool        HiggsMassPointInfo::doShapeAnalysis = true;
bool        HiggsMassPointInfo::do1D = true;
bool        HiggsMassPointInfo::doSS = true;
std::string HiggsMassPointInfo::treeFolder = "/home/avartak/CMS/Higgs/HZZ4L/CMSSW_4_2_8_patch7/src/WWAnalysis/AnalysisStep/trees/";

void doHZZAnalysis() {

    HiggsMassPointInfo h115(115, 100, 130, 200, 250);
    HiggsMassPointInfo h120(120, 100, 135, 201, 251);
    HiggsMassPointInfo h130(130, 110, 145, 202, 252);
    HiggsMassPointInfo h140(140, 120, 155, 203, 253);
    HiggsMassPointInfo h150(150, 130, 165, 204, 254);
    HiggsMassPointInfo h160(160, 140, 175, 205, 255);
    HiggsMassPointInfo h170(170, 150, 185, 206, 256);
    HiggsMassPointInfo h180(180, 160, 195, 207, 257);
    HiggsMassPointInfo h190(190, 170, 205, 208, 258);
    HiggsMassPointInfo h200(200, 180, 215, 209, 259);
    HiggsMassPointInfo h210(210, 190, 225, 210, 260);

    bool hinfodoshape = HiggsMassPointInfo::doShapeAnalysis;
    HiggsMassPointInfo::doShapeAnalysis = false;
    analysisEngine(h115);
    HiggsMassPointInfo::doShapeAnalysis = hinfodoshape;

    if (HiggsMassPointInfo::doShapeAnalysis) {
        //analysisEngine(h115);
        //analysisEngine(h120);
        //analysisEngine(h130);
        //analysisEngine(h140);
        //analysisEngine(h150);
        //analysisEngine(h160);
        //analysisEngine(h170);
        //analysisEngine(h180);
        analysisEngine(h190);
        analysisEngine(h200);
        analysisEngine(h210);
    }

    else {
        analysisEngine(h115);
        analysisEngine(h120);
        analysisEngine(h130);
        analysisEngine(h140);
        analysisEngine(h150);
        analysisEngine(h160);
        analysisEngine(h170);
        analysisEngine(h180);
        analysisEngine(h190);
        analysisEngine(h200);
        analysisEngine(h210);
    }
}
