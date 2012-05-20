#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH2D.h>
#include <TBranch.h>
#include <vector>
#include <iostream>
#include <utility>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "WWAnalysis/AnalysisStep/interface/scales.h"
#include "WWAnalysis/AnalysisStep/interface/findAndReplace.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"

#include <RooAddPdf.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooDataSet.h>
#include <RooExtendPdf.h>
#include <RooFitResult.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooRealVar.h> 
#include <RooWorkspace.h> 
#include <RooLandau.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h> 
#include <RooFFTConvPdf.h>

using namespace RooFit;


float mubarrel[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float muendcap[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float elbarrel[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float elendcap[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

float mubarrelerr[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float muendcaperr[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float elbarrelerr[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float elendcaperr[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

void domufakes(std::string path, float* barrel, float* endcap, float* barrelerr, float* endcaperr) {
    TFile* file = new TFile(path.c_str());
    TTree* tree = (TTree*)file->Get("zllmtree/probe_tree");

    Float_t bins[] = {0.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,80.0};
    Int_t binnum = 10;
    Float_t ybins[] = {0.0,1.479,2.5};
    Int_t ybinnum = 2;


    TH2D* hfake = new TH2D("hfake", "", binnum, bins, ybinnum, ybins);
    TH2D* hall  = new TH2D("hall",  "", binnum, bins, ybinnum, ybins);
    TH2D* hpass = new TH2D("hpass", "", binnum, bins, ybinnum, ybins);

    hall  ->Sumw2();
    hpass ->Sumw2();
    hfake ->Sumw2();

    TBranch *bmet     = tree->GetBranch("met");
    TBranch *bzmass   = tree->GetBranch("zmass");
    TBranch *bpt      = tree->GetBranch("pt");
    TBranch *beta     = tree->GetBranch("eta");
    TBranch *bid      = tree->GetBranch("pfid");
    TBranch *bbdtiso  = tree->GetBranch("bdtIso");
    TBranch *bchiso   = tree->GetBranch("pfIsoChHad04");
    TBranch *bphiso   = tree->GetBranch("pfIsoPhoton04_NoEA");
    TBranch *bnhiso   = tree->GetBranch("pfIsoNHad04_NoEA");
    TBranch *brho     = tree->GetBranch("rhoAA");
    TBranch *brun     = tree->GetBranch("run");
    TBranch *blum     = tree->GetBranch("lumi");
    TBranch *bevt     = tree->GetBranch("event");

    float met     = 0.0;
    float zmass   = 0.0;
    float pt      = 0.0;
    float eta     = 0.0;
    int   id      = 0;
    float bdtiso  = 0.0;
    float chiso   = 0.0;
    float phiso   = 0.0;
    float nhiso   = 0.0;
    float rho     = 0.0;
    int   run     = 0;
    int   lum     = 0;
    int   evt     = 0;

    bmet      ->SetAddress(&met);
    bzmass    ->SetAddress(&zmass);
    bpt       ->SetAddress(&pt);
    beta      ->SetAddress(&eta);
    bid       ->SetAddress(&id);
    bbdtiso   ->SetAddress(&bdtiso);
    bchiso    ->SetAddress(&chiso);
    bphiso    ->SetAddress(&phiso);
    bnhiso    ->SetAddress(&nhiso);
    brho      ->SetAddress(&rho);
    brun      ->SetAddress(&run);
    blum      ->SetAddress(&lum);
    bevt      ->SetAddress(&evt);

    MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_ = MuonEffectiveArea::kMuEAData2011;
    MuonEffectiveArea::MuonEffectiveAreaType   effArea_ = MuonEffectiveArea::kMuGammaAndNeutralHadronIso04;

    for (int i = 0; i < tree->GetEntries(); i++) {
        bmet     ->GetEvent(i);
        bzmass   ->GetEvent(i);
        bpt      ->GetEvent(i);
        beta     ->GetEvent(i);
        bid      ->GetEvent(i);
        bbdtiso  ->GetEvent(i);
        bchiso   ->GetEvent(i);
        bphiso   ->GetEvent(i);
        bnhiso   ->GetEvent(i);
        brho     ->GetEvent(i);
        brun     ->GetEvent(i);
        blum     ->GetEvent(i);
        bevt     ->GetEvent(i);

        float eff_area = MuonEffectiveArea::GetMuonEffectiveArea(effArea_, fabs(eta), effAreaTarget_);

        float iso = (nhiso + phiso - (eff_area)*rho);
        if (iso<0.0) iso = 0.0;
        iso += chiso;
        iso /= pt;

        if (zmass>40 && zmass<120) {
            if (fabs(eta)<1.479) {
                hall->Fill(pt,fabs(eta));
                if (id>0 && iso<0.4) {hpass->Fill(pt, fabs(eta));}
            }
            else {
                hall->Fill(pt,fabs(eta));
                if (id>0 && iso<0.4) {hpass->Fill(pt, fabs(eta));}
            }
        }

    }

    hfake ->Divide(hpass,  hall);
    for (std::size_t k = 0; k < 10; k++) {
        barrel[k] = hfake->GetBinContent(k+1, 1);
        endcap[k] = hfake->GetBinContent(k+1, 2);
        barrelerr[k] = hfake->GetBinError(k+1, 1);
        endcaperr[k] = hfake->GetBinError(k+1, 2);
    }

    std::cout << "---- Computing fake rates for muons ----" << std::endl;
    std::cout << "---- Muon Barrel ----" << std::endl;
    std::cout << "pT [ 0 -  5] : " << mubarrel[0] << " +/- " << mubarrelerr[0] << endl;
    std::cout << "pT [ 7 -  7] : " << mubarrel[1] << " +/- " << mubarrelerr[1] << endl;
    std::cout << "pT [ 7 - 10] : " << mubarrel[2] << " +/- " << mubarrelerr[2] << endl;
    std::cout << "pT [10 - 15] : " << mubarrel[3] << " +/- " << mubarrelerr[3] << endl;
    std::cout << "pT [15 - 20] : " << mubarrel[4] << " +/- " << mubarrelerr[4] << endl;
    std::cout << "pT [20 - 25] : " << mubarrel[5] << " +/- " << mubarrelerr[5] << endl;
    std::cout << "pT [25 - 30] : " << mubarrel[6] << " +/- " << mubarrelerr[6] << endl;
    std::cout << "pT [30 - 40] : " << mubarrel[7] << " +/- " << mubarrelerr[7] << endl;
    std::cout << "pT [40 - 50] : " << mubarrel[8] << " +/- " << mubarrelerr[8] << endl;
    std::cout << "pT [50 - 80] : " << mubarrel[9] << " +/- " << mubarrelerr[9] << endl;
    std::cout << std::endl;
    std::cout << "---- Muon Endcap ----" << std::endl;
    std::cout << "pT [ 0 -  5] : " << muendcap[0] << " +/- " << muendcaperr[0] << endl;
    std::cout << "pT [ 7 -  7] : " << muendcap[1] << " +/- " << muendcaperr[1] << endl;
    std::cout << "pT [ 7 - 10] : " << muendcap[2] << " +/- " << muendcaperr[2] << endl;
    std::cout << "pT [10 - 15] : " << muendcap[3] << " +/- " << muendcaperr[3] << endl;
    std::cout << "pT [15 - 20] : " << muendcap[4] << " +/- " << muendcaperr[4] << endl;
    std::cout << "pT [20 - 25] : " << muendcap[5] << " +/- " << muendcaperr[5] << endl;
    std::cout << "pT [25 - 30] : " << muendcap[6] << " +/- " << muendcaperr[6] << endl;
    std::cout << "pT [30 - 40] : " << muendcap[7] << " +/- " << muendcaperr[7] << endl;
    std::cout << "pT [40 - 50] : " << muendcap[8] << " +/- " << muendcaperr[8] << endl;
    std::cout << "pT [50 - 80] : " << muendcap[9] << " +/- " << muendcaperr[9] << endl;
    std::cout << std::endl;

    file->Close();
    delete file;
}


void doelfakes(std::string path, float* barrel, float* endcap, float* barrelerr, float*endcaperr) {
    TFile* file = new TFile(path.c_str());
    TTree* tree = (TTree*)file->Get("zlletree/probe_tree");
    
    Float_t bins[] = {0.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,80.0};
    Int_t binnum = 10;
    Float_t ybins[] = {0.0,1.479,2.5};
    Int_t ybinnum = 2;
    

    TH2D* hfake = new TH2D("hfake", "", binnum, bins, ybinnum, ybins);
    TH2D* hall  = new TH2D("hall",  "", binnum, bins, ybinnum, ybins);
    TH2D* hpass = new TH2D("hpass", "", binnum, bins, ybinnum, ybins);


    hall  ->Sumw2();
    hpass ->Sumw2();
    hfake ->Sumw2();

    TBranch *bzmass   = tree->GetBranch("zmass");
    TBranch *bpt      = tree->GetBranch("pt");
    TBranch *beta     = tree->GetBranch("eta");
    TBranch *bid      = tree->GetBranch("bdtIdNtDz");
    TBranch *bbdtiso  = tree->GetBranch("bdtIsoNtDz");
    TBranch *bmhits   = tree->GetBranch("nmisshits");
    TBranch *bchiso   = tree->GetBranch("pfIsoChHad04");
    TBranch *bphiso   = tree->GetBranch("pfIsoPhoton04_NoEA");
    TBranch *bnhiso   = tree->GetBranch("pfIsoNHad04_NoEA");
    TBranch *brho     = tree->GetBranch("rhoAA");
    TBranch *brun     = tree->GetBranch("run");
    TBranch *blum     = tree->GetBranch("lumi");
    TBranch *bevt     = tree->GetBranch("event");

    float zmass   = 0.0;
    float pt      = 0.0;
    float eta     = 0.0;
    float id      = 0.0;
    float bdtiso  = 0.0;
    float mhits   = 0.0;
    float chiso   = 0.0;
    float phiso   = 0.0;
    float nhiso   = 0.0;
    float rho     = 0.0;
    int   run     = 0;
    int   lum     = 0;
    int   evt     = 0;

    bzmass    ->SetAddress(&zmass);
    bpt       ->SetAddress(&pt);
    beta      ->SetAddress(&eta);
    bid       ->SetAddress(&id);
    bbdtiso   ->SetAddress(&bdtiso);
    bmhits    ->SetAddress(&mhits);
    bchiso    ->SetAddress(&chiso);
    bphiso    ->SetAddress(&phiso);
    bnhiso    ->SetAddress(&nhiso);
    brho      ->SetAddress(&rho);
    brun      ->SetAddress(&run);
    blum      ->SetAddress(&lum);
    bevt      ->SetAddress(&evt);

    ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_ = ElectronEffectiveArea::kEleEAData2011;
    ElectronEffectiveArea::ElectronEffectiveAreaType   effArea_ = ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04;

    for (int i = 0; i < tree->GetEntries(); i++) {
        bzmass   ->GetEvent(i);
        bpt      ->GetEvent(i);
        beta     ->GetEvent(i);
        bid      ->GetEvent(i);
        bbdtiso  ->GetEvent(i);
        bmhits   ->GetEvent(i);
        bchiso   ->GetEvent(i);
        bphiso   ->GetEvent(i);
        bnhiso   ->GetEvent(i);
        brho     ->GetEvent(i);
        brun     ->GetEvent(i);
        blum     ->GetEvent(i);
        bevt     ->GetEvent(i);

        float eff_area = ElectronEffectiveArea::GetElectronEffectiveArea(effArea_, fabs(eta), effAreaTarget_);

        float iso = (nhiso + phiso - (eff_area)*rho);
        if (iso<0.0) iso = 0.0;
        iso += chiso;
        iso /= pt;

        if (mhits <= 1 && zmass>40 && zmass<120) {
            if (fabs(eta)<1.479) {
                hall->Fill(pt,fabs(eta));
                if (pt <= 10.0 && ((fabs(eta)< 0.8 && id > 0.47) || (fabs(eta)>=0.8 && id>0.04)) && iso<0.4) {hpass->Fill(pt, fabs(eta));}
                if (pt >  10.0 && ((fabs(eta)< 0.8 && id > 0.50) || (fabs(eta)>=0.8 && id>0.12)) && iso<0.4) {hpass->Fill(pt, fabs(eta));}
            }
            else {
                hall->Fill(pt,fabs(eta));
                if (pt <= 10.0 && id>0.295 && iso<0.4) {hpass->Fill(pt, fabs(eta));}
                if (pt >  10.0 && id>0.600 && iso<0.4) {hpass->Fill(pt, fabs(eta));}
            }
        }

    }

    hfake ->Divide(hpass , hall);
    for (std::size_t k = 0; k < 10; k++) {
        barrel[k] = hfake->GetBinContent(k+1, 1);
        endcap[k] = hfake->GetBinContent(k+1, 2);
        barrelerr[k] = hfake->GetBinError(k+1, 1);
        endcaperr[k] = hfake->GetBinError(k+1, 2);
    }

    std::cout << "---- Computing fake rates for electrons ----" << std::endl;
    std::cout << "---- Electron Barrel ----" << std::endl;
    std::cout << "pT [ 0 -  5] : " << elbarrel[0] << " +/- " << elbarrelerr[0] << endl;
    std::cout << "pT [ 7 -  7] : " << elbarrel[1] << " +/- " << elbarrelerr[1] << endl;
    std::cout << "pT [ 7 - 10] : " << elbarrel[2] << " +/- " << elbarrelerr[2] << endl;
    std::cout << "pT [10 - 15] : " << elbarrel[3] << " +/- " << elbarrelerr[3] << endl;
    std::cout << "pT [15 - 20] : " << elbarrel[4] << " +/- " << elbarrelerr[4] << endl;
    std::cout << "pT [20 - 25] : " << elbarrel[5] << " +/- " << elbarrelerr[5] << endl;
    std::cout << "pT [25 - 30] : " << elbarrel[6] << " +/- " << elbarrelerr[6] << endl;
    std::cout << "pT [30 - 40] : " << elbarrel[7] << " +/- " << elbarrelerr[7] << endl;
    std::cout << "pT [40 - 50] : " << elbarrel[8] << " +/- " << elbarrelerr[8] << endl;
    std::cout << "pT [50 - 80] : " << elbarrel[9] << " +/- " << elbarrelerr[9] << endl;
    std::cout << std::endl;
    std::cout << "---- Electron Endcap ----" << std::endl;
    std::cout << "pT [ 0 -  5] : " << elendcap[0] << " +/- " << elendcaperr[0] << endl;
    std::cout << "pT [ 7 -  7] : " << elendcap[1] << " +/- " << elendcaperr[1] << endl;
    std::cout << "pT [ 7 - 10] : " << elendcap[2] << " +/- " << elendcaperr[2] << endl;
    std::cout << "pT [10 - 15] : " << elendcap[3] << " +/- " << elendcaperr[3] << endl;
    std::cout << "pT [15 - 20] : " << elendcap[4] << " +/- " << elendcaperr[4] << endl;
    std::cout << "pT [20 - 25] : " << elendcap[5] << " +/- " << elendcaperr[5] << endl;
    std::cout << "pT [25 - 30] : " << elendcap[6] << " +/- " << elendcaperr[6] << endl;
    std::cout << "pT [30 - 40] : " << elendcap[7] << " +/- " << elendcaperr[7] << endl;
    std::cout << "pT [40 - 50] : " << elendcap[8] << " +/- " << elendcaperr[8] << endl;
    std::cout << "pT [50 - 80] : " << elendcap[9] << " +/- " << elendcaperr[9] << endl;
    std::cout << std::endl;

   
    file->Close();
    delete file;
}

float getFakeRate(float pt, float eta, float id) {

    int bin = 0.0;
    if (pt>0  && pt<=5)  bin = 0;    
    if (pt>5  && pt<=7)  bin = 1;    
    if (pt>7  && pt<=10) bin = 2;    
    if (pt>10 && pt<=15) bin = 3;    
    if (pt>15 && pt<=20) bin = 4;    
    if (pt>20 && pt<=25) bin = 5;    
    if (pt>25 && pt<=30) bin = 6;    
    if (pt>30 && pt<=40) bin = 7;    
    if (pt>40 && pt<=50) bin = 8;    
    if (pt>50 && pt<=80) bin = 9;    

    if      (fabs(id) == 13 && fabs(eta)<1.479) return mubarrel[bin];
    else if (fabs(id) == 13 && fabs(eta)>1.479) return muendcap[bin];
    else if (fabs(id) == 11 && fabs(eta)<1.479) return elbarrel[bin];
    else if (fabs(id) == 11 && fabs(eta)>1.479) return elendcap[bin];
    else return 0.0;

}

float getFakeRateErr(float pt, float eta, float id) {

    int bin = 0.0;
    if (pt>0  && pt<=5)  bin = 0;    
    if (pt>5  && pt<=7)  bin = 1;    
    if (pt>7  && pt<=10) bin = 2;    
    if (pt>10 && pt<=15) bin = 3;    
    if (pt>15 && pt<=20) bin = 4;    
    if (pt>20 && pt<=25) bin = 5;    
    if (pt>25 && pt<=30) bin = 6;    
    if (pt>30 && pt<=40) bin = 7;    
    if (pt>40 && pt<=50) bin = 8;    
    if (pt>50 && pt<=80) bin = 9;    

    if      (fabs(id) == 13 && fabs(eta)<1.479) return mubarrelerr[bin];
    else if (fabs(id) == 13 && fabs(eta)>1.479) return muendcaperr[bin];
    else if (fabs(id) == 11 && fabs(eta)<1.479) return elbarrelerr[bin];
    else if (fabs(id) == 11 && fabs(eta)>1.479) return elendcaperr[bin];
    else return 0.0;

}

std::vector<float> getYield(std::string path, int ch, float scale, bool isMC, float z1min, float z2min, float min4l, float max4l, RooDataSet& dataset, RooArgSet& argset, int zxMode) {

    float yield = 0.0;
    float yield_count = 0.0;
    float yield_err = 0.0;
    float yield_up = 0.0;
    float yield_dn = 0.0;

    TFile* file = new TFile(path.c_str());
    TTree* tree = isMC ? (TTree*)file->Get("zz4lTree/probe_tree") : (TTree*)file->Get("zxTree/probe_tree");
    
    TBranch *bchannel   = tree->GetBranch("channel");
    TBranch *bz1mass    = tree->GetBranch("z1mass");
    TBranch *bz2mass    = tree->GetBranch("z2mass");
    TBranch *bmass      = tree->GetBranch("mass");
    TBranch *bnumsim    = tree->GetBranch("intimeSimVertices");
    TBranch *bl3pt      = tree->GetBranch("l3pt");    
    TBranch *bl3eta     = tree->GetBranch("l3eta");    
    TBranch *bl3id      = tree->GetBranch("l3idNew");    
    TBranch *bl3iso     = tree->GetBranch("l3pfIsoComb04EACorr");    
    TBranch *bl3pdgId   = tree->GetBranch("l3pdgId");    
    TBranch *bl4pt      = tree->GetBranch("l4pt");    
    TBranch *bl4eta     = tree->GetBranch("l4eta");    
    TBranch *bl4id      = tree->GetBranch("l4idNew");    
    TBranch *bl4iso     = tree->GetBranch("l4pfIsoComb04EACorr");    
    TBranch *bl4pdgId   = tree->GetBranch("l4pdgId");    
    TBranch *bevent     = tree->GetBranch("event");    
    TBranch *brun       = tree->GetBranch("run");    

    float channel   = 0.0;
    float z1mass    = 0.0;
    float z2mass    = 0.0;
    float mass      = 0.0;
    float numsim    = 0.0;
    float l3pt      = 0.0;
    float l3eta     = 0.0;
    float l3id      = 0.0;
    float l3iso     = 0.0;
    float l3pdgId   = 0.0;
    float l4pt      = 0.0;
    float l4eta     = 0.0;
    float l4id      = 0.0;
    float l4iso     = 0.0;
    float l4pdgId   = 0.0;
    int   event     = 0;
    int   run       = 0;
    
    bchannel   ->SetAddress(&channel);
    bz1mass    ->SetAddress(&z1mass);
    bz2mass    ->SetAddress(&z2mass);
    bmass      ->SetAddress(&mass);
    bnumsim    ->SetAddress(&numsim);
    bl3pt      ->SetAddress(&l3pt);
    bl3eta     ->SetAddress(&l3eta);
    bl3id      ->SetAddress(&l3id);
    bl3iso     ->SetAddress(&l3iso);
    bl3pdgId   ->SetAddress(&l3pdgId);
    bl4pt      ->SetAddress(&l4pt);
    bl4eta     ->SetAddress(&l4eta);
    bl4id      ->SetAddress(&l4id);
    bl4iso     ->SetAddress(&l4iso);
    bl4pdgId   ->SetAddress(&l4pdgId);
    bevent     ->SetAddress(&event);
    brun       ->SetAddress(&run);
    
    
    std::vector<std::pair<int, int> > runeventinfo;

    for (int i = 0; i < tree->GetEntries(); i++) {
        bchannel   ->GetEvent(i);
        bz1mass    ->GetEvent(i);
        bz2mass    ->GetEvent(i);
        bmass      ->GetEvent(i);
        bnumsim    ->GetEvent(i);
        bl3pt      ->GetEvent(i);    
        bl3eta     ->GetEvent(i);    
        bl3id      ->GetEvent(i);    
        bl3iso     ->GetEvent(i);    
        bl3pdgId   ->GetEvent(i);    
        bl4pt      ->GetEvent(i);    
        bl4eta     ->GetEvent(i);    
        bl4id      ->GetEvent(i);    
        bl4iso     ->GetEvent(i);    
        bl4pdgId   ->GetEvent(i);    
        bevent     ->GetEvent(i);
        brun       ->GetEvent(i);

        if (mass>min4l && mass<max4l && channel == (int)ch && z1mass>z1min && z2mass>z2min) {
    
            if (isMC) {
                yield += scale * getPUWeight((int)numsim);
                yield_err += (scale * getPUWeight((int)numsim))*(scale * getPUWeight((int)numsim));
                yield_count += 1.0;
                argset.setRealValue("mass", mass);
                argset.setRealValue("weight", scale * getPUWeight((int)numsim));
                dataset.add(argset,scale * getPUWeight((int)numsim));
            }

            else {
                bool existsAlready = false;
                for (std::size_t k = 0; k < runeventinfo.size(); k++) {
                    if (run == runeventinfo[k].first && event == runeventinfo[k].second) existsAlready = true;
                }

                if (!existsAlready) {
                    if (l3pdgId == -l4pdgId && (l3id==0 || l3iso/l3pt>0.4) && (l4id==0 || l4iso/l4pt>0.4)) {
                        runeventinfo.push_back(std::pair<int, int>(run, event));
                        float f1    = getFakeRate(l3pt, l3eta, l3pdgId);
                        float f2    = getFakeRate(l4pt, l4eta, l4pdgId);
                        float f1_up = getFakeRate(l3pt, l3eta, l3pdgId) + getFakeRateErr(l3pt, l3eta, l3pdgId);
                        float f2_up = getFakeRate(l4pt, l4eta, l4pdgId) + getFakeRateErr(l4pt, l4eta, l4pdgId);
                        float f1_dn = getFakeRate(l3pt, l3eta, l3pdgId) - getFakeRateErr(l3pt, l3eta, l3pdgId);
                        float f2_dn = getFakeRate(l4pt, l4eta, l4pdgId) - getFakeRateErr(l4pt, l4eta, l4pdgId);
                        yield += (f1/(1-f1))*(f2/(1-f2));
                        yield_up += (f1_up/(1-f1_up))*(f2_up/(1-f2_up));
                        yield_dn += (f1_dn/(1-f1_dn))*(f2_dn/(1-f2_dn));
                        yield_count += 1.0;
                        if (zxMode == 0) {
                            argset.setRealValue("mass", mass);
                            argset.setRealValue("weight", (f1/(1-f1))*(f2/(1-f2)));
                            dataset.add(argset, (f1/(1-f1))*(f2/(1-f2)));
                        }
                        else if (zxMode == 1) {
                            argset.setRealValue("mass", mass);
                            argset.setRealValue("weight", (f1_up/(1-f1_up))*(f2/(1-f2)));
                            dataset.add(argset, (f1_up/(1-f1_up))*(f2/(1-f2)));
                        }
                        else if (zxMode == 2) {
                            argset.setRealValue("mass", mass);
                            argset.setRealValue("weight", (f1_dn/(1-f1_dn))*(f2/(1-f2)));
                            dataset.add(argset, (f1_dn/(1-f1_dn))*(f2/(1-f2)));
                        }
                    }
                }
            }    

        }
    
    }

    if (isMC) yield_err = sqrt(yield_err);
    else yield_err = std::max<float>(fabs(yield_up-yield), fabs(yield_dn-yield));

    file->Close();
    delete file;


    std::vector<float> results;
    results.push_back(yield);
    results.push_back(yield_count);
    results.push_back(yield_err);
    return results;

}


std::string createCardTemplate(int channel) {

    std::string card = "";
        card += "imax 1\n";
        card += "jmax *\n";
        card += "kmax *\n";
        card += "------------\n";
        card += "bin         BIN\n";
        card += "observation OBS\n";
        card += "------------\n";
        card += "bin     BIN      BIN      BIN       BIN       BIN       BIN\n";
        card += "process sig_ggH sig_VBF sig_WH   bkg_qqzz bkg_ggzz bkg_zjets\n";
        card += "process -2     -1       0        1        2        3\n";
        card += "rate    SIG_GGH SIG_VBF SIG_WH   BKG_QQZZ BKG_GGZZ BKG_ZJETS\n";
        card += "------------\n";
        card += "lumi                      lnN        1.045  1.045  1.045  1.045   1.045  -\n";
        card += "pdf_gg                    lnN        1.0755 -      -      -       1.0708 -\n";
        card += "pdf_qqbar                 lnN         -     1.023  1.035  1.0341  -      -\n";
        card += "pdf_H4l_accept            lnN        1.02   1.02   1.02   -       -      -\n";
        card += "QCDscale_ggH              lnN        1.0975 -      -      -       -      -\n";
        card += "QCDscale_qqH              lnN        -      1.003  -      -       -      -\n";
        card += "QCDscale_VH               lnN        -      -      1.015  -       -      -\n";
        card += "theoryUncXS_HighMH        lnN        1      1      1      -       -      -\n";
        card += "QCDscale_ggVV             lnN        -      -      -      -       1.2431 -\n";
        card += "QCDscale_VV               lnN        -      -      -      1.0284  -      -\n";
        card += "BRhiggs_ZZ4l              lnN        1.02   1.02   1.02   -       -      -\n";
    if (channel == 0 || channel == 2) {
        card +=  "CMS_eff_m                 lnN        1.0224 1.0224 1.0224 1.0224  1.0224 -\n";
        card +=  "CMS_scale_m               lnN        1.01   1.01   1.01   1.01    1.01   -\n";
        card +=  "CMS_trigger_m             lnN        1.015  1.015  1.015  1.015   1.015  -\n";
    }
    if (channel == 1 || channel == 2) {
        card += "CMS_eff_e                 lnN        1.0447 1.0447 1.0447 1.0447  1.0447 -\n";
        card += "CMS_scale_e               lnN        1.02   1.02   1.02   1.02    1.02   -\n";
        card += "CMS_trigger_e             lnN        1.015  1.015  1.015  1.015   1.015  -\n";
    }
    if (channel == 0) {
        card += "CMS_hzz4mu_Zjets          gmN ZJEVT  -      -      -      -       -      ZJALPHA\n";
        card += "CMS_hzz4mu_Zjets_EXTRAP   lnN        -      -      -      -       -      ZJSYST\n";
    }
    if (channel == 1) {
        card += "CMS_hzz4e_Zjets           gmN ZJEVT  -      -      -      -       -      ZJALPHA\n";
        card += "CMS_hzz4e_Zjets_EXTRAP    lnN        -      -      -      -       -      ZJSYST\n";
    }
    if (channel == 2) {
        card += "CMS_hzz2e2mu_Zjets        gmN ZJEVT  -      -      -      -       -      ZJALPHA\n";
        card += "CMS_hzz2e2mu_Zjets_EXTRAP lnN        -      -      -      -       -      ZJSYST\n";
    }

    return card;
}


std::string createShapeCardTemplate(int channel, std::string outfilename) {

    std::string card = "";
        card += "imax 1\n";
        card += "jmax *\n";
        card += "kmax *\n";
        card += "------------\n";
    if (channel == 0) {
        card += std::string("shapes * * ") + outfilename + " w_4mu:$PROCESS\n";
    }
    if (channel == 1) {
        card += std::string("shapes * * ") + outfilename + " w_4e:$PROCESS\n";
    }
    if (channel == 2) {
        card += std::string("shapes * * ") + outfilename + " w_2e2mu:$PROCESS\n";
    }
        card += "------------\n";
        card += "bin         BIN\n";
        card += "observation OBS\n";
        card += "------------\n";
        card += "bin     BIN      BIN      BIN       BIN       BIN       BIN\n";
        card += "process sig_ggH sig_VBF sig_WH   bkg_qqzz bkg_ggzz bkg_zjets\n";
        card += "process -2     -1       0        1        2        3\n";
        card += "rate    SIG_GGH SIG_VBF SIG_WH   BKG_QQZZ BKG_GGZZ BKG_ZJETS\n";
        card += "------------\n";
        card += "lumi                      lnN        1.045  1.045  1.045  1.045   1.045  -\n";
        card += "pdf_gg                    lnN        1.0755 -      -      -       1.0708 -\n";
        card += "pdf_qqbar                 lnN         -     1.023  1.035  1.0341  -      -\n";
        card += "pdf_H4l_accept            lnN        1.02   1.02   1.02   -       -      -\n";
        card += "QCDscale_ggH              lnN        1.0975 -      -      -       -      -\n";
        card += "QCDscale_qqH              lnN        -      1.003  -      -       -      -\n";
        card += "QCDscale_VH               lnN        -      -      1.015  -       -      -\n";
        card += "theoryUncXS_HighMH        lnN        1      1      1      -       -      -\n";
        card += "QCDscale_ggVV             lnN        -      -      -      -       1.2431 -\n";
        card += "QCDscale_VV               lnN        -      -      -      1.0284  -      -\n";
        card += "BRhiggs_ZZ4l              lnN        1.02   1.02   1.02   -       -      -\n";
    if (channel == 0 || channel == 2) {
        card +=  "CMS_eff_m                 lnN        1.0224 1.0224 1.0224 1.0224  1.0224 -\n";
        card +=  "CMS_scale_m               lnN        1.01   1.01   1.01   1.01    1.01   -\n";
        card +=  "CMS_trigger_m             lnN        1.015  1.015  1.015  1.015   1.015  -\n";
    }
    if (channel == 1 || channel == 2) {
        card += "CMS_eff_e                 lnN        1.0447 1.0447 1.0447 1.0447  1.0447 -\n";
        card += "CMS_scale_e               lnN        1.02   1.02   1.02   1.02    1.02   -\n";
        card += "CMS_trigger_e             lnN        1.015  1.015  1.015  1.015   1.015  -\n";
    }
    if (channel == 0) {
        card += "CMS_hzz4mu_Zjets          gmN ZJEVT  -      -      -      -       -      ZJALPHA\n";
        card += "CMS_hzz4mu_Zjets_EXTRAP   lnN        -      -      -      -       -      ZJSYST\n";
    }
    if (channel == 1) {
        card += "CMS_hzz4e_Zjets           gmN ZJEVT  -      -      -      -       -      ZJALPHA\n";
        card += "CMS_hzz4e_Zjets_EXTRAP    lnN        -      -      -      -       -      ZJSYST\n";
    }
    if (channel == 2) {
        card += "CMS_hzz2e2mu_Zjets        gmN ZJEVT  -      -      -      -       -      ZJALPHA\n";
        card += "CMS_hzz2e2mu_Zjets_EXTRAP lnN        -      -      -      -       -      ZJSYST\n";
    }

    return card;
}


void getBkgFit(RooWorkspace& workspace, std::string varname, RooDataSet& data_bkg_red, RooRealVar& CMS_zz4l_mass, std::string outfilename, std::string pdfname) {

    RooRealVar CMS_zz4l_a1((varname+"_a1").c_str(),(varname+"_a1").c_str(),120.,100.,200.);
    RooRealVar CMS_zz4l_a2((varname+"_a2").c_str(),(varname+"_a2").c_str(),25.,10.,50.);
    RooRealVar CMS_zz4l_a3((varname+"_a3").c_str(),(varname+"_a3").c_str(),50.,20.,200.);
    RooGenericPdf background_low((varname+"_background_low").c_str(),(varname+"_background_low").c_str(),"(0.5 + 0.5 * sign(@0-@1) * TMath::Erf(TMath::Abs(@0-@1)/@2)) * exp(-1.*@0/@3)",RooArgList(CMS_zz4l_mass,CMS_zz4l_a1,CMS_zz4l_a2,CMS_zz4l_a3));

    RooRealVar CMS_zz4l_b1((varname+"_b1").c_str(),(varname+"_b1").c_str(),150.,130.,200.);
    RooRealVar CMS_zz4l_b2((varname+"_b2").c_str(),(varname+"_b2").c_str(),12.,1.,50.);
    RooRealVar CMS_zz4l_b3((varname+"_b3").c_str(),(varname+"_b3").c_str(),60.,20.,150.);
    RooGenericPdf background_base((varname+"_background_base").c_str(),(varname+"_background_base").c_str(),"(0.5 + 0.5 * sign(@0-@1) * TMath::Erf(TMath::Abs(@0-@1)/@2)) * exp(-1.*@0/@3)",RooArgList(CMS_zz4l_mass,CMS_zz4l_b1,CMS_zz4l_b2,CMS_zz4l_b3));

    RooRealVar CMS_zz4l_c1((varname+"_c1").c_str(),(varname+"_c1").c_str(),200.,170.,250.);
    RooRealVar CMS_zz4l_c2((varname+"_c2").c_str(),(varname+"_c2").c_str(),12.,1.,50.);
    RooRealVar CMS_zz4l_c3((varname+"_c3").c_str(),(varname+"_c3").c_str(),60.,20.,150.);
    RooGenericPdf background_high((varname+"_background_high").c_str(),(varname+"_background_high").c_str(),"(0.5 + 0.5 * sign(@0-@1) * TMath::Erf(TMath::Abs(@0-@1)/@2)) * exp(-1.*@0/@3)",RooArgList(CMS_zz4l_mass,CMS_zz4l_c1,CMS_zz4l_c2,CMS_zz4l_c3));

    RooRealVar CMS_zz4l_bkgfrac1((varname+"_bkgfrac1").c_str(),(varname+"_bkgfrac1").c_str(),0.182,0.,1.);
    RooAddPdf *background_nonorm1 = new RooAddPdf((varname+"_background_nonorm1").c_str(),(varname+"_background_nonorm1").c_str(),RooArgList(background_low,background_base),RooArgList(CMS_zz4l_bkgfrac1));

    RooRealVar CMS_zz4l_bkgfrac2((varname+"_bkgfrac2").c_str(),(varname+"_bkgfrac2").c_str(),0.182,0.,1.);
    RooAddPdf *background_nonorm = new RooAddPdf(pdfname.c_str(),pdfname.c_str(),RooArgList(*background_nonorm1,background_high),RooArgList(CMS_zz4l_bkgfrac2));

    RooRealVar CMS_zz4l_Nzz((varname+"_Nzz").c_str(),(varname+"_Nzz").c_str(),1.,-1.,1000.);
    CMS_zz4l_Nzz.setVal(data_bkg_red.sumEntries());
    RooExtendPdf background("background","background",*background_nonorm,CMS_zz4l_Nzz);

    background_nonorm->fitTo(data_bkg_red,RooFit::Save(1),RooFit::SumW2Error(kTRUE));

    CMS_zz4l_Nzz.setConstant(kTRUE);
    CMS_zz4l_a1.setConstant(kTRUE);
    CMS_zz4l_a2.setConstant(kTRUE);
    CMS_zz4l_a3.setConstant(kTRUE);
    CMS_zz4l_b1.setConstant(kTRUE);
    CMS_zz4l_b2.setConstant(kTRUE);
    CMS_zz4l_b3.setConstant(kTRUE);
    CMS_zz4l_bkgfrac1.setConstant(kTRUE);
    CMS_zz4l_bkgfrac2.setConstant(kTRUE);

    workspace.import(data_bkg_red);
    workspace.import(*background_nonorm);
    //workspace.import(CMS_zz4l_a1);
    //workspace.import(CMS_zz4l_a2);
    //workspace.import(CMS_zz4l_a3);
    //workspace.import(CMS_zz4l_b1);
    //workspace.import(CMS_zz4l_b2);
    //workspace.import(CMS_zz4l_b3);
    //workspace.import(CMS_zz4l_bkgfrac1);
    //workspace.import(CMS_zz4l_bkgfrac2);
    workspace.import(CMS_zz4l_Nzz);

    RooPlot *frame = CMS_zz4l_mass.frame(50);
    data_bkg_red.plotOn(frame,RooFit::DataError(RooAbsData::SumW2));
    background_nonorm->plotOn(frame);

    TCanvas canvas;
    canvas.cd();
    frame->Draw();
    canvas.SaveAs((outfilename+".gif").c_str());


}

void getZXFit(RooWorkspace& workspace, std::string varname, RooDataSet& data_bkg_red, RooRealVar& CMS_zz4l_mass, std::string outfilename, std::string pdfname) {
    RooRealVar CMS_zz4l_mean_zx((varname+"_mean_zx").c_str(),(varname+"_mean_zx").c_str(),130.,100.,200.);
    RooRealVar CMS_zz4l_sigma_zx((varname+"_sigma_zx").c_str(),(varname+"_sigma_zx").c_str(),10.,0.,100.);
    RooLandau zxLandau_nonorm(pdfname.c_str(),pdfname.c_str(),CMS_zz4l_mass,CMS_zz4l_mean_zx,CMS_zz4l_sigma_zx);

    RooRealVar CMS_zz4l_Nzx((varname+"_Nzx").c_str(),(varname+"_Nzx").c_str(),0.5,-1.,10.);
    CMS_zz4l_Nzx.setVal(data_bkg_red.sumEntries());
    RooExtendPdf zxLandau("zx","zx",zxLandau_nonorm,CMS_zz4l_Nzx);

    zxLandau_nonorm.fitTo(data_bkg_red,RooFit::Save(1),RooFit::SumW2Error(kFALSE));

    CMS_zz4l_mean_zx.setConstant(kTRUE);
    CMS_zz4l_sigma_zx.setConstant(kTRUE);
    CMS_zz4l_Nzx.setConstant(kTRUE);

    RooPlot *frame = CMS_zz4l_mass.frame(25);
    data_bkg_red.plotOn(frame,RooFit::DataError(RooAbsData::SumW2));
    zxLandau_nonorm.plotOn(frame);

    workspace.import(data_bkg_red);
    workspace.import(zxLandau_nonorm);
    //workspace.import(CMS_zz4l_mean_zx);
    //workspace.import(CMS_zz4l_sigma_zx);
    workspace.import(CMS_zz4l_Nzx);

    TCanvas c1;
    c1.cd();
    frame->Draw();
    c1.SaveAs((outfilename+".gif").c_str());

}



void getSigFit(RooWorkspace& workspace, std::string varname, RooDataSet& data_sig_red, RooRealVar& CMS_zz4l_mass, std::string outfilename, std::string pdfname, float mass, float low_M, float high_M, float gamma) {

    RooRealVar CMS_zz4l_mean_sig((varname+"_mean_sig").c_str(), (varname+"_mean_sig").c_str(),0.,-10.,10.);
    RooRealVar CMS_zz4l_sigma_sig((varname+"_sigma_sig").c_str(), (varname+"_sigma_sig").c_str(),3.,0.,30.);
    RooRealVar CMS_zz4l_alpha((varname+"_alpha").c_str(),(varname+"_alpha").c_str(),2.,0.,10.);
    RooRealVar CMS_zz4l_n((varname+"_n").c_str(),(varname+"_n").c_str(),1.,-10.,10.);
    RooCBShape signalCB((varname+"_signalCB").c_str(),(varname+"_signalCB").c_str(),CMS_zz4l_mass,CMS_zz4l_mean_sig,CMS_zz4l_sigma_sig,CMS_zz4l_alpha,CMS_zz4l_n);

    RooRealVar CMS_zz4l_mean_BW((varname+"_mean_BW").c_str(),(varname+"_mean_BW").c_str(), mass,low_M,high_M);
    RooRealVar CMS_zz4l_gamma_BW((varname+"_gamma_BW").c_str(),(varname+"_gamma_BW").c_str(),gamma);
    RooBreitWigner signalBW((varname+"_signalBW").c_str(), (varname+"_signalBW").c_str(),CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_gamma_BW);

    RooFFTConvPdf signal_nonorm(pdfname.c_str(),"BW (X) CB",CMS_zz4l_mass,signalBW,signalCB);
    signal_nonorm.setBufferFraction(0.5);

    RooRealVar CMS_zz4l_Nsig((varname+"_Nsig").c_str(),(varname+"_Nsig").c_str(),0.5,-1.,10.);
    CMS_zz4l_Nsig.setVal(data_sig_red.sumEntries());

    signal_nonorm.fitTo(data_sig_red,RooFit::Save(1),RooFit::SumW2Error(kFALSE));

    CMS_zz4l_mean_sig.setConstant(kTRUE);
    CMS_zz4l_sigma_sig.setConstant(kTRUE);
    CMS_zz4l_alpha.setConstant(kTRUE);
    CMS_zz4l_n.setConstant(kTRUE);
    CMS_zz4l_Nsig.setConstant(kTRUE);
    CMS_zz4l_mean_BW.setConstant(kTRUE);
    CMS_zz4l_gamma_BW.setConstant(kTRUE);

    workspace.import(data_sig_red);
    workspace.import(signal_nonorm);
    //workspace.import(CMS_zz4l_mean_sig);
    //workspace.import(CMS_zz4l_sigma_sig);
    //workspace.import(CMS_zz4l_alpha);
    //workspace.import(CMS_zz4l_n);
    //workspace.import(CMS_zz4l_mean_BW);
    //workspace.import(CMS_zz4l_gamma_BW);
    workspace.import(CMS_zz4l_Nsig);


    RooPlot *frame = CMS_zz4l_mass.frame(low_M,high_M,30);
    data_sig_red.plotOn(frame,RooFit::DataError(RooAbsData::SumW2));
    signal_nonorm.plotOn(frame);

    TCanvas c1;
    c1.cd();
    frame->Draw();
    c1.SaveAs((outfilename+".gif").c_str());

}


void doHZZAnalysis() {

    std::string base_folder = "/home/avartak/CMS/Higgs/CMSSW_4_2_8_patch7/src/WWAnalysis/AnalysisStep/trees/conf1/";
    std::string card_name   = "hzz_m120";
    float min4l = 100;
    float max4l = 600;
    float z1min = 40;
    float z2min = 12;
    float higgsMass = 120.;
    float higgsMassLow = 100.;
    float higgsMassHigh = 140.;
    float gamma = 3.47*0.001;
    bool doShapeAnalysis = true;
    std::string workspace_2e2mu = "hzz2e2mu_workspace.root";
    std::string workspace_4e    = "hzz4e_workspace.root";
    std::string workspace_4mu   = "hzz4mu_workspace.root";

    init();

    domufakes(base_folder+"hzzTree.root", mubarrel, muendcap, mubarrelerr, muendcaperr);
    doelfakes(base_folder+"hzzTree.root", elbarrel, elendcap, elbarrelerr, elendcaperr);

    float yield_qq_mm = 0.0;
    float yield_qq_ee = 0.0;
    float yield_qq_em = 0.0;

    float yield_gg_mm = 0.0;
    float yield_gg_ee = 0.0;
    float yield_gg_em = 0.0;

    float yield_ggh_mm= 0.0;
    float yield_ggh_ee= 0.0;
    float yield_ggh_em= 0.0;

    float yield_vbh_mm= 0.0;
    float yield_vbh_ee= 0.0;
    float yield_vbh_em= 0.0;

    float yield_wzh_mm= 0.0;
    float yield_wzh_ee= 0.0;
    float yield_wzh_em= 0.0;

    float yield_zj_mm = 0.0;
    float yield_zj_ee = 0.0;
    float yield_zj_em = 0.0;

    float yield_zj_mm_count = 0.0;
    float yield_zj_ee_count = 0.0;
    float yield_zj_em_count = 0.0;

    float yield_zj_mm_error = 0.0;
    float yield_zj_ee_error = 0.0;
    float yield_zj_em_error = 0.0;

    RooRealVar mass_2e2mu("mass", "mass", 100, 100, 600);
    RooRealVar mass_4e   ("mass", "mass", 100, 100, 600);
    RooRealVar mass_4mu  ("mass", "mass", 100, 100, 600);
    
    RooRealVar mass_gz2e2mu("mass", "mass", 100, 100, 600);
    RooRealVar mass_gz4e   ("mass", "mass", 100, 100, 600);
    RooRealVar mass_gz4mu  ("mass", "mass", 100, 100, 600);
    
    RooRealVar weight_2e2mu("weight", "weight", 0., 0., 10.);
    RooRealVar weight_4e   ("weight", "weight", 0., 0., 10.);
    RooRealVar weight_4mu  ("weight", "weight", 0., 0., 10.);

    RooRealVar weight_gz2e2mu("weight", "weight", 0., 0., 10.);
    RooRealVar weight_gz4e   ("weight", "weight", 0., 0., 10.);
    RooRealVar weight_gz4mu  ("weight", "weight", 0., 0., 10.);

    RooArgSet argset_2e2mu(mass_2e2mu, weight_2e2mu, "argset_2e2mu");
    RooArgSet argset_4e   (mass_4e,    weight_4e,    "argset_4e");
    RooArgSet argset_4mu  (mass_4mu,   weight_4mu,   "argset_4mu");

    RooArgSet argset_gz2e2mu(mass_gz2e2mu, weight_gz2e2mu, "argset_gz2e2mu");
    RooArgSet argset_gz4e   (mass_gz4e,    weight_gz4e,    "argset_gz4e");
    RooArgSet argset_gz4mu  (mass_gz4mu,   weight_gz4mu,   "argset_gz4mu");

    RooDataSet dataset_2e2mu("dataset_2e2mu", "", argset_2e2mu, RooFit::WeightVar("weight"));
    RooDataSet dataset_4e   ("dataset_4e"   , "", argset_4e,    RooFit::WeightVar("weight"));
    RooDataSet dataset_4mu  ("dataset_4mu"  , "", argset_4mu,   RooFit::WeightVar("weight"));

    RooDataSet dataset_gz2e2mu("dataset_gz2e2mu", "", argset_gz2e2mu, RooFit::WeightVar("weight"));
    RooDataSet dataset_gz4e   ("dataset_gz4e"   , "", argset_gz4e,    RooFit::WeightVar("weight"));
    RooDataSet dataset_gz4mu  ("dataset_gz4mu"  , "", argset_gz4mu,   RooFit::WeightVar("weight"));

    RooRealVar mass_g2e2mu("mass", "mass", higgsMassLow, higgsMassHigh);
    RooRealVar mass_g4e   ("mass", "mass", higgsMassLow, higgsMassHigh);
    RooRealVar mass_g4mu  ("mass", "mass", higgsMassLow, higgsMassHigh);

    RooRealVar mass_q2e2mu("mass", "mass", higgsMassLow, higgsMassHigh);
    RooRealVar mass_q4e   ("mass", "mass", higgsMassLow, higgsMassHigh);
    RooRealVar mass_q4mu  ("mass", "mass", higgsMassLow, higgsMassHigh);

    RooRealVar mass_v2e2mu("mass", "mass", higgsMassLow, higgsMassHigh);
    RooRealVar mass_v4e   ("mass", "mass", higgsMassLow, higgsMassHigh);
    RooRealVar mass_v4mu  ("mass", "mass", higgsMassLow, higgsMassHigh);

    mass_g2e2mu.setBins(100000,"fft");
    mass_g4e.setBins   (100000,"fft");
    mass_g4mu.setBins  (100000,"fft");

    mass_q2e2mu.setBins(100000,"fft");
    mass_q4e.setBins   (100000,"fft");
    mass_q4mu.setBins  (100000,"fft");

    mass_v2e2mu.setBins(100000,"fft");
    mass_v4e.setBins   (100000,"fft");
    mass_v4mu.setBins  (100000,"fft");

    RooRealVar weight_g2e2mu("weight", "weight", 0., 0., 10.);
    RooRealVar weight_g4e   ("weight", "weight", 0., 0., 10.);
    RooRealVar weight_g4mu  ("weight", "weight", 0., 0., 10.);

    RooRealVar weight_q2e2mu("weight", "weight", 0., 0., 10.);
    RooRealVar weight_q4e   ("weight", "weight", 0., 0., 10.);
    RooRealVar weight_q4mu  ("weight", "weight", 0., 0., 10.);

    RooRealVar weight_v2e2mu("weight", "weight", 0., 0., 10.);
    RooRealVar weight_v4e   ("weight", "weight", 0., 0., 10.);
    RooRealVar weight_v4mu  ("weight", "weight", 0., 0., 10.);

    RooArgSet argset_g2e2mu(mass_g2e2mu, weight_g2e2mu, "argset_g2e2mu");
    RooArgSet argset_g4e   (mass_g4e,    weight_g4e,    "argset_g4e");
    RooArgSet argset_g4mu  (mass_g4mu,   weight_g4mu,   "argset_g4mu");

    RooArgSet argset_q2e2mu(mass_q2e2mu, weight_q2e2mu, "argset_q2e2mu");
    RooArgSet argset_q4e   (mass_q4e,    weight_q4e,    "argset_q4e");
    RooArgSet argset_q4mu  (mass_q4mu,   weight_q4mu,   "argset_q4mu");

    RooArgSet argset_v2e2mu(mass_v2e2mu, weight_v2e2mu, "argset_v2e2mu");
    RooArgSet argset_v4e   (mass_v4e,    weight_v4e,    "argset_v4e");
    RooArgSet argset_v4mu  (mass_v4mu,   weight_v4mu,   "argset_v4mu");

    RooDataSet dataset_g2e2mu("dataset_g2e2mu", "", argset_g2e2mu, RooFit::WeightVar("weight"));
    RooDataSet dataset_g4e   ("dataset_g4e"   , "", argset_g4e,    RooFit::WeightVar("weight"));
    RooDataSet dataset_g4mu  ("dataset_g4mu"  , "", argset_g4mu,   RooFit::WeightVar("weight"));

    RooDataSet dataset_q2e2mu("dataset_q2e2mu", "", argset_q2e2mu, RooFit::WeightVar("weight"));
    RooDataSet dataset_q4e   ("dataset_q4e"   , "", argset_q4e,    RooFit::WeightVar("weight"));
    RooDataSet dataset_q4mu  ("dataset_q4mu"  , "", argset_q4mu,   RooFit::WeightVar("weight"));

    RooDataSet dataset_v2e2mu("dataset_v2e2mu", "", argset_v2e2mu, RooFit::WeightVar("weight"));
    RooDataSet dataset_v4e   ("dataset_v4e"   , "", argset_v4e,    RooFit::WeightVar("weight"));
    RooDataSet dataset_v4mu  ("dataset_v4mu"  , "", argset_v4mu,   RooFit::WeightVar("weight"));

    RooRealVar mass_zx2e2mu("mass", "mass", 100, 600);
    RooRealVar mass_zx4e   ("mass", "mass", 100, 600);
    RooRealVar mass_zx4mu  ("mass", "mass", 100, 600);
    
    RooRealVar weight_zx2e2mu("weight", "weight", 0., 0., 10.);
    RooRealVar weight_zx4e   ("weight", "weight", 0., 0., 10.);
    RooRealVar weight_zx4mu  ("weight", "weight", 0., 0., 10.);

    RooArgSet argset_zx2e2mu(mass_zx2e2mu, weight_zx2e2mu, "argset_zx2e2mu");
    RooArgSet argset_zx4e   (mass_zx4e,    weight_zx4e,    "argset_zx4e");
    RooArgSet argset_zx4mu  (mass_zx4mu,   weight_zx4mu,   "argset_zx4mu");

    RooDataSet dataset_zx2e2mu("dataset_zx2e2mu", "", argset_zx2e2mu, RooFit::WeightVar("weight"));
    RooDataSet dataset_zx4e   ("dataset_zx4e"   , "", argset_zx4e,    RooFit::WeightVar("weight"));
    RooDataSet dataset_zx4mu  ("dataset_zx4mu"  , "", argset_zx4mu,   RooFit::WeightVar("weight"));


    yield_zj_mm_count += getYield((base_folder+"hzzTree.root"),     0, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx4mu  , argset_zx4mu,   4)[1];
    yield_zj_ee_count += getYield((base_folder+"hzzTree.root"),     1, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx4e   , argset_zx4e,    4)[1];
    yield_zj_em_count += getYield((base_folder+"hzzTree.root"),     2, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx2e2mu, argset_zx2e2mu, 4)[1];
    yield_zj_em_count += getYield((base_folder+"hzzTree.root"),     3, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx2e2mu, argset_zx2e2mu, 4)[1];

    yield_zj_mm_error += getYield((base_folder+"hzzTree.root"),     0, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx4mu  , argset_zx4mu,   4)[2];
    yield_zj_ee_error += getYield((base_folder+"hzzTree.root"),     1, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx4e   , argset_zx4e,    4)[2];
    yield_zj_em_error += getYield((base_folder+"hzzTree.root"),     2, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx2e2mu, argset_zx2e2mu, 4)[2];
    yield_zj_em_error += getYield((base_folder+"hzzTree.root"),     3, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx2e2mu, argset_zx2e2mu, 4)[2];

    yield_zj_mm       += getYield(base_folder+"hzzTree.root",       0, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx4mu  , argset_zx4mu  , 0)[0];
    yield_zj_ee       += getYield(base_folder+"hzzTree.root",       1, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx4e   , argset_zx4e   , 0)[0];
    yield_zj_em       += getYield(base_folder+"hzzTree.root",       2, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx2e2mu, argset_zx2e2mu, 0)[0];
    yield_zj_em       += getYield(base_folder+"hzzTree.root",       3, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx2e2mu, argset_zx2e2mu, 0)[0];

    yield_ggh_mm      += getYield(base_folder+"hzzTree_id108.root", 0, (xsec_gghzz120 /evt_gghzz120 )    *5.0, 1, z1min, z2min, min4l, max4l, dataset_g4mu   , argset_g4mu   , 0)[0];
    yield_vbh_mm      += getYield(base_folder+"hzzTree_id109.root", 0, (xsec_vbfhzz120/evt_vbfhzz120)    *5.0, 1, z1min, z2min, min4l, max4l, dataset_q4mu   , argset_q4mu   , 0)[0];
    yield_wzh_mm      += getYield(base_folder+"hzzTree_id110.root", 0, (xsec_vhzz120  /evt_vhzz120  )    *5.0, 1, z1min, z2min, min4l, max4l, dataset_v4mu   , argset_v4mu   , 0)[0];

    yield_ggh_ee      += getYield(base_folder+"hzzTree_id108.root", 1, (xsec_gghzz120 /evt_gghzz120 )    *5.0, 1, z1min, z2min, min4l, max4l, dataset_g4e    , argset_g4e    , 0)[0];
    yield_vbh_ee      += getYield(base_folder+"hzzTree_id109.root", 1, (xsec_vbfhzz120/evt_vbfhzz120)    *5.0, 1, z1min, z2min, min4l, max4l, dataset_q4e    , argset_q4e    , 0)[0];
    yield_wzh_ee      += getYield(base_folder+"hzzTree_id110.root", 1, (xsec_vhzz120  /evt_vhzz120  )    *5.0, 1, z1min, z2min, min4l, max4l, dataset_v4e    , argset_v4e    , 0)[0];

    yield_ggh_em      += getYield(base_folder+"hzzTree_id108.root", 2, (xsec_gghzz120 /evt_gghzz120 )    *5.0, 1, z1min, z2min, min4l, max4l, dataset_g2e2mu , argset_g2e2mu , 0)[0];
    yield_vbh_em      += getYield(base_folder+"hzzTree_id109.root", 2, (xsec_vbfhzz120/evt_vbfhzz120)    *5.0, 1, z1min, z2min, min4l, max4l, dataset_q2e2mu , argset_q2e2mu , 0)[0];
    yield_wzh_em      += getYield(base_folder+"hzzTree_id110.root", 2, (xsec_vhzz120  /evt_vhzz120  )    *5.0, 1, z1min, z2min, min4l, max4l, dataset_v2e2mu , argset_v2e2mu , 0)[0];

    yield_ggh_em      += getYield(base_folder+"hzzTree_id108.root", 3, (xsec_gghzz120 /evt_gghzz120 )    *5.0, 1, z1min, z2min, min4l, max4l, dataset_g2e2mu , argset_g2e2mu , 0)[0];
    yield_vbh_em      += getYield(base_folder+"hzzTree_id109.root", 3, (xsec_vbfhzz120/evt_vbfhzz120)    *5.0, 1, z1min, z2min, min4l, max4l, dataset_q2e2mu , argset_q2e2mu , 0)[0];
    yield_wzh_em      += getYield(base_folder+"hzzTree_id110.root", 3, (xsec_vhzz120  /evt_vhzz120  )    *5.0, 1, z1min, z2min, min4l, max4l, dataset_v2e2mu , argset_v2e2mu , 0)[0];

    yield_qq_mm       += getYield(base_folder+"hzzTree_id103.root", 0, (xsec_qqzz4mu    /evt_qqzz4mu)    *5.0, 1, z1min, z2min, min4l, max4l, dataset_4mu    , argset_4mu    , 0)[0];
    yield_qq_mm       += getYield(base_folder+"hzzTree_id102.root", 0, (xsec_qqzz4e     /evt_qqzz4e)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_4mu    , argset_4mu    , 0)[0];
    yield_qq_mm       += getYield(base_folder+"hzzTree_id104.root", 0, (xsec_qqzz4tau   /evt_qqzz4tau)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_4mu    , argset_4mu    , 0)[0];
    yield_qq_mm       += getYield(base_folder+"hzzTree_id105.root", 0, (xsec_qqzz2e2mu  /evt_qqzz2e2mu)  *5.0, 1, z1min, z2min, min4l, max4l, dataset_4mu    , argset_4mu    , 0)[0];
    yield_qq_mm       += getYield(base_folder+"hzzTree_id106.root", 0, (xsec_qqzz2e2tau /evt_qqzz2e2tau) *5.0, 1, z1min, z2min, min4l, max4l, dataset_4mu    , argset_4mu    , 0)[0];
    yield_qq_mm       += getYield(base_folder+"hzzTree_id107.root", 0, (xsec_qqzz2mu2tau/evt_qqzz2mu2tau)*5.0, 1, z1min, z2min, min4l, max4l, dataset_4mu    , argset_4mu    , 0)[0];
    yield_gg_mm       += getYield(base_folder+"hzzTree_id100.root", 0, (xsec_ggzz2l2l   /evt_ggzz2l2l)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz4mu  , argset_gz4mu  , 0)[0];
    yield_gg_mm       += getYield(base_folder+"hzzTree_id101.root", 0, (xsec_ggzz4l     /evt_ggzz4l)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz4mu  , argset_gz4mu  , 0)[0];
    
    yield_qq_ee       += getYield(base_folder+"hzzTree_id103.root", 1, (xsec_qqzz4mu    /evt_qqzz4mu)    *5.0, 1, z1min, z2min, min4l, max4l, dataset_4e     , argset_4e     , 0)[0];
    yield_qq_ee       += getYield(base_folder+"hzzTree_id102.root", 1, (xsec_qqzz4e     /evt_qqzz4e)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_4e     , argset_4e     , 0)[0];
    yield_qq_ee       += getYield(base_folder+"hzzTree_id104.root", 1, (xsec_qqzz4tau   /evt_qqzz4tau)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_4e     , argset_4e     , 0)[0];
    yield_qq_ee       += getYield(base_folder+"hzzTree_id105.root", 1, (xsec_qqzz2e2mu  /evt_qqzz2e2mu)  *5.0, 1, z1min, z2min, min4l, max4l, dataset_4e     , argset_4e     , 0)[0];
    yield_qq_ee       += getYield(base_folder+"hzzTree_id106.root", 1, (xsec_qqzz2e2tau /evt_qqzz2e2tau) *5.0, 1, z1min, z2min, min4l, max4l, dataset_4e     , argset_4e     , 0)[0];
    yield_qq_ee       += getYield(base_folder+"hzzTree_id107.root", 1, (xsec_qqzz2mu2tau/evt_qqzz2mu2tau)*5.0, 1, z1min, z2min, min4l, max4l, dataset_4e     , argset_4e     , 0)[0];
    yield_gg_ee       += getYield(base_folder+"hzzTree_id100.root", 1, (xsec_ggzz2l2l   /evt_ggzz2l2l)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz4e   , argset_gz4e   , 0)[0];
    yield_gg_ee       += getYield(base_folder+"hzzTree_id101.root", 1, (xsec_ggzz4l     /evt_ggzz4l)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz4e   , argset_gz4e   , 0)[0];

    yield_qq_em       += getYield(base_folder+"hzzTree_id103.root", 2, (xsec_qqzz4mu    /evt_qqzz4mu)    *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset_2e2mu  , 0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id102.root", 2, (xsec_qqzz4e     /evt_qqzz4e)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset_2e2mu  , 0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id104.root", 2, (xsec_qqzz4tau   /evt_qqzz4tau)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset_2e2mu  , 0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id105.root", 2, (xsec_qqzz2e2mu  /evt_qqzz2e2mu)  *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset_2e2mu  , 0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id106.root", 2, (xsec_qqzz2e2tau /evt_qqzz2e2tau) *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset_2e2mu  , 0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id107.root", 2, (xsec_qqzz2mu2tau/evt_qqzz2mu2tau)*5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset_2e2mu  , 0)[0];
    yield_gg_em       += getYield(base_folder+"hzzTree_id100.root", 2, (xsec_ggzz2l2l   /evt_ggzz2l2l)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz2e2mu, argset_gz2e2mu, 0)[0];
    yield_gg_em       += getYield(base_folder+"hzzTree_id101.root", 2, (xsec_ggzz4l     /evt_ggzz4l)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz2e2mu, argset_gz2e2mu, 0)[0];

    yield_qq_em       += getYield(base_folder+"hzzTree_id103.root", 3, (xsec_qqzz4mu    /evt_qqzz4mu)    *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset_2e2mu  , 0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id102.root", 3, (xsec_qqzz4e     /evt_qqzz4e)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset_2e2mu  , 0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id104.root", 3, (xsec_qqzz4tau   /evt_qqzz4tau)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset_2e2mu  , 0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id105.root", 3, (xsec_qqzz2e2mu  /evt_qqzz2e2mu)  *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset_2e2mu  , 0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id106.root", 3, (xsec_qqzz2e2tau /evt_qqzz2e2tau) *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset_2e2mu  , 0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id107.root", 3, (xsec_qqzz2mu2tau/evt_qqzz2mu2tau)*5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset_2e2mu  , 0)[0];
    yield_gg_em       += getYield(base_folder+"hzzTree_id100.root", 3, (xsec_ggzz2l2l   /evt_ggzz2l2l)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz2e2mu, argset_gz2e2mu, 0)[0];
    yield_gg_em       += getYield(base_folder+"hzzTree_id101.root", 3, (xsec_ggzz4l     /evt_ggzz4l)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz2e2mu, argset_gz2e2mu, 0)[0];

    
    std::string card_4mu   = doShapeAnalysis == false? createCardTemplate(0) : createShapeCardTemplate(0, workspace_4mu.c_str());
    std::string card_4e    = doShapeAnalysis == false? createCardTemplate(1) : createShapeCardTemplate(1, workspace_4e.c_str());
    std::string card_2e2mu = doShapeAnalysis == false? createCardTemplate(2) : createShapeCardTemplate(2, workspace_2e2mu.c_str());

    card_2e2mu = findAndReplace(card_2e2mu, "SIG_GGH",  yield_ggh_em);
    card_4e    = findAndReplace(card_4e,    "SIG_GGH",  yield_ggh_ee);
    card_4mu   = findAndReplace(card_4mu  , "SIG_GGH",  yield_ggh_mm);
     
    card_2e2mu = findAndReplace(card_2e2mu, "SIG_VBF",  yield_vbh_em);
    card_4e    = findAndReplace(card_4e,    "SIG_VBF",  yield_vbh_ee);
    card_4mu   = findAndReplace(card_4mu,   "SIG_VBF",  yield_vbh_mm);
     
    card_2e2mu = findAndReplace(card_2e2mu, "SIG_WH",   yield_wzh_em);
    card_4e    = findAndReplace(card_4e,    "SIG_WH",   yield_wzh_ee);
    card_4mu   = findAndReplace(card_4mu,   "SIG_WH",   yield_wzh_mm);
     
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
     
    card_2e2mu = findAndReplace(card_2e2mu, "ZJSYST",   yield_zj_em_error / yield_zj_em);
    card_4e    = findAndReplace(card_4e,    "ZJSYST",   yield_zj_ee_error / yield_zj_ee);
    card_4mu   = findAndReplace(card_4mu,   "ZJSYST",   yield_zj_mm_error / yield_zj_mm);

    card_2e2mu = findAndReplace(card_2e2mu, "BIN",   "a3");
    card_4e    = findAndReplace(card_4e,    "BIN",   "a2");
    card_4mu   = findAndReplace(card_4mu,   "BIN",   "a1");

    if (doShapeAnalysis) {
        card_2e2mu = findAndReplace(card_2e2mu, "OBS",      yield_qq_em);
        card_4e    = findAndReplace(card_4e,    "OBS",      yield_qq_ee);
        card_4mu   = findAndReplace(card_4mu,   "OBS",      yield_qq_mm);
    }     
    else {
        card_2e2mu = findAndReplace(card_2e2mu, "OBS",      yield_zj_em + yield_gg_em + yield_qq_em);
        card_4e    = findAndReplace(card_4e,    "OBS",      yield_zj_ee + yield_gg_ee + yield_qq_ee);
        card_4mu   = findAndReplace(card_4mu,   "OBS",      yield_zj_mm + yield_gg_mm + yield_qq_mm);
    }
     
    RooWorkspace w_2e2mu("w_2e2mu", "");
    RooWorkspace w_4e("w_4e", "");
    RooWorkspace w_4mu("w_4mu", "");

    getSigFit(w_2e2mu, "CMS_gh2e2mu", dataset_g2e2mu  , mass_g2e2mu  , "fit_g2e2mu" , "sig_ggH", higgsMass, higgsMassLow, higgsMassHigh, gamma);
    getSigFit(w_4e,    "CMS_gh4e",    dataset_g4e     , mass_g4e     , "fit_g4e"    , "sig_ggH", higgsMass, higgsMassLow, higgsMassHigh, gamma);
    getSigFit(w_4mu,   "CMS_gh4mu",   dataset_g4mu    , mass_g4mu    , "fit_g4mu"   , "sig_ggH", higgsMass, higgsMassLow, higgsMassHigh, gamma);
    
    getSigFit(w_2e2mu, "CMS_vb2e2mu", dataset_q2e2mu  , mass_q2e2mu  , "fit_q2e2mu" , "sig_VBF", higgsMass, higgsMassLow, higgsMassHigh, gamma);
    getSigFit(w_4e,    "CMS_vb4e",    dataset_q4e     , mass_q4e     , "fit_q4e"    , "sig_VBF", higgsMass, higgsMassLow, higgsMassHigh, gamma);
    getSigFit(w_4mu,   "CMS_vb4mu",   dataset_q4mu    , mass_q4mu    , "fit_q4mu"   , "sig_VBF", higgsMass, higgsMassLow, higgsMassHigh, gamma);
    
    getSigFit(w_2e2mu, "CMS_vh2e2mu", dataset_v2e2mu  , mass_v2e2mu  , "fit_v2e2mu" , "sig_WH",  higgsMass, higgsMassLow, higgsMassHigh, gamma);
    getSigFit(w_4e,    "CMS_vh4e",    dataset_v4e     , mass_v4e     , "fit_v4e"    , "sig_WH",  higgsMass, higgsMassLow, higgsMassHigh, gamma);
    getSigFit(w_4mu,   "CMS_vh4mu",   dataset_v4mu    , mass_v4mu    , "fit_v4mu"   , "sig_WH",  higgsMass, higgsMassLow, higgsMassHigh, gamma);

    getZXFit (w_2e2mu, "CMS_zx2e2mu", dataset_zx2e2mu , mass_zx2e2mu , "fit_zx2e2mu", "bkg_zjets");
    getZXFit (w_4e,    "CMS_zx4e",    dataset_zx4e    , mass_zx4e    , "fit_zx4e"   , "bkg_zjets");
    getZXFit (w_4mu,   "CMS_zx4mu",   dataset_zx4mu   , mass_zx4mu   , "fit_zx4mu"  , "bkg_zjets");

    getBkgFit(w_2e2mu, "CMS_gz2e2mu", dataset_gz2e2mu , mass_gz2e2mu , "fit_gz2e2mu", "bkg_ggzz");
    getBkgFit(w_4e,    "CMS_gz4e",    dataset_gz4e    , mass_gz4e    , "fit_gz4e"   , "bkg_ggzz");
    getBkgFit(w_4mu,   "CMS_gz4mu",   dataset_gz4mu   , mass_gz4mu   , "fit_gz4mu"  , "bkg_ggzz");
    
    getBkgFit(w_2e2mu, "CMS_zz2e2mu", dataset_2e2mu   , mass_2e2mu   , "fit_2e2mu"  , "bkg_qqzz");
    getBkgFit(w_4e,    "CMS_zz4e",    dataset_4e      , mass_4e      , "fit_4e"     , "bkg_qqzz");
    getBkgFit(w_4mu,   "CMS_zz4mu",   dataset_4mu     , mass_4mu     , "fit_4mu"    , "bkg_qqzz");
    

    w_2e2mu.import(*((RooDataSet*)dataset_2e2mu.Clone("data_obs")));
    w_4mu.import  (*((RooDataSet*)dataset_4mu.Clone  ("data_obs")));
    w_4e.import   (*((RooDataSet*)dataset_4e.Clone   ("data_obs")));

    w_2e2mu.writeToFile(workspace_2e2mu.c_str());
    w_4e.writeToFile(workspace_4e.c_str());
    w_4mu.writeToFile(workspace_4mu.c_str());

    std::cout << "---- ZX Yields ----" << std::endl;
    std::cout << "4e     : " << yield_zj_ee << " +/- " << yield_zj_ee_error << " (" << yield_zj_ee_count << " events in CR)" << std::endl;
    std::cout << "4mu    : " << yield_zj_mm << " +/- " << yield_zj_mm_error << " (" << yield_zj_mm_count << " events in CR)" << std::endl;
    std::cout << "2e2mu  : " << yield_zj_em << " +/- " << yield_zj_em_error << " (" << yield_zj_em_count << " events in CR)" << std::endl;

    std::cout << std::endl;

    std::cout << "---- ZZ Yields ----" << std::endl;

    std::cout << "qq->ZZ->4e     : " << yield_qq_ee << std::endl;
    std::cout << "qq->ZZ->4mu    : " << yield_qq_mm << std::endl;
    std::cout << "qq->ZZ->2e2mu  : " << yield_qq_em << std::endl;

    std::cout << std::endl;

    std::cout << "gg->ZZ->4e     : " << yield_gg_ee << std::endl;
    std::cout << "gg->ZZ->4mu    : " << yield_gg_mm << std::endl;
    std::cout << "gg->ZZ->2e2mu  : " << yield_gg_em << std::endl;

    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "---- Higgs Yields ----" << std::endl;

    std::cout << "H->ZZ->4e      : " << yield_ggh_ee << std::endl;
    std::cout << "H->ZZ->4mu     : " << yield_ggh_mm << std::endl;
    std::cout << "H->ZZ->2e2mu   : " << yield_ggh_em << std::endl;

    std::cout << std::endl;

    std::cout << "qqH->ZZ->4e    : " << yield_vbh_ee << std::endl;
    std::cout << "qqH->ZZ->4mu   : " << yield_vbh_mm << std::endl;
    std::cout << "qqH->ZZ->2e2mu : " << yield_vbh_em << std::endl;

    std::cout << std::endl;

    std::cout << "VH->ZZ->4e     : " << yield_wzh_ee << std::endl;
    std::cout << "VH->ZZ->4mu    : " << yield_wzh_mm << std::endl;
    std::cout << "VH->ZZ->2e2mu  : " << yield_wzh_em << std::endl;

    std::cout << std::endl;


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


}
