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
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"

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

using namespace RooFit;


float mubarrel[10]    = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float muendcap[10]    = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float elbarrel[10]    = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float elendcap[10]    = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

float mubarrelerr[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float muendcaperr[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float elbarrelerr[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float elendcaperr[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

void domufakes(std::string path, float* barrel, float* endcap, float* barrelerr, float* endcaperr, std::ofstream& outfile) {
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

    outfile << "---- Computing fake rates for muons ----" << std::endl;
    outfile << "---- Muon Barrel ----" << std::endl;
    outfile << "pT [ 0 -  5] : " << mubarrel[0] << " +/- " << mubarrelerr[0] << endl;
    outfile << "pT [ 7 -  7] : " << mubarrel[1] << " +/- " << mubarrelerr[1] << endl;
    outfile << "pT [ 7 - 10] : " << mubarrel[2] << " +/- " << mubarrelerr[2] << endl;
    outfile << "pT [10 - 15] : " << mubarrel[3] << " +/- " << mubarrelerr[3] << endl;
    outfile << "pT [15 - 20] : " << mubarrel[4] << " +/- " << mubarrelerr[4] << endl;
    outfile << "pT [20 - 25] : " << mubarrel[5] << " +/- " << mubarrelerr[5] << endl;
    outfile << "pT [25 - 30] : " << mubarrel[6] << " +/- " << mubarrelerr[6] << endl;
    outfile << "pT [30 - 40] : " << mubarrel[7] << " +/- " << mubarrelerr[7] << endl;
    outfile << "pT [40 - 50] : " << mubarrel[8] << " +/- " << mubarrelerr[8] << endl;
    outfile << "pT [50 - 80] : " << mubarrel[9] << " +/- " << mubarrelerr[9] << endl;
    outfile << std::endl;
    outfile << "---- Muon Endcap ----" << std::endl;
    outfile << "pT [ 0 -  5] : " << muendcap[0] << " +/- " << muendcaperr[0] << endl;
    outfile << "pT [ 7 -  7] : " << muendcap[1] << " +/- " << muendcaperr[1] << endl;
    outfile << "pT [ 7 - 10] : " << muendcap[2] << " +/- " << muendcaperr[2] << endl;
    outfile << "pT [10 - 15] : " << muendcap[3] << " +/- " << muendcaperr[3] << endl;
    outfile << "pT [15 - 20] : " << muendcap[4] << " +/- " << muendcaperr[4] << endl;
    outfile << "pT [20 - 25] : " << muendcap[5] << " +/- " << muendcaperr[5] << endl;
    outfile << "pT [25 - 30] : " << muendcap[6] << " +/- " << muendcaperr[6] << endl;
    outfile << "pT [30 - 40] : " << muendcap[7] << " +/- " << muendcaperr[7] << endl;
    outfile << "pT [40 - 50] : " << muendcap[8] << " +/- " << muendcaperr[8] << endl;
    outfile << "pT [50 - 80] : " << muendcap[9] << " +/- " << muendcaperr[9] << endl;
    outfile << std::endl;

    file->Close();
    delete file;
}


void doelfakes(std::string path, float* barrel, float* endcap, float* barrelerr, float*endcaperr, std::ofstream& outfile) {
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

    outfile << "---- Computing fake rates for electrons ----" << std::endl;
    outfile << "---- Electron Barrel ----" << std::endl;
    outfile << "pT [ 0 -  5] : " << elbarrel[0] << " +/- " << elbarrelerr[0] << endl;
    outfile << "pT [ 7 -  7] : " << elbarrel[1] << " +/- " << elbarrelerr[1] << endl;
    outfile << "pT [ 7 - 10] : " << elbarrel[2] << " +/- " << elbarrelerr[2] << endl;
    outfile << "pT [10 - 15] : " << elbarrel[3] << " +/- " << elbarrelerr[3] << endl;
    outfile << "pT [15 - 20] : " << elbarrel[4] << " +/- " << elbarrelerr[4] << endl;
    outfile << "pT [20 - 25] : " << elbarrel[5] << " +/- " << elbarrelerr[5] << endl;
    outfile << "pT [25 - 30] : " << elbarrel[6] << " +/- " << elbarrelerr[6] << endl;
    outfile << "pT [30 - 40] : " << elbarrel[7] << " +/- " << elbarrelerr[7] << endl;
    outfile << "pT [40 - 50] : " << elbarrel[8] << " +/- " << elbarrelerr[8] << endl;
    outfile << "pT [50 - 80] : " << elbarrel[9] << " +/- " << elbarrelerr[9] << endl;
    outfile << std::endl;
    outfile << "---- Electron Endcap ----" << std::endl;
    outfile << "pT [ 0 -  5] : " << elendcap[0] << " +/- " << elendcaperr[0] << endl;
    outfile << "pT [ 7 -  7] : " << elendcap[1] << " +/- " << elendcaperr[1] << endl;
    outfile << "pT [ 7 - 10] : " << elendcap[2] << " +/- " << elendcaperr[2] << endl;
    outfile << "pT [10 - 15] : " << elendcap[3] << " +/- " << elendcaperr[3] << endl;
    outfile << "pT [15 - 20] : " << elendcap[4] << " +/- " << elendcaperr[4] << endl;
    outfile << "pT [20 - 25] : " << elendcap[5] << " +/- " << elendcaperr[5] << endl;
    outfile << "pT [25 - 30] : " << elendcap[6] << " +/- " << elendcaperr[6] << endl;
    outfile << "pT [30 - 40] : " << elendcap[7] << " +/- " << elendcaperr[7] << endl;
    outfile << "pT [40 - 50] : " << elendcap[8] << " +/- " << elendcaperr[8] << endl;
    outfile << "pT [50 - 80] : " << elendcap[9] << " +/- " << elendcaperr[9] << endl;
    outfile << std::endl;

   
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

        if (mass>100. && mass<600. && channel == (int)ch && z1mass>z1min && z2mass>z2min) {
    
            if (isMC) {
                if (mass>min4l && mass<max4l) {
                    yield += scale * getPUWeight((int)numsim);
                    yield_err += (scale * getPUWeight((int)numsim))*(scale * getPUWeight((int)numsim));
                    yield_count += 1.0;
                }
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
                        if (mass>min4l && mass<max4l) {
                            yield += (f1/(1-f1))*(f2/(1-f2));
                            yield_up += (f1_up/(1-f1_up))*(f2_up/(1-f2_up));
                            yield_dn += (f1_dn/(1-f1_dn))*(f2_dn/(1-f2_dn));
                            yield_count += 1.0;
                        }
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

std::string createCardTemplate(int channel, bool doShapeAnalysis, std::string outfilename) {

    std::string card = "";
        card += "imax 1\n";
        card += "jmax *\n";
        card += "kmax *\n";
        card += "------------\n";
    if (doShapeAnalysis && channel == 0) {
        card += std::string("shapes * * ") + outfilename + " w_4mu:$PROCESS\n";
    }
    if (doShapeAnalysis && channel == 1) {
        card += std::string("shapes * * ") + outfilename + " w_4e:$PROCESS\n";
    }
    if (doShapeAnalysis && channel == 2) {
        card += std::string("shapes * * ") + outfilename + " w_2e2mu:$PROCESS\n";
    }
        card += "------------\n";
        card += "bin         BIN\n";
        card += "observation OBS\n";
        card += "------------\n";
        card += "bin     BIN      BIN        BIN       BIN       BIN\n";
        card += "process sig_ggH sig_VBF  bkg_qqzz bkg_ggzz bkg_zjets\n";
        card += "process -2     -1        1        2        3\n";
        card += "rate    SIG_GGH SIG_VBF  BKG_QQZZ BKG_GGZZ BKG_ZJETS\n";
        card += "------------\n";
        card += "lumi                      lnN        1.045  1.045  1.045   1.045  -\n";
        card += "pdf_gg                    lnN        1.0755 -      -       1.0708 -\n";
        card += "pdf_qqbar                 lnN         -     1.023  1.0341  -      -\n";
        card += "pdf_H4l_accept            lnN        1.02   1.02   -       -      -\n";
        card += "QCDscale_ggH              lnN        1.0975 -      -       -      -\n";
        card += "QCDscale_qqH              lnN        -      1.003  -       -      -\n";
        card += "QCDscale_VH               lnN        -      -      -       -      -\n";
        card += "theoryUncXS_HighMH        lnN        1      1      -       -      -\n";
        card += "QCDscale_ggVV             lnN        -      -      -       1.2431 -\n";
        card += "QCDscale_VV               lnN        -      -      1.0284  -      -\n";
        card += "BRhiggs_ZZ4l              lnN        1.02   1.02   -       -      -\n";
    if (channel == 0 || channel == 2) {
        card +=  "CMS_eff_m                lnN        1.0224 1.0224 1.0224  1.0224 -\n";
        card +=  "CMS_scale_m              lnN        1.01   1.01   1.01    1.01   -\n";
        card +=  "CMS_trigger_m            lnN        1.015  1.015  1.015   1.015  -\n";
    }
    if (channel == 1 || channel == 2) {
        card += "CMS_eff_e                 lnN        1.0447 1.0447 1.0447  1.0447 -\n";
        card += "CMS_scale_e               lnN        1.02   1.02   1.02    1.02   -\n";
        card += "CMS_trigger_e             lnN        1.015  1.015  1.015   1.015  -\n";
    }
    if (channel == 0) {
        card += "CMS_hzz4mu_Zjets          gmN ZJEVT  -      -      -       -      ZJALPHA\n";
        card += "CMS_hzz4mu_Zjets_EXTRAP   lnN        -      -      -       -      ZJSYST\n";
    }
    if (channel == 1) {
        card += "CMS_hzz4e_Zjets           gmN ZJEVT  -      -      -       -      ZJALPHA\n";
        card += "CMS_hzz4e_Zjets_EXTRAP    lnN        -      -      -       -      ZJSYST\n";
    }
    if (channel == 2) {
        card += "CMS_hzz2e2mu_Zjets        gmN ZJEVT  -      -      -       -      ZJALPHA\n";
        card += "CMS_hzz2e2mu_Zjets_EXTRAP lnN        -      -      -       -      ZJSYST\n";
    }

    return card;
}


void getBkgFit(RooWorkspace& workspace, std::string varname, RooDataSet& data_bkg_red, RooRealVar& CMS_zz4l_mass, std::string outfilename, std::string pdfname, float low_M, float high_M, bool isgg) {
    RooRealVar CMS_zz4l_a0 ((varname+"_a0" ).c_str(),(varname+"_a0" ).c_str(),120.,100.,200.);
    RooRealVar CMS_zz4l_a1 ((varname+"_a1" ).c_str(),(varname+"_a1" ).c_str(),25.,10.,50.);
    RooRealVar CMS_zz4l_a2 ((varname+"_a2" ).c_str(),(varname+"_a2" ).c_str(),50.,20.,200.);
    RooRealVar CMS_zz4l_a3 ((varname+"_a3" ).c_str(),(varname+"_a3" ).c_str(),0.01,0.,1.);
    RooRealVar CMS_zz4l_a4 ((varname+"_a4" ).c_str(),(varname+"_a4" ).c_str(),200.,100.,400.);
    RooRealVar CMS_zz4l_a5 ((varname+"_a5" ).c_str(),(varname+"_a5" ).c_str(),10.,0.,150.);
    RooRealVar CMS_zz4l_a6 ((varname+"_a6" ).c_str(),(varname+"_a6" ).c_str(),10.,0.,100.);
    RooRealVar CMS_zz4l_a7 ((varname+"_a7" ).c_str(),(varname+"_a7" ).c_str(),0.1,0.,1.);
    RooRealVar CMS_zz4l_a8 ((varname+"_a8" ).c_str(),(varname+"_a8" ).c_str(),50.,0.,150.);
    RooRealVar CMS_zz4l_a9 ((varname+"_a9" ).c_str(),(varname+"_a9" ).c_str(),0.01,0.,1.);
    RooRealVar CMS_zz4l_a10((varname+"_a10").c_str(),(varname+"_a10").c_str(),60.,20.,200.);
    RooRealVar CMS_zz4l_a11((varname+"_a11").c_str(),(varname+"_a11").c_str(),-0.5,-1.,1.);
    RooRealVar CMS_zz4l_a12((varname+"_a12").c_str(),(varname+"_a12").c_str(),5000.,0.,10000.);
    RooRealVar CMS_zz4l_a13((varname+"_a13").c_str(),(varname+"_a13").c_str(),0.2,0.,1.);


    if (isgg) {
        RooggZZPdf_v2 background(pdfname.c_str(),pdfname.c_str(),CMS_zz4l_mass,CMS_zz4l_a0,CMS_zz4l_a1,CMS_zz4l_a2,CMS_zz4l_a3,CMS_zz4l_a4,CMS_zz4l_a5,CMS_zz4l_a6,CMS_zz4l_a7,CMS_zz4l_a8,CMS_zz4l_a9);
        background.fitTo(data_bkg_red,RooFit::SumW2Error(kTRUE));
        
        CMS_zz4l_a0.setConstant(kTRUE);
        CMS_zz4l_a1.setConstant(kTRUE);
        CMS_zz4l_a2.setConstant(kTRUE);
        CMS_zz4l_a3.setConstant(kTRUE);
        CMS_zz4l_a4.setConstant(kTRUE);
        CMS_zz4l_a5.setConstant(kTRUE);
        CMS_zz4l_a6.setConstant(kTRUE);
        CMS_zz4l_a7.setConstant(kTRUE);
        CMS_zz4l_a8.setConstant(kTRUE);
        CMS_zz4l_a9.setConstant(kTRUE);
        CMS_zz4l_a10.setConstant(kTRUE);
        CMS_zz4l_a11.setConstant(kTRUE);
        CMS_zz4l_a12.setConstant(kTRUE);
        CMS_zz4l_a13.setConstant(kTRUE);
        
        
        RooPlot *frame = CMS_zz4l_mass.frame(50);
        data_bkg_red.plotOn(frame,RooFit::DataError(RooAbsData::SumW2));
        background.plotOn(frame);

        workspace.import(background);

        TCanvas canvas;
        canvas.cd();
        frame->Draw();
        canvas.SaveAs((outfilename+".gif").c_str());

    }

    else {
        RooqqZZPdf_v2 background(pdfname.c_str(),pdfname.c_str(),CMS_zz4l_mass,CMS_zz4l_a0,CMS_zz4l_a1,CMS_zz4l_a2,CMS_zz4l_a3,CMS_zz4l_a4,CMS_zz4l_a5,CMS_zz4l_a6,CMS_zz4l_a7,CMS_zz4l_a8,CMS_zz4l_a9,CMS_zz4l_a10,CMS_zz4l_a11,CMS_zz4l_a12,CMS_zz4l_a13);

        background.fitTo(data_bkg_red,RooFit::SumW2Error(kTRUE));

        CMS_zz4l_a0.setConstant(kTRUE);
        CMS_zz4l_a1.setConstant(kTRUE);
        CMS_zz4l_a2.setConstant(kTRUE);
        CMS_zz4l_a3.setConstant(kTRUE);
        CMS_zz4l_a4.setConstant(kTRUE);
        CMS_zz4l_a5.setConstant(kTRUE);
        CMS_zz4l_a6.setConstant(kTRUE);
        CMS_zz4l_a7.setConstant(kTRUE);
        CMS_zz4l_a8.setConstant(kTRUE);
        CMS_zz4l_a9.setConstant(kTRUE);
        CMS_zz4l_a10.setConstant(kTRUE);
        CMS_zz4l_a11.setConstant(kTRUE);
        CMS_zz4l_a12.setConstant(kTRUE);
        CMS_zz4l_a13.setConstant(kTRUE);

        RooPlot *frame = CMS_zz4l_mass.frame(50);
        data_bkg_red.plotOn(frame,RooFit::DataError(RooAbsData::SumW2));
        background.plotOn(frame);

        workspace.import(background);

        TCanvas canvas;
        canvas.cd();
        frame->Draw();
        canvas.SaveAs((outfilename+".gif").c_str());
    }
    

}

void getZXFit(RooWorkspace& workspace, std::string varname, RooDataSet& data_bkg_red, RooRealVar& CMS_zz4l_mass, std::string outfilename, std::string pdfname, float low_M, float high_M) {
    RooRealVar CMS_zz4l_mean_zx((varname+"_mean_zx").c_str(),(varname+"_mean_zx").c_str(),130.,100.,200.);
    RooRealVar CMS_zz4l_sigma_zx((varname+"_sigma_zx").c_str(),(varname+"_sigma_zx").c_str(),10.,0.,100.);
    RooLandau zxLandau(pdfname.c_str(),pdfname.c_str(),CMS_zz4l_mass,CMS_zz4l_mean_zx,CMS_zz4l_sigma_zx);

    zxLandau.fitTo(data_bkg_red,RooFit::SumW2Error(kTRUE));

    CMS_zz4l_mean_zx.setConstant(kTRUE);
    CMS_zz4l_sigma_zx.setConstant(kTRUE);

    RooPlot *frame = CMS_zz4l_mass.frame(25);
    data_bkg_red.plotOn(frame,RooFit::DataError(RooAbsData::SumW2));
    zxLandau.plotOn(frame);

    workspace.import(zxLandau);

    TCanvas c1;
    c1.cd();
    frame->Draw();
    c1.SaveAs((outfilename+".gif").c_str());

}



void getSigFit(RooWorkspace& workspace, std::string varname, RooDataSet& data_sig_red, RooRealVar& CMS_zz4l_mass, std::string outfilename, std::string pdfname, float mass, float low_M, float high_M) {
    CMS_zz4l_mass.setRange(low_M, high_M);
    CMS_zz4l_mass.setBins(100000, "fft");
    //CMS_zz4l_mass.setBins(10000);


    RooRealVar CMS_zz4l_mean_sig((varname+"_mean_sig").c_str(), (varname+"_mean_sig").c_str(),0.,-10.,10.);
    RooRealVar CMS_zz4l_sigma_sig((varname+"_sigma_sig").c_str(), (varname+"_sigma_sig").c_str(),3.,0.,30.);
    RooRealVar CMS_zz4l_alpha((varname+"_alpha").c_str(),(varname+"_alpha").c_str(),2.,0.,10.);
    RooRealVar CMS_zz4l_n((varname+"_n").c_str(),(varname+"_n").c_str(),1.,-10.,10.);
    RooCBShape signalCB((varname+"_signalCB").c_str(),(varname+"_signalCB").c_str(),CMS_zz4l_mass,CMS_zz4l_mean_sig,CMS_zz4l_sigma_sig,CMS_zz4l_alpha,CMS_zz4l_n);

    RooRealVar CMS_zz4l_mean_BW((varname+"_mean_BW").c_str(),(varname+"_mean_BW").c_str(), mass);
    RooRealVar CMS_zz4l_gamma_BW((varname+"_gamma_BW").c_str(),(varname+"_gamma_BW").c_str(),1.0);
    RooRelBWUFParam signalBW((varname+"_signalBW").c_str(), (varname+"_signalBW").c_str(),CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_gamma_BW);

    RooFFTConvPdf signal(pdfname.c_str(),"BW (X) CB",CMS_zz4l_mass,signalBW,signalCB,2);
    signal.setBufferFraction(0.2);

    signal.fitTo(data_sig_red,RooFit::SumW2Error(kTRUE));

    CMS_zz4l_mean_sig.setConstant(kTRUE);
    CMS_zz4l_sigma_sig.setConstant(kTRUE);
    CMS_zz4l_alpha.setConstant(kTRUE);
    CMS_zz4l_n.setConstant(kTRUE);
    CMS_zz4l_mean_BW.setConstant(kTRUE);
    CMS_zz4l_gamma_BW.setConstant(kTRUE);

    workspace.import(signal);

    RooPlot *frame = CMS_zz4l_mass.frame(30);
    data_sig_red.plotOn(frame,RooFit::DataError(RooAbsData::SumW2));
    signal.plotOn(frame);

    TCanvas c1;
    c1.cd();
    frame->Draw();
    c1.SaveAs((outfilename+".gif").c_str());

}

void makeBkgWorkspace(RooWorkspace& w_tmp, std::string varname, std::string pdfname, std::string outfilename, RooWorkspace& w, RooRealVar& CMS_zz4l_mass, float low_M, float high_M, bool isgg) {
    //RooRealVar CMS_zz4l_mass("mass", "mass", low_M, high_M);

    if (!isgg) {
        RooRealVar CMS_zz4l_a0  = *w_tmp.var((varname+"_a0" ).c_str());
        RooRealVar CMS_zz4l_a1  = *w_tmp.var((varname+"_a1" ).c_str());
        RooRealVar CMS_zz4l_a2  = *w_tmp.var((varname+"_a2" ).c_str());
        RooRealVar CMS_zz4l_a3  = *w_tmp.var((varname+"_a3" ).c_str());
        RooRealVar CMS_zz4l_a4  = *w_tmp.var((varname+"_a4" ).c_str());
        RooRealVar CMS_zz4l_a5  = *w_tmp.var((varname+"_a5" ).c_str());
        RooRealVar CMS_zz4l_a6  = *w_tmp.var((varname+"_a6" ).c_str());
        RooRealVar CMS_zz4l_a7  = *w_tmp.var((varname+"_a7" ).c_str());
        RooRealVar CMS_zz4l_a8  = *w_tmp.var((varname+"_a8" ).c_str());
        RooRealVar CMS_zz4l_a9  = *w_tmp.var((varname+"_a9" ).c_str());
        RooRealVar CMS_zz4l_a10 = *w_tmp.var((varname+"_a10").c_str());
        RooRealVar CMS_zz4l_a11 = *w_tmp.var((varname+"_a11").c_str());
        RooRealVar CMS_zz4l_a12 = *w_tmp.var((varname+"_a12").c_str());
        RooRealVar CMS_zz4l_a13 = *w_tmp.var((varname+"_a13").c_str());

        RooqqZZPdf_v2 background(pdfname.c_str(),pdfname.c_str(),CMS_zz4l_mass,CMS_zz4l_a0,CMS_zz4l_a1,CMS_zz4l_a2,CMS_zz4l_a3,CMS_zz4l_a4,CMS_zz4l_a5,CMS_zz4l_a6,CMS_zz4l_a7,CMS_zz4l_a8,CMS_zz4l_a9,CMS_zz4l_a10,CMS_zz4l_a11,CMS_zz4l_a12,CMS_zz4l_a13);

        CMS_zz4l_a0.setConstant(kTRUE);
        CMS_zz4l_a1.setConstant(kTRUE);
        CMS_zz4l_a2.setConstant(kTRUE);
        CMS_zz4l_a3.setConstant(kTRUE); 
        CMS_zz4l_a4.setConstant(kTRUE);
        CMS_zz4l_a5.setConstant(kTRUE);
        CMS_zz4l_a6.setConstant(kTRUE);
        CMS_zz4l_a7.setConstant(kTRUE);
        CMS_zz4l_a8.setConstant(kTRUE);
        CMS_zz4l_a9.setConstant(kTRUE);
        CMS_zz4l_a10.setConstant(kTRUE);
        CMS_zz4l_a11.setConstant(kTRUE);
        CMS_zz4l_a12.setConstant(kTRUE);
        CMS_zz4l_a13.setConstant(kTRUE);

        w.import(background);

        RooPlot *frame = CMS_zz4l_mass.frame(25);
        background.plotOn(frame);
        
        TCanvas c1;
        c1.cd();
        frame->Draw();
        c1.SaveAs((outfilename+"_sigregion.gif").c_str());

    }

    else {
        RooRealVar CMS_zz4l_a0  = *w_tmp.var((varname+"_a0" ).c_str());
        RooRealVar CMS_zz4l_a1  = *w_tmp.var((varname+"_a1" ).c_str());
        RooRealVar CMS_zz4l_a2  = *w_tmp.var((varname+"_a2" ).c_str());
        RooRealVar CMS_zz4l_a3  = *w_tmp.var((varname+"_a3" ).c_str());
        RooRealVar CMS_zz4l_a4  = *w_tmp.var((varname+"_a4" ).c_str());
        RooRealVar CMS_zz4l_a5  = *w_tmp.var((varname+"_a5" ).c_str());
        RooRealVar CMS_zz4l_a6  = *w_tmp.var((varname+"_a6" ).c_str());
        RooRealVar CMS_zz4l_a7  = *w_tmp.var((varname+"_a7" ).c_str());
        RooRealVar CMS_zz4l_a8  = *w_tmp.var((varname+"_a8" ).c_str());
        RooRealVar CMS_zz4l_a9  = *w_tmp.var((varname+"_a9" ).c_str());

        RooggZZPdf_v2 background(pdfname.c_str(),pdfname.c_str(),CMS_zz4l_mass,CMS_zz4l_a0,CMS_zz4l_a1,CMS_zz4l_a2,CMS_zz4l_a3,CMS_zz4l_a4,CMS_zz4l_a5,CMS_zz4l_a6,CMS_zz4l_a7,CMS_zz4l_a8,CMS_zz4l_a9);

        CMS_zz4l_a0.setConstant(kTRUE);
        CMS_zz4l_a1.setConstant(kTRUE);
        CMS_zz4l_a2.setConstant(kTRUE);
        CMS_zz4l_a3.setConstant(kTRUE);
        CMS_zz4l_a4.setConstant(kTRUE);
        CMS_zz4l_a5.setConstant(kTRUE);
        CMS_zz4l_a6.setConstant(kTRUE);
        CMS_zz4l_a7.setConstant(kTRUE);
        CMS_zz4l_a8.setConstant(kTRUE);
        CMS_zz4l_a9.setConstant(kTRUE);

        w.import(background);

        RooPlot *frame = CMS_zz4l_mass.frame(25);
        background.plotOn(frame);
        
        TCanvas c1;
        c1.cd();
        frame->Draw();
        c1.SaveAs((outfilename+"_sigregion.gif").c_str());

    }

}

void makeZXWorkspace(RooWorkspace& w_tmp, std::string varname, std::string pdfname, std::string outfilename, RooWorkspace& w, RooRealVar& CMS_zz4l_mass, float low_M, float high_M) {

    //RooRealVar CMS_zz4l_mass("mass", "mass", low_M, high_M);

    RooRealVar CMS_zz4l_mean_zx  = *w_tmp.var((varname+"_mean_zx").c_str());
    RooRealVar CMS_zz4l_sigma_zx = *w_tmp.var((varname+"_sigma_zx").c_str());
    RooLandau zxLandau(pdfname.c_str(),pdfname.c_str(),CMS_zz4l_mass,CMS_zz4l_mean_zx,CMS_zz4l_sigma_zx);

    CMS_zz4l_mean_zx.setConstant(kTRUE);
    CMS_zz4l_sigma_zx.setConstant(kTRUE);

    w.import(zxLandau);

    RooPlot *frame = CMS_zz4l_mass.frame(25);
    zxLandau.plotOn(frame);

    TCanvas c1;
    c1.cd();
    frame->Draw();
    c1.SaveAs((outfilename+"_sigregion.gif").c_str());

}

void makeSigWorkspace(RooWorkspace& w_tmp, std::string varname, std::string pdfname, RooWorkspace& w, RooRealVar& CMS_zz4l_mass, float low_M, float high_M) {

    //RooRealVar CMS_zz4l_mass("mass", "mass", low_M, high_M);
    CMS_zz4l_mass.setBins(100000, "fft");

    RooRealVar CMS_zz4l_mean_sig  = *w_tmp.var((varname+"_mean_sig").c_str());
    RooRealVar CMS_zz4l_sigma_sig = *w_tmp.var((varname+"_sigma_sig").c_str());
    RooRealVar CMS_zz4l_alpha     = *w_tmp.var((varname+"_alpha").c_str());
    RooRealVar CMS_zz4l_n         = *w_tmp.var((varname+"_n").c_str());
    RooCBShape signalCB((varname+"_signalCB").c_str(),(varname+"_signalCB").c_str(),CMS_zz4l_mass,CMS_zz4l_mean_sig,CMS_zz4l_sigma_sig,CMS_zz4l_alpha,CMS_zz4l_n);

    RooRealVar CMS_zz4l_mean_BW   = *w_tmp.var((varname+"_mean_BW").c_str());
    RooRealVar CMS_zz4l_gamma_BW  = *w_tmp.var((varname+"_gamma_BW").c_str());
    RooRelBWUFParam signalBW((varname+"_signalBW").c_str(), (varname+"_signalBW").c_str(),CMS_zz4l_mass,CMS_zz4l_mean_BW,CMS_zz4l_gamma_BW);

    RooFFTConvPdf signal(pdfname.c_str(),"BW (X) CB",CMS_zz4l_mass,signalBW,signalCB,2);
    signal.setBufferFraction(0.2);

    CMS_zz4l_mean_sig.setConstant(kTRUE);
    CMS_zz4l_sigma_sig.setConstant(kTRUE);
    CMS_zz4l_alpha.setConstant(kTRUE);
    CMS_zz4l_n.setConstant(kTRUE);
    CMS_zz4l_mean_BW.setConstant(kTRUE);
    CMS_zz4l_gamma_BW.setConstant(kTRUE);

    w.import(signal);
}

void analysisEngine(std::string mass_str, std::string gghid, std::string vbfid, float higgsMass, float higgsMassLow, float higgsMassHigh, bool doShapeAnalysis, RooWorkspace& w_2e2mu, RooWorkspace& w_4e, RooWorkspace& w_4mu, bool doBkgFit, bool doSigFit) {

    std::cout << "Analyzing " << mass_str << " GeV mass point ... " << std::endl;

    std::string base_folder = "/home/avartak/CMS/Higgs/CMSSW_4_2_8_patch7/src/WWAnalysis/AnalysisStep/trees/conf1/";
    std::string card_name   = std::string("hzz_m")+mass_str;
    float min4l = higgsMassLow;
    float max4l = higgsMassHigh;
    float z1min = 40;
    float z2min = 12;
    std::string workspace_2e2mu = std::string("hzz2e2mu_m")+mass_str+"_workspace.root";
    std::string workspace_4e    = std::string("hzz4e_m")   +mass_str+"_workspace.root";
    std::string workspace_4mu   = std::string("hzz4mu_m")  +mass_str+"_workspace.root";

    init();

    ofstream result_out;
    result_out.open ((std::string("m")+mass_str +"_yields.txt").c_str());

    domufakes(base_folder+"hzzTree.root", mubarrel, muendcap, mubarrelerr, muendcaperr, result_out);
    doelfakes(base_folder+"hzzTree.root", elbarrel, elendcap, elbarrelerr, elendcaperr, result_out);


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

    float yield_zj_mm = 0.0;
    float yield_zj_ee = 0.0;
    float yield_zj_em = 0.0;

    float yield_zj_mm_count = 0.0;
    float yield_zj_ee_count = 0.0;
    float yield_zj_em_count = 0.0;

    float yield_zj_mm_error = 0.0;
    float yield_zj_ee_error = 0.0;
    float yield_zj_em_error = 0.0;

    std::string ggh_rootfile = "hzzTree_id";
    ggh_rootfile += gghid;
    ggh_rootfile += ".root";

    std::string vbf_rootfile = "hzzTree_id";
    vbf_rootfile += vbfid;
    vbf_rootfile += ".root";

    RooRealVar mass("mass",   "mass",   100, 100, 600);
    RooRealVar weight("weight", "weight", 0.,  0.,  10.);

    RooArgSet argset     (mass,     weight, "argset");

    RooDataSet dataset_2e2mu  ("dataset_2e2mu"  , "", argset    , RooFit::WeightVar("weight"));
    RooDataSet dataset_4e     ("dataset_4e"     , "", argset    , RooFit::WeightVar("weight"));
    RooDataSet dataset_4mu    ("dataset_4mu"    , "", argset    , RooFit::WeightVar("weight"));

    RooDataSet dataset_gz2e2mu("dataset_gz2e2mu", "", argset    , RooFit::WeightVar("weight"));
    RooDataSet dataset_gz4e   ("dataset_gz4e"   , "", argset    , RooFit::WeightVar("weight"));
    RooDataSet dataset_gz4mu  ("dataset_gz4mu"  , "", argset    , RooFit::WeightVar("weight"));

    RooDataSet dataset_g2e2mu ("dataset_g2e2mu" , "", argset    , RooFit::WeightVar("weight"));
    RooDataSet dataset_g4e    ("dataset_g4e"    , "", argset    , RooFit::WeightVar("weight"));
    RooDataSet dataset_g4mu   ("dataset_g4mu"   , "", argset    , RooFit::WeightVar("weight"));

    RooDataSet dataset_q2e2mu ("dataset_q2e2mu" , "", argset    , RooFit::WeightVar("weight"));
    RooDataSet dataset_q4e    ("dataset_q4e"    , "", argset    , RooFit::WeightVar("weight"));
    RooDataSet dataset_q4mu   ("dataset_q4mu"   , "", argset    , RooFit::WeightVar("weight"));

    RooDataSet dataset_zx2e2mu("dataset_zx2e2mu", "", argset    , RooFit::WeightVar("weight"));
    RooDataSet dataset_zx4e   ("dataset_zx4e"   , "", argset    , RooFit::WeightVar("weight"));
    RooDataSet dataset_zx4mu  ("dataset_zx4mu"  , "", argset    , RooFit::WeightVar("weight"));

    yield_zj_mm_count += getYield((base_folder+"hzzTree.root"),     0, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx4mu  , argset    ,  4)[1];
    yield_zj_ee_count += getYield((base_folder+"hzzTree.root"),     1, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx4e   , argset    ,  4)[1];
    yield_zj_em_count += getYield((base_folder+"hzzTree.root"),     2, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx2e2mu, argset    ,  4)[1];
    yield_zj_em_count += getYield((base_folder+"hzzTree.root"),     3, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx2e2mu, argset    ,  4)[1];

    yield_zj_mm_error += getYield((base_folder+"hzzTree.root"),     0, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx4mu  , argset    ,  4)[2];
    yield_zj_ee_error += getYield((base_folder+"hzzTree.root"),     1, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx4e   , argset    ,  4)[2];
    yield_zj_em_error += getYield((base_folder+"hzzTree.root"),     2, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx2e2mu, argset    ,  4)[2];
    yield_zj_em_error += getYield((base_folder+"hzzTree.root"),     3, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx2e2mu, argset    ,  4)[2];

    yield_zj_mm       += getYield(base_folder+"hzzTree.root",       0, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx4mu  , argset    ,  0)[0];
    yield_zj_ee       += getYield(base_folder+"hzzTree.root",       1, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx4e   , argset    ,  0)[0];
    yield_zj_em       += getYield(base_folder+"hzzTree.root",       2, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx2e2mu, argset    ,  0)[0];
    yield_zj_em       += getYield(base_folder+"hzzTree.root",       3, 1.0,                                    0, z1min, z2min, min4l, max4l, dataset_zx2e2mu, argset    ,  0)[0];

    yield_qq_mm       += getYield(base_folder+"hzzTree_id103b.root",0, (xsec_qqzz4mu    /evt_qqzz4mu)    *5.0, 1, z1min, z2min, min4l, max4l, dataset_4mu    , argset    ,  0)[0];
    yield_qq_mm       += getYield(base_folder+"hzzTree_id102b.root",0, (xsec_qqzz4e     /evt_qqzz4e)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_4mu    , argset    ,  0)[0];
    yield_qq_mm       += getYield(base_folder+"hzzTree_id104b.root",0, (xsec_qqzz4tau   /evt_qqzz4tau)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_4mu    , argset    ,  0)[0];
    yield_qq_mm       += getYield(base_folder+"hzzTree_id105b.root",0, (xsec_qqzz2e2mu  /evt_qqzz2e2mu)  *5.0, 1, z1min, z2min, min4l, max4l, dataset_4mu    , argset    ,  0)[0];
    yield_qq_mm       += getYield(base_folder+"hzzTree_id106b.root",0, (xsec_qqzz2e2tau /evt_qqzz2e2tau) *5.0, 1, z1min, z2min, min4l, max4l, dataset_4mu    , argset    ,  0)[0];
    yield_qq_mm       += getYield(base_folder+"hzzTree_id107b.root",0, (xsec_qqzz2mu2tau/evt_qqzz2mu2tau)*5.0, 1, z1min, z2min, min4l, max4l, dataset_4mu    , argset    ,  0)[0];
    yield_gg_mm       += getYield(base_folder+"hzzTree_id100.root" ,0, (xsec_ggzz2l2l   /evt_ggzz2l2l)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz4mu  , argset    ,  0)[0];
    yield_gg_mm       += getYield(base_folder+"hzzTree_id101.root" ,0, (xsec_ggzz4l     /evt_ggzz4l)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz4mu  , argset    ,  0)[0];
    
    yield_qq_ee       += getYield(base_folder+"hzzTree_id103b.root",1, (xsec_qqzz4mu    /evt_qqzz4mu)    *5.0, 1, z1min, z2min, min4l, max4l, dataset_4e     , argset    ,  0)[0];
    yield_qq_ee       += getYield(base_folder+"hzzTree_id102b.root",1, (xsec_qqzz4e     /evt_qqzz4e)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_4e     , argset    ,  0)[0];
    yield_qq_ee       += getYield(base_folder+"hzzTree_id104b.root",1, (xsec_qqzz4tau   /evt_qqzz4tau)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_4e     , argset    ,  0)[0];
    yield_qq_ee       += getYield(base_folder+"hzzTree_id105b.root",1, (xsec_qqzz2e2mu  /evt_qqzz2e2mu)  *5.0, 1, z1min, z2min, min4l, max4l, dataset_4e     , argset    ,  0)[0];
    yield_qq_ee       += getYield(base_folder+"hzzTree_id106b.root",1, (xsec_qqzz2e2tau /evt_qqzz2e2tau) *5.0, 1, z1min, z2min, min4l, max4l, dataset_4e     , argset    ,  0)[0];
    yield_qq_ee       += getYield(base_folder+"hzzTree_id107b.root",1, (xsec_qqzz2mu2tau/evt_qqzz2mu2tau)*5.0, 1, z1min, z2min, min4l, max4l, dataset_4e     , argset    ,  0)[0];
    yield_gg_ee       += getYield(base_folder+"hzzTree_id100.root" ,1, (xsec_ggzz2l2l   /evt_ggzz2l2l)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz4e   , argset    ,  0)[0];
    yield_gg_ee       += getYield(base_folder+"hzzTree_id101.root" ,1, (xsec_ggzz4l     /evt_ggzz4l)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz4e   , argset    ,  0)[0];

    yield_qq_em       += getYield(base_folder+"hzzTree_id103b.root",2, (xsec_qqzz4mu    /evt_qqzz4mu)    *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset    ,  0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id102b.root",2, (xsec_qqzz4e     /evt_qqzz4e)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset    ,  0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id104b.root",2, (xsec_qqzz4tau   /evt_qqzz4tau)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset    ,  0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id105b.root",2, (xsec_qqzz2e2mu  /evt_qqzz2e2mu)  *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset    ,  0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id106b.root",2, (xsec_qqzz2e2tau /evt_qqzz2e2tau) *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset    ,  0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id107b.root",2, (xsec_qqzz2mu2tau/evt_qqzz2mu2tau)*5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset    ,  0)[0];
    yield_gg_em       += getYield(base_folder+"hzzTree_id100.root" ,2, (xsec_ggzz2l2l   /evt_ggzz2l2l)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz2e2mu, argset    ,  0)[0];
    yield_gg_em       += getYield(base_folder+"hzzTree_id101.root" ,2, (xsec_ggzz4l     /evt_ggzz4l)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz2e2mu, argset    ,  0)[0];

    yield_qq_em       += getYield(base_folder+"hzzTree_id103b.root",3, (xsec_qqzz4mu    /evt_qqzz4mu)    *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset    ,  0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id102b.root",3, (xsec_qqzz4e     /evt_qqzz4e)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset    ,  0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id104b.root",3, (xsec_qqzz4tau   /evt_qqzz4tau)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset    ,  0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id105b.root",3, (xsec_qqzz2e2mu  /evt_qqzz2e2mu)  *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset    ,  0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id106b.root",3, (xsec_qqzz2e2tau /evt_qqzz2e2tau) *5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset    ,  0)[0];
    yield_qq_em       += getYield(base_folder+"hzzTree_id107b.root",3, (xsec_qqzz2mu2tau/evt_qqzz2mu2tau)*5.0, 1, z1min, z2min, min4l, max4l, dataset_2e2mu  , argset    ,  0)[0];
    yield_gg_em       += getYield(base_folder+"hzzTree_id100.root" ,3, (xsec_ggzz2l2l   /evt_ggzz2l2l)   *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz2e2mu, argset    ,  0)[0];
    yield_gg_em       += getYield(base_folder+"hzzTree_id101.root" ,3, (xsec_ggzz4l     /evt_ggzz4l)     *5.0, 1, z1min, z2min, min4l, max4l, dataset_gz2e2mu, argset    ,  0)[0];

    yield_ggh_mm      += getYield(base_folder+ggh_rootfile        , 0, (xggh[mass_str]/egg[mass_str])    *5.0, 1, z1min, z2min, min4l, max4l, dataset_g4mu   , argset    ,  0)[0];
    yield_vbh_mm      += getYield(base_folder+vbf_rootfile        , 0, (xvbf[mass_str]/evb[mass_str])    *5.0, 1, z1min, z2min, min4l, max4l, dataset_q4mu   , argset    ,  0)[0];

    yield_ggh_ee      += getYield(base_folder+ggh_rootfile        , 1, (xggh[mass_str]/egg[mass_str])    *5.0, 1, z1min, z2min, min4l, max4l, dataset_g4e    , argset    ,  0)[0];
    yield_vbh_ee      += getYield(base_folder+vbf_rootfile        , 1, (xvbf[mass_str]/evb[mass_str])    *5.0, 1, z1min, z2min, min4l, max4l, dataset_q4e    , argset    ,  0)[0];

    yield_ggh_em      += getYield(base_folder+ggh_rootfile        , 2, (xggh[mass_str]/egg[mass_str])    *5.0, 1, z1min, z2min, min4l, max4l, dataset_g2e2mu , argset    ,  0)[0];
    yield_vbh_em      += getYield(base_folder+vbf_rootfile        , 2, (xvbf[mass_str]/evb[mass_str])    *5.0, 1, z1min, z2min, min4l, max4l, dataset_q2e2mu , argset    ,  0)[0];

    yield_ggh_em      += getYield(base_folder+ggh_rootfile        , 3, (xggh[mass_str]/egg[mass_str])    *5.0, 1, z1min, z2min, min4l, max4l, dataset_g2e2mu , argset    ,  0)[0];
    yield_vbh_em      += getYield(base_folder+vbf_rootfile        , 3, (xvbf[mass_str]/evb[mass_str])    *5.0, 1, z1min, z2min, min4l, max4l, dataset_q2e2mu , argset    ,  0)[0];

    
    std::string card_4mu   = createCardTemplate(0, doShapeAnalysis, workspace_4mu.c_str());
    std::string card_4e    = createCardTemplate(1, doShapeAnalysis, workspace_4e.c_str());
    std::string card_2e2mu = createCardTemplate(2, doShapeAnalysis, workspace_2e2mu.c_str());

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
     
    card_2e2mu = findAndReplace(card_2e2mu, "ZJSYST",   yield_zj_em_error / yield_zj_em);
    card_4e    = findAndReplace(card_4e,    "ZJSYST",   yield_zj_ee_error / yield_zj_ee);
    card_4mu   = findAndReplace(card_4mu,   "ZJSYST",   yield_zj_mm_error / yield_zj_mm);

    card_2e2mu = findAndReplace(card_2e2mu, "BIN",   "a3");
    card_4e    = findAndReplace(card_4e,    "BIN",   "a2");
    card_4mu   = findAndReplace(card_4mu,   "BIN",   "a1");

    if (doShapeAnalysis) {
        card_2e2mu = findAndReplace(card_2e2mu, "OBS",      2.0);
        card_4e    = findAndReplace(card_4e,    "OBS",      2.0);
        card_4mu   = findAndReplace(card_4mu,   "OBS",      2.0);
    }     
    else {
        card_2e2mu = findAndReplace(card_2e2mu, "OBS",      yield_zj_em + yield_gg_em + yield_qq_em);
        card_4e    = findAndReplace(card_4e,    "OBS",      yield_zj_ee + yield_gg_ee + yield_qq_ee);
        card_4mu   = findAndReplace(card_4mu,   "OBS",      yield_zj_mm + yield_gg_mm + yield_qq_mm);
    }
  
    RooWorkspace w_2e2mu_tmp("w_2e2mu_tmp", "");
    RooWorkspace w_4e_tmp   ("w_4e_tmp"   , "");
    RooWorkspace w_4mu_tmp  ("w_4mu_tmp"  , "");
 
    RooRealVar CMS_zz4l_mass("mass", "mass", higgsMassLow, higgsMassHigh);
 
   if (doShapeAnalysis) {  
        if (doSigFit) { 
            getSigFit(w_2e2mu_tmp, "CMS_gh2e2mu", dataset_g2e2mu  , mass     , std::string("m")+mass_str+"_fit_g2e2mu" , "sig_ggH", higgsMass, higgsMassLow, higgsMassHigh);
            getSigFit(w_4e_tmp,    "CMS_gh4e",    dataset_g4e     , mass     , std::string("m")+mass_str+"_fit_g4e"    , "sig_ggH", higgsMass, higgsMassLow, higgsMassHigh);
            getSigFit(w_4mu_tmp,   "CMS_gh4mu",   dataset_g4mu    , mass     , std::string("m")+mass_str+"_fit_g4mu"   , "sig_ggH", higgsMass, higgsMassLow, higgsMassHigh);
        
            makeSigWorkspace(w_2e2mu_tmp, "CMS_gh2e2mu","sig_ggH", w_2e2mu, CMS_zz4l_mass, higgsMassLow, higgsMassHigh);
            makeSigWorkspace(w_4e_tmp,    "CMS_gh4e",   "sig_ggH", w_4e   , CMS_zz4l_mass, higgsMassLow, higgsMassHigh);
            makeSigWorkspace(w_4mu_tmp,   "CMS_gh4mu",  "sig_ggH", w_4mu  , CMS_zz4l_mass, higgsMassLow, higgsMassHigh);
        
            getSigFit(w_2e2mu_tmp, "CMS_vb2e2mu", dataset_q2e2mu  , mass     , std::string("m")+mass_str+"_fit_q2e2mu" , "sig_VBF", higgsMass, higgsMassLow, higgsMassHigh);
            getSigFit(w_4e_tmp,    "CMS_vb4e",    dataset_q4e     , mass     , std::string("m")+mass_str+"_fit_q4e"    , "sig_VBF", higgsMass, higgsMassLow, higgsMassHigh);
            getSigFit(w_4mu_tmp,   "CMS_vb4mu",   dataset_q4mu    , mass     , std::string("m")+mass_str+"_fit_q4mu"   , "sig_VBF", higgsMass, higgsMassLow, higgsMassHigh);

            makeSigWorkspace(w_2e2mu_tmp, "CMS_vb2e2mu","sig_VBF", w_2e2mu, CMS_zz4l_mass, higgsMassLow, higgsMassHigh);
            makeSigWorkspace(w_4e_tmp,    "CMS_vb4e",   "sig_VBF", w_4e   , CMS_zz4l_mass, higgsMassLow, higgsMassHigh);
            makeSigWorkspace(w_4mu_tmp,   "CMS_vb4mu",  "sig_VBF", w_4mu  , CMS_zz4l_mass, higgsMassLow, higgsMassHigh);
        }        

        if (doBkgFit) {
            getZXFit (w_2e2mu_tmp, "CMS_zx2e2mu", dataset_zx2e2mu , mass     , std::string("m")+mass_str+"_fit_zx2e2mu", "bkg_zjets", higgsMassLow, higgsMassHigh);
            getZXFit (w_4e_tmp,    "CMS_zx4e",    dataset_zx4e    , mass     , std::string("m")+mass_str+"_fit_zx4e"   , "bkg_zjets", higgsMassLow, higgsMassHigh);
            getZXFit (w_4mu_tmp,   "CMS_zx4mu",   dataset_zx4mu   , mass     , std::string("m")+mass_str+"_fit_zx4mu"  , "bkg_zjets", higgsMassLow, higgsMassHigh);
            
            makeZXWorkspace(w_2e2mu_tmp, "CMS_zx2e2mu", "bkg_zjets", "shape_zjets2e2mu", w_2e2mu, CMS_zz4l_mass, higgsMassLow, higgsMassHigh);
            makeZXWorkspace(w_4e_tmp,    "CMS_zx4e",    "bkg_zjets", "shape_zjets4e",    w_4e   , CMS_zz4l_mass, higgsMassLow, higgsMassHigh);
            makeZXWorkspace(w_4mu_tmp,   "CMS_zx4mu",   "bkg_zjets", "shape_zjets4mu",   w_4mu  , CMS_zz4l_mass, higgsMassLow, higgsMassHigh);
            
            getBkgFit(w_2e2mu_tmp, "CMS_gz2e2mu", dataset_gz2e2mu , mass     , std::string("m")+mass_str+"_fit_gz2e2mu", "bkg_ggzz" , higgsMassLow, higgsMassHigh, true);
            getBkgFit(w_4e_tmp,    "CMS_gz4e",    dataset_gz4e    , mass     , std::string("m")+mass_str+"_fit_gz4e"   , "bkg_ggzz" , higgsMassLow, higgsMassHigh, true);
            getBkgFit(w_4mu_tmp,   "CMS_gz4mu",   dataset_gz4mu   , mass     , std::string("m")+mass_str+"_fit_gz4mu"  , "bkg_ggzz" , higgsMassLow, higgsMassHigh, true);
            
            makeBkgWorkspace(w_2e2mu_tmp, "CMS_gz2e2mu", "bkg_ggzz", "shape_gz2e2mu", w_2e2mu   , CMS_zz4l_mass, higgsMassLow, higgsMassHigh, true);
            makeBkgWorkspace(w_4e_tmp,    "CMS_gz4e",    "bkg_ggzz", "shape_gz4e",    w_4e      , CMS_zz4l_mass, higgsMassLow, higgsMassHigh, true);
            makeBkgWorkspace(w_4mu_tmp,   "CMS_gz4mu",   "bkg_ggzz", "shape_gz4mu",   w_4mu     , CMS_zz4l_mass, higgsMassLow, higgsMassHigh, true);
            
            getBkgFit(w_2e2mu_tmp, "CMS_zz2e2mu", dataset_2e2mu   , mass     , std::string("m")+mass_str+"_fit_2e2mu"  , "bkg_qqzz" , higgsMassLow, higgsMassHigh, false);
            getBkgFit(w_4e_tmp,    "CMS_zz4e",    dataset_4e      , mass     , std::string("m")+mass_str+"_fit_4e"     , "bkg_qqzz" , higgsMassLow, higgsMassHigh, false);
            getBkgFit(w_4mu_tmp,   "CMS_zz4mu",   dataset_4mu     , mass     , std::string("m")+mass_str+"_fit_4mu"    , "bkg_qqzz" , higgsMassLow, higgsMassHigh, false);

            makeBkgWorkspace(w_2e2mu_tmp, "CMS_zz2e2mu", "bkg_qqzz", "shape_2e2mu",   w_2e2mu   , CMS_zz4l_mass, higgsMassLow, higgsMassHigh, false);
            makeBkgWorkspace(w_4e_tmp,    "CMS_zz4e",    "bkg_qqzz", "shape_4e",      w_4e      , CMS_zz4l_mass, higgsMassLow, higgsMassHigh, false);
            makeBkgWorkspace(w_4mu_tmp,   "CMS_zz4mu",   "bkg_qqzz", "shape_4mu",     w_4mu     , CMS_zz4l_mass, higgsMassLow, higgsMassHigh, false);


        }
        
        RooArgSet argset_obs(CMS_zz4l_mass, weight, "argset_obs");
        RooDataSet data_obs("data_obs", "data_obs", argset_obs, RooFit::WeightVar("weight"));
        
        argset_obs.setRealValue("mass", higgsMassLow+1.0);
        argset_obs.setRealValue("weight", 1.0);
        data_obs.add(argset_obs, 1.0);
        
        argset_obs.setRealValue("mass", higgsMassHigh-1.0);
        argset_obs.setRealValue("weight", 1.0);
        data_obs.add(argset_obs, 1.0);
        
        w_2e2mu.import(data_obs);
        w_4e.import   (data_obs);
        w_4mu.import  (data_obs);

        
        if (doSigFit) {
            w_2e2mu.writeToFile(workspace_2e2mu.c_str());
            w_4e.writeToFile(workspace_4e.c_str());
            w_4mu.writeToFile(workspace_4mu.c_str());
        }
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
    result_out << std::endl;

    result_out << "---- Higgs Yields ----" << std::endl;

    result_out << "H->ZZ->4e      : " << yield_ggh_ee << std::endl;
    result_out << "H->ZZ->4mu     : " << yield_ggh_mm << std::endl;
    result_out << "H->ZZ->2e2mu   : " << yield_ggh_em << std::endl;

    result_out << std::endl;

    result_out << "qqH->ZZ->4e    : " << yield_vbh_ee << std::endl;
    result_out << "qqH->ZZ->4mu   : " << yield_vbh_mm << std::endl;
    result_out << "qqH->ZZ->2e2mu : " << yield_vbh_em << std::endl;

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


}


void doHZZAnalysis() {

    bool doShapeAnalysis = false;

    RooWorkspace w_2e2mu_tmp("w_2e2mu", "");
    RooWorkspace w_4e_tmp   ("w_4e"   , "");
    RooWorkspace w_4mu_tmp  ("w_4mu"  , "");
           
    if (doShapeAnalysis) {
           
        analysisEngine("115", "200", "250", 115.,  100., 130., true, w_2e2mu_tmp, w_4e_tmp, w_4mu_tmp, true, false);
        
        RooWorkspace w_2e2mu_115 = w_2e2mu_tmp;
        RooWorkspace w_2e2mu_120 = w_2e2mu_tmp;
        RooWorkspace w_2e2mu_130 = w_2e2mu_tmp;
        RooWorkspace w_2e2mu_140 = w_2e2mu_tmp;
        RooWorkspace w_2e2mu_150 = w_2e2mu_tmp;
        RooWorkspace w_2e2mu_160 = w_2e2mu_tmp;
        
        RooWorkspace w_4e_115    = w_4e_tmp;
        RooWorkspace w_4e_120    = w_4e_tmp;
        RooWorkspace w_4e_130    = w_4e_tmp;
        RooWorkspace w_4e_140    = w_4e_tmp;
        RooWorkspace w_4e_150    = w_4e_tmp;
        RooWorkspace w_4e_160    = w_4e_tmp;
        
        RooWorkspace w_4mu_115   = w_4mu_tmp;
        RooWorkspace w_4mu_120   = w_4mu_tmp;
        RooWorkspace w_4mu_130   = w_4mu_tmp;
        RooWorkspace w_4mu_140   = w_4mu_tmp;
        RooWorkspace w_4mu_150   = w_4mu_tmp;
        RooWorkspace w_4mu_160   = w_4mu_tmp;
        
        analysisEngine("115", "200", "250", 115., 100., 130., true, w_2e2mu_115, w_4e_115, w_4mu_115, false, true);
        analysisEngine("120", "201", "251", 120., 100., 135., true, w_2e2mu_120, w_4e_120, w_4mu_120, false, true);
        analysisEngine("130", "202", "252", 130., 110., 145., true, w_2e2mu_130, w_4e_130, w_4mu_130, false, true);
        analysisEngine("140", "203", "253", 140., 120., 155., true, w_2e2mu_140, w_4e_140, w_4mu_140, false, true);
        analysisEngine("150", "204", "254", 150., 130., 165., true, w_2e2mu_150, w_4e_150, w_4mu_150, false, true);
        analysisEngine("160", "205", "255", 160., 140., 175., true, w_2e2mu_160, w_4e_160, w_4mu_160, false, true);
    }

    else {
        analysisEngine("115", "200", "250", 115., 100., 130., false, w_2e2mu_tmp, w_4e_tmp, w_4mu_tmp, false, false);
        analysisEngine("120", "201", "251", 120., 100., 135., false, w_2e2mu_tmp, w_4e_tmp, w_4mu_tmp, false, false);
        analysisEngine("130", "202", "252", 130., 110., 145., false, w_2e2mu_tmp, w_4e_tmp, w_4mu_tmp, false, false);
        analysisEngine("140", "203", "253", 140., 120., 155., false, w_2e2mu_tmp, w_4e_tmp, w_4mu_tmp, false, false);
        analysisEngine("150", "204", "254", 150., 130., 165., false, w_2e2mu_tmp, w_4e_tmp, w_4mu_tmp, false, false);
        analysisEngine("160", "205", "255", 160., 140., 175., false, w_2e2mu_tmp, w_4e_tmp, w_4mu_tmp, false, false);
    }

}
