#include <TFile.h>
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


float mubarrel[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float muendcap[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float elbarrel[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
float elendcap[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

void domufakes(std::string path, float* barrel, float* endcap) {
    TFile* file = new TFile(path.c_str());
    TTree* tree = (TTree*)file->Get("zllmtree/probe_tree");

    Float_t bins[] = {0.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,80.0};
    Int_t binnum = 10;
    Float_t ybins[] = {0.0,1.0,2.0};
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

        if (zmass>60 && zmass<120) {
            if (fabs(eta)<1.479) {
                hall->Fill(pt,0.5);
                if (id>0 && iso<0.4) {hpass->Fill(pt, 0.5);}
            }
            else {
                hall->Fill(pt,1.5);
                if (id>0 && iso<0.4) {hpass->Fill(pt, 1.5);}
            }
        }

    }

    hfake ->Divide(hpass,  hall);
    for (std::size_t k = 0; k < 10; k++) {
        barrel[k] = hfake->GetBinContent(k+1, 1);
        endcap[k] = hfake->GetBinContent(k+1, 2);
    }
    /*
    file->Close();
    delete file;
    delete hfake;
    delete hall;
    delete hpass;
    */
}


void doelfakes(std::string path, float* barrel, float* endcap) {
    TFile* file = new TFile(path.c_str());
    TTree* tree = (TTree*)file->Get("zlletree/probe_tree");
    
    Float_t bins[] = {0.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,80.0};
    Int_t binnum = 10;
    Float_t ybins[] = {0.0,1.0,2.0};
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

        if (mhits <= 1 && zmass>60 && zmass<120) {
            if (fabs(eta)<1.479) {
                hall->Fill(pt,0.5);
                if (pt <= 10.0 && ((fabs(eta)< 0.8 && id > 0.47) || (fabs(eta)>=0.8 && id>0.04)) && iso<0.4) {hpass->Fill(pt, 0.5);}
                if (pt >  10.0 && ((fabs(eta)< 0.8 && id > 0.50) || (fabs(eta)>=0.8 && id>0.12)) && iso<0.4) {hpass->Fill(pt, 0.5);}
            }
            else {
                hall->Fill(pt,1.5);
                if (pt <= 10.0 && id>0.295 && iso<0.4) {hpass->Fill(pt, 1.5);}
                if (pt >  10.0 && id>0.600 && iso<0.4) {hpass->Fill(pt, 1.5);}
            }
        }

    }

    hfake ->Divide(hpass , hall);
    for (std::size_t k = 0; k < 10; k++) {
        barrel[k] = hfake->GetBinContent(k+1, 1);
        endcap[k] = hfake->GetBinContent(k+1, 2);
    }
   
    /* 
    file->Close();
    delete file;
    delete hfake;
    delete hall;
    delete hpass;     
    */
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


std::pair<float, float> getYield(std::string path, int ch, float scale, bool isMC, float z1min, float z2min, float min4l, float max4l) {

    float yield = 0.0;
    float yield_count = 0.0;


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
                yield_count += 1.0;
            }

            else {
                bool existsAlready = false;
                for (std::size_t k = 0; k < runeventinfo.size(); k++) {
                    if (run == runeventinfo[k].first && event == runeventinfo[k].second) existsAlready = true;
                }

                if (!existsAlready) {
                    if (l3pdgId == -l4pdgId && (l3id==0 || l3iso>0.4) && (l4id==0 || l4iso>0.4)) {
                        runeventinfo.push_back(std::pair<int, int>(run, event));
                        float f1 = getFakeRate(l3pt, l3eta, l3pdgId);
                        float f2 = getFakeRate(l4pt, l4eta, l4pdgId);
                        yield += (f1/(1-f1))*(f2/(1-f2));
                        yield_count += 1.0;
                    }
                }
            }    

        }
    
    }
    file->Close();
    delete file;


    return std::pair<float, float>(yield, yield_count);

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
        card += "## mass window [100,140]\n";
        card += "bin     BIN      BIN      BIN       BIN       BIN       BIN\n";
        card += "process sig_ggH sig_VBF sig_WH   bkg_qqzz bkg_ggzz bkg_zjets\n";
        card += "process -4     -3       -2       1        2        3\n";
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
    if (channel == 0 || channel == 3) {
        card +=  "CMS_eff_m                 lnN        1.0224 1.0224 1.0224 1.0224  1.0224 -\n";
        card +=  "CMS_scale_m               lnN        1.01   1.01   1.01   1.01    1.01   -\n";
        card +=  "CMS_trigger_m             lnN        1.015  1.015  1.015  1.015   1.015  -\n";
    }
    if (channel == 1 || channel == 3) {
        card += "CMS_eff_e                 lnN        1.0447 1.0447 1.0447 1.0447  1.0447 -\n";
        card += "CMS_scale_e               lnN        1.02   1.02   1.02   1.02    1.02   -\n";
        card += "CMS_trigger_e             lnN        1.015  1.015  1.015  1.015   1.015  -\n";
    }
        card += "CMS_hzz2e2mu_Zjets        gmN ZJEVT  -      -      -      -       -      ZJALPHA\n";
        card += "CMS_hzz2e2mu_Zjets_EXTRAP lnN        -      -      -      -       -      ZJSYST\n";


    return card;
}

void doHZZCutnCountAnalysis() {

    std::string base_folder = "/home/avartak/CMS/Higgs/CMSSW_4_2_8_patch7/src/WWAnalysis/AnalysisStep/trees/conf2/";
    float min4l = 100;
    float max4l = 140;
    float z1min = 40;
    float z2min = 12;

    init();

    std::cout << "---- Computing fake rates for muons ----" << std::endl;
    domufakes(base_folder+"hzzTree.root", mubarrel, muendcap);
    std::cout << "---- Muon Barrel ----" << std::endl; 
    std::cout << "pT [ 0 -  5] : " << mubarrel[0] << endl;
    std::cout << "pT [ 7 -  7] : " << mubarrel[1] << endl;
    std::cout << "pT [ 7 - 10] : " << mubarrel[2] << endl;
    std::cout << "pT [10 - 15] : " << mubarrel[3] << endl;
    std::cout << "pT [15 - 20] : " << mubarrel[4] << endl;
    std::cout << "pT [20 - 25] : " << mubarrel[5] << endl;
    std::cout << "pT [25 - 30] : " << mubarrel[6] << endl;
    std::cout << "pT [30 - 40] : " << mubarrel[7] << endl;
    std::cout << "pT [40 - 50] : " << mubarrel[8] << endl;
    std::cout << "pT [50 - 80] : " << mubarrel[9] << endl;
    std::cout << std::endl;
    std::cout << "---- Muon Endcap ----" << std::endl;
    std::cout << "pT [ 0 -  5] : " << muendcap[0] << endl;
    std::cout << "pT [ 7 -  7] : " << muendcap[1] << endl;
    std::cout << "pT [ 7 - 10] : " << muendcap[2] << endl;
    std::cout << "pT [10 - 15] : " << muendcap[3] << endl;
    std::cout << "pT [15 - 20] : " << muendcap[4] << endl;
    std::cout << "pT [20 - 25] : " << muendcap[5] << endl;
    std::cout << "pT [25 - 30] : " << muendcap[6] << endl;
    std::cout << "pT [30 - 40] : " << muendcap[7] << endl;
    std::cout << "pT [40 - 50] : " << muendcap[8] << endl;
    std::cout << "pT [50 - 80] : " << muendcap[9] << endl;
    std::cout << std::endl;

    std::cout << "---- Computing fake rates for electrons ----" << std::endl;
    doelfakes(base_folder+"hzzTree.root", elbarrel, elendcap);
    std::cout << "---- Electron Barrel ----" << std::endl;
    std::cout << "pT [ 0 -  5] : " << elbarrel[0] << endl;
    std::cout << "pT [ 7 -  7] : " << elbarrel[1] << endl;
    std::cout << "pT [ 7 - 10] : " << elbarrel[2] << endl;
    std::cout << "pT [10 - 15] : " << elbarrel[3] << endl;
    std::cout << "pT [15 - 20] : " << elbarrel[4] << endl;
    std::cout << "pT [20 - 25] : " << elbarrel[5] << endl;
    std::cout << "pT [25 - 30] : " << elbarrel[6] << endl;
    std::cout << "pT [30 - 40] : " << elbarrel[7] << endl;
    std::cout << "pT [40 - 50] : " << elbarrel[8] << endl;
    std::cout << "pT [50 - 80] : " << elbarrel[9] << endl;
    std::cout << std::endl;
    std::cout << "---- Electron Endcap ----" << std::endl;
    std::cout << "pT [ 0 -  5] : " << elendcap[0] << endl;
    std::cout << "pT [ 7 -  7] : " << elendcap[1] << endl;
    std::cout << "pT [ 7 - 10] : " << elendcap[2] << endl;
    std::cout << "pT [10 - 15] : " << elendcap[3] << endl;
    std::cout << "pT [15 - 20] : " << elendcap[4] << endl;
    std::cout << "pT [20 - 25] : " << elendcap[5] << endl;
    std::cout << "pT [25 - 30] : " << elendcap[6] << endl;
    std::cout << "pT [30 - 40] : " << elendcap[7] << endl;
    std::cout << "pT [40 - 50] : " << elendcap[8] << endl;
    std::cout << "pT [50 - 80] : " << elendcap[9] << endl;
    std::cout << std::endl;

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

    yield_zj_mm_count += getYield((base_folder+"hzzTree.root"),     0, 1.0,                                    0, z1min, z2min, min4l, max4l).second;
    yield_zj_ee_count += getYield((base_folder+"hzzTree.root"),     1, 1.0,                                    0, z1min, z2min, min4l, max4l).second;
    yield_zj_em_count += getYield((base_folder+"hzzTree.root"),     2, 1.0,                                    0, z1min, z2min, min4l, max4l).second;
    yield_zj_em_count += getYield((base_folder+"hzzTree.root"),     3, 1.0,                                    0, z1min, z2min, min4l, max4l).second;

    yield_zj_mm       += getYield(base_folder+"hzzTree.root",       0, 1.0,                                    0, z1min, z2min, min4l, max4l).first;
    yield_zj_ee       += getYield(base_folder+"hzzTree.root",       1, 1.0,                                    0, z1min, z2min, min4l, max4l).first;
    yield_zj_em       += getYield(base_folder+"hzzTree.root",       2, 1.0,                                    0, z1min, z2min, min4l, max4l).first;
    yield_zj_em       += getYield(base_folder+"hzzTree.root",       3, 1.0,                                    0, z1min, z2min, min4l, max4l).first;

    yield_ggh_mm      += getYield(base_folder+"hzzTree_id108.root", 0, (xsec_gghzz120 /evt_gghzz120 )    *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_vbh_mm      += getYield(base_folder+"hzzTree_id109.root", 0, (xsec_vbfhzz120/evt_vbfhzz120)    *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_wzh_mm      += getYield(base_folder+"hzzTree_id110.root", 0, (xsec_vhzz120  /evt_vhzz120  )    *5.0, 1, z1min, z2min, min4l, max4l).first;

    yield_ggh_ee      += getYield(base_folder+"hzzTree_id108.root", 1, (xsec_gghzz120 /evt_gghzz120 )    *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_vbh_ee      += getYield(base_folder+"hzzTree_id109.root", 1, (xsec_vbfhzz120/evt_vbfhzz120)    *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_wzh_ee      += getYield(base_folder+"hzzTree_id110.root", 1, (xsec_vhzz120  /evt_vhzz120  )    *5.0, 1, z1min, z2min, min4l, max4l).first;

    yield_ggh_em      += getYield(base_folder+"hzzTree_id108.root", 2, (xsec_gghzz120 /evt_gghzz120 )    *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_vbh_em      += getYield(base_folder+"hzzTree_id109.root", 2, (xsec_vbfhzz120/evt_vbfhzz120)    *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_wzh_em      += getYield(base_folder+"hzzTree_id110.root", 2, (xsec_vhzz120  /evt_vhzz120  )    *5.0, 1, z1min, z2min, min4l, max4l).first;

    yield_ggh_em      += getYield(base_folder+"hzzTree_id108.root", 3, (xsec_gghzz120 /evt_gghzz120 )    *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_vbh_em      += getYield(base_folder+"hzzTree_id109.root", 3, (xsec_vbfhzz120/evt_vbfhzz120)    *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_wzh_em      += getYield(base_folder+"hzzTree_id110.root", 3, (xsec_vhzz120  /evt_vhzz120  )    *5.0, 1, z1min, z2min, min4l, max4l).first;

    yield_qq_mm       += getYield(base_folder+"hzzTree_id103.root", 0, (xsec_qqzz4mu    /evt_qqzz4mu)    *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_mm       += getYield(base_folder+"hzzTree_id102.root", 0, (xsec_qqzz4e     /evt_qqzz4e)     *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_mm       += getYield(base_folder+"hzzTree_id104.root", 0, (xsec_qqzz4tau   /evt_qqzz4tau)   *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_mm       += getYield(base_folder+"hzzTree_id105.root", 0, (xsec_qqzz2e2mu  /evt_qqzz2e2mu)  *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_mm       += getYield(base_folder+"hzzTree_id106.root", 0, (xsec_qqzz2e2tau /evt_qqzz2e2tau) *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_mm       += getYield(base_folder+"hzzTree_id107.root", 0, (xsec_qqzz2mu2tau/evt_qqzz2mu2tau)*5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_gg_mm       += getYield(base_folder+"hzzTree_id100.root", 0, (xsec_ggzz2l2l   /evt_ggzz2l2l)   *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_gg_mm       += getYield(base_folder+"hzzTree_id101.root", 0, (xsec_ggzz4l     /evt_ggzz4l)     *5.0, 1, z1min, z2min, min4l, max4l).first;
    
    yield_qq_ee       += getYield(base_folder+"hzzTree_id103.root", 1, (xsec_qqzz4mu    /evt_qqzz4mu)    *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_ee       += getYield(base_folder+"hzzTree_id102.root", 1, (xsec_qqzz4e     /evt_qqzz4e)     *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_ee       += getYield(base_folder+"hzzTree_id104.root", 1, (xsec_qqzz4tau   /evt_qqzz4tau)   *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_ee       += getYield(base_folder+"hzzTree_id105.root", 1, (xsec_qqzz2e2mu  /evt_qqzz2e2mu)  *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_ee       += getYield(base_folder+"hzzTree_id106.root", 1, (xsec_qqzz2e2tau /evt_qqzz2e2tau) *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_ee       += getYield(base_folder+"hzzTree_id107.root", 1, (xsec_qqzz2mu2tau/evt_qqzz2mu2tau)*5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_gg_ee       += getYield(base_folder+"hzzTree_id100.root", 1, (xsec_ggzz2l2l   /evt_ggzz2l2l)   *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_gg_ee       += getYield(base_folder+"hzzTree_id101.root", 1, (xsec_ggzz4l     /evt_ggzz4l)     *5.0, 1, z1min, z2min, min4l, max4l).first;

    yield_qq_em       += getYield(base_folder+"hzzTree_id103.root", 2, (xsec_qqzz4mu    /evt_qqzz4mu)    *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_em       += getYield(base_folder+"hzzTree_id102.root", 2, (xsec_qqzz4e     /evt_qqzz4e)     *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_em       += getYield(base_folder+"hzzTree_id104.root", 2, (xsec_qqzz4tau   /evt_qqzz4tau)   *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_em       += getYield(base_folder+"hzzTree_id105.root", 2, (xsec_qqzz2e2mu  /evt_qqzz2e2mu)  *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_em       += getYield(base_folder+"hzzTree_id106.root", 2, (xsec_qqzz2e2tau /evt_qqzz2e2tau) *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_em       += getYield(base_folder+"hzzTree_id107.root", 2, (xsec_qqzz2mu2tau/evt_qqzz2mu2tau)*5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_gg_em       += getYield(base_folder+"hzzTree_id100.root", 2, (xsec_ggzz2l2l   /evt_ggzz2l2l)   *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_gg_em       += getYield(base_folder+"hzzTree_id101.root", 2, (xsec_ggzz4l     /evt_ggzz4l)     *5.0, 1, z1min, z2min, min4l, max4l).first;

    yield_qq_em       += getYield(base_folder+"hzzTree_id103.root", 3, (xsec_qqzz4mu    /evt_qqzz4mu)    *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_em       += getYield(base_folder+"hzzTree_id102.root", 3, (xsec_qqzz4e     /evt_qqzz4e)     *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_em       += getYield(base_folder+"hzzTree_id104.root", 3, (xsec_qqzz4tau   /evt_qqzz4tau)   *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_em       += getYield(base_folder+"hzzTree_id105.root", 3, (xsec_qqzz2e2mu  /evt_qqzz2e2mu)  *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_em       += getYield(base_folder+"hzzTree_id106.root", 3, (xsec_qqzz2e2tau /evt_qqzz2e2tau) *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_qq_em       += getYield(base_folder+"hzzTree_id107.root", 3, (xsec_qqzz2mu2tau/evt_qqzz2mu2tau)*5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_gg_em       += getYield(base_folder+"hzzTree_id100.root", 3, (xsec_ggzz2l2l   /evt_ggzz2l2l)   *5.0, 1, z1min, z2min, min4l, max4l).first;
    yield_gg_em       += getYield(base_folder+"hzzTree_id101.root", 3, (xsec_ggzz4l     /evt_ggzz4l)     *5.0, 1, z1min, z2min, min4l, max4l).first;

    std::string card_4mu   = createCardTemplate(0);
    std::string card_4e    = createCardTemplate(1);
    std::string card_2e2mu = createCardTemplate(2);

    card_2e2mu = findAndReplaceFloat(card_2e2mu, "SIG_GGH",  yield_ggh_em);
    card_4e    = findAndReplaceFloat(card_4e,    "SIG_GGH",  yield_ggh_ee);
    card_4mu   = findAndReplaceFloat(card_4mu  , "SIG_GGH",  yield_ggh_mm);
     
    card_2e2mu = findAndReplaceFloat(card_2e2mu, "SIG_VBF",  yield_vbh_em);
    card_4e    = findAndReplaceFloat(card_4e,    "SIG_VBF",  yield_vbh_ee);
    card_4mu   = findAndReplaceFloat(card_4mu,   "SIG_VBF",  yield_vbh_mm);
     
    card_2e2mu = findAndReplaceFloat(card_2e2mu, "SIG_WH",   yield_wzh_em);
    card_4e    = findAndReplaceFloat(card_4e,    "SIG_WH",   yield_wzh_ee);
    card_4mu   = findAndReplaceFloat(card_4mu,   "SIG_WH",   yield_wzh_mm);
     
    card_2e2mu = findAndReplaceFloat(card_2e2mu, "BKG_QQZZ", yield_qq_em);
    card_4e    = findAndReplaceFloat(card_4e,    "BKG_QQZZ", yield_qq_ee);
    card_4mu   = findAndReplaceFloat(card_4mu,   "BKG_QQZZ", yield_qq_mm);
     
    card_2e2mu = findAndReplaceFloat(card_2e2mu, "BKG_GGZZ", yield_gg_em);
    card_4e    = findAndReplaceFloat(card_4e,    "BKG_GGZZ", yield_gg_ee);
    card_4mu   = findAndReplaceFloat(card_4mu,   "BKG_GGZZ", yield_gg_mm);
     
    card_2e2mu = findAndReplaceFloat(card_2e2mu, "BKG_ZJETS",yield_zj_em);
    card_4e    = findAndReplaceFloat(card_4e,    "BKG_ZJETS",yield_zj_ee);
    card_4mu   = findAndReplaceFloat(card_4mu,   "BKG_ZJETS",yield_zj_mm);
     
    card_2e2mu = findAndReplaceFloat(card_2e2mu, "OBS",      yield_zj_em + yield_gg_em + yield_qq_em);
    card_4e    = findAndReplaceFloat(card_4e,    "OBS",      yield_zj_ee + yield_gg_ee + yield_qq_ee);
    card_4mu   = findAndReplaceFloat(card_4mu,   "OBS",      yield_zj_mm + yield_gg_mm + yield_qq_mm);
     
    card_2e2mu = findAndReplaceFloat(card_2e2mu, "ZJEVT",    yield_zj_em_count);
    card_4e    = findAndReplaceFloat(card_4e,    "ZJEVT",    yield_zj_ee_count);
    card_4mu   = findAndReplaceFloat(card_4mu,   "ZJEVT",    yield_zj_mm_count);
     
    card_2e2mu = findAndReplaceFloat(card_2e2mu, "ZJALPHA",  yield_zj_em / yield_zj_em_count);
    card_4e    = findAndReplaceFloat(card_4e,    "ZJALPHA",  yield_zj_ee / yield_zj_ee_count);
    card_4mu   = findAndReplaceFloat(card_4mu,   "ZJALPHA",  yield_zj_mm / yield_zj_mm_count);
     
    card_2e2mu = findAndReplaceFloat(card_2e2mu, "ZJSYST",   1.33);
    card_4e    = findAndReplaceFloat(card_4e,    "ZJSYST",   1.33);
    card_4mu   = findAndReplaceFloat(card_4mu,   "ZJSYST",   1.33);

    card_2e2mu = findAndReplace(card_2e2mu, "BIN",   "a3");
    card_4e    = findAndReplace(card_4e,    "BIN",   "a2");
    card_4mu   = findAndReplace(card_4mu,   "BIN",   "a1");

    std::cout << "---- ZX Yields ----" << std::endl;
    std::cout << "4e     : " << yield_zj_ee << " (" << yield_zj_ee_count << ")" << std::endl;
    std::cout << "4mu    : " << yield_zj_mm << " (" << yield_zj_mm_count << ")" << std::endl;
    std::cout << "2e2mu  : " << yield_zj_em << " (" << yield_zj_em_count << ")" << std::endl;

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
    file_2e2mu.open ("hzz4l_2e2mu.txt");
    file_2e2mu << card_2e2mu;
    file_2e2mu.close();

    ofstream file_4e;
    file_4e.open ("hzz4l_4e.txt");
    file_4e << card_4e;
    file_4e.close();

    ofstream file_4mu;
    file_4mu.open ("hzz4l_4mu.txt");
    file_4mu << card_4mu;
    file_4mu.close();

}
