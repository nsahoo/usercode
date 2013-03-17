#ifndef FAKERATECALCULATOR_H
#define FAKERATECALCULATOR_H

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "scales2.h"

class FakeRateCalculator {

    private:
        std::vector<std::pair<int, int> > runeventinfo;    
        std::vector<std::pair<int, int> > runeventinfo2;    

        std::vector<float> mubarrel;
        std::vector<float> muendcap;
        std::vector<float> elbarrel;
        std::vector<float> elendcap;

        std::vector<float> mubarrelerr;
        std::vector<float> muendcaperr;
        std::vector<float> elbarrelerr;
        std::vector<float> elendcaperr;

        std::vector<std::pair<float, float> > elbarrelnew;
        std::vector<std::pair<float, float> > elendcapnew;

        std::vector<float> elbarrelavg;
        std::vector<float> elendcapavg;

        bool dofsr;

        std::vector<TGraphErrors*> elbarrelgr;
        std::vector<TGraphErrors*> elendcapgr;

        TH2D  *havg1, *havg2, *havg3, *havg4;

        float getMass(float l1pt, float l1eta, float l1phi, float l2pt, float l2eta, float l2phi) {
            float px1 = l1pt*cos(l1phi);
            float py1 = l1pt*sin(l1phi);
            float pz1 = l1pt*sinh(l1eta);
            float p1  = sqrt(px1*px1 + py1*py1 + pz1*pz1);
            
            float px2 = l2pt*cos(l2phi);
            float py2 = l2pt*sin(l2phi);
            float pz2 = l2pt*sinh(l2eta);
            float p2  = sqrt(px2*px2 + py2*py2 + pz2*pz2);
            
            TLorentzVector v1(px1, py1, pz1, p1);
            TLorentzVector v2(px2, py2, pz2, p2);
            
            v1 += v2;
            
            return v1.M();
        }

        float getMass2(float l1pt, float l1eta, float l1phi, float l2pt, float l2eta, float l2phi, float l3pt, float l3eta, float l3phi) {
            float px1 = l1pt*cos(l1phi);
            float py1 = l1pt*sin(l1phi);
            float pz1 = l1pt*sinh(l1eta);
            float p1  = sqrt(px1*px1 + py1*py1 + pz1*pz1);

            float px2 = l2pt*cos(l2phi);
            float py2 = l2pt*sin(l2phi);
            float pz2 = l2pt*sinh(l2eta);
            float p2  = sqrt(px2*px2 + py2*py2 + pz2*pz2);

            float px3 = l3pt*cos(l3phi);
            float py3 = l3pt*sin(l3phi);
            float pz3 = l3pt*sinh(l3eta);
            float p3  = sqrt(px3*px3 + py3*py3 + pz3*pz3);

            TLorentzVector v1(px1, py1, pz1, p1);
            TLorentzVector v2(px2, py2, pz2, p2);
            TLorentzVector v3(px3, py3, pz3, p3);

            v1 += v2;
            v1 += v3;

            return v1.M();
        }

        void domufakes(std::string path, bool do7TeV, float zmin, float zmax, float drcut, float mcut) {
            runeventinfo.clear();
            TFile* file = TFile::Open(path.c_str());
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
            TBranch *bl1pt    = tree->GetBranch("l1pt");
            TBranch *bl1eta   = tree->GetBranch("l1eta");
            TBranch *bl1phi   = tree->GetBranch("l1phi");
            TBranch *bl2pt    = tree->GetBranch("l2pt");
            TBranch *bl2eta   = tree->GetBranch("l2eta");
            TBranch *bl2phi   = tree->GetBranch("l2phi");
            TBranch *bpt      = tree->GetBranch("pt");
            TBranch *beta     = tree->GetBranch("eta");
            TBranch *bphi     = tree->GetBranch("phi");
            TBranch *bid      = tree->GetBranch("newID");
            TBranch *biso     = dofsr ? tree->GetBranch("pfCombRelIso04EACorr_WithFSR") : tree->GetBranch("pfCombRelIso04EACorr");
            TBranch *brun     = tree->GetBranch("run");
            TBranch *blum     = tree->GetBranch("lumi");
            TBranch *bevt     = tree->GetBranch("event");
            
            float met     = 0.0;
            float l1pt    = 0.0;
            float l1eta   = 0.0;
            float l1phi   = 0.0;
            float l2pt    = 0.0;
            float l2eta   = 0.0;
            float l2phi   = 0.0;
            float zmass   = 0.0;
            float pt      = 0.0;
            float eta     = 0.0;
            float phi     = 0.0;
            float iso     = 0.0;
            int   id      = 0;
            int   run     = 0;
            int   lum     = 0;
            int   evt     = 0;
            
            bmet      ->SetAddress(&met);
            bzmass    ->SetAddress(&zmass);
            bl1pt     ->SetAddress(&l1pt);
            bl1eta    ->SetAddress(&l1eta);
            bl1phi    ->SetAddress(&l1phi);
            bl2pt     ->SetAddress(&l2pt);
            bl2eta    ->SetAddress(&l2eta);
            bl2phi    ->SetAddress(&l2phi);
            bpt       ->SetAddress(&pt);
            beta      ->SetAddress(&eta);
            bphi      ->SetAddress(&phi);
            bid       ->SetAddress(&id);
            biso      ->SetAddress(&iso);
            brun      ->SetAddress(&run);
            blum      ->SetAddress(&lum);
            bevt      ->SetAddress(&evt);
           
            runeventinfo.clear(); 
            for (int i = 0; i < tree->GetEntries(); i++) {
                bmet     ->GetEvent(i);
                bzmass   ->GetEvent(i);
                bl1pt    ->GetEvent(i);
                bl1eta   ->GetEvent(i);
                bl1phi   ->GetEvent(i);
                bl2pt    ->GetEvent(i);
                bl2eta   ->GetEvent(i);
                bl2phi   ->GetEvent(i);
                bpt      ->GetEvent(i);
                beta     ->GetEvent(i);
                bphi     ->GetEvent(i);
                bid      ->GetEvent(i);
                biso     ->GetEvent(i);
                brun     ->GetEvent(i);
                blum     ->GetEvent(i);
                bevt     ->GetEvent(i);
            
                bool existsAlready = false;
                for (std::size_t k = 0; k < runeventinfo.size(); k++) {
                    if (run == runeventinfo[k].first && evt == runeventinfo[k].second) existsAlready = true;
                }
       
                if (existsAlready) continue;
                runeventinfo.push_back(std::pair<int, int>(run, evt));
 
                float dR1 = reco::deltaR(l1eta, l1phi, eta, phi); 
                float dR2 = reco::deltaR(l2eta, l2phi, eta, phi); 
            
                float m1 = getMass(l1pt, l1eta, l1phi, pt, eta, phi); 
                float m2 = getMass(l2pt, l2eta, l2phi, pt, eta, phi); 
            
                if (zmass>zmin && zmass<zmax && dR1 > drcut && dR2 > drcut && m1>mcut && m2>mcut) {
                    hall->Fill(pt,fabs(eta));
                    if (id>0 && iso<0.4) {
                        hpass->Fill(pt, fabs(eta));
                    }
                }
            
            }
            
            hfake ->Divide(hpass,  hall);
            
            for (std::size_t k = 0; k < 10; k++) {
                mubarrel[k] = hfake->GetBinContent(k+1, 1);
                muendcap[k] = hfake->GetBinContent(k+1, 2);
                mubarrelerr[k] = hfake->GetBinError(k+1, 1);
                muendcaperr[k] = hfake->GetBinError(k+1, 2);
            }
            
            file->Close();
            delete file;
        }

        void doelavg(std::string path, float z1min, float z2min, float massmin) {
            runeventinfo2.clear();
            TFile* file = TFile::Open(path.c_str());
            TTree* tree = (TTree*)file->Get("zxTreeSS/probe_tree");
            
            Float_t bins[] = {0.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,80.0};
            Int_t binnum = 10;
            Float_t ybins[] = {0.0,1.479,2.5};
            Int_t ybinnum = 2;
           
            TH2D* havgnum = new TH2D("havgnum", "", binnum, bins, ybinnum, ybins);
            TH2D* havgden = new TH2D("havgden", "", binnum, bins, ybinnum, ybins);
            
            TBranch *bchannel   = tree->GetBranch("channel");
            TBranch *bz1mass    = tree->GetBranch("z1mass");
            TBranch *bz2mass    = tree->GetBranch("z2mass");
            TBranch *bmass      = tree->GetBranch("mass");
            TBranch *bl3pt      = tree->GetBranch("l3pt");
            TBranch *bl3eta     = tree->GetBranch("l3eta");
            TBranch *bl3phi     = tree->GetBranch("l3phi");
            TBranch *bl3id      = tree->GetBranch("l3idNew");
            TBranch *bl3iso     = tree->GetBranch("l3pfIsoComb04EACorr");
            TBranch *bl3pdgId   = tree->GetBranch("l3pdgId");
            TBranch *bl3nmhits  = tree->GetBranch("l3nmisshits");
            TBranch *bl4pt      = tree->GetBranch("l4pt");
            TBranch *bl4eta     = tree->GetBranch("l4eta");
            TBranch *bl4phi     = tree->GetBranch("l4phi");
            TBranch *bl4id      = tree->GetBranch("l4idNew");
            TBranch *bl4iso     = tree->GetBranch("l4pfIsoComb04EACorr");
            TBranch *bl4pdgId   = tree->GetBranch("l4pdgId");
            TBranch *bl4nmhits  = tree->GetBranch("l4nmisshits");
            TBranch *bl1pt      = tree->GetBranch("l1pt");
            TBranch *bl1eta     = tree->GetBranch("l1eta");
            TBranch *bl1phi     = tree->GetBranch("l1phi");
            TBranch *bl1pdgId   = tree->GetBranch("l1pdgId");
            TBranch *bl2pt      = tree->GetBranch("l2pt");
            TBranch *bl2eta     = tree->GetBranch("l2eta");
            TBranch *bl2phi     = tree->GetBranch("l2phi");
            TBranch *bl2pdgId   = tree->GetBranch("l2pdgId");
            TBranch *bevent     = tree->GetBranch("event");
            TBranch *brun       = tree->GetBranch("run");

            float channel   = 0.0;
            float z1mass    = 0.0;
            float z2mass    = 0.0;
            float mass      = 0.0;
            float l3pt      = 0.0;
            float l3eta     = 0.0;
            float l3phi     = 0.0;
            float l3id      = 0.0;
            float l3iso     = 0.0;
            float l3pdgId   = 0.0;
            float l3nmhits  = 0.0;
            float l4pt      = 0.0;
            float l4eta     = 0.0;
            float l4phi     = 0.0;
            float l4id      = 0.0;
            float l4iso     = 0.0;
            float l4pdgId   = 0.0;
            float l4nmhits  = 0.0;
            float l1pt      = 0.0;
            float l1eta     = 0.0;
            float l1phi     = 0.0;
            float l1pdgId   = 0.0;
            float l2pt      = 0.0;
            float l2eta     = 0.0;
            float l2phi     = 0.0;
            float l2pdgId   = 0.0;
            int   event     = 0;
            int   run       = 0;

            bchannel   ->SetAddress(&channel);
            bz1mass    ->SetAddress(&z1mass);
            bz2mass    ->SetAddress(&z2mass);
            bmass      ->SetAddress(&mass);
            bl3pt      ->SetAddress(&l3pt);
            bl3eta     ->SetAddress(&l3eta);
            bl3phi     ->SetAddress(&l3phi);
            bl3id      ->SetAddress(&l3id);
            bl3iso     ->SetAddress(&l3iso);
            bl3pdgId   ->SetAddress(&l3pdgId);
            bl3nmhits  ->SetAddress(&l3nmhits);
            bl4pt      ->SetAddress(&l4pt);
            bl4eta     ->SetAddress(&l4eta);
            bl4phi     ->SetAddress(&l4phi);
            bl4id      ->SetAddress(&l4id);
            bl4iso     ->SetAddress(&l4iso);
            bl4pdgId   ->SetAddress(&l4pdgId);
            bl4nmhits  ->SetAddress(&l4nmhits);
            bl1pt      ->SetAddress(&l1pt);
            bl1eta     ->SetAddress(&l1eta);
            bl1phi     ->SetAddress(&l1phi);
            bl1pdgId   ->SetAddress(&l1pdgId);
            bl2pt      ->SetAddress(&l2pt);
            bl2eta     ->SetAddress(&l2eta);
            bl2phi     ->SetAddress(&l2phi);
            bl2pdgId   ->SetAddress(&l2pdgId);
            bevent     ->SetAddress(&event);
            brun       ->SetAddress(&run);

            for (int i = 0; i < tree->GetEntries(); i++) {
                bchannel   ->GetEvent(i);
                bz1mass    ->GetEvent(i);
                bz2mass    ->GetEvent(i);
                bmass      ->GetEvent(i);
                bl3pt      ->GetEvent(i);
                bl3eta     ->GetEvent(i);
                bl3phi     ->GetEvent(i);
                bl3id      ->GetEvent(i);
                bl3iso     ->GetEvent(i);
                bl3pdgId   ->GetEvent(i);
                bl3nmhits  ->GetEvent(i);
                bl4pt      ->GetEvent(i);
                bl4eta     ->GetEvent(i);
                bl4phi     ->GetEvent(i);
                bl4id      ->GetEvent(i);
                bl4iso     ->GetEvent(i);
                bl4pdgId   ->GetEvent(i);
                bl4nmhits  ->GetEvent(i);
                bevent     ->GetEvent(i);
                brun       ->GetEvent(i);
                bl1pt      ->GetEvent(i);
                bl1eta     ->GetEvent(i);
                bl1phi     ->GetEvent(i);
                bl1pdgId   ->GetEvent(i);
                bl2pt      ->GetEvent(i);
                bl2eta     ->GetEvent(i);
                bl2phi     ->GetEvent(i);
                bl2pdgId   ->GetEvent(i);

                bool existsAlready = false; 
                for (std::size_t k = 0; k < runeventinfo2.size(); k++) {
                    if (run == runeventinfo2[k].first && event == runeventinfo2[k].second) existsAlready = true;
                }

                if (!existsAlready && z1mass>z1min && z1mass<120. && z2mass>z2min && z2mass<120. && mass>massmin) {
                    runeventinfo2.push_back(std::pair<int, int>(run, event));
                    if (l3nmhits >= 0.0) {
                        havgden->Fill(l3pt,fabs(l3eta));
                        havgnum->Fill(l3pt,fabs(l3eta), l3nmhits);
                    }
                    if (l4nmhits >= 0.0) {
                        havgden->Fill(l4pt,fabs(l4eta));
                        havgnum->Fill(l4pt,fabs(l4eta), l4nmhits);
                    }
                }
            }

            havgnum->Divide(havgden);
            for (std::size_t k = 0; k < 10; k++) {
                elbarrelavg[k] = havgnum->GetBinContent(k+1, 1);
                elendcapavg[k] = havgnum->GetBinContent(k+1, 2);
            }    
        }

        void doelfakes(std::string path, bool do7TeV, float zmin, float zmax, float drcut, float mcut, float zgammacut=5.0) {
            runeventinfo.clear();
            TFile* file = TFile::Open(path.c_str());
            TTree* tree = (TTree*)file->Get("zlletree/probe_tree");
            
            
            Float_t bins[] = {0.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,80.0};
            Int_t binnum = 10;
            Float_t ybins[] = {0.0,1.479,2.5};
            Int_t ybinnum = 2;
           
            TH2D* hfake = new TH2D("hfake", "", binnum, bins, ybinnum, ybins);
            TH2D* hall  = new TH2D("hall",  "", binnum, bins, ybinnum, ybins);
            TH2D* hpass = new TH2D("hpass", "", binnum, bins, ybinnum, ybins);
            
            TH2D* hfake1 = new TH2D("hfake1", "", binnum, bins, ybinnum, ybins);
            TH2D* hall1  = new TH2D("hall1",  "", binnum, bins, ybinnum, ybins);
            TH2D* hpass1 = new TH2D("hpass1", "", binnum, bins, ybinnum, ybins);
            
            TH2D* hfake2 = new TH2D("hfake2", "", binnum, bins, ybinnum, ybins);
            TH2D* hall2  = new TH2D("hall2",  "", binnum, bins, ybinnum, ybins);
            TH2D* hpass2 = new TH2D("hpass2", "", binnum, bins, ybinnum, ybins);
            
            TH2D* hfake3 = new TH2D("hfake3", "", binnum, bins, ybinnum, ybins);
            TH2D* hall3  = new TH2D("hall3",  "", binnum, bins, ybinnum, ybins);
            TH2D* hpass3 = new TH2D("hpass3", "", binnum, bins, ybinnum, ybins);
            
            TH2D* hfake4 = new TH2D("hfake4", "", binnum, bins, ybinnum, ybins);
            TH2D* hall4  = new TH2D("hall4",  "", binnum, bins, ybinnum, ybins);
            TH2D* hpass4 = new TH2D("hpass4", "", binnum, bins, ybinnum, ybins);
            
            TH2D* hden1   = new TH2D("hden1", "", binnum, bins, ybinnum, ybins);
            TH2D* hden2   = new TH2D("hden2", "", binnum, bins, ybinnum, ybins);
            TH2D* hden3   = new TH2D("hden3", "", binnum, bins, ybinnum, ybins);
            TH2D* hden4   = new TH2D("hden4", "", binnum, bins, ybinnum, ybins);

            hall  ->Sumw2();
            hpass ->Sumw2();
            hfake ->Sumw2();
            
            hall1  ->Sumw2();
            hpass1 ->Sumw2();
            hfake1 ->Sumw2();
            
            hall2  ->Sumw2();
            hpass2 ->Sumw2();
            hfake2 ->Sumw2();
            
            hall3  ->Sumw2();
            hpass3 ->Sumw2();
            hfake3 ->Sumw2();
            
            hall4  ->Sumw2();
            hpass4 ->Sumw2();
            hfake4 ->Sumw2();
            
            TBranch *bmet     = tree->GetBranch("met");
            TBranch *bzmass   = tree->GetBranch("zmass");
            TBranch *bl1pt    = tree->GetBranch("l1pt");
            TBranch *bl1eta   = tree->GetBranch("l1eta");
            TBranch *bl1phi   = tree->GetBranch("l1phi");
            TBranch *bl2pt    = tree->GetBranch("l2pt");
            TBranch *bl2eta   = tree->GetBranch("l2eta");
            TBranch *bl2phi   = tree->GetBranch("l2phi");
            TBranch *bpt      = tree->GetBranch("pt");
            TBranch *beta     = tree->GetBranch("eta");
            TBranch *bphi     = tree->GetBranch("phi");
            TBranch *bid      = tree->GetBranch("newID");
            TBranch *biso     = dofsr ? tree->GetBranch("pfCombRelIso04EACorr_WithFSR") : tree->GetBranch("pfCombRelIso04EACorr");
            TBranch *brun     = tree->GetBranch("run");
            TBranch *blum     = tree->GetBranch("lumi");
            TBranch *bevt     = tree->GetBranch("event");
            TBranch *bnmiss   = tree->GetBranch("nmisshits");
           
            float met     = 0.0;
            float l1pt    = 0.0;
            float l1eta   = 0.0;
            float l1phi   = 0.0;
            float l2pt    = 0.0;
            float l2eta   = 0.0;
            float l2phi   = 0.0;
            float zmass   = 0.0;
            float pt      = 0.0;
            float eta     = 0.0;
            float phi     = 0.0;
            float iso     = 0.0;
            float nmiss   = 0.0;
            int   id      = 0;
            int   run     = 0;
            int   lum     = 0;
            int   evt     = 0;

            bmet      ->SetAddress(&met);
            bzmass    ->SetAddress(&zmass);
            bl1pt     ->SetAddress(&l1pt);
            bl1eta    ->SetAddress(&l1eta);
            bl1phi    ->SetAddress(&l1phi);
            bl2pt     ->SetAddress(&l2pt);
            bl2eta    ->SetAddress(&l2eta);
            bl2phi    ->SetAddress(&l2phi);
            bpt       ->SetAddress(&pt);
            beta      ->SetAddress(&eta);
            bphi      ->SetAddress(&phi);
            bid       ->SetAddress(&id);
            biso      ->SetAddress(&iso);
            bnmiss    ->SetAddress(&nmiss);
            brun      ->SetAddress(&run);
            blum      ->SetAddress(&lum);
            bevt      ->SetAddress(&evt);
 
            runeventinfo.clear(); 
            for (int i = 0; i < tree->GetEntries(); i++) {
                bmet     ->GetEvent(i);
                bzmass   ->GetEvent(i);
                bl1pt    ->GetEvent(i);
                bl1eta   ->GetEvent(i);
                bl1phi   ->GetEvent(i);
                bl2pt    ->GetEvent(i);
                bl2eta   ->GetEvent(i);
                bl2phi   ->GetEvent(i);
                bpt      ->GetEvent(i);
                beta     ->GetEvent(i);
                bphi     ->GetEvent(i);
                bid      ->GetEvent(i);
                biso     ->GetEvent(i);
                bnmiss   ->GetEvent(i);
                brun     ->GetEvent(i);
                blum     ->GetEvent(i);
                bevt     ->GetEvent(i);
            
                bool existsAlready = false;
                for (std::size_t k = 0; k < runeventinfo.size(); k++) {
                    if (run == runeventinfo[k].first && evt == runeventinfo[k].second) existsAlready = true;
                }
      
                if (existsAlready) continue;            
                runeventinfo.push_back(std::pair<int, int>(run, evt));
 
                float dR1 = reco::deltaR(l1eta, l1phi, eta, phi); 
                float dR2 = reco::deltaR(l2eta, l2phi, eta, phi); 
            
                float m1 = getMass(l1pt, l1eta, l1phi, pt, eta, phi); 
                float m2 = getMass(l2pt, l2eta, l2phi, pt, eta, phi); 
            
                if (zmass>zmin && zmass<zmax && dR1 > drcut && dR2 > drcut && m1>mcut && m2>mcut) {
                    hall->Fill(pt,fabs(eta));
                    if (id>0 && iso<0.4) {
                        hpass->Fill(pt, fabs(eta));
                    }
                }

                if (dR1 > drcut && dR2 > drcut && m1>mcut && m2>mcut) {
                    if (zmass>40 && zmass<60) {
                        hall1->Fill(pt, abs(eta));
                        havg1->Fill(pt, fabs(eta), nmiss);
                        hden1->Fill(pt, fabs(eta));
                        if (id>0 && iso<0.4) hpass1->Fill(pt, fabs(eta));
                    }
                    else if (zmass>60 && zmass<80) {
                        hall2->Fill(pt, abs(eta));
                        havg2->Fill(pt, fabs(eta), nmiss);
                        hden2->Fill(pt, fabs(eta));
                        if (id>0 && iso<0.4) hpass2->Fill(pt, fabs(eta));
                    }
                    else if ((zmass>80 && zmass<85) || zmass>95) {
                        hall3->Fill(pt, abs(eta));
                        havg3->Fill(pt, fabs(eta), nmiss);
                        hden3->Fill(pt, fabs(eta));
                        if (id>0 && iso<0.4) hpass3->Fill(pt, fabs(eta));
                    }
                    else {
                        hall4->Fill(pt, abs(eta));
                        havg4->Fill(pt, fabs(eta), nmiss);
                        hden4->Fill(pt, fabs(eta));
                        if (id>0 && iso<0.4) hpass4->Fill(pt, fabs(eta));
                    }

                }
           
                hden1->Fill(1, 0.0);
                hden2->Fill(1, 0.0);
                hden3->Fill(1, 0.0);
                hden4->Fill(1, 0.0);
                hden1->Fill(1, 1.5);
                hden2->Fill(1, 1.5);
                hden3->Fill(1, 1.5);
                hden4->Fill(1, 0.0);
                hden1->Fill(6, 0.0);
                hden2->Fill(6, 0.0);
                hden3->Fill(6, 0.0);
                hden4->Fill(6, 0.0);
                hden1->Fill(6, 1.5);
                hden2->Fill(6, 1.5);
                hden3->Fill(6, 1.5);
                hden4->Fill(6, 0.0);
 
            }
            
            havg1 ->Divide(hden1);
            havg2 ->Divide(hden2);
            havg3 ->Divide(hden3);
            havg4 ->Divide(hden4);
            hfake1->Divide(hpass1  , hall1);
            hfake2->Divide(hpass2  , hall2);
            hfake3->Divide(hpass3  , hall3);
            hfake4->Divide(hpass4  , hall4);
            hfake ->Divide(hpass   , hall);
            
            TF1* f = new TF1("f", "[0]*x+[1]", 0, 1);
            for (std::size_t k = 0; k < 10; k++) {
                elbarrel[k] = hfake->GetBinContent(k+1, 1);
                elendcap[k] = hfake->GetBinContent(k+1, 2);
                elbarrelerr[k] = hfake->GetBinError(k+1, 1);
                elendcaperr[k] = hfake->GetBinError(k+1, 2);

                float bx[4] = {havg1->GetBinContent(k+1, 1), havg2->GetBinContent(k+1, 1), havg3->GetBinContent(k+1, 1), havg4->GetBinContent(k+1, 1)};
                float bxerr[4] = {0.0, 0.0, 0.0, 0.0};
                float by[4] = {hfake1->GetBinContent(k+1, 1), hfake2->GetBinContent(k+1, 1), hfake3->GetBinContent(k+1, 1), hfake4->GetBinContent(k+1, 1)};
                float byerr[4] = {hfake1->GetBinError(k+1, 1), hfake2->GetBinError(k+1, 1), hfake3->GetBinError(k+1, 1), hfake4->GetBinError(k+1, 1)};

                float ex[4] = {havg1->GetBinContent(k+1, 2), havg2->GetBinContent(k+1, 2), havg3->GetBinContent(k+1, 2), havg4->GetBinContent(k+1, 2)};
                float exerr[4] = {0.0, 0.0, 0.0, 0.0};
                float ey[4] = {hfake1->GetBinContent(k+1, 2), hfake2->GetBinContent(k+1, 2), hfake3->GetBinContent(k+1, 2), hfake4->GetBinContent(k+1, 2)};
                float eyerr[4] = {hfake1->GetBinError(k+1, 2), hfake2->GetBinError(k+1, 2), hfake3->GetBinError(k+1, 2), hfake4->GetBinError(k+1, 2)};
               
                std::stringstream bss;
                std::stringstream ess;
                bss << "bgr" <<k;
                ess << "egr" <<k;
 
                TGraphErrors* bgr = new TGraphErrors(4, bx, by, bxerr, byerr);
                bgr->SetName(bss.str().c_str());
                bgr->Fit("f", "Q");
                float b1 = f->GetParameter(0);
                float b2 = f->GetParameter(1);

                TGraphErrors* egr = new TGraphErrors(4, ex, ey, exerr, eyerr);
                egr->SetName(ess.str().c_str());
                egr->Fit("f", "Q");
                float e1 = f->GetParameter(0);
                float e2 = f->GetParameter(1);                

                elbarrelnew[k] = std::pair<float, float>(b1, b2);
                elendcapnew[k] = std::pair<float, float>(e1, e2);

                elbarrelgr.push_back(bgr);
                elendcapgr.push_back(egr);
            }
            
            file->Close();
            delete file;
        }

    public:
        FakeRateCalculator():
            mubarrel(std::vector<float>(11, 0.0)),
            muendcap(std::vector<float>(11, 0.0)),
            elbarrel(std::vector<float>(11, 0.0)),
            elendcap(std::vector<float>(11, 0.0)),

            mubarrelerr(std::vector<float>(11, 0.0)),
            muendcaperr(std::vector<float>(11, 0.0)),
            elbarrelerr(std::vector<float>(11, 0.0)),
            elendcaperr(std::vector<float>(11, 0.0)),

            elbarrelnew(std::vector<std::pair<float, float> >(11, std::pair<float, float>(0.0, 0.0))),
            elendcapnew(std::vector<std::pair<float, float> >(11, std::pair<float, float>(0.0, 0.0))),

            elbarrelavg(std::vector<float>(11, 0.0)),
            elendcapavg(std::vector<float>(11, 0.0)),

            dofsr(false) 
        {
            Float_t bins[] = {0.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,80.0};
            Int_t binnum = 10;
            Float_t ybins[] = {0.0,1.479,2.5};
            Int_t ybinnum = 2;
            havg1   = new TH2D("havg1", "", binnum, bins, ybinnum, ybins);
            havg2   = new TH2D("havg2", "", binnum, bins, ybinnum, ybins);
            havg3   = new TH2D("havg3", "", binnum, bins, ybinnum, ybins);
            havg4   = new TH2D("havg4", "", binnum, bins, ybinnum, ybins);
        }

        FakeRateCalculator(std::string path, bool do7TeV, float zmin, float zmax, float drcut, float mcut, bool df):
            mubarrel(std::vector<float>(11, 0.0)),
            muendcap(std::vector<float>(11, 0.0)),
            elbarrel(std::vector<float>(11, 0.0)),
            elendcap(std::vector<float>(11, 0.0)),

            mubarrelerr(std::vector<float>(11, 0.0)),
            muendcaperr(std::vector<float>(11, 0.0)),
            elbarrelerr(std::vector<float>(11, 0.0)),
            elendcaperr(std::vector<float>(11, 0.0)),

            elbarrelnew(std::vector<std::pair<float, float> >(11, std::pair<float, float>(0.0, 0.0))),
            elendcapnew(std::vector<std::pair<float, float> >(11, std::pair<float, float>(0.0, 0.0))),

            elbarrelavg(std::vector<float>(11, 0.0)),
            elendcapavg(std::vector<float>(11, 0.0)),

            dofsr(df) 
        {
            init(do7TeV);
            Float_t bins[] = {0.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,80.0};
            Int_t binnum = 10;
            Float_t ybins[] = {0.0,1.479,2.5};
            Int_t ybinnum = 2;
            havg1   = new TH2D("havg1", "", binnum, bins, ybinnum, ybins);
            havg2   = new TH2D("havg2", "", binnum, bins, ybinnum, ybins);
            havg3   = new TH2D("havg3", "", binnum, bins, ybinnum, ybins);
            havg4   = new TH2D("havg4", "", binnum, bins, ybinnum, ybins);
            doelavg(path, 40., 12., 100.);
            domufakes(path, do7TeV, zmin, zmax, drcut, mcut);
            doelfakes(path, do7TeV, zmin, zmax, drcut, mcut);
        }

        float getPromptRate(float pt, float eta, float id) {
            return getPR(pt, eta, id);
        }

        float getFakeRate(float pt, float eta, float id) {
            int bin = 0;
          
            if (pt>0  && pt<= 5) bin = 0;
            if (pt>5  && pt<= 7) bin = 1;
            if (pt>7  && pt<=10) bin = 2;
            if (pt>10 && pt<=15) bin = 3;
            if (pt>15 && pt<=20) bin = 4;
            if (pt>20 && pt<=25) bin = 5;
            if (pt>25 && pt<=30) bin = 6;
            if (pt>30 && pt<=40) bin = 7;
            if (pt>40 && pt<=50) bin = 8;
            if (pt>50 && pt<=80) bin = 9;
            if (pt>80)           bin = 9;

            if      (fabs(id) == 13 && fabs(eta)<1.479) return mubarrel[bin];
            else if (fabs(id) == 13 && fabs(eta)>1.479) return muendcap[bin];
            else if (fabs(id) == 11 && fabs(eta)<1.479) return elbarrel[bin];
            else if (fabs(id) == 11 && fabs(eta)>1.479) return elendcap[bin];
            else return 0.0;

        }

        float getFakeRateNew(float pt, float eta, float id) {
            int bin = 0;

            if (pt>0  && pt<= 5) bin = 0;
            if (pt>5  && pt<= 7) bin = 1;
            if (pt>7  && pt<=10) bin = 2;
            if (pt>10 && pt<=15) bin = 3;
            if (pt>15 && pt<=20) bin = 4;
            if (pt>20 && pt<=25) bin = 5;
            if (pt>25 && pt<=30) bin = 6;
            if (pt>30 && pt<=40) bin = 7;
            if (pt>40 && pt<=50) bin = 8;
            if (pt>50 && pt<=80) bin = 9;
            if (pt>80)           bin = 9;

            if      (fabs(id) == 13 && fabs(eta)<1.479) return mubarrel[bin];
            else if (fabs(id) == 13 && fabs(eta)>1.479) return muendcap[bin];
            else if (fabs(id) == 11 && fabs(eta)<1.479) return elbarrelnew[bin].second + elbarrelnew[bin].first*elbarrelavg[bin];
            else if (fabs(id) == 11 && fabs(eta)>1.479) return elendcapnew[bin].second + elendcapnew[bin].first*elendcapavg[bin];
            else return 0.0;

        }

        float getFakeRateErr(float pt, float eta, float id) {
            int bin = 0;

            if (pt>0  && pt<= 5) bin = 0;
            if (pt>5  && pt<= 7) bin = 1;
            if (pt>7  && pt<=10) bin = 2;
            if (pt>10 && pt<=15) bin = 3;
            if (pt>15 && pt<=20) bin = 4;
            if (pt>20 && pt<=25) bin = 5;
            if (pt>25 && pt<=30) bin = 6;
            if (pt>30 && pt<=40) bin = 7;
            if (pt>40 && pt<=50) bin = 8;
            if (pt>50 && pt<=80) bin = 9;
            if (pt>80)           bin = 9;

            if      (fabs(id) == 13 && fabs(eta)<1.479) return mubarrelerr[bin];
            else if (fabs(id) == 13 && fabs(eta)>1.479) return muendcaperr[bin];
            else if (fabs(id) == 11 && fabs(eta)<1.479) return elbarrelerr[bin];
            else if (fabs(id) == 11 && fabs(eta)>1.479) return elendcaperr[bin];
            else return 0.0;
        }

        void printMuonFakeRate(std::ostream& out) {

            out << "---- Muon Barrel ----" << std::endl;
            out << "pT [ 0 -  5] : " << mubarrel[0]  << " +/- " << mubarrelerr[0]  << endl;
            out << "pT [ 5 -  7] : " << mubarrel[1]  << " +/- " << mubarrelerr[1]  << endl;
            out << "pT [ 7 - 10] : " << mubarrel[2]  << " +/- " << mubarrelerr[2]  << endl;
            out << "pT [10 - 15] : " << mubarrel[3]  << " +/- " << mubarrelerr[3]  << endl;
            out << "pT [15 - 20] : " << mubarrel[4]  << " +/- " << mubarrelerr[4]  << endl;
            out << "pT [20 - 25] : " << mubarrel[5]  << " +/- " << mubarrelerr[5]  << endl;
            out << "pT [25 - 30] : " << mubarrel[6]  << " +/- " << mubarrelerr[6]  << endl;
            out << "pT [30 - 40] : " << mubarrel[7]  << " +/- " << mubarrelerr[7]  << endl;
            out << "pT [40 - 50] : " << mubarrel[8]  << " +/- " << mubarrelerr[8]  << endl;
            out << "pT [50 - 80] : " << mubarrel[9]  << " +/- " << mubarrelerr[9]  << endl;
            out << std::endl;
            out << "---- Muon Endcap ----" << std::endl;
            out << "pT [ 0 -  5] : " << muendcap[0]  << " +/- " << muendcaperr[0]  << endl;
            out << "pT [ 5 -  7] : " << muendcap[1]  << " +/- " << muendcaperr[1]  << endl;
            out << "pT [ 7 - 10] : " << muendcap[2]  << " +/- " << muendcaperr[2]  << endl;
            out << "pT [10 - 15] : " << muendcap[3]  << " +/- " << muendcaperr[3]  << endl;
            out << "pT [15 - 20] : " << muendcap[4]  << " +/- " << muendcaperr[4]  << endl;
            out << "pT [20 - 25] : " << muendcap[5]  << " +/- " << muendcaperr[5]  << endl;
            out << "pT [25 - 30] : " << muendcap[6]  << " +/- " << muendcaperr[6]  << endl;
            out << "pT [30 - 40] : " << muendcap[7]  << " +/- " << muendcaperr[7]  << endl;
            out << "pT [40 - 50] : " << muendcap[8]  << " +/- " << muendcaperr[8]  << endl;
            out << "pT [50 - 80] : " << muendcap[9]  << " +/- " << muendcaperr[9]  << endl;
            out << std::endl;


        }

        void printElectronFakeRate(std::ostream& out) {

            out << "---- Electron Barrel ----" << std::endl;
            out << "pT [ 0 -  5] : " << elbarrel[0]  << " +/- " << elbarrelerr[0]  << endl;
            out << "pT [ 5 -  7] : " << elbarrel[1]  << " +/- " << elbarrelerr[1]  << endl;
            out << "pT [ 7 - 10] : " << elbarrel[2]  << " +/- " << elbarrelerr[2]  << endl;
            out << "pT [10 - 15] : " << elbarrel[3]  << " +/- " << elbarrelerr[3]  << endl;
            out << "pT [15 - 20] : " << elbarrel[4]  << " +/- " << elbarrelerr[4]  << endl;
            out << "pT [20 - 25] : " << elbarrel[5]  << " +/- " << elbarrelerr[5]  << endl;
            out << "pT [25 - 30] : " << elbarrel[6]  << " +/- " << elbarrelerr[6]  << endl;
            out << "pT [30 - 40] : " << elbarrel[7]  << " +/- " << elbarrelerr[7]  << endl;
            out << "pT [40 - 50] : " << elbarrel[8]  << " +/- " << elbarrelerr[8]  << endl;
            out << "pT [50 - 80] : " << elbarrel[9]  << " +/- " << elbarrelerr[9]  << endl;
            out << std::endl;
            out << "---- Electron Endcap ----" << std::endl;
            out << "pT [ 0 -  5] : " << elendcap[0]  << " +/- " << elendcaperr[0]  << endl;
            out << "pT [ 5 -  7] : " << elendcap[1]  << " +/- " << elendcaperr[1]  << endl;
            out << "pT [ 7 - 10] : " << elendcap[2]  << " +/- " << elendcaperr[2]  << endl;
            out << "pT [10 - 15] : " << elendcap[3]  << " +/- " << elendcaperr[3]  << endl;
            out << "pT [15 - 20] : " << elendcap[4]  << " +/- " << elendcaperr[4]  << endl;
            out << "pT [20 - 25] : " << elendcap[5]  << " +/- " << elendcaperr[5]  << endl;
            out << "pT [25 - 30] : " << elendcap[6]  << " +/- " << elendcaperr[6]  << endl;
            out << "pT [30 - 40] : " << elendcap[7]  << " +/- " << elendcaperr[7]  << endl;
            out << "pT [40 - 50] : " << elendcap[8]  << " +/- " << elendcaperr[8]  << endl;
            out << "pT [50 - 80] : " << elendcap[9]  << " +/- " << elendcaperr[9]  << endl;
            out << std::endl;


        }

        void printGraphs() { 

            TCanvas* c = new TCanvas("c", "");
            for (std::size_t i = 2; i < elbarrelgr.size(); i++) {
                elbarrelgr[i]->Draw("AP");
                c->Print((std::string(elbarrelgr[i]->GetName())+".pdf").c_str());
            }
            for (std::size_t i = 2; i < elendcapgr.size(); i++) {
                elendcapgr[i]->Draw("AP");
                c->Print((std::string(elendcapgr[i]->GetName())+".pdf").c_str());
            }
        }


};

#endif 
