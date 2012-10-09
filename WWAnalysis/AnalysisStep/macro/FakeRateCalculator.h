#ifndef FAKERATECALCULATOR_H
#define FAKERATECALCULATOR_H

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
#include "DataFormats/Math/interface/deltaR.h"


class FakeRateCalculator {

    private:
        std::vector<float> mubarrel;
        std::vector<float> muendcap;
        std::vector<float> elbarrel;
        std::vector<float> elendcap;

        std::vector<float> mubarrelerr;
        std::vector<float> muendcaperr;
        std::vector<float> elbarrelerr;
        std::vector<float> elendcaperr;

        bool dofsr;

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

        void domufakes(std::string path, bool is2011, float zmin, float zmax, float drcut, float mcut) {
            TFile* file = new TFile(path.c_str());
            //TTree* tree = (TTree*)file->Get("zllmtreeNoOR/probe_tree");
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


        void doelfakes(std::string path, bool is2011, float zmin, float zmax, float drcut, float mcut) {
            TFile* file = new TFile(path.c_str());
            //TTree* tree = (TTree*)file->Get("zlletreeNoOR/probe_tree");
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
            
            hfake ->Divide(hpass , hall);
            
            for (std::size_t k = 0; k < 10; k++) {
                elbarrel[k] = hfake->GetBinContent(k+1, 1);
                elendcap[k] = hfake->GetBinContent(k+1, 2);
                elbarrelerr[k] = hfake->GetBinError(k+1, 1);
                elendcaperr[k] = hfake->GetBinError(k+1, 2);
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

            dofsr(false) 
        {}

        FakeRateCalculator(std::string path, bool is2011, float zmin, float zmax, float drcut, float mcut, bool df):
            mubarrel(std::vector<float>(11, 0.0)),
            muendcap(std::vector<float>(11, 0.0)),
            elbarrel(std::vector<float>(11, 0.0)),
            elendcap(std::vector<float>(11, 0.0)),

            mubarrelerr(std::vector<float>(11, 0.0)),
            muendcaperr(std::vector<float>(11, 0.0)),
            elbarrelerr(std::vector<float>(11, 0.0)),
            elendcaperr(std::vector<float>(11, 0.0)),

            dofsr(df) 
        {
            domufakes(path, is2011, zmin, zmax, drcut, mcut);
            doelfakes(path, is2011, zmin, zmax, drcut, mcut);
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

};

#endif 
