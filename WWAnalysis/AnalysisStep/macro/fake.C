#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLorentzVector.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>
#include <TF1.h>
#include <string>
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "DataFormats/Math/interface/deltaR.h"

void prep(TH1* hist1, TH1* hist2, const char* filename) {
    hist1->GetXaxis()->SetLabelSize(0.035);
    hist1->GetYaxis()->SetLabelSize(0.035);
    hist1->GetXaxis()->SetTitleSize(0.035);
    hist1->GetYaxis()->SetTitleSize(0.035);
    hist1->GetXaxis()->SetTitleOffset(1.2);
    hist1->GetYaxis()->SetTitleOffset(1.8);
    hist1->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    hist1->GetYaxis()->SetTitle("Fake Rate (Reco#rightarrow ID+Isolation)");
    hist1->SetMarkerStyle(20);
    hist1->SetMarkerColor(4);
    hist1->SetMarkerSize(0.7);
    hist1->SetLineColor(4);
    hist1->SetLineWidth(1);
    hist1->SetStats(0);
    hist1->GetYaxis()->SetRangeUser(0.0,0.4);

    hist2->GetXaxis()->SetLabelSize(0.035);
    hist2->GetYaxis()->SetLabelSize(0.035);
    hist2->GetXaxis()->SetTitleSize(0.035);
    hist2->GetYaxis()->SetTitleSize(0.035);
    hist2->GetXaxis()->SetTitleOffset(1.2);
    hist2->GetYaxis()->SetTitleOffset(1.8);
    hist2->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    hist2->GetYaxis()->SetTitle("Fake Rate (Reco#rightarrow ID+Isolation)");
    hist2->SetMarkerStyle(20);
    hist2->SetMarkerColor(2);
    hist2->SetMarkerSize(0.7);
    hist2->SetLineColor(2);
    hist2->SetLineWidth(1);
    hist2->SetStats(0);
    hist2->GetYaxis()->SetRangeUser(0.0,0.4);
    
    TCanvas* c = new TCanvas("canvas", "", 500, 500);

    TLegend* leg = new TLegend(0.6, 0.8, 0.7, 0.9);
    leg->SetFillColor(0);
    leg->AddEntry(hist1, "Barrel");
    leg->AddEntry(hist2, "Endcap");

    hist1->Draw("PE1");
    hist2->Draw("PE1 SAME");
    leg->Draw("SAME");

    c->Print(filename);
}

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


std::pair<TH1*, TH1*> dofakes(std::string path, bool isMu, bool dofsr, float zmin, float zmax) {
    TFile* file = new TFile(path.c_str());
    TTree* tree = isMu ? (TTree*)file->Get("zllmtree/probe_tree") : (TTree*)file->Get("zlletree/probe_tree") ;

    Float_t bins[] = {0.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,80.0};
    Int_t binnum = 10;

    TH1D* hfake1 = new TH1D("hfake1", "", binnum, bins);
    TH1D* hall1  = new TH1D("hall1",  "", binnum, bins);
    TH1D* hpass1 = new TH1D("hpass1", "", binnum, bins);
    TH1D* hfake2 = new TH1D("hfake2", "", binnum, bins);
    TH1D* hall2  = new TH1D("hall2",  "", binnum, bins);
    TH1D* hpass2 = new TH1D("hpass2", "", binnum, bins);

    hall1  ->Sumw2();
    hpass1 ->Sumw2();
    hfake1 ->Sumw2();
    hall2  ->Sumw2();
    hpass2 ->Sumw2();
    hfake2 ->Sumw2();

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

        if (zmass>zmin && zmass<zmax) {
            if (fabs(eta)<1.479) hall1->Fill(pt);
            else                 hall2->Fill(pt);
            if (id>0 && iso<0.4) {
                if (fabs(eta)<1.479) hpass1->Fill(pt);
                else                 hpass2->Fill(pt);
            }
        }

    }

    hfake1 ->Divide(hpass1,  hall1);
    hfake2 ->Divide(hpass2,  hall2);
    return std::pair<TH1*, TH1*>(hfake1, hfake2);
}


void fake() {

    std::pair<TH1*, TH1*> hists = dofakes("/home/avartak/CMS/Higgs/HZZ4L/CMSSW_5_2_4_patch4/src/WWAnalysis/AnalysisStep/trees/7TeV/hzzTree.root",  false, true,  40., 120.);

    prep(hists.first, hists.second, "elfakes7TeV.pdf");    

}
