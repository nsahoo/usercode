#define estimateElectronFakeRateHzz4lTree_cxx
#include "include/estimateElectronFakeRateHzz4lTree.h"
#include "include/EfficiencyEvaluator.hh"
#include "include/eIDHzzCuts.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>

using namespace std;

int estimateElectronFakeRateHzz4lTree::etabin(float eta) {
  if(fabs(eta)<0.8) return 0;
  else if(fabs(eta)<1.479) return 1;
  else return 2;
}

int estimateElectronFakeRateHzz4lTree::ptbinTrg(float pt) {
  if(pt<20) return 0;
  else return 1;
}

int estimateElectronFakeRateHzz4lTree::ptbinNoTrg(float pt) {
  if(pt<10) return 0;
  else return 1;
}

Bool_t estimateElectronFakeRateHzz4lTree::passRefEleSel() {
  
  float abseta=fabs(eta);

  ElectronEffectiveArea::ElectronEffectiveAreaTarget effAreaTarget_ = ElectronEffectiveArea::kEleEAData2011;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaGamma_   = ElectronEffectiveArea::kEleGammaIso04;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaNeutralHad_ = ElectronEffectiveArea::kEleNeutralHadronIso04;
  ElectronEffectiveArea::ElectronEffectiveAreaType effAreaGammaAndNeutralHad_ =  ElectronEffectiveArea::kEleGammaAndNeutralHadronIso04;

  float eff_area_ga  = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGamma_, abseta, effAreaTarget_);
  float eff_area_nh  = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaNeutralHad_, abseta, effAreaTarget_);
  float eff_area_ganh = ElectronEffectiveArea::GetElectronEffectiveArea(effAreaGammaAndNeutralHad_, abseta, effAreaTarget_);

  float iso = pfIsoChHad04;
  iso += max<float>(0.,pfIsoNHad04_NoEA+pfIsoPhoton04_NoEA - eff_area_ganh*rhoAA);

  int etab=etabin(eta);
  int ptb=ptbinNoTrg(pt);
  float bdtcut=cutrefbdt_hzz[ptb][etab];
  return (bdtIdNtDz>bdtcut && iso/pt<0.25 && fabs(sip3d)<4.0 && nmisshits<=1);
} 

Bool_t estimateElectronFakeRateHzz4lTree::passOptimizedEleSel(int wp) {

  int etab=etabin(eta);
  int ptb=ptbinNoTrg(pt);
  float bdtcut=999.;
  float isocut=-999.;

  // optimized in 2D space
  bdtcut=cutbdt_hzz[ptb][etab][wp];
  isocut=cutmvaiso_hzz[ptb][etab][wp];
  
  return (bdtIdNtDz>bdtcut && bdtIsoNtDz>isocut && fabs(sip3d)<4.0 && nmisshits<=1);
} 

void estimateElectronFakeRateHzz4lTree::Loop(const char *outname)
{
  if (fChain == 0) return;

  std::vector<TString> NoTrgElectronID;
  NoTrgElectronID.push_back("hzzPfIso");  // reference pf iso
  NoTrgElectronID.push_back("hzzMvaLoose");  // duncan's optimization for HZZ
  NoTrgElectronID.push_back("hzzMvaTight");  // duncan's optimization for HZZ

  // -----------------------------------------------------------------------
  // study vs eta

  Float_t LowerEta[5];
  LowerEta[0]=0.0;
  LowerEta[1]=1.0;
  LowerEta[2]=1.479;  
  LowerEta[3]=2.0;
  LowerEta[4]=2.5;
  TH1F *RecoEtaHighPt   = new TH1F( "RecoEtaHighPt",    "reconstructed #eta", 4, LowerEta);
  TH1F *RecoEtaLowPt    = new TH1F( "RecoEtaLowPt",     "reconstructed #eta", 4, LowerEta);
  
  // eta, high pT
  std::vector<TH1F*> NoTrgElectronEtaHighPt;
  for (int i=0;i<(int)NoTrgElectronID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgElectron"+TString(NoTrgElectronID[i])+"EtaHighPt",   "HZZ BDT ID #eta", 4, LowerEta);     NoTrgElectronEtaHighPt.push_back(aHisto);
  }

  // eta, low pT
  std::vector<TH1F*> NoTrgElectronEtaLowPt;
  for (int i=0;i<(int)NoTrgElectronID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgElectron"+TString(NoTrgElectronID[i])+"EtaLowPt",   "HZZ BDT ID #eta", 4, LowerEta);
    NoTrgElectronEtaLowPt.push_back(aHisto);
  }

  // -----------------------------------------------------------------------
  // study vs pT
  Float_t LowerPt[11] = {0.0,5.0,7.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,80.0};
  TH1F *RecoPtBarrel   = new TH1F( "RecoPtBarrel",    "reconstructed p_{T} (GeV)", 10, LowerPt);
  TH1F *RecoPtEndcap   = new TH1F( "RecoPtEndcap",    "reconstructed p_{T} (GeV)", 10, LowerPt);

  // to have the full picture in the barrel
  std::vector<TH1F*> NoTrgElectronPtBarrel;
  for (int i=0;i<(int)NoTrgElectronID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgElectron"+TString(NoTrgElectronID[i])+"PtBarrel", "HZZ BDT ID #eta",   10, LowerPt );
    NoTrgElectronPtBarrel.push_back(aHisto);
  }

  // to have the full picture in the endcap
  std::vector<TH1F*> NoTrgElectronPtEndcap;
  for (int i=0;i<(int)NoTrgElectronID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgElectron"+TString(NoTrgElectronID[i])+"PtEndcap", "HZZ BDT ID #eta",   10, LowerPt);
    NoTrgElectronPtEndcap.push_back(aHisto);
  }

  // -----------------------------------------------------------------------
  // study vs PU
  Float_t LowerPU[11];
  LowerPU[0] = 1;
  LowerPU[1] = 4;
  for(int i=2;i<6;++i) LowerPU[i]=i+3;
  LowerPU[6] = 10;
  LowerPU[7] = 15;
  LowerPU[8] = 25;
  LowerPU[9] = 30;
  LowerPU[10] = 35;

  TH1F *RecoPUBarrel   = new TH1F( "RecoPUBarrel",   "reconstructed nPU", 10, LowerPU);
  TH1F *RecoPUEndcap   = new TH1F( "RecoPUEndcap",   "reconstructed nPU", 10, LowerPU);

  //  barrel
  std::vector<TH1F*> NoTrgElectronPUBarrel;
  for (int i=0;i<(int)NoTrgElectronID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgElectron"+TString(NoTrgElectronID[i])+"PUBarrel", "HZZ BDT ID #eta",   10, LowerPU );
    NoTrgElectronPUBarrel.push_back(aHisto);
  }

  // endcap
  std::vector<TH1F*> NoTrgElectronPUEndcap;
  for (int i=0;i<(int)NoTrgElectronID.size();++i) {
    TH1F* aHisto = new TH1F( "NoTrgElectron"+TString(NoTrgElectronID[i])+"PUEndcap", "HZZ BDT ID #eta",   10, LowerPU );
    NoTrgElectronPUEndcap.push_back(aHisto);
  }

  // loop on events
  Long64_t nbytes = 0, nb = 0;
  Long64_t nentries = fChain->GetEntries();
  std::cout << "Number of entries = " << nentries << std::endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (jentry%100000 == 0) std::cout << ">>> Processing event # " << jentry << std::endl;

    // fill the denominator: take only the highest pT denominator candidate
    float etaFake = fabs(eta);
    float etFake  = pt;
    bool isInEB   = fabs(eta)<1.479;
    bool isInEE   = !isInEB;
    bool highPt   = (pt>10.);
    bool lowPt    = (pt<=10.);

    // pass the denominator object
    if(etFake<5 || etaFake>2.5 || fabs(sip3d)>4 || nmisshits>1) continue;

    // filling
    if (highPt) RecoEtaHighPt -> Fill(etaFake);  //, theWeight); 
    if (lowPt)  RecoEtaLowPt  -> Fill(etaFake);  // , theWeight);

    if (isInEB) {
      RecoPtBarrel -> Fill(etFake);  //, theWeight);
      RecoPUBarrel -> Fill(numvertices);  //, theWeight);
    }
    if (isInEE) {
      RecoPtEndcap -> Fill(etFake);  //, theWeight);
      RecoPUEndcap -> Fill(numvertices);  //, theWeight);
    }


    // fill the numerator(s)
    // === MVA ele-ID, combined Iso ===
    if(passRefEleSel()) {
      if (highPt) NoTrgElectronEtaHighPt[0]->Fill(etaFake);
      if (lowPt)  NoTrgElectronEtaLowPt[0] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElectronPtBarrel[0] ->Fill(etFake);
	NoTrgElectronPUBarrel[0] ->Fill(numvertices);
      }
      if (isInEE) {
	NoTrgElectronPtEndcap[0] ->Fill(etFake);
	NoTrgElectronPUEndcap[0] ->Fill(numvertices);
      }
    }

    // === MVA ele-ID, MVA Iso, Loose ===
    if(passOptimizedEleSel(kLoose)) {
      if (highPt) NoTrgElectronEtaHighPt[1]->Fill(etaFake);
      if (lowPt)  NoTrgElectronEtaLowPt[1] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElectronPtBarrel[1] ->Fill(etFake);
	NoTrgElectronPUBarrel[1] ->Fill(numvertices);
      }
      if (isInEE) {
	NoTrgElectronPtEndcap[1] ->Fill(etFake);
	NoTrgElectronPUEndcap[1] ->Fill(numvertices);
      }
    }

    // === MVA ele-ID, MVA Iso, Tight ===
    if(passOptimizedEleSel(kTight)) {
      if (highPt) NoTrgElectronEtaHighPt[2]->Fill(etaFake);
      if (lowPt)  NoTrgElectronEtaLowPt[2] ->Fill(etaFake);
      if (isInEB) { 
	NoTrgElectronPtBarrel[2] ->Fill(etFake);
	NoTrgElectronPUBarrel[2] ->Fill(numvertices);
      }
      if (isInEE) {
	NoTrgElectronPtEndcap[2] ->Fill(etFake);
	NoTrgElectronPUEndcap[2] ->Fill(numvertices);
      }
    }

  } // loop over events

    // saving efficiency histos
  // === as a function of eta ===
  char filename[500];
  sprintf(filename,"%s-ElectronMisidEtaHighPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaHighPt(filename);
  ElectronEffEtaHighPt.AddNumerator(RecoEtaHighPt);
  for (int icut=0;icut<(int)NoTrgElectronID.size();++icut){
    ElectronEffEtaHighPt.AddNumerator(NoTrgElectronEtaHighPt[icut]);
  }

  ElectronEffEtaHighPt.SetDenominator(RecoEtaHighPt);
  ElectronEffEtaHighPt.ComputeEfficiencies();
  ElectronEffEtaHighPt.SetTitle("fake rate vs pT");
  ElectronEffEtaHighPt.SetXaxisTitle("electron p_{T}");
  ElectronEffEtaHighPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaHighPt.SetYaxisMin(0.0);
  ElectronEffEtaHighPt.Write();

  sprintf(filename,"%s-ElectronMisidEtaLowPt.root",outname);
  EfficiencyEvaluator ElectronEffEtaLowPt(filename);
  ElectronEffEtaLowPt.AddNumerator(RecoEtaLowPt);
  for (int icut=0;icut<(int)NoTrgElectronID.size();++icut){
    ElectronEffEtaLowPt.AddNumerator(NoTrgElectronEtaLowPt[icut]);
  }

  ElectronEffEtaLowPt.SetDenominator(RecoEtaLowPt);
  ElectronEffEtaLowPt.ComputeEfficiencies();
  ElectronEffEtaLowPt.SetTitle("fake rate vs eta");
  ElectronEffEtaLowPt.SetXaxisTitle("electron #eta");
  ElectronEffEtaLowPt.SetYaxisTitle("Fake rate");
  ElectronEffEtaLowPt.SetYaxisMin(0.0);
  ElectronEffEtaLowPt.Write();

  // === as a function of pt ===
  sprintf(filename,"%s-ElectronMisidPtBarrel.root",outname);
  EfficiencyEvaluator ElectronEffPtBarrel(filename);
  ElectronEffPtBarrel.AddNumerator(RecoPtBarrel);
  for (int icut=0;icut<(int)NoTrgElectronID.size();++icut){
    ElectronEffPtBarrel.AddNumerator(NoTrgElectronPtBarrel[icut]);
  }

  ElectronEffPtBarrel.SetDenominator(RecoPtBarrel);
  ElectronEffPtBarrel.ComputeEfficiencies();
  ElectronEffPtBarrel.SetTitle("fake rate vs pT");
  ElectronEffPtBarrel.SetXaxisTitle("electron pT");
  ElectronEffPtBarrel.SetYaxisTitle("Fake rate");
  ElectronEffPtBarrel.SetYaxisMin(0.0);
  ElectronEffPtBarrel.Write();

  sprintf(filename,"%s-ElectronMisidPtEndcap.root",outname);
  EfficiencyEvaluator ElectronEffPtEndcap(filename);
  ElectronEffPtEndcap.AddNumerator(RecoPtEndcap);
  for (int icut=0;icut<(int)NoTrgElectronID.size();++icut){
    ElectronEffPtEndcap.AddNumerator(NoTrgElectronPtEndcap[icut]);
  }

  ElectronEffPtEndcap.SetDenominator(RecoPtEndcap);
  ElectronEffPtEndcap.ComputeEfficiencies();
  ElectronEffPtEndcap.SetTitle("fake rate vs pT");
  ElectronEffPtEndcap.SetXaxisTitle("electron pT");
  ElectronEffPtEndcap.SetYaxisTitle("Fake rate");
  ElectronEffPtEndcap.SetYaxisMin(0.0);
  ElectronEffPtEndcap.Write();

  // === as a function of PU ===
  sprintf(filename,"%s-ElectronMisidPUBarrel.root",outname);
  EfficiencyEvaluator ElectronEffPUBarrel(filename);
  ElectronEffPUBarrel.AddNumerator(RecoPUBarrel);
  for (int icut=0;icut<(int)NoTrgElectronID.size();++icut){
    ElectronEffPUBarrel.AddNumerator(NoTrgElectronPUBarrel[icut]);
  }

  ElectronEffPUBarrel.SetDenominator(RecoPUBarrel);
  ElectronEffPUBarrel.ComputeEfficiencies();
  ElectronEffPUBarrel.SetTitle("fake rate vs vertices");
  ElectronEffPUBarrel.SetXaxisTitle("# vertices");
  ElectronEffPUBarrel.SetYaxisTitle("Fake rate");
  ElectronEffPUBarrel.SetYaxisMin(0.0);
  ElectronEffPUBarrel.Write();

  sprintf(filename,"%s-ElectronMisidPUEndcap.root",outname);
  EfficiencyEvaluator ElectronEffPUEndcap(filename);
  ElectronEffPUEndcap.AddNumerator(RecoPUEndcap);
  for (int icut=0;icut<(int)NoTrgElectronID.size();++icut){
    ElectronEffPUEndcap.AddNumerator(NoTrgElectronPUEndcap[icut]);
  }

  ElectronEffPUEndcap.SetDenominator(RecoPUEndcap);
  ElectronEffPUEndcap.ComputeEfficiencies();
  ElectronEffPUEndcap.SetTitle("fake rate vs vertices");
  ElectronEffPUEndcap.SetXaxisTitle("# vertices");
  ElectronEffPUEndcap.SetYaxisTitle("Fake rate");
  ElectronEffPUEndcap.SetYaxisMin(0.0);
  ElectronEffPUEndcap.Write();

}
