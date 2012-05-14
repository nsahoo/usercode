#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TPaveText.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <TF1.h>

#define YEAR 2011

enum {
  kMuRef = 0,
  kMuMvaIso
};

enum {
  kEleRef = 0,
  kEleMvaLoose,
  kEleMvaTight
};

float computeSF(float pt, float abseta, bool isMu, int elwp, int muwp) {

  // 2012 fake rates parameterizations (PromptReco)
  float mu2012_bar_par0[2] = { 0.088146, 0.172263};
  float mu2012_bar_par1[2] = { 0.813482, 0.267368};
  float mu2012_bar_par2[2] = { 0.100546, 0.099457};
  float mu2012_bar_par3[2] = { 0.078343, 0.026675};
    
  float mu2012_end_par0[2] = { 0.097366, 0.233072 };
  float mu2012_end_par1[2] = { 0.829411, 0.332927 };
  float mu2012_end_par2[2] = { 0.10532,  0.126052 };
  float mu2012_end_par3[2] = { 0.081078, 0.020154 };
    
  float el2012_bar_par0[3] = { 0.0545295, 0.0736869, 0.0560685 };
  float el2012_bar_par1[3] = { 1.04825,   1.17273,   0.77083   };
  float el2012_bar_par2[3] = { 0.214458,  0.18919,   0.149952  };
  float el2012_bar_par3[3] = { 0.194721,  0.165607,  0.131008  };
  
  float el2012_end_par0[3] = { 0.172428,  0.192236,  0.18395  };
  float el2012_end_par1[3] = { 2.59677,   1.6312,    0.592073 };
  float el2012_end_par2[3] = { 0.26902,   0.252914,  0.236726 };
  float el2012_end_par3[3] = { 0.217469,  0.175229,  0.100146 };

  // 2011 fake rates parameterizations (Jan16 rereco)
  float mu2011_bar_par0[2] = { 0.0442711, 0.0924824 };
  float mu2011_bar_par1[2] = { 1.6023,    0.0 };
  float mu2011_bar_par2[2] = { 0.398586,  0.101534 };
  float mu2011_bar_par3[2] = { 0.386322,  0.0943444 };
  
  float mu2011_end_par0[2] = { 0.0364204, 0.13092 };
  float mu2011_end_par1[2] = { 0.562334,  1.40327 };
  float mu2011_end_par2[2] = { 0.171554,  0.0345104 };
  float mu2011_end_par3[2] = { 0.167615,  0.0291992 };
  
  float el2011_bar_par0[3] = { 0.0786844, 0.0868894, 0.072289 };
  float el2011_bar_par1[3] = { 1.84035,   1.3331,    1.74931  };
  float el2011_bar_par2[3] = { 0.251348,  0.25139,   0.208941 };
  float el2011_bar_par3[3] = { 0.224511,  0.211936,  0.189535 };
 
  float el2011_end_par0[3] = { 0.207228,  0.221917, 0.173104 };
  float el2011_end_par1[3] = { 3.79389,   0.75421,  1.95527  };
  float el2011_end_par2[3] = { 0.326206,  1.76081,  0.229775 };
  float el2011_end_par3[3] = { 0.261076,  0.155266, 0.176875 };

  if (YEAR == 2012) {
    if (!isMu && abseta<1.479)   return el2012_bar_par0[elwp] + el2012_bar_par1[elwp]*(exp(-el2012_bar_par2[elwp]*pt) - exp(-el2012_bar_par3[elwp]*pt));
    if (!isMu && abseta>1.479)   return el2012_end_par0[elwp] + el2012_end_par1[elwp]*(exp(-el2012_end_par2[elwp]*pt) - exp(-el2012_end_par3[elwp]*pt));
    
    if (isMu && abseta<1.479)    return mu2012_bar_par0[muwp] + mu2012_bar_par1[muwp]*(exp(-mu2012_bar_par2[muwp]*pt) - exp(-mu2012_bar_par3[muwp]*pt));
    if (isMu && abseta>1.479)    return mu2012_end_par0[muwp] + mu2012_end_par1[muwp]*(exp(-mu2012_end_par2[muwp]*pt) - exp(-mu2012_end_par3[muwp]*pt));
  } else { 
    if (!isMu && abseta<1.479)   return el2011_bar_par0[elwp] + el2011_bar_par1[elwp]*(exp(-el2011_bar_par2[elwp]*pt) - exp(-el2011_bar_par3[elwp]*pt));
    if (!isMu && abseta>1.479)   return el2011_end_par0[elwp] + el2011_end_par1[elwp]*(exp(-el2011_end_par2[elwp]*pt) - exp(-el2011_end_par3[elwp]*pt));
    
    if (isMu && abseta<1.479)    return mu2011_bar_par0[muwp] + mu2011_bar_par1[muwp]*(exp(-mu2011_bar_par2[muwp]*pt) - exp(-mu2011_bar_par3[muwp]*pt));
    if (isMu && abseta>1.479)    return mu2011_end_par0[muwp] + mu2011_end_par1[muwp]*(exp(-mu2011_end_par2[muwp]*pt) - exp(-mu2011_end_par3[muwp]*pt));
  }

  return 0;
}

float dozx(const char* path, int ch, int elwp, int muwp) {
    TFile* file = new TFile(path);
    TTree* tree = (TTree*)file->Get("zxTree/probe_tree");

    TBranch *bl1pt      = tree->GetBranch("l3pt");
    TBranch *bl1eta     = tree->GetBranch("l3eta");
    TBranch *bl2pt      = tree->GetBranch("l4pt");
    TBranch *bl2eta     = tree->GetBranch("l4eta");

    TBranch *bmass      = tree->GetBranch("mass");
    TBranch *bz1mass    = tree->GetBranch("z1mass");
    TBranch *bz2mass    = tree->GetBranch("z2mass");
    TBranch *bchannel   = tree->GetBranch("channel");
    TBranch *bl3pdgId   = tree->GetBranch("l3pdgId");
    TBranch *bl4pdgId   = tree->GetBranch("l4pdgId");
    TBranch *bl3idNew   = tree->GetBranch("l3idNew");
    TBranch *bl4idNew   = tree->GetBranch("l4idNew");

    TBranch *bevent     = tree->GetBranch("event");
    TBranch *brun       = tree->GetBranch("run");

    float l1pt      = 0.0;
    float l1eta     = 0.0;
    float l2pt      = 0.0;
    float l2eta     = 0.0;

    float mass      = 0.0;
    float z1mass    = 0.0;
    float z2mass    = 0.0;
    float channel   = 0.0;
    float l3pdgId   = 0.0;
    float l4pdgId   = 0.0;
    float l3idNew   = 0.0;
    float l4idNew   = 0.0;

    int   event     = 0.0;
    int   run       = 0.0;

    bl1pt   ->SetAddress(&l1pt);
    bl1eta  ->SetAddress(&l1eta);
    bl2pt   ->SetAddress(&l2pt);
    bl2eta  ->SetAddress(&l2eta);

    bmass   ->SetAddress(&mass);
    bz1mass ->SetAddress(&z1mass);
    bz2mass ->SetAddress(&z2mass);
    bchannel->SetAddress(&channel);
    bl3pdgId ->SetAddress(&l3pdgId);
    bl4pdgId ->SetAddress(&l4pdgId);
    bl3idNew ->SetAddress(&l3idNew);
    bl4idNew ->SetAddress(&l4idNew);

    bevent  ->SetAddress(&event);
    brun    ->SetAddress(&run);

    float yield = 0.0;
    for (int i = 0; i < tree->GetEntries(); i++) {
        bl1pt   ->GetEvent(i);
        bl1eta  ->GetEvent(i);
        bl2pt   ->GetEvent(i);
        bl2eta  ->GetEvent(i);

        bmass   ->GetEvent(i);
        bz1mass ->GetEvent(i);
        bz2mass ->GetEvent(i);
        bchannel->GetEvent(i);
	bl3pdgId ->GetEvent(i);
	bl4pdgId ->GetEvent(i);
	bl3idNew ->GetEvent(i);
	bl4idNew ->GetEvent(i);

        bevent  ->GetEvent(i);
        brun    ->GetEvent(i);

        // float z1min = 40;
        // float z2min = 4;

	float z1min = 50;
	float z2min = 12;

        if (z1mass<z1min || z2mass<z2min) continue;
        if (channel != ch) continue;
        if (mass<100 || mass>600) continue;
	//if (l3pdgId!=l4pdgId) continue; // using SS leptons
	bool osfail = l3pdgId!=l4pdgId && fabs(l3pdgId)==fabs(l4pdgId) && l3idNew==0 && l4idNew==0;
	if (!osfail) continue;
        //if (z1mass>120 || z2mass>120) continue;

        float sf1 = 0.0; 
        float sf2 = 0.0; 

        bool isMu = false;
        if (channel == 0 || channel == 2) isMu = true;
        
        float osss = 0.0;
        // if (channel == 0) osss = 1.28;    
        // if (channel == 1) osss = 0.93;    
        // if (channel == 2 || channel == 3) osss = 0.94;    
        if (channel == 0) osss = 1.0;    
	if (channel == 1) osss = 1.0;    
	if (channel == 2 || channel == 3) osss = 1.0;    
	
        sf1 = computeSF(l1pt, fabs(l1eta), isMu, elwp, muwp);
        sf2 = computeSF(l2pt, fabs(l2eta), isMu, elwp, muwp);

        //yield += (sf11*sf21 + 0.5*sf11*sf22 + 0.5*sf21*sf12)*osss;  // PRL triangle
	//yield += sf1*sf2*osss;  // ss
	//yield += sf1/(1-sf1)*sf2/(1-sf2); // osfail
	yield += 1.0;

        //cout << run  << " " << event << endl;
    
    }

    return yield;

}

void zx(int elwp, int muwp) {

    float yield_ee = 0.0;
    float yield_mm = 0.0;
    float yield_em = 0.0;

    yield_ee += dozx("results_data/hzzTree.root",   1, elwp, muwp);
    yield_mm += dozx("results_data/hzzTree.root",   0, elwp, muwp);
    yield_em += dozx("results_data/hzzTree.root",   2, elwp, muwp);
    yield_em += dozx("results_data/hzzTree.root",   3, elwp, muwp);

    cout << "yield(4e)   = " << yield_ee << std::endl;
    cout << "yield(4mu)  = " << yield_mm << std::endl;
    cout << "yeld(2e2mu) = " << yield_em << std::endl;
}
