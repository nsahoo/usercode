#define event_cxx
#include "event.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void event::Loop()
{
	//   In a ROOT session, you can do:
	//      Root > .L event.C
	//      Root > event t
	//      Root > t.GetEntry(12); // Fill t data members with entry number 12
	//      Root > t.Show();       // Show values of entry 12
	//      Root > t.Show(16);     // Read and show values of entry 16
	//      Root > t.Loop();       // Loop on all entries
	//
	
	//     This is the loop skeleton where:
	//    jentry is the global entry number in the chain
	//    ientry is the entry number in the current Tree
	//  Note that the argument to GetEntry must be:
	//    jentry for TChain::GetEntry
	//    ientry for TTree::GetEntry and TBranch::GetEntry
	//
	//       To read only selected branches, Insert statements like:
	// METHOD1:
	//    fChain->SetBranchStatus("*",0);  // disable all branches
	//    fChain->SetBranchStatus("branchname",1);  // activate branchname
	// METHOD2: replace line
	//    fChain->GetEntry(jentry);       //read all branches
	//by  b_branchname->GetEntry(ientry); //read only this branch
	if (fChain == 0) return;
	
	Long64_t nentries = fChain->GetEntriesFast();
	
	// book some histos:
	
	TH1F* scall = new TH1F("scall","DeltaR: MC - ECAL S.C. (all)",50,0.0,0.21);
	TH1F* scfid = new TH1F("scfid","DeltaR: MC - ECAL S.C. (fid)",50,0.0,0.21);

	TH1F* qh = new TH1F("qh","Reco M(Z1,Z2), loose ele cuts (CMS)",50,0.0,200.);
	TH1F* rh = new TH1F("rh","Reco M(Z1,Z2), loose ele cuts (S&T)",50,0.0,200.);

	
	int fid, qel, rel;
	
	Long64_t nbytes = 0, nb = 0;
	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = LoadTree(jentry);
		if (ientry < 0) break;
		nb = fChain->GetEntry(jentry);   nbytes += nb;
		
		// Set cuts:
	    fid = 1;
	  	for (int j=0; j<4; j++) { if ( fabs(mc_eta[j])>2.5 ) fid = 0;	}
	    qel = 1 ;
	  	for (int j=0; j<4; j++) { if ( el_dr[j] >=0.1 || el_eopout[j] <0.5 || fabs(el_detain[j]) >0.02 || el_hoe[j] >0.2 || el_dphiin[j]>0.1 ) qel = 0; }
	    rel = 1;
	  	for (int j=0; j<4; j++) { if ( el1_dr[j] >=0.1 || el1_eopout[j] <0.5 || fabs(el1_detain[j]) >0.02 || el1_hoe[j] >0.2 || el1_dphiin[j]>0.1 ) rel = 0; }
		
		
		// fill histos:
		if (qel==1) {qh->Fill(qmh);}
		if (rel==1) {rh->Fill(rmh);}
		
		for (int j=0; j<4; j++) {
			scall->Fill(sc_dr[j]);
			if (fid==1) {scfid->Fill(sc_dr[j]);}
		}
		
		// if (Cut(ientry) < 0) continue;
	}
}
