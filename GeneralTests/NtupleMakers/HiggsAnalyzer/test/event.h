//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 29 13:01:56 2007 by ROOT version 5.14/00
// from TTree event/Event data
// found on file: /Users/yagil/data/eid/checks/hzz4e_150.root
//////////////////////////////////////////////////////////

#ifndef event_h
#define event_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class event {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Int_t           run;
   Int_t           id;
   Int_t           mc_n;
   Float_t         mc_pt[4];   //[mc_n]
   Float_t         mc_eta[4];   //[mc_n]
   Float_t         mc_phi[4];   //[mc_n]
   Int_t           mc_id[4];   //[mc_n]
   Float_t         mc_e[4];   //[mc_n]
   Int_t           mc_mother[4];   //[mc_n]
   Int_t           mc_crack[4];   //[mc_n]
   Float_t         sc_e[4];   //[mc_n]
   Float_t         sc_rawe[4];   //[mc_n]
   Float_t         sc_et[4];   //[mc_n]
   Float_t         sc_eta[4];   //[mc_n]
   Float_t         sc_phi[4];   //[mc_n]
   Float_t         sc_dr[4];   //[mc_n]
   Int_t           sc_type[4];   //[mc_n]
   Float_t         el_pt[4];   //[mc_n]
   Float_t         el_e[4];   //[mc_n]
   Float_t         el_eta[4];   //[mc_n]
   Float_t         el_phi[4];   //[mc_n]
   Float_t         el_dr[4];   //[mc_n]
   Float_t         el_eopin[4];   //[mc_n]
   Float_t         el_eopout[4];   //[mc_n]
   Float_t         el_pout[4];   //[mc_n]
   Float_t         el_fbrem[4];   //[mc_n]
   Float_t         el_hoe[4];   //[mc_n]
   Float_t         el_detain[4];   //[mc_n]
   Float_t         el_dphiin[4];   //[mc_n]
   Float_t         el_detaout[4];   //[mc_n]
   Float_t         el_dphiout[4];   //[mc_n]
   Float_t         el_e3x3[4];   //[mc_n]
   Float_t         el_e5x5[4];   //[mc_n]
   Float_t         el_eseed[4];   //[mc_n]
   Float_t         el_spp[4];   //[mc_n]
   Float_t         el_see[4];   //[mc_n]
   Int_t           el_class[4];   //[mc_n]
   Float_t         el1_pt[4];   //[mc_n]
   Float_t         el1_e[4];   //[mc_n]
   Float_t         el1_eta[4];   //[mc_n]
   Float_t         el1_phi[4];   //[mc_n]
   Float_t         el1_dr[4];   //[mc_n]
   Float_t         el1_eopin[4];   //[mc_n]
   Float_t         el1_eopout[4];   //[mc_n]
   Float_t         el1_pout[4];   //[mc_n]
   Float_t         el1_fbrem[4];   //[mc_n]
   Float_t         el1_hoe[4];   //[mc_n]
   Float_t         el1_detain[4];   //[mc_n]
   Float_t         el1_dphiin[4];   //[mc_n]
   Float_t         el1_detaout[4];   //[mc_n]
   Float_t         el1_dphiout[4];   //[mc_n]
   Float_t         el1_e3x3[4];   //[mc_n]
   Float_t         el1_e5x5[4];   //[mc_n]
   Float_t         el1_eseed[4];   //[mc_n]
   Float_t         el1_spp[4];   //[mc_n]
   Float_t         el1_see[4];   //[mc_n]
   Int_t           el1_class[4];   //[mc_n]
   Int_t           nz;
   Float_t         mz1;
   Float_t         mz2;
   Float_t         mzh;
   Float_t         mzl;
   Float_t         mh;
   Float_t         rmz1;
   Float_t         rmz2;
   Float_t         rmh;
   Float_t         qmz1;
   Float_t         qmz2;
   Float_t         qmh;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_id;   //!
   TBranch        *b_mc_n;   //!
   TBranch        *b_mc_pt;   //!
   TBranch        *b_mc_eta;   //!
   TBranch        *b_mc_phi;   //!
   TBranch        *b_mc_id;   //!
   TBranch        *b_mc_e;   //!
   TBranch        *b_mc_mother;   //!
   TBranch        *b_mc_crack;   //!
   TBranch        *b_sc_e;   //!
   TBranch        *b_sc_rawe;   //!
   TBranch        *b_sc_et;   //!
   TBranch        *b_sc_eta;   //!
   TBranch        *b_sc_phi;   //!
   TBranch        *b_sc_dr;   //!
   TBranch        *b_sc_type;   //!
   TBranch        *b_el_pt;   //!
   TBranch        *b_el_e;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_dr;   //!
   TBranch        *b_el_eopin;   //!
   TBranch        *b_el_eopout;   //!
   TBranch        *b_el_pout;   //!
   TBranch        *b_el_fbrem;   //!
   TBranch        *b_el_hoe;   //!
   TBranch        *b_el_detain;   //!
   TBranch        *b_el_dphiin;   //!
   TBranch        *b_el_detaout;   //!
   TBranch        *b_el_dphiout;   //!
   TBranch        *b_el_e3x3;   //!
   TBranch        *b_el_e5x5;   //!
   TBranch        *b_el_eseed;   //!
   TBranch        *b_el_spp;   //!
   TBranch        *b_el_see;   //!
   TBranch        *b_el_class;   //!
   TBranch        *b_el1_pt;   //!
   TBranch        *b_el1_e;   //!
   TBranch        *b_el1_eta;   //!
   TBranch        *b_el1_phi;   //!
   TBranch        *b_el1_dr;   //!
   TBranch        *b_el1_eopin;   //!
   TBranch        *b_el1_eopout;   //!
   TBranch        *b_el1_pout;   //!
   TBranch        *b_el1_fbrem;   //!
   TBranch        *b_el1_hoe;   //!
   TBranch        *b_el1_detain;   //!
   TBranch        *b_el1_dphiin;   //!
   TBranch        *b_el1_detaout;   //!
   TBranch        *b_el1_dphiout;   //!
   TBranch        *b_el1_e3x3;   //!
   TBranch        *b_el1_e5x5;   //!
   TBranch        *b_el1_eseed;   //!
   TBranch        *b_el1_spp;   //!
   TBranch        *b_el1_see;   //!
   TBranch        *b_el1_class;   //!
   TBranch        *b_nz;   //!
   TBranch        *b_mz1;   //!
   TBranch        *b_mz2;   //!
   TBranch        *b_mzh;   //!
   TBranch        *b_mzl;   //!
   TBranch        *b_mh;   //!
   TBranch        *b_rmz1;   //!
   TBranch        *b_rmz2;   //!
   TBranch        *b_rmh;   //!
   TBranch        *b_qmz1;   //!
   TBranch        *b_qmz2;   //!
   TBranch        *b_qmh;   //!

   event(TTree *tree=0);
   virtual ~event();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef event_cxx
event::event(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Users/yagil/data/eid/checks/hzz4e_150.root");
      if (!f) {
         f = new TFile("/Users/yagil/data/eid/checks/hzz4e_150.root");
      }
      tree = (TTree*)gDirectory->Get("event");

   }
   Init(tree);
}

event::~event()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t event::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t event::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void event::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("id", &id, &b_id);
   fChain->SetBranchAddress("mc_n", &mc_n, &b_mc_n);
   fChain->SetBranchAddress("mc_pt", mc_pt, &b_mc_pt);
   fChain->SetBranchAddress("mc_eta", mc_eta, &b_mc_eta);
   fChain->SetBranchAddress("mc_phi", mc_phi, &b_mc_phi);
   fChain->SetBranchAddress("mc_id", mc_id, &b_mc_id);
   fChain->SetBranchAddress("mc_e", mc_e, &b_mc_e);
   fChain->SetBranchAddress("mc_mother", mc_mother, &b_mc_mother);
   fChain->SetBranchAddress("mc_crack", mc_crack, &b_mc_crack);
   fChain->SetBranchAddress("sc_e", sc_e, &b_sc_e);
   fChain->SetBranchAddress("sc_rawe", sc_rawe, &b_sc_rawe);
   fChain->SetBranchAddress("sc_et", sc_et, &b_sc_et);
   fChain->SetBranchAddress("sc_eta", sc_eta, &b_sc_eta);
   fChain->SetBranchAddress("sc_phi", sc_phi, &b_sc_phi);
   fChain->SetBranchAddress("sc_dr", sc_dr, &b_sc_dr);
   fChain->SetBranchAddress("sc_type", sc_type, &b_sc_type);
   fChain->SetBranchAddress("el_pt", el_pt, &b_el_pt);
   fChain->SetBranchAddress("el_e", el_e, &b_el_e);
   fChain->SetBranchAddress("el_eta", el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_phi", el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_dr", el_dr, &b_el_dr);
   fChain->SetBranchAddress("el_eopin", el_eopin, &b_el_eopin);
   fChain->SetBranchAddress("el_eopout", el_eopout, &b_el_eopout);
   fChain->SetBranchAddress("el_pout", el_pout, &b_el_pout);
   fChain->SetBranchAddress("el_fbrem", el_fbrem, &b_el_fbrem);
   fChain->SetBranchAddress("el_hoe", el_hoe, &b_el_hoe);
   fChain->SetBranchAddress("el_detain", el_detain, &b_el_detain);
   fChain->SetBranchAddress("el_dphiin", el_dphiin, &b_el_dphiin);
   fChain->SetBranchAddress("el_detaout", el_detaout, &b_el_detaout);
   fChain->SetBranchAddress("el_dphiout", el_dphiout, &b_el_dphiout);
   fChain->SetBranchAddress("el_e3x3", el_e3x3, &b_el_e3x3);
   fChain->SetBranchAddress("el_e5x5", el_e5x5, &b_el_e5x5);
   fChain->SetBranchAddress("el_eseed", el_eseed, &b_el_eseed);
   fChain->SetBranchAddress("el_spp", el_spp, &b_el_spp);
   fChain->SetBranchAddress("el_see", el_see, &b_el_see);
   fChain->SetBranchAddress("el_class", el_class, &b_el_class);
   fChain->SetBranchAddress("el1_pt", el1_pt, &b_el1_pt);
   fChain->SetBranchAddress("el1_e", el1_e, &b_el1_e);
   fChain->SetBranchAddress("el1_eta", el1_eta, &b_el1_eta);
   fChain->SetBranchAddress("el1_phi", el1_phi, &b_el1_phi);
   fChain->SetBranchAddress("el1_dr", el1_dr, &b_el1_dr);
   fChain->SetBranchAddress("el1_eopin", el1_eopin, &b_el1_eopin);
   fChain->SetBranchAddress("el1_eopout", el1_eopout, &b_el1_eopout);
   fChain->SetBranchAddress("el1_pout", el1_pout, &b_el1_pout);
   fChain->SetBranchAddress("el1_fbrem", el1_fbrem, &b_el1_fbrem);
   fChain->SetBranchAddress("el1_hoe", el1_hoe, &b_el1_hoe);
   fChain->SetBranchAddress("el1_detain", el1_detain, &b_el1_detain);
   fChain->SetBranchAddress("el1_dphiin", el1_dphiin, &b_el1_dphiin);
   fChain->SetBranchAddress("el1_detaout", el1_detaout, &b_el1_detaout);
   fChain->SetBranchAddress("el1_dphiout", el1_dphiout, &b_el1_dphiout);
   fChain->SetBranchAddress("el1_e3x3", el1_e3x3, &b_el1_e3x3);
   fChain->SetBranchAddress("el1_e5x5", el1_e5x5, &b_el1_e5x5);
   fChain->SetBranchAddress("el1_eseed", el1_eseed, &b_el1_eseed);
   fChain->SetBranchAddress("el1_spp", el1_spp, &b_el1_spp);
   fChain->SetBranchAddress("el1_see", el1_see, &b_el1_see);
   fChain->SetBranchAddress("el1_class", el1_class, &b_el1_class);
   fChain->SetBranchAddress("nz", &nz, &b_nz);
   fChain->SetBranchAddress("mz1", &mz1, &b_mz1);
   fChain->SetBranchAddress("mz2", &mz2, &b_mz2);
   fChain->SetBranchAddress("mzh", &mzh, &b_mzh);
   fChain->SetBranchAddress("mzl", &mzl, &b_mzl);
   fChain->SetBranchAddress("mh", &mh, &b_mh);
   fChain->SetBranchAddress("rmz1", &rmz1, &b_rmz1);
   fChain->SetBranchAddress("rmz2", &rmz2, &b_rmz2);
   fChain->SetBranchAddress("rmh", &rmh, &b_rmh);
   fChain->SetBranchAddress("qmz1", &qmz1, &b_qmz1);
   fChain->SetBranchAddress("qmz2", &qmz2, &b_qmz2);
   fChain->SetBranchAddress("qmh", &qmh, &b_qmh);
   Notify();
}

Bool_t event::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void event::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t event::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef event_cxx
