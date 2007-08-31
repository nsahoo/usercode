{

  gROOT->ProcessLine(".L ~yagil/ana/rootTools/binerr.C");
  gROOT->ProcessLine(".L ~yagil/ana/rootTools/mk.C");   
  gROOT->ProcessLine("gStyle->SetOptStat(010011)");

//TFile *feid = TFile::Open("~yagil/data/eid/checks/hzz4e_150.root");
TFile *feid = TFile::Open("merged.root");
TTree *eid  = (TTree*)(gROOT->FindObject("event")); 

TCut acc(" fabs(mc_eta[0])<2.5 && fabs(mc_eta[1])<2.5 && fabs(mc_eta[2])<2.5 && fabs(mc_eta[3])<2.5 "); 

// simple loose ele Q-cuts
// CMS collection
TCut qel0  ("el_dr[0] <0.1&&el_eopout[0] >0.5 &&fabs(el_detain[0]) <0.02&&el_hoe[0] <0.2 && el_dphiin[0]<0.1");
TCut qel1  ("el_dr[1] <0.1&&el_eopout[1] >0.5 &&fabs(el_detain[1]) <0.02&&el_hoe[1] <0.2 && el_dphiin[1]<0.1");
TCut qel2  ("el_dr[2] <0.1&&el_eopout[2] >0.5 &&fabs(el_detain[2]) <0.02&&el_hoe[2] <0.2 && el_dphiin[2]<0.1");
TCut qel3  ("el_dr[3] <0.1&&el_eopout[3] >0.5 &&fabs(el_detain[3]) <0.02&&el_hoe[3] <0.2 && el_dphiin[3]<0.1");
TCut qel = qel0+qel1+qel2+qel3;
//SNT collection
TCut rel0  ("el1_dr[0] <0.1&&el1_eopout[0] >0.5 &&fabs(el1_detain[0]) <0.02&&el1_hoe[0] <0.2 && el1_dphiin[0]<0.1");
TCut rel1  ("el1_dr[1] <0.1&&el1_eopout[1] >0.5 &&fabs(el1_detain[1]) <0.02&&el1_hoe[1] <0.2 && el1_dphiin[1]<0.1");
TCut rel2  ("el1_dr[2] <0.1&&el1_eopout[2] >0.5 &&fabs(el1_detain[2]) <0.02&&el1_hoe[2] <0.2 && el1_dphiin[2]<0.1");
TCut rel3  ("el1_dr[3] <0.1&&el1_eopout[3] >0.5 &&fabs(el1_detain[3]) <0.02&&el1_hoe[3] <0.2 && el1_dphiin[3]<0.1");
TCut rel = rel0+rel1+rel2+rel3;



gStyle->SetOptStat("eou");
TCanvas *b = new TCanvas("b", "MC info",  50,50,800,800);
  b->Divide(1,2);
  TH1F* h1 = new TH1F("h1","h1",26,-6,6);
  h1->SetTitle("MC electron: Eta ");
  TH1F* h2 = new TH1F("h2","h2",20,0.0,100.);
  h2->SetTitle("MC electron: PT ");
  b->cd(1); eid->Draw("mc_eta>>h1");
  b->cd(2); eid->Draw("mc_pt>>h2");
  b->cd(1); h1->Draw();
  b->cd(2); h2->Draw();
//  b->Print("~yagil/data/eid/checks/hzz/mc.gif");    

TCanvas *z = new TCanvas("z", "Z/Z* MC info",  50,50,800,800);
  z->Divide(1,2);
  TH1F* h6 = new TH1F("h6","h6",50,0,200);
  h6->SetTitle("MC M(Z1)");
  TH1F* h7 = new TH1F("h7","h7",50,0.0,200.);
  h7->SetTitle("MC M(Z2) ");
  z->cd(1); eid->Draw("mz1>>h6");
  z->cd(2); eid->Draw("mz2>>h7");
  z->cd(1); h6->Draw();
  z->cd(2); h7->Draw();
    
TCanvas *b1 = new TCanvas("b1", "SC matching to MC electrons",  50,50,800,800);
  b1->Divide(1,2);
  TH1F* h3 = new TH1F("h3","h3",50,0.0,0.21);
  h3->SetTitle("DeltaR: MC - ECAL S.C. (all)");
  TH1F* h4 = new TH1F("h4","h4",50,0.0,0.21);
  h4->SetTitle("DeltaR: MC - ECAL S.C. (fid)");
  b1->cd(1); eid->Draw("sc_dr>>h3");    
  b1->cd(2); eid->Draw("sc_dr>>h4",acc);
  b1->cd(1); h3->Draw();
  b1->cd(2); h4->Draw();
// b1->Print("~yagil/data/eid/checks/pics/sc.gif");    


TCanvas *x = new TCanvas("x", "H MC info",  50,50,800,800);
  x->Divide(1,2);
  TH1F* h8 = new TH1F("h8","h8",50,0,200);
  h8->SetTitle("MC M(Z1,Z2), all");
  TH1F* h9 = new TH1F("h9","h9",50,0.0,200.);
  h9->SetTitle("MC M(Z1,Z2), |eta|<2.5 ");
  x->cd(1); eid->Draw("mh>>h8");
  x->cd(2); eid->Draw("mh>>h9",acc);
  x->cd(1); h8->Draw();
  x->cd(2); h9->Draw();
  
  
  TCanvas *y = new TCanvas("y"," Reco - Mz1 Vs. Mz2 ",  50,50,800,800);
  y->Divide(1,2);
  TH1F* qz1 = new TH1F("qz1","qz1",50,0,200);
  qz1->SetTitle("Reco M(Z1,Z2), loose ele cuts (CMS)");
  TH1F* rz1 = new TH1F("rz1","rz1",50,0,200);
  rz1->SetTitle("Reco M(Z1,Z2), loose ele cuts (S&T) ");
  y->cd(1); eid->Draw("qmh>>qz1",qel);
  y->cd(2); eid->Draw("rmh>>rz1",rel);
  y->cd(1); qz1->Draw();
  y->cd(2); rz1->Draw();

TCanvas *x = new TCanvas("x", "Higgs Reco ",  50,50,800,800);
  x->Divide(1,2);
  TH1F* q1 = new TH1F("q1","q1",50,0,200);
  q1->SetTitle("Reco M(Z1,Z2), loose ele cuts (CMS)");
  TH1F* r1 = new TH1F("r1","r1",50,0,200);
  r1->SetTitle("Reco M(Z1,Z2), loose ele cuts (S&T) ");
  x->cd(1); eid->Draw("qmh>>q1",qel);
  x->cd(2); eid->Draw("rmh>>r1",rel);
  x->cd(1); q1->Draw();
  x->cd(2); r1->Draw();

}
