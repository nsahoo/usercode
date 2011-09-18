#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCut.h"
#include "THStack.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <sstream>
#include <math.h>


int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
   }


void makeFakeRatePlot(TTree* treeZ,TTree* treeW,TTree* treeQCD,TTree* treeData,
		      double scale1,double scale2,double scale3,
		      float* ptBin,int ptBinsSize,
		      TCut selection, string variable,
		      string labelX,
		      TH1D* frPt1,TH1D* frPt2,TH1D* frPt3,TH1D* frPt4,TH1D* frPt5);
	
void makeStackPlot(TTree* tree1,TTree* tree2,TTree* tree3,TTree* tree4,
		   TTree* treeD1,TTree* treeD2,TTree* treeD3,
		   double scale1,double scale2,double scale3,
		   //float* bins,int binsSize,
		   TCut selection, string variable,
		   string labelX,
		   TCanvas* canv);

void makeStackPlotJet(TTree* tree1,TTree* tree2,TTree* tree3,TTree* tree4,
		      TTree* treeD1,TTree* treeD2,TTree* treeD3,
		      double scale1,double scale2,double scale3,
		      //float* bins,int binsSize,
		      TCut selection, string variable,
		      string labelX,
		      TCanvas* canv);


void setHistoAttributes(TH1D* h,
			int lCol,int lWidth,float max, float min,string xTitle, string yTitle,string Title);	      

void muFR() {

  // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  //gStyle->SetOptStat("iourme");
  //gStyle->SetOptStat("rme");
  gStyle->SetOptStat("");
  gStyle->SetOptFit(1111);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleOffset(0.80,"X");
  gStyle->SetTitleOffset(0.80,"Y");
  // ------------------------------ 

  /* NO MC triggers
  TFile* f1 = new TFile("output/tree_004.root");
  TFile* f2 = new TFile("output/tree_026.root");
  TFile* f3 = new TFile("output/tree_101.root");
  TFile* f4 = new TFile("output/tree_151.root");
  */

  // yes MC triggers
  TFile* f1 = new TFile("last/tree_004.root");
  TFile* f2 = new TFile("last/tree_026.root");
  TFile* f3 = new TFile("last/tree_101.root");
  TFile* f4 = new TFile("last/tree_151.root");
  TFile* fD1 = new TFile("last/tree_151_mu8.root");
  TFile* fD2 = new TFile("last/tree_151_mu15.root");
  TFile* fD3 = new TFile("last/tree_151_mu24.root");
  

  double scale1 = 0.8396525673;
  double scale2 = 2.0722688028;


  //double sigma3 = 296600000.*0.00118;
  //double scale3 = 1000.*sigma3/7892468;
  //double scale3 = 1000.*sigma3/7642468.;
  double scale3 = 43.5174874180;
    
  //scale to 1 /pb
  scale1 *= 1./1000.;
  scale2 *= 1./1000.;
  scale3 *= 1./1000.;
  
  

  TCut commonSelection = "dR>1.0 && jetPt >20 && mt < 20 && pfmet < 20";
  TCut zVeto = "(!(mass>0 && mass<12) && !(mass>0 && abs(mass-91.1876) < 15  ))";
  TCut selection = commonSelection && zVeto;

  TTree* tree1 = (TTree*)f1->Get("treeL/probe_tree");
  TTree* tree2 = (TTree*)f2->Get("treeL/probe_tree");
  TTree* tree3 = (TTree*)f3->Get("treeL/probe_tree");
  TTree* tree4 = (TTree*)f4->Get("treeL/probe_tree");

  TTree* treeD1 = (TTree*)fD1->Get("treeL/probe_tree");
  TTree* treeD2 = (TTree*)fD2->Get("treeL/probe_tree");
  TTree* treeD3 = (TTree*)fD3->Get("treeL/probe_tree");


  //====================================
  const int ptBinsSize(5);
  float ptBins[ptBinsSize+1];
  ptBins[0] = 10;
  ptBins[1] = 14;
  ptBins[2] = 20;
  ptBins[3] = 26;
  ptBins[4] = 50;
  ptBins[5] = 100;

  const int etaBinsSize(3);
  float etaBins[etaBinsSize+1];
  etaBins[0] = -2.5;
  etaBins[1] = -1.5;
  etaBins[2] =  1.5;
  etaBins[3] =  2.5;
  
  
  const int phiBinsSize(10);
  float phiBins[phiBinsSize+1];
  float rangePhi = 2*M_PI;
  float binSizePhi = rangePhi/10.;
  for(unsigned int i=0;i<=10;i++){
    phiBins[i] = -M_PI + i*binSizePhi;
  }
  
  
  //====================================
  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  c1->Divide(2,2); c1->cd(1);
  

  string ptLabelX("pT (GeV)");
  string etaLabelX("eta");
  string phiLabelX("phi");

  //------
  
  TH1D* frPt1 = new TH1D("frPt1","frPt1",ptBinsSize,ptBins);
  TH1D* frPt2 = new TH1D("frPt2","frPt2",ptBinsSize,ptBins);
  TH1D* frPt3 = new TH1D("frPt3","frPt3",ptBinsSize,ptBins);
  TH1D* frPt4 = new TH1D("frPt4","frPt4",ptBinsSize,ptBins);
  TH1D* frPt5 = new TH1D("frPt5","frPt5",ptBinsSize,ptBins);

  string varPt("muPt");
  makeFakeRatePlot(tree1,tree2,tree3,tree4,
		   scale1,scale2,scale3,
		   ptBins,ptBinsSize,
		   selection, varPt,
		   ptLabelX,
		   frPt1,frPt2,frPt3,frPt4,frPt5);

  c1->cd(1);  frPt1->Draw();
  c1->cd(2);  frPt2->Draw();
  c1->cd(3);  frPt3->Draw();
  c1->cd(4);  frPt5->Draw("hist");frPt4->Draw("same P");

  
  TLegend* leg = new TLegend(0.5,0.7,0.9,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(frPt5,"Z+W+QCD MCs","F");
  leg->AddEntry(frPt4,"May10ReReco data","L");
  leg->Draw();

  gPad->Update();
  c1->Print("fakeRatePt.pdf"); Wait();
  


  //------
  
  TH1D* frEta1 = new TH1D("frEta1","frEta1",etaBinsSize,etaBins);
  TH1D* frEta2 = new TH1D("frEta2","frEta2",etaBinsSize,etaBins);
  TH1D* frEta3 = new TH1D("frEta3","frEta3",etaBinsSize,etaBins);
  TH1D* frEta4 = new TH1D("frEta4","frEta4",etaBinsSize,etaBins);
  TH1D* frEta5 = new TH1D("frEta5","frEta5",etaBinsSize,etaBins);

  string varEta("muEta");
  makeFakeRatePlot(tree1,tree2,tree3,tree4,
		   scale1,scale2,scale3,
		   etaBins,etaBinsSize,
		   selection, varEta,
		   etaLabelX,
		   frEta1,frEta2,frEta3,frEta4,frEta5);
  
  c1->cd(1);  frEta1->Draw();
  c1->cd(2);  frEta2->Draw();
  c1->cd(3);  frEta3->Draw();
  c1->cd(4);  frEta5->Draw("hist");frEta4->Draw("same P");


  TLegend* leg2 = new TLegend(0.5,0.7,0.9,0.89);
  leg2->SetFillColor(0);
  leg2->AddEntry(frEta5,"Z+W+QCD MCs","F");
  leg2->AddEntry(frEta4,"May10ReReco data","L");
  leg2->Draw();

  gPad->Update();
  c1->Print("fakeRateEta.pdf"); Wait();
  

  //------
  /*
  TH1D* frPhi1 = new TH1D("frPhi1","frPhi1",phiBinsSize,phiBins);
  TH1D* frPhi2 = new TH1D("frPhi2","frPhi2",phiBinsSize,phiBins);
  TH1D* frPhi3 = new TH1D("frPhi3","frPhi3",phiBinsSize,phiBins);
  TH1D* frPhi4 = new TH1D("frPhi4","frPhi4",phiBinsSize,phiBins);
  TH1D* frPhi5 = new TH1D("frPhi5","frPhi5",phiBinsSize,phiBins);

  string varPhi("muPhi");
  makeFakeRatePlot(tree1,tree2,tree3,tree4,
		   scale1,scale2,scale3,
		   phiBins,phiBinsSize,
		   selection, varPhi,
		   phiLabelX,
		   frPhi1,frPhi2,frPhi3,frPhi4,frPhi5);

  c1->cd(1);  frPhi1->Draw();
  c1->cd(2);  frPhi2->Draw();
  c1->cd(3);  frPhi3->Draw();
  c1->cd(4);  frPhi5->Draw("hist");frPhi4->Draw("same P");


  TLegend* leg3 = new TLegend(0.5,0.7,0.9,0.89);
  leg3->SetFillColor(0);
  leg3->AddEntry(frPhi5,"Z+W+QCD MCs","F");
  leg3->AddEntry(frPhi4,"May10ReReco data","L");
  leg3->Draw();

  gPad->Update();
  c1->Print("fakeRatePhi.pdf"); Wait();
  */

  // ---------
  c1->Divide(1,1); c1->cd();

  //string varPt("muPt");


  makeStackPlot(tree1,tree2,tree3,tree4,
		treeD1,treeD2,treeD3,
		scale1,scale2,scale3,
		selection,varPt,
		ptLabelX,
		c1);


  string varJetPt("jetPt");
  makeStackPlotJet(tree1,tree2,tree3,tree4,
		treeD1,treeD2,treeD3,
		scale1,scale2,scale3,
		selection,varJetPt,
		ptLabelX,
		c1);

}


   
void makeFakeRatePlot(TTree* tree1,TTree* tree2,TTree* tree3,TTree* tree4,
		      double scale1,double scale2,double scale3,
		      float* bins,int binsSize,
		      TCut selection, string variable,
		      string labelX,
		      TH1D* hRatio1,TH1D* hRatio2,TH1D* hRatio3,TH1D* hRatio4,TH1D* hRatio5){

  TH1D* hD1 = new TH1D("hD1","histo den1",binsSize,bins);
  TH1D* hD2 = new TH1D("hD2","histo den1",binsSize,bins);
  TH1D* hD3 = new TH1D("hD3","histo den1",binsSize,bins);
  TH1D* hD4 = new TH1D("hD4","histo den1",binsSize,bins);

  TH1D* hN1 = new TH1D("hN1","histo num1",binsSize,bins);
  TH1D* hN2 = new TH1D("hN2","histo num2",binsSize,bins);
  TH1D* hN3 = new TH1D("hN3","histo num3",binsSize,bins);
  TH1D* hN4 = new TH1D("hN4","histo num4",binsSize,bins);
    
  stringstream streamD1; streamD1 << variable << ">>hD1";
  stringstream streamD2; streamD2 << variable << ">>hD2";
  stringstream streamD3; streamD3 << variable << ">>hD3";
  stringstream streamD4; streamD4 << variable << ">>hD4";

  stringstream streamN1; streamN1 << variable << ">>hN1";
  stringstream streamN2; streamN2 << variable << ">>hN2";
  stringstream streamN3; streamN3 << variable << ">>hN3";
  stringstream streamN4; streamN4 << variable << ">>hN4";

  tree1->Draw(streamD1.str().c_str(),selection && "isLoose");
  tree1->Draw(streamN1.str().c_str(),selection && "isLoose && isTight");

  tree2->Draw(streamD2.str().c_str(),selection && "isLoose");
  tree2->Draw(streamN2.str().c_str(),selection && "isLoose && isTight");

  tree3->Draw(streamD3.str().c_str(),selection && "isLoose");
  tree3->Draw(streamN3.str().c_str(),selection && "isLoose && isTight");

  tree4->Draw(streamD4.str().c_str(),selection && "isLoose");
  tree4->Draw(streamN4.str().c_str(),selection && "isLoose && isTight");

  hN1->Sumw2(); hD1->Sumw2();
  hN2->Sumw2(); hD2->Sumw2();
  hN3->Sumw2(); hD3->Sumw2();
  hN4->Sumw2(); hD4->Sumw2();

  TH1D* hN5 = (TH1D*) hN1->Clone(); hN5->SetName("hN5"); hN5->Reset();
  TH1D* hD5 = (TH1D*) hD1->Clone(); hD5->SetName("hD5"); hD5->Reset();

  hN1->Scale(scale1); 
  hN2->Scale(scale2); 
  hN3->Scale(scale3); 

  hN5->Add(hN1);
  hN5->Add(hN2);
  hN5->Add(hN3);


  hD1->Scale(scale1);   
  hD2->Scale(scale2); 
  hD3->Scale(scale3); 

  hD5->Add(hD1);
  hD5->Add(hD2);
  hD5->Add(hD3);

  hRatio1->Divide(hN1,hD1,1.,1.,"B");
  hRatio2->Divide(hN2,hD2,1.,1.,"B");
  hRatio3->Divide(hN3,hD3,1.,1.,"B");
  hRatio4->Divide(hN4,hD4,1.,1.,"B");
  hRatio5->Divide(hN5,hD5,1.,1.,"B");

  string yTitle("#LooseAndTight / #Loose");
  string xTitle(labelX);

  setHistoAttributes(hRatio1,2,2,1.01,0.0,xTitle,yTitle,"ZJets MC");
  setHistoAttributes(hRatio2,2,2,1.01,0.0,xTitle,yTitle,"WJets MC");
  setHistoAttributes(hRatio3,2,2,1.01,0.0,xTitle,yTitle,"QCD MC");
  setHistoAttributes(hRatio4,1,2,1.01,0.0,xTitle,yTitle,"");
  setHistoAttributes(hRatio5,2,2,1.01,0.0,xTitle,yTitle,"");

  hRatio5->SetFillColor(2);

}


void setHistoAttributes(TH1D* h,
			int lCol,int lWidth,float max, float min,string xTitle, string yTitle,string title){
  
  h->SetLineColor(lCol);
  h->SetLineWidth(lWidth);
  h->SetMaximum(max);
  h->SetMinimum(min);

  h->GetXaxis()->SetTitle(xTitle.c_str());
  h->GetYaxis()->SetTitle(yTitle.c_str());
  h->SetTitle(title.c_str());



}	      





void makeStackPlot(TTree* tree1,TTree* tree2,TTree* tree3,TTree* tree4,
		   TTree* treeD1,TTree* treeD2,TTree* treeD3,
		   double scale1,double scale2,double scale3,
		   //float* bins,int binsSize,
		   TCut selection, string variable,
		   string labelX,
		   TCanvas* canv){

  canv->SetLogy(1);
  canv->cd();

  THStack *hs2 = new THStack("hs2","Stacked 1D histograms");
  TH1F* hpt1 = new TH1F("hpt1","histo pt1",100,10.,100);
  TH1F* hpt2 = new TH1F("hpt2","histo pt2",100,10.,100);
  TH1F* hpt3 = new TH1F("hpt3","histo pt3",100,10.,100);
  TH1F* hpt4 = new TH1F("hpt4","histo pt4",100,10.,100);

  TH1F* hpt_mu8 = new TH1F("hpt_mu8","histo pt mu8",6,10.,16);
  TH1F* hpt_mu15 = new TH1F("hpt_mu15","histo pt mu15",11,16.,27);
  TH1F* hpt_mu24 = new TH1F("hpt_mu24","histo pt mu24",73,27.,100);



  stringstream stream1; stream1 << variable << ">>hpt1";
  stringstream stream2; stream2 << variable << ">>hpt2";
  stringstream stream3; stream3 << variable << ">>hpt3";
  stringstream stream4; stream4 << variable << ">>hpt4";

  stringstream streamMu8; streamMu8 << variable << ">>hpt_mu8";
  stringstream streamMu15; streamMu15 << variable << ">>hpt_mu15";
  stringstream streamMu24; streamMu24 << variable << ">>hpt_mu24";


  tree1->Draw(stream1.str().c_str(),selection && "isLoose"); 
  tree2->Draw(stream2.str().c_str(),selection && "isLoose");
  tree3->Draw(stream3.str().c_str(),selection && "isLoose");
  tree4->Draw(stream4.str().c_str(),selection && "isLoose");

  //hpt4->Draw(); gPad->Update(); Wait();

  treeD1->Draw(streamMu8.str().c_str(),selection && "isLoose");
  treeD2->Draw(streamMu15.str().c_str(),selection && "isLoose");
  treeD3->Draw(streamMu24.str().c_str(),selection && "isLoose");

  hpt1->SetFillColor(kBlue);
  hpt2->SetFillColor(kRed);
  hpt3->SetFillColor(kGreen);

  hpt1->Scale(scale1);
  hpt2->Scale(scale2);
  hpt3->Scale(scale3);


  hs2->Add(hpt1);
  hs2->Add(hpt2);
  hs2->Add(hpt3);



  hpt_mu8->Sumw2();
  hpt_mu15->Sumw2();
  hpt_mu24->Sumw2();

  hpt_mu8->Scale(1./0.61337);
  hpt_mu15->Scale(1./48.32);
  hpt_mu24->Scale(1./210.42);


  hpt_mu8->SetLineWidth(2);
  hpt_mu15->SetLineWidth(2);
  hpt_mu24->SetLineWidth(2);
  
  float maximum = hpt_mu8->GetMaximum();

  hs2->SetMaximum(maximum);
  hs2->Draw("");
  
  hpt_mu8->Draw("same p");
  hpt_mu15->Draw("same p");
  hpt_mu24->Draw("same p");


  hs2->GetXaxis()->SetTitle(labelX.c_str());
  
  TLegend* leg = new TLegend(0.5,0.7,0.9,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(hpt1,"ZJets","F");
  leg->AddEntry(hpt2,"WJets","F");
  leg->AddEntry(hpt3,"QCD","F");
  leg->Draw();
  //hpt4->Scale(0.5);
  //hpt4->Draw("same");
  gPad->Update();
  gPad->Print("stackPt.pdf");Wait();

  
  
}

void makeStackPlotJet(TTree* tree1,TTree* tree2,TTree* tree3,TTree* tree4,
		      TTree* treeD1,TTree* treeD2,TTree* treeD3,
		      double scale1,double scale2,double scale3,
		      //float* bins,int binsSize,
		      TCut selection, string variable,
		      string labelX,
		      TCanvas* canv){

  canv->SetLogy(1);
  canv->cd();

  THStack *hs2 = new THStack("hs2","Stacked 1D histograms");
  TH1F* hpt1 = new TH1F("hpt1","histo pt1",10,10.,160);
  TH1F* hpt2 = new TH1F("hpt2","histo pt2",10,10.,160);
  TH1F* hpt3 = new TH1F("hpt3","histo pt3",10,10.,160);
  TH1F* hpt4 = new TH1F("hpt4","histo pt4",10,10.,160);

  TH1F* hpt_mu8 = new TH1F("hpt_mu8","histo pt mu8",10,10.,160);
  TH1F* hpt_mu15 = new TH1F("hpt_mu15","histo pt mu15",10,10.,160);
  TH1F* hpt_mu24 = new TH1F("hpt_mu24","histo pt mu24",10,10.,160);



  stringstream stream1; stream1 << variable << ">>hpt1";
  stringstream stream2; stream2 << variable << ">>hpt2";
  stringstream stream3; stream3 << variable << ">>hpt3";
  stringstream stream4; stream4 << variable << ">>hpt4";

  stringstream streamMu8; streamMu8 << variable << ">>hpt_mu8";
  stringstream streamMu15; streamMu15 << variable << ">>hpt_mu15";
  stringstream streamMu24; streamMu24 << variable << ">>hpt_mu24";


  tree1->Draw(stream1.str().c_str(),selection && "isLoose"); 
  tree2->Draw(stream2.str().c_str(),selection && "isLoose");
  tree3->Draw(stream3.str().c_str(),selection && "isLoose");
  tree4->Draw(stream4.str().c_str(),selection && "isLoose");

  //hpt4->Draw(); gPad->Update(); Wait();

  treeD1->Draw(streamMu8.str().c_str(),selection && "isLoose");
  treeD2->Draw(streamMu15.str().c_str(),selection && "isLoose");
  treeD3->Draw(streamMu24.str().c_str(),selection && "isLoose");

  hpt1->SetFillColor(kBlue);
  hpt2->SetFillColor(kRed);
  hpt3->SetFillColor(kGreen);

  hpt1->Scale(scale1);
  hpt2->Scale(scale2);
  hpt3->Scale(scale3);


  hs2->Add(hpt1);
  hs2->Add(hpt2);
  hs2->Add(hpt3);



  hpt_mu8->Sumw2();
  hpt_mu15->Sumw2();
  hpt_mu24->Sumw2();

  hpt_mu8->Scale(1./0.61337);
  hpt_mu15->Scale(1./48.32);
  hpt_mu24->Scale(1./210.42);


  hpt_mu8->SetLineWidth(2);
  hpt_mu15->SetLineWidth(2);
  hpt_mu24->SetLineWidth(2);
  
  hpt_mu8->Add(hpt_mu15);
  hpt_mu8->Add(hpt_mu24);
  
  //float maximum = hpt_mu8->GetMaximum();

  //hs2->SetMaximum(maximum);
  hs2->Draw("");
  
  hpt_mu8->Draw("same p");
  //hpt_mu15->Draw("same p");
  //hpt_mu24->Draw("same p");


  hs2->GetXaxis()->SetTitle(labelX.c_str());
  
  TLegend* leg = new TLegend(0.5,0.7,0.9,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(hpt1,"ZJets","F");
  leg->AddEntry(hpt2,"WJets","F");
  leg->AddEntry(hpt3,"QCD","F");
  leg->Draw();
  //hpt4->Scale(0.5);
  //hpt4->Draw("same");
  gPad->Update();
  gPad->Print("stackPt.pdf");Wait();

  
  
}






  //====================================
  /*
  TCanvas* c0 = new TCanvas("c0","c0",400,400);
  TH1F* hmet1 = new TH1F("hmet1","histo met1",100,0,25);
  TH1F* hmet2 = new TH1F("hmet2","histo met2",100,0,25);
  TH1F* hmet3 = new TH1F("hmet3","histo met3",100,0,25);
  TH1F* hmet4 = new TH1F("hmet4","histo met4",100,0,25);

  THStack *hs = new THStack("hs","Stacked 1D histograms");


  tree1->Draw("pfmet>>hmet1",commonSelection && "isLoose"); 
  tree2->Draw("pfmet>>hmet2",commonSelection && "isLoose");
  tree3->Draw("pfmet>>hmet3",commonSelection && "isLoose");
  tree4->Draw("pfmet>>hmet4",commonSelection && "isLoose");

  hmet1->Scale(scale1);
  hmet2->Scale(scale2);
  hmet3->Scale(scale3);

  c0->cd();
  hmet1->Draw();gPad->Update();Wait();
  hmet2->Draw();gPad->Update();Wait();
  hmet3->Draw();gPad->Update();Wait();

  hmet1->SetFillColor(kBlue);
  hmet2->SetFillColor(kRed);
  hmet3->SetFillColor(kGreen);

  hs->Add(hmet1);
  hs->Add(hmet2);
  hs->Add(hmet3);
  c0->cd(); hs->Draw(); 
  //gPad->WaitPrimitive();
  
  TLegend* leg = new TLegend(0.5,0.7,0.9,0.89);
  leg->SetFillColor(0);
  leg->AddEntry(hmet1,"ZJets","F");
  leg->AddEntry(hmet2,"WJets","F");
  leg->AddEntry(hmet3,"QCD","F");
  leg->Draw();
  hmet4->Scale(0.5);
  //hmet4->Draw("same");
  gPad->Update();Wait();
  */
