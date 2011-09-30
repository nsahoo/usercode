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

using namespace std;

int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
   }

double fake(float pt,float eta,TH2F* fakeH2);
double prompt(float pt,float eta,TH2F* promptH2);

double eps(double fake){
  return fake/(1-fake);
}

double eta(double prompt){
  return (1-prompt)/prompt;
}

double wFailFail(double eps1,double eps2,
		 double eta1,double eta2){
  double result = -2*eps1*eps2;
  return result;
}

double wPassPass(double eps1,double eps2,
		 double eta1,double eta2){

  double result = -(eps1*eta1+eps2*eta2);
  return result;
}

double wPassFail(double epsPass, double etaPass,
		 double epsFail, double etaFail){
  double result = epsFail+epsPass*etaPass*epsFail;
  return result;		 
}


void setHisto2D(TH2F* fakeHisto,TH2F* promptHisto);

int addWJetsWeights(string input="latino_ll.test.root",bool addWeight=true) {

  // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  //gStyle->SetOptStat("iourme");
  gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(1111);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"Y");
  gStyle->SetTitleOffset(0.80,"X");
  gStyle->SetTitleOffset(0.80,"Y");
  // ------------------------------ 


  string openMode;
  if(addWeight)
    openMode="update";//xAlicia
  else{
    openMode="read";
  }

  TH1F* h_mll = new TH1F("h_mll","h_mll",25,0,200);
  TH1F* h_ptMax = new TH1F("h_ptMax","h_ptMax",85,15,100);
  TH1F* h_ptMin = new TH1F("h_ptMin","h_ptMin",95,5,70);


  //
  TFile* file = new TFile(input.c_str(),openMode.c_str());
  TTree* tree = (TTree*)file->Get("latino");


  const int ptBinsSize(5);
  float ptBins[ptBinsSize+1];
  ptBins[0] = 10;
  ptBins[1] = 15;
  ptBins[2] = 20;
  ptBins[3] = 25;
  ptBins[4] = 50;
  ptBins[5] = 1000;

  const int etaBinsSize(2);
  float etaBins[etaBinsSize+1];
  etaBins[0] = 0;
  etaBins[1] = 1.5;
  etaBins[2] = 2.5;

  TH2F* fakeH2 = new TH2F("fakeH2","fakeH2",ptBinsSize,ptBins,etaBinsSize,etaBins);
  TH2F* promptH2 = new TH2F("promptH2","promptH2",ptBinsSize,ptBins,etaBinsSize,etaBins);


  setHisto2D(fakeH2,promptH2);

  enum ch {mm=0,ee=1,em=2,me=3};
  
  float weight;
  float channel,pt1,pt2,eta1,eta2;
  int passTight1,passTight2;
  
  float trigger,pfmet,ppfmet,mll,mpmet,njet,nextra;
  int zveto,dphiveto,bveto_mu,bveto_ip,sameflav;

  TBranch *newBranch;
  if(addWeight)
    newBranch = tree->Branch("fakeW",&weight,"fakeW/F");//xAlicia

  tree->SetBranchAddress("channel",&channel);
  tree->SetBranchAddress("pt1",&pt1);
  tree->SetBranchAddress("pt2",&pt2);
  tree->SetBranchAddress("eta1",&eta1);
  tree->SetBranchAddress("eta2",&eta2);
  tree->SetBranchAddress("passTight1",&passTight1);
  tree->SetBranchAddress("passTight2",&passTight2);

  tree->SetBranchAddress("trigger",&trigger);
  tree->SetBranchAddress("pfmet",&pfmet);
  tree->SetBranchAddress("ppfmet",&ppfmet);
  tree->SetBranchAddress("mll",&mll);
  tree->SetBranchAddress("zveto",&zveto);
  tree->SetBranchAddress("mpmet",&mpmet);
  tree->SetBranchAddress("sameflav",&sameflav);
  tree->SetBranchAddress("njet",&njet);
  tree->SetBranchAddress("dphiveto",&dphiveto);
  tree->SetBranchAddress("bveto_mu",&bveto_mu);
  tree->SetBranchAddress("nextra",&nextra);
  tree->SetBranchAddress("bveto_ip",&bveto_ip);
  //tree->SetBranchAddress("",&);
  
  double sumOfWeights(0.);
  double sum00(0),sum10(0.),sum11(0.);
  double checksum(0.);

  unsigned nentries = tree->GetEntries();
  for(unsigned int i = 0; i<nentries; i++){
    weight = 1.;//xAlicia
    tree->GetEntry(i);
    bool ww0JetSel = trigger && pfmet>20 && mll >12 && zveto && mpmet > (20+20*sameflav);
    ww0JetSel = ww0JetSel && njet==0 && (dphiveto || !sameflav) && bveto_mu && nextra == 0 && bveto_ip;
    //if(!ww0JetSel) continue;

    double w(999999999.);
    if( channel == mm){
      checksum +=1;

      double fake1 = fake(pt1,eta1,fakeH2);
      double fake2 = fake(pt2,eta2,fakeH2);
      double prompt1 = prompt(pt1,eta1,promptH2);
      double prompt2 = prompt(pt2,eta2,promptH2);

      double Eps1 = eps(fake1);
      double Eps2 = eps(fake2);

      double Eta1 = eta(prompt1);
      double Eta2 = eta(prompt2);

      if(passTight1 && passTight2){
	w = wPassPass(Eps1,Eps2,Eta1,Eta2);
	sum11 += w/(1-Eps1*Eta1)/(1-Eps2*Eta2);
      }else if(!passTight1 && !passTight2){
	w = wFailFail(Eps1,Eps2,Eta1,Eta2);
	sum00 += w/(1-Eps1*Eta1)/(1-Eps2*Eta2);
      }else if(passTight1 && !passTight2){
	w = wPassFail(Eps1,Eta1,
		      Eps2,Eta2);
	sum10 += w/(1-Eps1*Eta1)/(1-Eps2*Eta2);
      }else if(!passTight1 && passTight2){ 
	w = wPassFail(Eps2,Eta2,
		      Eps1,Eta1);
	sum10 += w/(1-Eps1*Eta1)/(1-Eps2*Eta2);
      }else{
	cout << "ERROR: logical problem" << endl;
	return 1;
      }
      w = w/(1-Eps1*Eta1)/(1-Eps2*Eta2);
      weight = w; //converting to float for saving into tree //xAlicia

      //------ summing the weights only for those events passing the full selection
      if(ww0JetSel) {
	h_mll->Fill(mll,w);
	h_ptMax->Fill(std::max(pt1,pt2),w);
	h_ptMin->Fill(std::min(pt1,pt2),w);
	sumOfWeights += w; 
      }
    }else{
      w = 1.0;
    }
    if(addWeight) newBranch->Fill();//xAlicia
  }

  
  cout << "sum00 mm: " << sum00 << endl;
  cout << "sum10 mm: " << sum10 << endl;
  cout << "sum11 mm: " << sum11 << endl;

  cout << "sumOfWeights mm: " << sumOfWeights << endl;
  cout << "checksum mm: " << checksum << endl;

  //====================================
  /*
  cout << "mll integral: " << h_mll->Integral() << endl;
  h_mll->Draw();gPad->Update(); Wait();
  cout << "ptMax integral: " << h_ptMax->Integral() << endl;
  h_ptMax->Draw();gPad->Update(); Wait();
  cout << "ptMin integral: " << h_ptMin->Integral() << endl;
  h_ptMin->Draw();gPad->Update(); Wait();
  */

  tree->Write("",TObject::kOverwrite);//xAlicia
  delete file;//xAlicia
  return 0;
}	      



double fake(float pt,float eta,TH2F* h2){
  int bin = h2->FindFixBin(pt,fabs(eta));
  double result = h2->GetBinContent(bin);
  //cout << "fake: " << result << endl;
  return result;
}

double prompt(float pt,float eta,TH2F* h2){
  int bin = h2->FindFixBin(pt,fabs(eta));
  double result = h2->GetBinContent(bin);
  //cout << "prompt: " << result << endl;
  return result;
}



void setHisto2D(TH2F* fakeH2,TH2F* promptH2){
  // Alicia's numbers from 
  //https://espace.cern.ch/calderon/Documents/HWW%20WJets%20results%202011/lumi_1.5fb_LP.pdf
  fakeH2->SetBinContent(1,1,0.235);
  fakeH2->SetBinContent(2,1,0.118);
  fakeH2->SetBinContent(3,1,0.251);
  fakeH2->SetBinContent(4,1,0.215);
  fakeH2->SetBinContent(5,1,0.215);
  //
  fakeH2->SetBinContent(1,2,0.286);
  fakeH2->SetBinContent(2,2,0.189);
  fakeH2->SetBinContent(3,2,0.198);
  fakeH2->SetBinContent(4,2,0.199);
  fakeH2->SetBinContent(5,2,0.199);
  //
  //
  promptH2->SetBinContent(1,1,0.757);
  promptH2->SetBinContent(2,1,0.790);
  promptH2->SetBinContent(3,1,0.938);
  promptH2->SetBinContent(4,1,0.972);
  promptH2->SetBinContent(5,1,0.988);
  //
  promptH2->SetBinContent(1,2,0.784);
  promptH2->SetBinContent(2,2,0.808);
  promptH2->SetBinContent(3,2,0.927);
  promptH2->SetBinContent(4,2,0.960);
  promptH2->SetBinContent(5,2,0.985);


}
