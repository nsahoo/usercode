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


void setHisto2D(TH2F* fakeHistoMu,TH2F* promptHistoMu,
		TH2F* fakeHistoMu,TH2F* promptHistoMu);

int addWJetsWeights(string input="latino_ll_LP_test.root",bool addWeight=false) {

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


  const int ptBinsSizeMu(5);
  float ptBinsMu[ptBinsSizeMu+1];
  ptBinsMu[0] = 10;
  ptBinsMu[1] = 15;
  ptBinsMu[2] = 20;
  ptBinsMu[3] = 25;
  ptBinsMu[4] = 50;
  ptBinsMu[5] = 1000;

  const int etaBinsSizeMu(2);
  float etaBinsMu[etaBinsSizeMu+1];
  etaBinsMu[0] = 0;
  etaBinsMu[1] = 1.5;
  etaBinsMu[2] = 2.5;

  const int ptBinsSizeEl(5);
  float ptBinsEl[ptBinsSizeEl+1];
  ptBinsEl[0] = 10;
  ptBinsEl[1] = 15;
  ptBinsEl[2] = 20;
  ptBinsEl[3] = 25;
  ptBinsEl[4] = 50;
  ptBinsEl[5] = 1000;

  const int etaBinsSizeElFake(2);
  float etaBinsElFake[etaBinsSizeElFake+1];
  etaBinsElFake[0] = 0;
  etaBinsElFake[1] = 1.5;
  etaBinsElFake[2] = 2.5;

  const int etaBinsSizeElPrompt(2);
  float etaBinsElPrompt[etaBinsSizeElPrompt+1];
  etaBinsElPrompt[0] = 0;
  etaBinsElPrompt[1] = 1.5;
  etaBinsElPrompt[2] = 2.5;


  TH2F* fakeMuonH2 = new TH2F("fakeMuonH2","fakeMuonH2",ptBinsSizeMu,ptBinsMu,etaBinsSizeMu,etaBinsMu);
  TH2F* promptMuonH2 = new TH2F("promptMuonH2","promptMuonH2",ptBinsSizeMu,ptBinsMu,etaBinsSizeMu,etaBinsMu);

  TH2F* fakeElH2 = new TH2F("fakeElH2","fakeElH2",ptBinsSizeEl,ptBinsEl,etaBinsSizeElFake,etaBinsElFake);
  TH2F* promptElH2 = new TH2F("promptElH2","promptElH2",ptBinsSizeEl,ptBinsEl,etaBinsSizeElPrompt,etaBinsElPrompt);


  setHisto2D(fakeMuonH2,promptMuonH2,
	     fakeElH2,promptElH2);

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
  
  int mmPP(0),mmPF(0),mmFF(0);
  int mePP(0),mePF(0),meFF(0);
  int emPP(0),emPF(0),emFF(0);
  int eePP(0),eePF(0),eeFF(0);

  double sumOfWeightsMuMu(0.);
  double sumOfWeightsMuEl(0.);
  double sumOfWeightsElMu(0.);
  double sumOfWeightsElEl(0.);
  double sumOfWeightsTrueElFakeMu(0.);
  double sumOfWeightsTrueMuFakeEl(0.);
  double sum00(0),sum10(0.),sum11(0.);
  double checksum(0.);

  unsigned nentries = tree->GetEntries();
  for(unsigned int i = 0; i<nentries; i++){
    weight = 1.;//xAlicia
    tree->GetEntry(i);
    bool ww0JetSel = trigger && pfmet>20 && mll >12 && zveto && mpmet > (20+20*sameflav);
    ww0JetSel = ww0JetSel && njet==0 && (dphiveto || !sameflav) && bveto_mu && nextra == 0 && bveto_ip;
    if(!ww0JetSel) continue;

    double wSumUp(999999999.);
    double wTrueElFakeMu(999999999.);
    double wTrueMuFakeEl(999999999.);
    //======== HERE IS THE PART COMMON TO ALL CHANNELS
    if( 1){  
      double fake1 = fake(pt1,eta1,(channel == mm || channel == me) ? fakeMuonH2 : fakeElH2);
      double fake2 = fake(pt2,eta2,(channel == ee || channel == me) ? fakeElH2 : fakeMuonH2);
      double prompt1 = prompt(pt1,eta1,(channel == mm || channel == me) ? promptMuonH2 : promptElH2);
      double prompt2 = prompt(pt2,eta2,(channel == ee || channel == me) ? promptElH2 : promptMuonH2);

      double Eps1 = eps(fake1);
      double Eps2 = eps(fake2);

      double Eta1 = eta(prompt1);
      double Eta2 = eta(prompt2);

      if(passTight1 && passTight2){
	wSumUp = -(Eps1*Eta1+Eps2*Eta2);
      }else if(!passTight1 && !passTight2){
	wSumUp = -2*Eps1*Eps2;
      }else if(passTight1 && !passTight2){
	wSumUp = Eps2+Eps1*Eta1*Eps2;
      }else if(!passTight1 && passTight2){ 
	wSumUp = Eps1+Eps2*Eta2*Eps1;
      }else{
	cout << "ERROR: logical problem" << endl;
	return 1;
      }
      double norm = 1./(1-Eps1*Eta1)/(1-Eps2*Eta2);
      wSumUp = wSumUp*norm;
      weight = wSumUp; //converting to float for saving into tree //xAlicia

      //------ summing the weights only for those events passing the full selection
      if(ww0JetSel) {
	if(channel == mm) sumOfWeightsMuMu += wSumUp; 
	if(channel == me) sumOfWeightsMuEl += wSumUp; 
	if(channel == em) sumOfWeightsElMu += wSumUp; 
	if(channel == ee) sumOfWeightsElEl += wSumUp; 

	if(passTight1 && passTight2){
	  if(channel == mm) mmPP++;
	  if(channel == me) mePP++;
	  if(channel == em) emPP++;
	  if(channel == ee) eePP++;
	}
 	if((passTight1 && !passTight2) || (!passTight1 && passTight2)){
	  if(channel == mm) mmPF++;
	  if(channel == me) mePF++;
	  if(channel == em) emPF++;
	  if(channel == ee) eePF++;
	}
 	if(!passTight1 && !passTight2){
	  if(channel == mm) mmFF++;
	  if(channel == me) meFF++;
	  if(channel == em) emFF++;
	  if(channel == ee) eeFF++;
	}
      }




    }

    //======== HERE IS THE PART FOR MIXED CHANNEL
    if( channel == em || channel == me){
      double ptEl = (channel == em) ? pt1 : pt2;
      double ptMu = (channel == em) ? pt2 : pt1;

      double psrapEl = (channel == em) ? eta1 : eta2;
      double psrapMu = (channel == em) ? eta2 : eta1;

      double fakeEl = fake(ptEl,psrapEl, fakeElH2);
      double fakeMu = fake(ptMu,psrapMu, fakeMuonH2);
      double promptEl = prompt(ptEl,psrapEl,promptElH2);
      double promptMu = prompt(ptMu,psrapMu, promptMuonH2);

      double epsEl = eps(fakeEl);
      double epsMu = eps(fakeMu);

      double etaEl = eta(promptEl);
      double etaMu = eta(promptMu);

      bool passEl,passMu;
      if(channel == em)
	{passEl = passTight1;passMu = passTight2;}
      else
	{passEl = passTight2; passMu = passTight1;}
      

      if(passEl && passMu){
	wTrueElFakeMu = -epsEl*epsMu;
	wTrueMuFakeEl = -epsMu*epsEl;
      }else if(!passEl && !passMu){
	wTrueElFakeMu = -epsMu*etaMu;
	wTrueMuFakeEl = -epsEl*etaEl;
      }else if(passEl && !passMu){
	wTrueElFakeMu = epsMu;
	wTrueMuFakeEl = epsEl*etaEl*epsMu;
      }else if(!passEl && passMu){ 
	wTrueElFakeMu = epsMu*epsEl*etaMu;
	wTrueMuFakeEl = epsEl;
      }else{
	cout << "ERROR: logical problem" << endl;
	return 1;
      }
      double norm = 1./(1-epsEl*etaEl)/(1-epsMu*etaMu);
      wTrueElFakeMu = wTrueElFakeMu*norm;
      wTrueMuFakeEl = wTrueMuFakeEl*norm;
      //weight = wSF; //TO BE FIXED

      //------ summing the weights only for those events passing the full selection
      if(ww0JetSel) {
	sumOfWeightsTrueElFakeMu += wTrueElFakeMu; 
	sumOfWeightsTrueMuFakeEl += wTrueMuFakeEl; 
      }
    }

    if(addWeight) newBranch->Fill();//xAlicia
  }

  

  cout << "sumOfWeights mm: " << sumOfWeightsMuMu << endl;
  cout << "sumOfWeights me: " << sumOfWeightsMuEl << endl;
  cout << "sumOfWeights em: " << sumOfWeightsElMu << endl;
  cout << "sumOfWeights ee: " << sumOfWeightsElEl << endl;

  cout << endl << "sumOfWeights em+me: " << sumOfWeightsElMu+sumOfWeightsMuEl << endl;

  cout << "sumOfWeights m-FakeEl: " << sumOfWeightsTrueMuFakeEl << endl;
  cout << "sumOfWeights e-FakeMu: " << sumOfWeightsTrueElFakeMu << endl;
  cout << "sumOfWeights e-FakeMu+m-FakeEl: " << (sumOfWeightsTrueElFakeMu+sumOfWeightsTrueMuFakeEl) << endl;

  float scaleToLP(7.38);
  cout << "=== Extrapolation to LP ===" << endl;
  cout << "mm to LP: " << sumOfWeightsMuMu*scaleToLP << endl;
  cout << "me to LP: " << sumOfWeightsMuEl*scaleToLP << endl;
  cout << "em to LP: " << sumOfWeightsElMu*scaleToLP << endl;
  cout << "ee to LP: " << sumOfWeightsElEl*scaleToLP << endl;

  cout << "m-FakeEl to LP: " << sumOfWeightsTrueMuFakeEl*scaleToLP << endl;
  cout << "e-FakeMu to LP: " << sumOfWeightsTrueElFakeMu*scaleToLP << endl;  

  cout << "checksum mm: " << checksum << endl;


  cout << "mm PP,PF,FF: " << mmPP << " , " << mmPF << " , " << mmFF << endl;
  cout << "me PP,PF,FF: " << mePP << " , " << mePF << " , " << meFF << endl;
  cout << "em PP,PF,FF: " << emPP << " , " << emPF << " , " << emFF << endl;
  cout << "ee PP,PF,FF: " << eePP << " , " << eePF << " , " << eeFF << endl;

  //====================================
  /*
  cout << "mll integral: " << h_mll->Integral() << endl;
  h_mll->Draw();gPad->Update(); Wait();
  cout << "ptMax integral: " << h_ptMax->Integral() << endl;
  h_ptMax->Draw();gPad->Update(); Wait();
  cout << "ptMin integral: " << h_ptMin->Integral() << endl;
  h_ptMin->Draw();gPad->Update(); Wait();
  */

  if(addWeight) tree->Write("",TObject::kOverwrite);//xAlicia
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



void setHisto2D(TH2F* fakeH2Mu,TH2F* promptH2Mu,
		TH2F* fakeH2El,TH2F* promptH2El){
  /*
  //numbers from AN2011_364_v3, Tab.8 and Tab.9
  fakeH2Mu->SetBinContent(1,1,0.235);
  fakeH2Mu->SetBinContent(2,1,0.118);
  fakeH2Mu->SetBinContent(3,1,0.251);
  fakeH2Mu->SetBinContent(4,1,0.215);
  fakeH2Mu->SetBinContent(5,1,0.215);
  //
  fakeH2Mu->SetBinContent(1,2,0.286);
  fakeH2Mu->SetBinContent(2,2,0.189);
  fakeH2Mu->SetBinContent(3,2,0.198);
  fakeH2Mu->SetBinContent(4,2,0.199);
  fakeH2Mu->SetBinContent(5,2,0.199);
  //
  //
  promptH2Mu->SetBinContent(1,1,0.757);
  promptH2Mu->SetBinContent(2,1,0.790);
  promptH2Mu->SetBinContent(3,1,0.938);
  promptH2Mu->SetBinContent(4,1,0.972);
  promptH2Mu->SetBinContent(5,1,0.988);
  //
  promptH2Mu->SetBinContent(1,2,0.784);
  promptH2Mu->SetBinContent(2,2,0.808);
  promptH2Mu->SetBinContent(3,2,0.927);
  promptH2Mu->SetBinContent(4,2,0.960);
  promptH2Mu->SetBinContent(5,2,0.985);

  //Electrons 
  fakeH2El->SetBinContent(1,1,0.071);
  fakeH2El->SetBinContent(2,1,0.060);
  fakeH2El->SetBinContent(3,1,0.097);
  fakeH2El->SetBinContent(4,1,0.076);
  fakeH2El->SetBinContent(5,1,0.076);
  //
  fakeH2El->SetBinContent(1,2,0.017);
  fakeH2El->SetBinContent(2,2,0.017);
  fakeH2El->SetBinContent(3,2,0.045);
  fakeH2El->SetBinContent(4,2,0.040);
  fakeH2El->SetBinContent(5,2,0.040);
  //
  //
  promptH2El->SetBinContent(1,1,0.69);
  promptH2El->SetBinContent(2,1,0.77);
  promptH2El->SetBinContent(3,1,0.83);
  promptH2El->SetBinContent(4,1,0.89);
  promptH2El->SetBinContent(5,1,0.91);
  //
  promptH2El->SetBinContent(1,2,0.51);
  promptH2El->SetBinContent(2,2,0.48);
  promptH2El->SetBinContent(3,2,0.62);
  promptH2El->SetBinContent(4,2,0.69);
  promptH2El->SetBinContent(5,2,0.75);
  //
  promptH2El->SetBinContent(1,3,0.40);
  promptH2El->SetBinContent(2,3,0.55);
  promptH2El->SetBinContent(3,3,0.65);
  promptH2El->SetBinContent(4,3,0.73);
  promptH2El->SetBinContent(5,3,0.81);
  */

  fakeH2Mu->SetBinContent(1,1,0.234618);
  fakeH2Mu->SetBinContent(2,1,0.117557);
  fakeH2Mu->SetBinContent(3,1,0.251346);
  fakeH2Mu->SetBinContent(4,1,0.215158);
  fakeH2Mu->SetBinContent(5,1,0.215158);
  //
  fakeH2Mu->SetBinContent(1,2,0.285527);
  fakeH2Mu->SetBinContent(2,2,0.189152);
  fakeH2Mu->SetBinContent(3,2,0.197514);
  fakeH2Mu->SetBinContent(4,2,0.199441);
  fakeH2Mu->SetBinContent(5,2,0.199441);
  //
  //
  promptH2Mu->SetBinContent(1,1,0.756528);
  promptH2Mu->SetBinContent(2,1,0.789526);
  promptH2Mu->SetBinContent(3,1,0.93848);
  promptH2Mu->SetBinContent(4,1,0.972116);
  promptH2Mu->SetBinContent(5,1,0.988179);
  //
  promptH2Mu->SetBinContent(1,2,0.783703);
  promptH2Mu->SetBinContent(2,2,0.808175);
  promptH2Mu->SetBinContent(3,2,0.926781);
  promptH2Mu->SetBinContent(4,2,0.95974);
  promptH2Mu->SetBinContent(5,2,0.985484);

  //Electrons 
  fakeH2El->SetBinContent(1,1,0.0706892);
  fakeH2El->SetBinContent(2,1,0.0600548);
  fakeH2El->SetBinContent(3,1,0.0967004);
  fakeH2El->SetBinContent(4,1,0.0759246);
  fakeH2El->SetBinContent(5,1,0.0757698);
  //
  fakeH2El->SetBinContent(1,2,0.0169693);
  fakeH2El->SetBinContent(2,2,0.0166916);
  fakeH2El->SetBinContent(3,2,0.0455474);
  fakeH2El->SetBinContent(4,2,0.0400856);
  fakeH2El->SetBinContent(5,2,0.0268137);
  //
  //
  promptH2El->SetBinContent(1,1,0.69);
  promptH2El->SetBinContent(2,1,0.77);
  promptH2El->SetBinContent(3,1,0.83);
  promptH2El->SetBinContent(4,1,0.89);
  promptH2El->SetBinContent(5,1,0.91);
  //
  promptH2El->SetBinContent(1,2,0.51);
  promptH2El->SetBinContent(2,2,0.48);
  promptH2El->SetBinContent(3,2,0.62);
  promptH2El->SetBinContent(4,2,0.69);
  promptH2El->SetBinContent(5,2,0.75);
  //
  promptH2El->SetBinContent(1,3,0.40);
  promptH2El->SetBinContent(2,3,0.55);
  promptH2El->SetBinContent(3,3,0.65);
  promptH2El->SetBinContent(4,3,0.73);
  promptH2El->SetBinContent(5,3,0.81);
}
