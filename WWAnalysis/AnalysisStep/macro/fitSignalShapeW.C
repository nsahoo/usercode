#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "THStack.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "./scales2.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ4LRooPdfs.h"


#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "Math/MinimizerOptions.h"

using namespace RooFit;



#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
}

float getFitEdge(float mass, float width, bool low);

void fitSignalShapeW(int massBin=120, int id=1120, int ch=0, int year = 2011,
		     float lumi=10, bool doSfLepton=false,double rangeLow=105., double rangeHigh=130.,
		     double bwSigma=1.,
		     double fitValues[5]=0, double fitErrors[5]=0);

void all(int channels=0,int year=2012, bool doSfLepton=true){
  /*
    channels = 0 --> 4mu
    channels = 1 --> 4el
    channels = 2 --> 2mu2el

   */
  //init();

  double bwSigma[40];
  int mass[40]; int id[40]; double xLow[40]; double xHigh[40];  
  int maxMassBin;

  XSecProvider xsecs;
  xsecs.initHiggs4lWidth();

  if(year==2011){
    init(true);
    mass[0] = 115; id[0]=200; xLow[0] = 100; xHigh[0] = 125; bwSigma[0] = 3.1/1000.;
    mass[1] = 120; id[1]=201; xLow[1] = 105; xHigh[1] = 130; bwSigma[1] = 3.5/1000.;
    mass[2] = 130; id[2]=202; xLow[2] = 115; xHigh[2] = 140; bwSigma[2] = 4.9/1000.;
    mass[3] = 140; id[3]=203; xLow[3] = 125; xHigh[3] = 150; bwSigma[3] = 8.1/1000.;
    mass[4] = 150; id[4]=204; xLow[4] = 135; xHigh[4] = 160; bwSigma[4] = 1.7/100.;
    mass[5] = 160; id[5]=205; xLow[5] = 145; xHigh[5] = 170; bwSigma[5] = 8.3/100.;
    mass[6] = 170; id[6]=206; xLow[6] = 155; xHigh[6] = 180; bwSigma[6] = 3.8/10.;
    mass[7] = 180; id[7]=207; xLow[7] = 165; xHigh[7] = 190; bwSigma[7] = 6.3/10.;
    mass[8] = 190; id[8]=208; xLow[8] = 170; xHigh[8] = 200; bwSigma[8] = 1.04;
    mass[9] = 200; id[9]=209; xLow[9] = 175; xHigh[9] = 210; bwSigma[9] = 1.43;
    mass[10] = 210; id[10]=210; xLow[10] = 185; xHigh[10] = 220; bwSigma[10] = 1.85;
    mass[11] = 220; id[11]=211; xLow[11] = 190; xHigh[11] = 230; bwSigma[11] = 2.31;
    maxMassBin = 12;
  }


  if(year==2012){
    init(false);
    float masses[32] = {115,117,119,120,121,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,220,250,300,350,400,450,500,550,650,700};
    for(int i=0;i<32;++i) {
      mass[i] = masses[i]; 
      if(masses[i]<1000) id[i]=1000+masses[i]; 
      else id[i]=11000;
      float width = xsecs.getHZZ4lWidth(masses[i]);
      xLow[i] = getFitEdge(masses[i],width,true); 
      xHigh[i] = getFitEdge(masses[i],width,false); 
      bwSigma[i] = width;
    }
    maxMassBin = 32;
  }
  // -----------------------



  double massV[40],massE[40];
  for(int i=0; i<maxMassBin;++i){
    massV[i]=mass[i];
    massE[i]=0;
  }

  double aVal[40],aErr[40];
  double nVal[40],nErr[40];
  double meanCBVal[40],meanCBErr[40];
  double sigmaCBVal[40],sigmaCBErr[40];
  double meanBWVal[40],meanBWErr[40];

  double fitValues[5];
  double fitErrors[5];

  double extendL(1),extendH(1);

  if(channels==0) {extendL=1.0;extendH=1.0;}
  if(channels==1) {extendL=0.90;extendH=1.05;}
  if(channels==2) {extendL=0.95;extendH=1.04;}



  //fitSignalShapeW(mass[0],channels,10.,doSfLepton,max(95.,xLow[0]*0.90),xHigh[0]*1.1,bwSigma[12],
  // 		  fitValues,fitErrors);  Wait();
  
  for(int i=0; i<maxMassBin;++i){

    fitSignalShapeW(mass[i],id[i],channels,year,10.,doSfLepton,max(95.,xLow[i]*extendL),xHigh[i]*extendH,bwSigma[i],
		    fitValues,fitErrors);  
  
    cout << "a value,error: " << fitValues[0] << " , " << fitErrors[0] << endl; 
    aVal[i]=fitValues[0]; aErr[i]=fitErrors[0];

    cout << "n value,error: " << fitValues[3] << " , " << fitErrors[3] << endl; 
    nVal[i]=fitValues[3]; nErr[i]=fitErrors[3];

    cout << "meanCB value,error: " << fitValues[1] << " , " << fitErrors[1] << endl;
    meanCBVal[i]=fitValues[1]; meanCBErr[i]=fitErrors[1];
    
    cout << "sigmaCB value,error: " << fitValues[4] << " , " << fitErrors[4] << endl; 
    sigmaCBVal[i]=fitValues[4]; sigmaCBErr[i]=fitErrors[4];

    cout << "meanBW value,error: " << fitValues[2] << " , " << fitErrors[2] << endl; 
    meanBWVal[i]=fitValues[2]; meanBWErr[i]=fitErrors[2];

    //Wait();
  }
  

  TGraphErrors* gA = new TGraphErrors(maxMassBin,massV,aVal,massE,aErr);
  TGraphErrors* gN = new TGraphErrors(maxMassBin,massV,nVal,massE,nErr);
  TGraphErrors* gMeanCB = new TGraphErrors(maxMassBin,massV,meanCBVal,massE,meanCBErr);
  TGraphErrors* gSigmaCB = new TGraphErrors(maxMassBin,massV,sigmaCBVal,massE,sigmaCBErr);
  TGraphErrors* gMeanBW = new TGraphErrors(maxMassBin,massV,meanBWVal,massE,meanBWErr);

  gA->SetMarkerStyle(20);   gA->SetMarkerSize(1);
  gN->SetMarkerStyle(20);   gN->SetMarkerSize(1);
  gMeanCB->SetMarkerStyle(20);   gMeanCB->SetMarkerSize(1);
  gSigmaCB->SetMarkerStyle(20);   gSigmaCB->SetMarkerSize(1);
  gMeanBW->SetMarkerStyle(20);   gMeanBW->SetMarkerSize(1);
  

  gA->SetTitle("");
  gA->GetXaxis()->SetTitle("mass (GeV)");
  gA->GetYaxis()->SetTitle("CB a-parameter");

  gN->SetTitle("");
  gN->GetXaxis()->SetTitle("mass (GeV)");
  gN->GetYaxis()->SetTitle("CB n-parameter");

  gMeanCB->SetTitle("");
  gMeanCB->GetXaxis()->SetTitle("mass (GeV)");
  gMeanCB->GetYaxis()->SetTitle("CB mean");

  gSigmaCB->SetTitle("");
  gSigmaCB->GetXaxis()->SetTitle("mass (GeV)");
  gSigmaCB->GetYaxis()->SetTitle("CB sigma");

  gMeanBW->SetTitle("");
  gMeanBW->GetXaxis()->SetTitle("mass (GeV)");
  gMeanBW->GetYaxis()->SetTitle("BW mean");



  gA->Fit("pol0"); gA->Draw("Ap"); gPad->Update(); gPad->Print("aFit.root"); Wait();
  gN->Fit("pol1"); gN->Draw("Ap"); gPad->Update(); gPad->Print("nFit.root"); Wait();
  gMeanCB->Fit("pol1"); gMeanCB->Draw("Ap"); gPad->Update(); gPad->Print("meanFit.root"); Wait();
  gSigmaCB->Fit("pol1"); gSigmaCB->Draw("Ap"); gPad->Update(); gPad->Print("sigmaFit.root"); Wait();
  //gMeanBW->Fit("pol1"); gMeanBW->Draw("Ap"); gPad->Update(); Wait();


}



void fitSignalShapeW(int massBin,int id, int channels, int year, 
		     float lumi, bool doSfLepton,double rangeLow, double rangeHigh,
		     double bwSigma,
		     double fitValues[5], double fitErrors[5]){
 // ------ root settings ---------
  gROOT->Reset();  
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  //gStyle->SetOptStat("kKsSiourRmMen");
  gStyle->SetOptStat("iourme");
  //gStyle->SetOptStat("rme");
  //gStyle->SetOptStat("");
  gStyle->SetOptFit(11);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------ 

  ROOT::Math::MinimizerOptions::SetDefaultTolerance( 1.E-7);

  TH1F* hPuWeights;


  stringstream ggFileName,vbfFileName;
  ggFileName << "root://pcmssd12//data/hzz4l/step2/HZZ4L_53X_S1_V10_S2_V06/MC/hzzTree_id" << id << ".root";
  //ggFileName << "/data/hzz4l/step2/HZZ4L_52X_S1_V05_S2_V00/hzzTree_id" << id << ".root"; 
  //ggFileName << "/data/hzz4l/step2/ichep2012/" << year << "/MC/hzzTree_id" << id << ".root"; 
  //ggFileName << "/data/hzz4l/step2/ichep2012/2011/MC/hzzTree_id207.root";


  TFile* ggFile = TFile::Open(ggFileName.str().c_str()); 

  TTree* ggTree = (TTree*) ggFile->Get("zz4lTree/probe_tree");


  TH1F* ggMass = new TH1F("ggMass","ggMass",200,100,300.);


  float nvtx,channel;
  //float ntrueint;
  float nTrueInt,nObsInt,rho,mass,m4l,trueInt;

  float pt1,eta1,id1,pt2,eta2,id2,pt3,eta3,id3,pt4,eta4,id4;
  float lsW;
  float sfLepton(1.);


int  nentries = ggTree->GetEntries();
 
 
  //--- ggTree part
  ggTree->SetBranchAddress("channel",&channel);
  ggTree->SetBranchAddress("numTrueInteractions",&nTrueInt);
  ggTree->SetBranchAddress("rho",&rho);
  ggTree->SetBranchAddress("mass",&mass);
  ggTree->SetBranchAddress("m4l",&m4l);
  ggTree->SetBranchAddress("l1pt",&pt1);ggTree->SetBranchAddress("l1eta",&eta1);ggTree->SetBranchAddress("l1pdgId",&id1);
  ggTree->SetBranchAddress("l2pt",&pt2);ggTree->SetBranchAddress("l2eta",&eta2);ggTree->SetBranchAddress("l2pdgId",&id2);
  ggTree->SetBranchAddress("l3pt",&pt3);ggTree->SetBranchAddress("l3eta",&eta3);ggTree->SetBranchAddress("l3pdgId",&id3);
  ggTree->SetBranchAddress("l4pt",&pt4);ggTree->SetBranchAddress("l4eta",&eta4);ggTree->SetBranchAddress("l4pdgId",&id4);
  ggTree->SetBranchAddress("genhiggsmassweight",&lsW);

  int massesWithLSW[13] = {1400,1450,1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,11000};
  bool hasLsW=false;
  for(int i=0;i<13;++i) {
    if(id==massesWithLSW[i]) {
      hasLsW=true;
      break;
    }
  }

  if(hasLsW) cout << "For ID = " << id << " using the Lineshape reweighting..." << endl;

   /*
  for(int k=0; k<nentries; k++){
    ggTree->GetEvent(k);
    if(doSfLepton) sfLepton = getSF(pt1, eta1, id1,typeErrSfLep)
		     *getSF(pt2, eta2, id2,typeErrSfLep)
		     *getSF(pt3, eta3, id3,typeErrSfLep)
		     *getSF(pt4, eta4, id4,typeErrSfLep);
    else sfLepton = 1.0;
    ggMass->Fill(mass,weightTrue2012(nTrueInt)*xsecweights[1000+massBin]*lumi*sfLepton);
  }
  */





  /*TCanvas* myC = new TCanvas("myCanvas","myCanvas",0,0,500,500); 
   myC->cd(1);

   vbfMass->Draw(); gPad->Update(); Wait();
  */

  //cout << "yield gg: " << ggMass->Integral() << endl;
  //cout << "yield vbf: " << vbfMass->Integral() << endl;
  //cout << "yield total for mass " << massBin << " is: " << ggMass->Integral()<< endl;
  //ggMass->Draw(); gPad->Update(); Wait();




  //--- rooFit part
  double xMin,xMax,xInit;
  xInit = (double) massBin;
  xMin = rangeLow;
  xMax = rangeHigh ;

  cout << "Summary: " << xMin << " , " << xMax << "  /  " << xInit << endl;
  //xMin = 114.;
  //xMax = 122 ;



  /*
  TCut cut = "channel == 0";
  RooRealVar x("mass","mass",xInit,xMin,xMax);
  RooRealVar ch("channel","channel",0,0,10);
  RooRealVar w("baseW","baseW",1.0,0.,1000.);
  RooArgSet ntupleVarSet(x,ch,w);
  RooDataSet dataset("mass4l","mass4l",ntupleVarSet,WeightVar(w),Import(*ggTree),Cut("channel==0"));
  */

  //---------  
  TCut cut = "channel == 0";
  RooRealVar x("mass","mass (GeV)",xInit,xMin,xMax);
  //RooRealVar x("m4l","m4l",xInit,xMin,xMax);
  RooRealVar w("myW","myW",1.0,0.,100000.);
  RooArgSet ntupleVarSet(x,w);
  RooDataSet dataset("mass4l","mass4l",ntupleVarSet,WeightVar("myW"));

  LeptSfProvider SfProvider;
  if(year==2011) SfProvider.init(true);
  else SfProvider.init(false);

  for(int k=0; k<nentries; k++){
    ggTree->GetEvent(k);
    if(doSfLepton) sfLepton = SfProvider.getSF(pt1, eta1, id1)
      *SfProvider.getSF(pt2, eta2, id2)
      *SfProvider.getSF(pt3, eta3, id3)
      *SfProvider.getSF(pt4, eta4, id4);
    else sfLepton = 1.0;

    //if(channel==ch1 || channel==ch2){      
    if(channels==0){if(channel != 0) continue;}
    if(channels==1){if(channel != 1) continue;}
    if(channels==2){if(channel != 2 && channel != 3) continue;}

    ntupleVarSet.setRealValue("mass",mass);
    ntupleVarSet.setRealValue("m4l",m4l);
    double localW(1);
    localW = getPUWeight(nTrueInt);
    localW = localW*sfLepton;
    if(hasLsW) localW = localW*lsW;

    ntupleVarSet.setRealValue("myW",localW);
    if(x.getVal()>xMin && x.getVal()<xMax)
      dataset.add(ntupleVarSet, localW);
    
  }
  //---------

  //cout << "dataset n entries: " << dataset.sumEntries() << endl;
  //cout << "reduced data n entries: " << reducedData->sumEntries() << endl;


  //--- simple CrystalBall
  float meanR = (xInit<400) ? 5.0 : 100;
  RooRealVar mean("mean","mean of gaussian",0,-5.,meanR) ;
  RooRealVar sigma("sigma","width of gaussian",1.5,0.,30.); 
  RooRealVar a("a","a",1.46,0.,4.);
  RooRealVar n("n","n",1.92,0.,25.);   
  RooCBShape CBall("CBall","Crystal ball",x, mean,sigma, a,n);
  

  //--- simple Gaussian
  RooRealVar mean2("mean2","mean2 of gaussian",xInit,xInit*0.80,xInit*1.2) ;
  RooRealVar sigma2("sigma2","width2 of gaussian",10.,5.,500.); 
  RooGaussian tailCatcher("tailCatcher","tailCatcher",x,mean2,sigma2);
  RooRealVar fsig("fsig","signal fraction",0.95,0.7,1.);


  //--- Breit-Wigner
  RooRealVar mean3("mean3","mean3",xInit) ;
  RooRealVar sigma3("sigma3","width3",bwSigma); 
  RooRealVar scale3("scale3","scale3 ",1.); 
  RooRelBWUFParam bw("bw","bw",x,mean3,scale3);

  //RooAddPdf model("model","model",RooArgList(CBall,tailCatcher),fsig);
  //RooCBShape model("model","model",x, mean,sigma, a,n);
  x.setBins(100000,"fft");
  RooFFTConvPdf model("model","model",x,bw,CBall);

  

  model.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8));

  stringstream frameTitle;
  if(channels==0){frameTitle << "4mu, mH = ";}
  if(channels==1){frameTitle << "4el, mH = ";}
  if(channels==2){frameTitle << "2mu2el, mH = ";}
  frameTitle << massBin << " GeV";

  RooPlot* xframe = x.frame(Title(frameTitle.str().c_str() )) ;
  dataset.plotOn(xframe,DataError(RooAbsData::SumW2) );
  model.plotOn(xframe);
  //model.paramOn(xframe);


  stringstream nameFile;
  nameFile << "fitM" << massBin << "_cha" << channels << ".pdf";
  xframe->Draw(); gPad->Update(); gPad->Print(nameFile.str().c_str());


  if(fitValues!=0){
    fitValues[0] = a.getVal();
    fitValues[1] = mean.getVal();
    fitValues[2] = mean3.getVal();
    fitValues[3] = n.getVal();
    fitValues[4] = sigma.getVal();
  }

  if(fitErrors!=0){
    fitErrors[0] = a.getError();
    fitErrors[1] = mean.getError();
    fitErrors[2] = mean3.getError();
    fitErrors[3] = n.getError();
    fitErrors[4] = sigma.getError();
  }

}

float getFitEdge(float mass, float width, bool low) {
  double windowVal = max(width, float(1.));
  double lowside = (mass >= 275) ? 180. : 100.;
  if (low) return std::max((mass - 20.*windowVal), lowside);
  else return std::min((mass + 15.*windowVal), 900.);
}
