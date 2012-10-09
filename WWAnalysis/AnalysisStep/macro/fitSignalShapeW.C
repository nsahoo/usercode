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
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "WWAnalysis/TreeModifiers/src/RooHighMassBW.h"

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
float getFitEdgeHighMass(float mass, float width, bool low);

void fitSignalShapeW(int massBin=120, int id=1120, int ch=0, int year = 2011,
		     float lumi=10, bool doSfLepton=false,double rangeLow=105., double rangeHigh=130.,
		     double bwSigma=1.,
		     double fitValues[7]=0, double fitErrors[7]=0);

void fitSignalHighMassShapeW(int massBin=120, int id=1120, int ch=0, int year = 2011,
                             float lumi=10, bool doSfLepton=false,double rangeLow=105., double rangeHigh=130.,
                             double bwSigma=1.,
                             double fitValues[8]=0, double fitErrors[8]=0);

void allHighMass(int channels=0,int year=2012, bool doSfLepton=false);

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
    init(false);     //  0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25 
    float masses[26] = {115,117,119,120,121,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,220,250,300,350};
    for(int i=0;i<26;++i) {
      mass[i] = masses[i]; 
      id[i]=1000+masses[i]; 
      float width = xsecs.getHZZ4lWidth(masses[i]);
      xLow[i] = getFitEdge(masses[i],width,true); 
      xHigh[i] = getFitEdge(masses[i],width,false); 
      //cout << "For mass = " << masses[i] << " width = " << width << "; => Fit Range = [" << xLow[i] << "," << xHigh[i] << "]" << endl;
      bwSigma[i] = width;
    }
    maxMassBin = 26;
  }
  // -----------------------



  double massV[40],massE[40];
  for(int i=0; i<maxMassBin;++i){
    massV[i]=mass[i];
    massE[i]=0;
  }

  double a1Val[40],a1Err[40];
  double a2Val[40],a2Err[40];
  double n1Val[40],n1Err[40];
  double n2Val[40],n2Err[40];
  double meanCBVal[40],meanCBErr[40];
  double sigmaCBVal[40],sigmaCBErr[40];
  double meanBWVal[40],meanBWErr[40];

  double fitValues[8];
  double fitErrors[8];


  //fitSignalShapeW(mass[0],channels,10.,doSfLepton,max(95.,xLow[0]*0.90),xHigh[0]*1.1,bwSigma[12],
  // 		  fitValues,fitErrors);  Wait();
  
  for(int i=0; i<maxMassBin;++i){
    // skip the gg masses not yet ready
    if(mass[i]==250 || mass[i]==350) continue;
    fitSignalShapeW(mass[i],id[i],channels,year,10.,doSfLepton,xLow[i],xHigh[i],bwSigma[i],
		    fitValues,fitErrors);  
  
    cout << "a1 value,error: " << fitValues[0] << " , " << fitErrors[0] << endl; 
    a1Val[i]=fitValues[0]; a1Err[i]=fitErrors[0];

    cout << "a2 value,error: " << fitValues[1] << " , " << fitErrors[1] << endl; 
    a2Val[i]=fitValues[1]; a2Err[i]=fitErrors[1];

    cout << "n1 value,error: " << fitValues[4] << " , " << fitErrors[4] << endl; 
    n1Val[i]=fitValues[4]; n1Err[i]=fitErrors[4];

    cout << "n2 value,error: " << fitValues[5] << " , " << fitErrors[5] << endl; 
    n2Val[i]=fitValues[5]; n2Err[i]=fitErrors[5];

    cout << "meanCB value,error: " << fitValues[2] << " , " << fitErrors[2] << endl;
    meanCBVal[i]=fitValues[2]; meanCBErr[i]=fitErrors[2];
    
    cout << "sigmaCB value,error: " << fitValues[6] << " , " << fitErrors[6] << endl; 
    sigmaCBVal[i]=fitValues[6]; sigmaCBErr[i]=fitErrors[6];

    cout << "meanBW value,error: " << fitValues[3] << " , " << fitErrors[3] << endl; 
    meanBWVal[i]=fitValues[3]; meanBWErr[i]=fitErrors[3];

    //Wait();
  }
  

  stringstream namefile;
  namefile << "parameters_channel" << channels << ".root";
  TFile *resultfile = TFile::Open(namefile.str().c_str(),"RECREATE");

  TGraphErrors* gA1 = new TGraphErrors(maxMassBin,massV,a1Val,massE,a1Err);
  TGraphErrors* gA2 = new TGraphErrors(maxMassBin,massV,a2Val,massE,a2Err);
  TGraphErrors* gN1 = new TGraphErrors(maxMassBin,massV,n1Val,massE,n1Err);
  TGraphErrors* gN2 = new TGraphErrors(maxMassBin,massV,n2Val,massE,n2Err);
  TGraphErrors* gMeanCB = new TGraphErrors(maxMassBin,massV,meanCBVal,massE,meanCBErr);
  TGraphErrors* gSigmaCB = new TGraphErrors(maxMassBin,massV,sigmaCBVal,massE,sigmaCBErr);
  TGraphErrors* gMeanBW = new TGraphErrors(maxMassBin,massV,meanBWVal,massE,meanBWErr);

  gA1->SetName("gA1");  gA1->SetMarkerStyle(20);   gA1->SetMarkerSize(1);
  gA2->SetName("gA2");  gA2->SetMarkerStyle(20);   gA2->SetMarkerSize(1);
  gN1->SetName("gN1");  gN1->SetMarkerStyle(20);   gN1->SetMarkerSize(1);
  gN2->SetName("gN2");  gN2->SetMarkerStyle(20);   gN2->SetMarkerSize(1);
  gMeanCB->SetName("gMeanCB"); gMeanCB->SetMarkerStyle(20);   gMeanCB->SetMarkerSize(1);
  gSigmaCB->SetName("gSigmaCB"); gSigmaCB->SetMarkerStyle(20);   gSigmaCB->SetMarkerSize(1);
  gMeanBW->SetName("gMeanBW"); gMeanBW->SetMarkerStyle(20);   gMeanBW->SetMarkerSize(1);
  
  gA1->SetTitle("");
  gA1->GetXaxis()->SetTitle("mass (GeV)");
  gA1->GetYaxis()->SetTitle("CB a-parameter");

  gA2->SetTitle("");
  gA2->GetXaxis()->SetTitle("mass (GeV)");
  gA2->GetYaxis()->SetTitle("CB a-parameter");

  gN1->SetTitle("");
  gN1->GetXaxis()->SetTitle("mass (GeV)");
  gN1->GetYaxis()->SetTitle("CB n-parameter");

  gN2->SetTitle("");
  gN2->GetXaxis()->SetTitle("mass (GeV)");
  gN2->GetYaxis()->SetTitle("CB n-parameter");

  gMeanCB->SetTitle("");
  gMeanCB->GetXaxis()->SetTitle("mass (GeV)");
  gMeanCB->GetYaxis()->SetTitle("CB mean");

  gSigmaCB->SetTitle("");
  gSigmaCB->GetXaxis()->SetTitle("mass (GeV)");
  gSigmaCB->GetYaxis()->SetTitle("CB sigma");

  gMeanBW->SetTitle("");
  gMeanBW->GetXaxis()->SetTitle("mass (GeV)");
  gMeanBW->GetYaxis()->SetTitle("BW mean");
  
  resultfile->cd();
  gA1->Fit("pol0"); gA1->Draw("Ap"); gA1->Write();
  gA2->Fit("pol0"); gA2->Draw("Ap"); gA2->Write();
  gN1->Fit("pol1"); gN1->Draw("Ap"); gN1->Write();
  gN2->Fit("pol1"); gN2->Draw("Ap"); gN2->Write();
  gMeanCB->Fit("pol1"); gMeanCB->Draw("Ap"); gMeanCB->Write();
  gSigmaCB->Fit("pol1"); gSigmaCB->Draw("Ap"); gSigmaCB->Write();
  
  resultfile->Close();

}

void allHighMass(int channels,int year, bool doSfLepton){
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
    maxMassBin = 0;
  }


  if(year==2012){
    init(false);   
                         //  0   1   2   3   4   5   6   7   8   9   10   11  
    float highmasses[12] = {400,450,500,600,650,700,750,800,850,900,950,1000};
    for(int i=0;i<12;++i) {
      mass[i] = highmasses[i]; 
      if(mass[i]<1000) id[i]=1000+mass[i]; 
      else id[i]=11000;
      float width = xsecs.getHZZ4lWidth(mass[i]);
      xLow[i] = getFitEdgeHighMass(mass[i],width,true); 
      xHigh[i] = getFitEdgeHighMass(mass[i],width,false); 
      //cout << "For mass = " << masses[i] << " width = " << width << "; => Fit Range = [" << xLow[i] << "," << xHigh[i] << "]" << endl;
      bwSigma[i] = width;
    }
    maxMassBin = 12;
  }
  // -----------------------



  double massV[40],massE[40];
  for(int i=0; i<maxMassBin;++i){
    massV[i]=mass[i];
    massE[i]=0;
  }

  double a1Val[40],a1Err[40];
  double a2Val[40],a2Err[40];
  double n1Val[40],n1Err[40];
  double n2Val[40],n2Err[40];
  double meanCBVal[40],meanCBErr[40];
  double sigmaCBVal[40],sigmaCBErr[40];
  double meanBWVal[40],meanBWErr[40];
  double sigmaBWVal[40],sigmaBWErr[40];

  double fitValues[8];
  double fitErrors[8];


  //fitSignalShapeW(mass[0],channels,10.,doSfLepton,max(95.,xLow[0]*0.90),xHigh[0]*1.1,bwSigma[12],
  // 		  fitValues,fitErrors);  Wait();
  
  for(int i=0; i<maxMassBin;++i){

    fitSignalHighMassShapeW(mass[i],id[i],channels,year,10.,doSfLepton,xLow[i],xHigh[i],bwSigma[i],
                            fitValues,fitErrors);  
  
    cout << "a1 value,error: " << fitValues[0] << " , " << fitErrors[0] << endl; 
    a1Val[i]=fitValues[0]; a1Err[i]=fitErrors[0];

    cout << "a2 value,error: " << fitValues[1] << " , " << fitErrors[1] << endl; 
    a2Val[i]=fitValues[1]; a2Err[i]=fitErrors[1];

    cout << "n1 value,error: " << fitValues[4] << " , " << fitErrors[4] << endl; 
    n1Val[i]=fitValues[4]; n1Err[i]=fitErrors[4];

    cout << "n2 value,error: " << fitValues[5] << " , " << fitErrors[5] << endl; 
    n2Val[i]=fitValues[5]; n2Err[i]=fitErrors[5];

    cout << "meanCB value,error: " << fitValues[2] << " , " << fitErrors[2] << endl;
    meanCBVal[i]=fitValues[2]; meanCBErr[i]=fitErrors[2];
    
    cout << "sigmaCB value,error: " << fitValues[6] << " , " << fitErrors[6] << endl; 
    sigmaCBVal[i]=fitValues[6]; sigmaCBErr[i]=fitErrors[6];

    cout << "meanBW value,error: " << fitValues[3] << " , " << fitErrors[3] << endl; 
    meanBWVal[i]=fitValues[3]; meanBWErr[i]=fitErrors[3];

    cout << "sigmaBW value,error: " << fitValues[7] << " , " << fitErrors[7] << endl; 
    sigmaBWVal[i]=fitValues[7]; sigmaBWErr[i]=fitErrors[7];

    //Wait();
  }
  

  stringstream namefile;
  namefile << "parameters_channel" << channels << ".root";
  TFile *resultfile = TFile::Open(namefile.str().c_str(),"RECREATE");

  TGraphErrors* gA1 = new TGraphErrors(maxMassBin,massV,a1Val,massE,a1Err);
  TGraphErrors* gA2 = new TGraphErrors(maxMassBin,massV,a2Val,massE,a2Err);
  TGraphErrors* gN1 = new TGraphErrors(maxMassBin,massV,n1Val,massE,n1Err);
  TGraphErrors* gN2 = new TGraphErrors(maxMassBin,massV,n2Val,massE,n2Err);
  TGraphErrors* gMeanCB = new TGraphErrors(maxMassBin,massV,meanCBVal,massE,meanCBErr);
  TGraphErrors* gSigmaCB = new TGraphErrors(maxMassBin,massV,sigmaCBVal,massE,sigmaCBErr);
  TGraphErrors* gMeanBW = new TGraphErrors(maxMassBin,massV,meanBWVal,massE,meanBWErr);
  TGraphErrors* gSigmaBW = new TGraphErrors(maxMassBin,massV,sigmaBWVal,massE,sigmaBWErr);

  gA1->SetName("gA1");  gA1->SetMarkerStyle(20);   gA1->SetMarkerSize(1);
  gA2->SetName("gA2");  gA2->SetMarkerStyle(20);   gA2->SetMarkerSize(1);
  gN1->SetName("gN1");  gN1->SetMarkerStyle(20);   gN1->SetMarkerSize(1);
  gN2->SetName("gN2");  gN2->SetMarkerStyle(20);   gN2->SetMarkerSize(1);
  gMeanCB->SetName("gMeanCB"); gMeanCB->SetMarkerStyle(20);   gMeanCB->SetMarkerSize(1);
  gSigmaCB->SetName("gSigmaCB"); gSigmaCB->SetMarkerStyle(20);   gSigmaCB->SetMarkerSize(1);
  gMeanBW->SetName("gMeanBW"); gMeanBW->SetMarkerStyle(20);   gMeanBW->SetMarkerSize(1);
  gSigmaBW->SetName("gSigmaBW"); gSigmaBW->SetMarkerStyle(20);   gSigmaBW->SetMarkerSize(1);
  
  gA1->SetTitle("");
  gA1->GetXaxis()->SetTitle("mass (GeV)");
  gA1->GetYaxis()->SetTitle("CB a-parameter");

  gA2->SetTitle("");
  gA2->GetXaxis()->SetTitle("mass (GeV)");
  gA2->GetYaxis()->SetTitle("CB a-parameter");

  gN1->SetTitle("");
  gN1->GetXaxis()->SetTitle("mass (GeV)");
  gN1->GetYaxis()->SetTitle("CB n-parameter");

  gN2->SetTitle("");
  gN2->GetXaxis()->SetTitle("mass (GeV)");
  gN2->GetYaxis()->SetTitle("CB n-parameter");

  gMeanCB->SetTitle("");
  gMeanCB->GetXaxis()->SetTitle("mass (GeV)");
  gMeanCB->GetYaxis()->SetTitle("CB mean");

  gSigmaCB->SetTitle("");
  gSigmaCB->GetXaxis()->SetTitle("mass (GeV)");
  gSigmaCB->GetYaxis()->SetTitle("CB sigma");

  gMeanBW->SetTitle("");
  gMeanBW->GetXaxis()->SetTitle("mass (GeV)");
  gMeanBW->GetYaxis()->SetTitle("BW mean");

  gSigmaBW->SetTitle("");
  gSigmaBW->GetXaxis()->SetTitle("mass (GeV)");
  gSigmaBW->GetYaxis()->SetTitle("BW #Gamma");
  
  resultfile->cd();
  gA1->Fit("pol0"); gA1->Draw("Ap"); gA1->Write();
  gA2->Fit("pol0"); gA2->Draw("Ap"); gA2->Write();
  gN1->Fit("pol1"); gN1->Draw("Ap"); gN1->Write();
  gN2->Fit("pol1"); gN2->Draw("Ap"); gN2->Write();
  gMeanCB->Fit("pol1"); gMeanCB->Draw("Ap"); gMeanCB->Write();
  gSigmaCB->Fit("pol1"); gSigmaCB->Draw("Ap"); gSigmaCB->Write(); 
  gSigmaBW->Fit("pol1"); gSigmaBW->Draw("Ap"); gSigmaBW->Write();
 
  resultfile->Close();

}

void interpolateAgain(int channel) {

  stringstream filename;
  filename << "parameters_channel" << channel << ".root";
  TFile *resultfile = TFile::Open(filename.str().c_str());
  TGraphErrors *gA1 = (TGraphErrors*)resultfile->Get("gA1");
  TGraphErrors *gA2 = (TGraphErrors*)resultfile->Get("gA2");
  TGraphErrors *gN1 = (TGraphErrors*)resultfile->Get("gN1");
  TGraphErrors *gN2 = (TGraphErrors*)resultfile->Get("gN2"); 
  TGraphErrors *gMeanCB = (TGraphErrors*)resultfile->Get("gMeanCB");
  TGraphErrors *gSigmaCB = (TGraphErrors*)resultfile->Get("gSigmaCB");

  // 
  gA1->GetYaxis()->SetRangeUser(-1,6);
  gA2->GetYaxis()->SetRangeUser(-1,6);
  gN1->GetYaxis()->SetRangeUser(0,50); 
  gN2->GetYaxis()->SetRangeUser(0,50); 
  gMeanCB->GetYaxis()->SetRangeUser(-2,2); 
  gSigmaCB->GetYaxis()->SetRangeUser(0,5); 

  // skip some mass points where we did not have MC
  gA1->SetPointError(23,0,1000);
  gA2->SetPointError(23,0,1000); 
  gN1->SetPointError(23,0,1000);
  gN2->SetPointError(23,0,1000);
  gMeanCB->SetPointError(23,0,1000);
  gSigmaCB->SetPointError(23,0,1000);
  
  gA1->SetPointError(25,0,1000);
  gA2->SetPointError(25,0,1000); 
  gN1->SetPointError(25,0,1000);
  gN2->SetPointError(25,0,1000);
  gMeanCB->SetPointError(25,0,1000);
  gSigmaCB->SetPointError(25,0,1000);
  
  TCanvas *c1 = new TCanvas("c1","c1");

  stringstream channame;
  channame << "_channel" << channel << ".pdf";
  gA1->Fit("pol2","","",115,300); gA1->Draw("Ap"); gPad->Update(); gPad->Print((string("gA1")+channame.str()).c_str());
  gA2->Fit("pol2","","",115,300); gA2->Draw("Ap"); gPad->Update(); gPad->Print((string("gA2")+channame.str()).c_str());
  gN1->Fit("pol2","","",115,300); gN1->Draw("Ap"); gPad->Update(); gPad->Print((string("gN1")+channame.str()).c_str());
  gN2->Fit("pol2","","",115,300); gN2->Draw("Ap"); gPad->Update(); gPad->Print((string("gN2")+channame.str()).c_str());
  gMeanCB->Fit("pol2","","",115,300); gMeanCB->Draw("Ap"); gPad->Update(); gPad->Print((string("gMeanCB")+channame.str()).c_str());
  gSigmaCB->Fit("pol1","","",115,300); gSigmaCB->Draw("Ap"); gPad->Update(); gPad->Print((string("gSigmaCB")+channame.str()).c_str());

}

void fitSignalShapeW(int massBin,int id, int channels, int year, 
		     float lumi, bool doSfLepton,double rangeLow, double rangeHigh,
		     double bwSigma,
		     double fitValues[7], double fitErrors[7]){
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
  ggFileName << "/cmsrm/pc21_2/emanuele/data/hzz4l/HZZ4L_53X_S1_V10_S2_V14/MC/hzzTree_id" << id << ".root"; 
  
  cout << "Using " << ggFileName.str() << endl;

  TFile* ggFile = new TFile(ggFileName.str().c_str()); 

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

  //--- rooFit part
  double xMin,xMax,xInit;
  xInit = (double) massBin;
  xMin = rangeLow;
  xMax = rangeHigh ;
  cout << "Fit range: [" << xMin << " , " << xMax << "]. Init value = " << xInit << endl;

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
  RooRealVar w("myW","myW",1.0,0.,1000.);
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


  //--- double CrystalBall
  float meanR = (xInit<400) ? 5.0 : 100;
  RooRealVar mean("mean","mean of gaussian",0,-5.,meanR) ;
  RooRealVar sigma("sigma","width of gaussian",1.5,0.,30.); 
  RooRealVar a1("a1","a1",1.46,0.,4.);
  RooRealVar n1("n1","n1",1.92,0.,50.);   
  RooRealVar a2("a2","a2",1.46,0.,4.);
  RooRealVar n2("n2","n2",1.92,0.,25.);   
  RooDoubleCB DCBall("DCBall","Double Crystal ball",x,mean,sigma,a1,n1,a2,n2);

  //--- Breit-Wigner
  float bwSigmaMax,bwSigmaMin;
  if(massBin<400) bwSigmaMin=bwSigmaMax=bwSigma;
  else { 
    bwSigmaMin=bwSigma-20.; 
    bwSigmaMax=bwSigma+20.; 
  }
  RooRealVar mean3("mean3","mean3",xInit) ;
  RooRealVar sigma3("sigma3","width3",bwSigma,bwSigmaMin,bwSigmaMax); 
  RooRealVar scale3("scale3","scale3 ",1.); 
  RooRelBWUFParam bw("bw","bw",x,mean3,scale3);

  x.setBins(10000,"fft");
  RooFFTConvPdf model("model","model",x,bw,DCBall);

  

  model.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8));

  stringstream frameTitle;
  if(channels==0){frameTitle << "4#mu, m_{H} = ";}
  if(channels==1){frameTitle << "4e, m_{H} = ";}
  if(channels==2){frameTitle << "2#mu 2e, m_{H} = ";}
  frameTitle << massBin << " GeV";

  RooPlot* xframe = x.frame(Title(frameTitle.str().c_str() )) ;
  dataset.plotOn(xframe,DataError(RooAbsData::SumW2) );
  int col;
  if(channels==0) col=kOrange+7;
  if(channels==1) col=kAzure+2;
  if(channels==2) col=kGreen+3;
  model.plotOn(xframe,LineColor(col));
  //model.paramOn(xframe);


  stringstream nameFile;
  nameFile << "fitM" << massBin << "_channel" << channels << ".pdf";
  xframe->Draw(); gPad->Update(); gPad->Print(nameFile.str().c_str());


  if(fitValues!=0){
    fitValues[0] = a1.getVal();
    fitValues[1] = a2.getVal();
    fitValues[2] = mean.getVal();
    fitValues[3] = mean3.getVal();
    fitValues[4] = n1.getVal();
    fitValues[5] = n2.getVal();
    fitValues[6] = sigma.getVal();
  }  

  if(fitErrors!=0){
    fitErrors[0] = a1.getError();
    fitErrors[1] = a2.getError();
    fitErrors[2] = mean.getError();
    fitErrors[3] = mean3.getError();
    fitErrors[4] = n1.getError();
    fitErrors[5] = n2.getError();
    fitErrors[6] = sigma.getError();
  }

}

void fitSignalHighMassShapeW(int massBin,int id, int channels, int year, 
                             float lumi, bool doSfLepton,double rangeLow, double rangeHigh,
                             double bwSigma,
                             double fitValues[8], double fitErrors[8]){
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
  ggFileName << "/cmsrm/pc21_2/emanuele/data/hzz4l/HZZ4L_53X_S1_V10_S2_V06/MC/hzzTree_id" << id << ".root"; 
  
  cout << "Using " << ggFileName.str() << endl;

  TFile* ggFile = new TFile(ggFileName.str().c_str()); 

  TTree* ggTree = (TTree*) ggFile->Get("zz4lTree/probe_tree");


  TH1F* ggMass = new TH1F("ggMass","ggMass",200,100,300.);


  float nvtx,channel;
  //float ntrueint;
  float nTrueInt,nObsInt,rho,mass,genhiggsmass,m4l,trueInt;

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
  ggTree->SetBranchAddress("genhiggsmass",&genhiggsmass);
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

  //--- rooFit part
  double xMin,xMax,xInit;
  xInit = (double) massBin;
  xMin = rangeLow;
  xMax = rangeHigh ;
  cout << "Fit range: [" << xMin << " , " << xMax << "]. Init value = " << xInit << endl;

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
  RooRealVar w("myW","myW",1.0,0.,1000.);
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
  RooRealVar mean("mean","mean of gaussian",0,-20.,20.) ;
  RooRealVar sigma("sigma","width of gaussian",1.5,0.,100.); 
  RooRealVar a1("a1","a1",1.46,0.,4.);
  RooRealVar n1("n1","n1",1.92,0.,50.);
  RooRealVar a2("a2","a2",1.46,0.,4.);
  RooRealVar n2("n2","n2",1.92,0.,25.);      
  RooDoubleCB DCBall("DCBall","Double Crystal ball",x,mean,sigma,a1,n1,a2,n2);

  //--- Breit-Wigner
  float bwSigmaMax,bwSigmaMin;
  if(massBin>=650) bwSigma=100;
  if(massBin<400) bwSigmaMin=bwSigmaMax=bwSigma;
  else { 
    bwSigmaMin=std::max(1.,double(bwSigma-100.)); 
    bwSigmaMax=bwSigma+100.; 
  }
  cout << "BW initial / range = " << bwSigma << " [" << bwSigmaMin << "," << bwSigmaMax << endl;
  RooRealVar mean3("mean3","mean3",xInit) ;
  RooRealVar sigma3("sigma3","width3",bwSigma,bwSigmaMin,bwSigmaMax); 
  RooRealVar scale3("scale3","scale3 ",1.); 
  RooHighMassBW bw("bw","bw",x,mean3,sigma3);

  x.setBins(10000,"fft");
  RooFFTConvPdf model("model","model",x,bw,DCBall);

  

  model.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8));

  stringstream frameTitle;
  if(channels==0){frameTitle << "4#mu, m_{H} = ";}
  if(channels==1){frameTitle << "4e, m_{H} = ";}
  if(channels==2){frameTitle << "2#mu 2e, m_{H} = ";}
  frameTitle << massBin << " GeV";

  RooPlot* xframe = x.frame(Title(frameTitle.str().c_str() )) ;
  dataset.plotOn(xframe,DataError(RooAbsData::SumW2) );
  int col;
  if(channels==0) col=kOrange+7;
  if(channels==1) col=kAzure+2;
  if(channels==2) col=kGreen+3;
  model.plotOn(xframe,LineColor(col));
  //model.paramOn(xframe);


  stringstream nameFile;
  nameFile << "fitM" << massBin << "_channel" << channels << ".pdf";
  xframe->Draw(); gPad->Update(); gPad->Print(nameFile.str().c_str());


  if(fitValues!=0){
    fitValues[0] = a1.getVal();
    fitValues[1] = a2.getVal();
    fitValues[2] = mean.getVal();
    fitValues[3] = mean3.getVal();
    fitValues[4] = n1.getVal();
    fitValues[5] = n2.getVal();
    fitValues[6] = sigma.getVal();
    fitValues[7] = sigma3.getVal();
  }  

  if(fitErrors!=0){
    fitErrors[0] = a1.getError();
    fitErrors[1] = a2.getError();
    fitErrors[2] = mean.getError();
    fitErrors[3] = mean3.getError();
    fitErrors[4] = n1.getError();
    fitErrors[5] = n2.getError();
    fitErrors[6] = sigma.getError();
    fitErrors[7] = sigma3.getError();
  }

}

float getFitEdge(float mass, float width, bool low) {
  double windowVal = max(width, float(1.));
  double lowside = (mass >= 275) ? 180. : 100.;
  double highside = (mass >= 650) ? 1500. : 800.;
  if (low) return std::max((mass - 20.*windowVal), lowside);
  else return std::min((mass + 15.*windowVal), highside);
}


float getFitEdgeHighMass(float mass, float width, bool low) {
  if (low) return std::max(200.,double(mass-5*width));
  else return std::min(1500.,double(mass+5*width));
}

