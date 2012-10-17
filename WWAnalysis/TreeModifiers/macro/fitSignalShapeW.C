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
#include "WWAnalysis/TreeModifiers/interface/HiggsMassWeightProvider.h"
#include "WWAnalysis/TreeModifiers/src/RooHighMassBW.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "Math/MinimizerOptions.h"
#include "SignalInterpolationStrings.h"

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

void fitSignalShapeW(int massBin, int id, int ch, int year,
		     float lumi, bool doSfLepton,double rangeLow, double rangeHigh,
		     double bwSigma,
		     double fitValues[7], double fitErrors[7], double covQual[1]);

void fitSignalHighMassShapeW(int massBin, int id, int ch, int year,
                             float lumi, bool doSfLepton,double rangeLow, double rangeHigh,
                             double bwSigma,
                             double fitValues[8], double fitErrors[8], double covQual[1], int nSigmaUp=0);

void allHighMass(int channels=0,int year=2012, bool doSfLepton=true, int nSigmaUp=0);


/// validation methods
void validateInterpolation(int massBin,int id, int channels, int year, 
                           float lumi, bool doSfLepton,double rangeLow, double rangeHigh,
                           double bwSigma);
void validateInterpolationHighMass(int massBin,int id, int channels, int year, 
                                   float lumi, bool doSfLepton,double rangeLow, double rangeHigh,
                                   double bwSigma, int nSigmaUp=0);
void validateall(int channels=0,int year=2012, bool doSfLepton=true);
void validateallHighMass(int channels=0,int year=2012, bool doSfLepton=true, int nSigmaUp=0);
////

bool do7TeV;
bool doFFT;

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
    init(true);      //  0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21
    float masses[22] = {115,120,124,125,126,130,140,150,160,170,180,190,200,210,220,230,250,275,300,325,350,375};
    for(int i=0;i<22;++i) {
      mass[i] = masses[i]; 
      id[i]=1000+masses[i]; 
      float width = xsecs.getHZZ4lWidth(masses[i]);
      xLow[i] = getFitEdge(masses[i],width,true); 
      xHigh[i] = getFitEdge(masses[i],width,false); 
      //cout << "For mass = " << masses[i] << " width = " << width << "; => Fit Range = [" << xLow[i] << "," << xHigh[i] << "]" << endl;
      bwSigma[i] = width;
    }
    maxMassBin = 22;
  }


  if(year==2012){
    init(false);     //  0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28
    float masses[29] = {115,117,119,120,121,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,220,250,275,300,325,350,375};
    for(int i=0;i<29;++i) {
      mass[i] = masses[i]; 
      id[i]=1000+masses[i]; 
      float width = xsecs.getHZZ4lWidth(masses[i]);
      xLow[i] = getFitEdge(masses[i],width,true); 
      xHigh[i] = getFitEdge(masses[i],width,false); 
      //cout << "For mass = " << masses[i] << " width = " << width << "; => Fit Range = [" << xLow[i] << "," << xHigh[i] << "]" << endl;
      bwSigma[i] = width;
    }
    maxMassBin = 29;
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
  double covQualVal[40];

  double fitValues[8];
  double fitErrors[8];
  double covQual[1];

  //fitSignalShapeW(mass[0],channels,10.,doSfLepton,max(95.,xLow[0]*0.90),xHigh[0]*1.1,bwSigma[12],
  // 		  fitValues,fitErrors);  Wait();
  
  for(int i=0; i<maxMassBin;++i){
    fitSignalShapeW(mass[i],id[i],channels,year,10.,doSfLepton,xLow[i],xHigh[i],bwSigma[i],
		    fitValues,fitErrors,covQual);  
  
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

    cout << "covQual of the fit: " << covQual[0] << endl;
    covQualVal[i] = covQual[0];

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
  TGraph* gCovQual = new TGraph(maxMassBin,massV,covQualVal);

  gA1->SetName("gA1");  gA1->SetMarkerStyle(20);   gA1->SetMarkerSize(1);
  gA2->SetName("gA2");  gA2->SetMarkerStyle(20);   gA2->SetMarkerSize(1);
  gN1->SetName("gN1");  gN1->SetMarkerStyle(20);   gN1->SetMarkerSize(1);
  gN2->SetName("gN2");  gN2->SetMarkerStyle(20);   gN2->SetMarkerSize(1);
  gMeanCB->SetName("gMeanCB"); gMeanCB->SetMarkerStyle(20);   gMeanCB->SetMarkerSize(1);
  gSigmaCB->SetName("gSigmaCB"); gSigmaCB->SetMarkerStyle(20);   gSigmaCB->SetMarkerSize(1);
  gMeanBW->SetName("gMeanBW"); gMeanBW->SetMarkerStyle(20);   gMeanBW->SetMarkerSize(1);
  gCovQual->SetName("gCovQual"); gCovQual->SetMarkerStyle(20);   gCovQual->SetMarkerSize(1);
  
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

  gCovQual->SetTitle("");
  gCovQual->GetXaxis()->SetTitle("mass (GeV)");
  gCovQual->GetYaxis()->SetTitle("cov. matrix qual.");
  
  resultfile->cd();
  gA1->Fit("pol0"); gA1->Draw("Ap"); gA1->Write();
  gA2->Fit("pol0"); gA2->Draw("Ap"); gA2->Write();
  gN1->Fit("pol1"); gN1->Draw("Ap"); gN1->Write();
  gN2->Fit("pol1"); gN2->Draw("Ap"); gN2->Write();
  gMeanCB->Fit("pol1"); gMeanCB->Draw("Ap"); gMeanCB->Write();
  gSigmaCB->Fit("pol1"); gSigmaCB->Draw("Ap"); gSigmaCB->Write();
  gCovQual->Write();

  resultfile->Close();

}

void allHighMass(int channels,int year, bool doSfLepton, int nSigmaUp){
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
    init(true);          //  0   1   2   3   4   5   6   7   8   9   10  11  12  13  14
    float highmasses[15] = {400,425,450,475,500,525,550,575,600,650,700,750,800,900,950};
    for(int i=0;i<15;++i) {
      mass[i] = highmasses[i]; 
      if(mass[i]<1000) id[i]=1000+mass[i]; 
      else id[i]=11000;
      float width = xsecs.getHZZ4lWidth(mass[i]);
      xLow[i] = getFitEdgeHighMass(mass[i],width,true); 
      xHigh[i] = getFitEdgeHighMass(mass[i],width,false); 
      bwSigma[i] = width;
    }
    maxMassBin = 15;
  }


  if(year==2012){
    init(false);   
                         //  0   1   2   3   4   5   6   7   8   9   10  11   12  
    float highmasses[13] = {400,450,500,550,600,650,700,750,800,850,900,950,1000};
    for(int i=0;i<13;++i) {
      mass[i] = highmasses[i]; 
      if(mass[i]<1000) id[i]=1000+mass[i]; 
      else id[i]=11000;
      float width = xsecs.getHZZ4lWidth(mass[i]);
      xLow[i] = getFitEdgeHighMass(mass[i],width,true); 
      xHigh[i] = getFitEdgeHighMass(mass[i],width,false); 
      //cout << "For mass = " << masses[i] << " width = " << width << "; => Fit Range = [" << xLow[i] << "," << xHigh[i] << "]" << endl;
      bwSigma[i] = width;
    }
    maxMassBin = 13;
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
  double covQualVal[40];

  double fitValues[8];
  double fitErrors[8];
  double covQual[1];

  //fitSignalShapeW(mass[0],channels,10.,doSfLepton,max(95.,xLow[0]*0.90),xHigh[0]*1.1,bwSigma[12],
  // 		  fitValues,fitErrors);  Wait();
  
  for(int i=0; i<maxMassBin;++i){

    fitSignalHighMassShapeW(mass[i],id[i],channels,year,10.,doSfLepton,xLow[i],xHigh[i],bwSigma[i],
                            fitValues,fitErrors,covQual,nSigmaUp);  
  
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

    cout << "covQual of the fit: " << covQual[0] << endl;
    covQualVal[i] = covQual[0];

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
  TGraph* gCovQual = new TGraph(maxMassBin,massV,covQualVal);

  gA1->SetName("gA1");  gA1->SetMarkerStyle(20);   gA1->SetMarkerSize(1);
  gA2->SetName("gA2");  gA2->SetMarkerStyle(20);   gA2->SetMarkerSize(1);
  gN1->SetName("gN1");  gN1->SetMarkerStyle(20);   gN1->SetMarkerSize(1);
  gN2->SetName("gN2");  gN2->SetMarkerStyle(20);   gN2->SetMarkerSize(1);
  gMeanCB->SetName("gMeanCB"); gMeanCB->SetMarkerStyle(20);   gMeanCB->SetMarkerSize(1);
  gSigmaCB->SetName("gSigmaCB"); gSigmaCB->SetMarkerStyle(20);   gSigmaCB->SetMarkerSize(1);
  gMeanBW->SetName("gMeanBW"); gMeanBW->SetMarkerStyle(20);   gMeanBW->SetMarkerSize(1);
  gSigmaBW->SetName("gSigmaBW"); gSigmaBW->SetMarkerStyle(20);   gSigmaBW->SetMarkerSize(1);
  gCovQual->SetName("gCovQual"); gCovQual->SetMarkerStyle(20);   gCovQual->SetMarkerSize(1);
  
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
  
  gCovQual->SetTitle("");
  gCovQual->GetXaxis()->SetTitle("mass (GeV)");
  gCovQual->GetYaxis()->SetTitle("cov. matrix qual.");

  resultfile->cd();
  gA1->Fit("pol0"); gA1->Draw("Ap"); gA1->Write();
  gA2->Fit("pol0"); gA2->Draw("Ap"); gA2->Write();
  gN1->Fit("pol1"); gN1->Draw("Ap"); gN1->Write();
  gN2->Fit("pol1"); gN2->Draw("Ap"); gN2->Write();
  gMeanCB->Fit("pol1"); gMeanCB->Draw("Ap"); gMeanCB->Write();
  gSigmaCB->Fit("pol1"); gSigmaCB->Draw("Ap"); gSigmaCB->Write(); 
  gSigmaBW->Fit("pol1"); gSigmaBW->Draw("Ap"); gSigmaBW->Write();
  gCovQual->Write();
 
  resultfile->Close();

}

void validateall(int channels,int year, bool doSfLepton){
  /*
    channels = 0 --> 4mu
    channels = 1 --> 4el
    channels = 2 --> 2mu2el

   */
  //init();

  do7TeV = (year==2012) ? false : true;
  doFFT = true;

  double bwSigma[40];
  int mass[40]; int id[40]; double xLow[40]; double xHigh[40];  
  int maxMassBin;

  XSecProvider xsecs;
  xsecs.initHiggs4lWidth();

  if(year==2011){
    init(true);
    float masses[22] = {115,120,124,125,126,130,140,150,160,170,180,190,200,210,220,230,250,275,300,325,350,375};
    for(int i=0;i<22;++i) {
      mass[i] = masses[i]; 
      id[i]=1000+masses[i]; 
      float width = xsecs.getHZZ4lWidth(masses[i]);
      xLow[i] = getFitEdge(masses[i],width,true); 
      xHigh[i] = getFitEdge(masses[i],width,false); 
      //cout << "For mass = " << masses[i] << " width = " << width << "; => Fit Range = [" << xLow[i] << "," << xHigh[i] << "]" << endl;
      bwSigma[i] = width;
    }
    maxMassBin = 22;
  }


  if(year==2012){
    init(false);     //  0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28 
    float masses[29] = {115,117,119,120,121,123,124,125,126,127,128,129,130,135,140,145,150,160,170,180,190,200,220,250,275,300,325,350,375};
    for(int i=0;i<29;++i) {
      mass[i] = masses[i]; 
      id[i]=1000+masses[i]; 
      float width = xsecs.getHZZ4lWidth(masses[i]);
      xLow[i] = getFitEdge(masses[i],width,true); 
      xHigh[i] = getFitEdge(masses[i],width,false); 
      //cout << "For mass = " << masses[i] << " width = " << width << "; => Fit Range = [" << xLow[i] << "," << xHigh[i] << "]" << endl;
      bwSigma[i] = width;
    }
    maxMassBin = 29;
  }
  // -----------------------



  for(int i=0; i<maxMassBin;++i){
    // skip the gg masses not yet ready
    validateInterpolation(mass[i],id[i],channels,year,10.,doSfLepton,xLow[i],xHigh[i],bwSigma[i]);  
  }
  
}

void validateallHighMass(int channels,int year, bool doSfLepton, int nSigmaUp){
  /*
    channels = 0 --> 4mu
    channels = 1 --> 4el
    channels = 2 --> 2mu2el

   */
  //init();

  do7TeV = (year==2012) ? false : true;
  doFFT = true;

  double bwSigma[40];
  int mass[40]; int id[40]; double xLow[40]; double xHigh[40];  
  int maxMassBin;

  XSecProvider xsecs;
  xsecs.initHiggs4lWidth();

  if(year==2011){
    init(true);          //  0   1   2   3   4   5   6   7   8   9   10  11  12  13  14
    float highmasses[15] = {400,425,450,475,500,525,550,575,600,650,700,750,800,900,950};
    for(int i=0;i<15;++i) {
      mass[i] = highmasses[i]; 
      if(mass[i]<1000) id[i]=1000+mass[i]; 
      else id[i]=11000;
      float width = xsecs.getHZZ4lWidth(mass[i]);
      xLow[i] = getFitEdgeHighMass(mass[i],width,true); 
      xHigh[i] = getFitEdgeHighMass(mass[i],width,false); 
      bwSigma[i] = width;
    }
    maxMassBin = 15;
  }


  if(year==2012){
    init(false);   
                         //  0   1   2   3   4   5   6   7   8   9   10  11  12  
    float highmasses[13] = {400,450,500,550,600,650,700,750,800,850,900,950,1000};
    for(int i=0;i<13;++i) {
      mass[i] = highmasses[i]; 
      if(mass[i]<1000) id[i]=1000+mass[i]; 
      else id[i]=11000;
      float width = xsecs.getHZZ4lWidth(mass[i]);
      xLow[i] = getFitEdgeHighMass(mass[i],width,true); 
      xHigh[i] = getFitEdgeHighMass(mass[i],width,false); 
      //cout << "For mass = " << masses[i] << " width = " << width << "; => Fit Range = [" << xLow[i] << "," << xHigh[i] << "]" << endl;
      bwSigma[i] = width;
    }
    maxMassBin = 13;
  }
  // -----------------------



  for(int i=0; i<maxMassBin;++i){
    // skip the gg masses not yet ready
    validateInterpolationHighMass(mass[i],id[i],channels,year,10.,doSfLepton,xLow[i],xHigh[i],bwSigma[i],nSigmaUp);  
  }
  
}


void fitSignalShapeW(int massBin,int id, int channels, int year, 
		     float lumi, bool doSfLepton,double rangeLow, double rangeHigh,
		     double bwSigma,
		     double fitValues[7], double fitErrors[7], double covQual[1]){
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

  stringstream ggFileName,vbfFileName;
  if(year==2012) ggFileName << "/cmsrm/pc21_2/emanuele/data/hzz4l/HZZ4L_53X_S1_V11_S2_V02/MC/hzzTree_id" << id << ".root"; 
  else if(year==2011) ggFileName << "/cmsrm/pc21_2/emanuele/data/hzz4l/HZZ4L_42X_S1_V11_S2_V02/MC/7TeV/yesRegrYesCalibYesMu/hzzTree_id" << id << ".root";
  else {
    cout << "Wrong year." << endl;
    return;
  }
    

  cout << "Using " << ggFileName.str() << endl;

  TFile* ggFile = new TFile(ggFileName.str().c_str()); 

  TTree* ggTree = (TTree*) ggFile->Get("zz4lTree/probe_tree");


  float channel;
  //float ntrueint;
  float nTrueInt,rho,mass,m4l;

  float pt1,eta1,id1,pt2,eta2,id2,pt3,eta3,id3,pt4,eta4,id4;
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

    ntupleVarSet.setRealValue("myW",localW);
    if(x.getVal()>xMin && x.getVal()<xMax)
      dataset.add(ntupleVarSet, localW);
    
  }
  //---------

  //cout << "dataset n entries: " << dataset.sumEntries() << endl;
  //cout << "reduced data n entries: " << reducedData->sumEntries() << endl;


  //--- double CrystalBall
  RooRealVar mean("bias","mean of gaussian",0,-5.,5.) ;
  RooRealVar sigma("sigma","width of gaussian",1.5,0.,30.); 
  RooRealVar a1("a1","a1",1.46,0.5,5.);
  RooRealVar n1("n1","n1",1.92,0.,10.);   
  RooRealVar a2("a2","a2",1.46,1.,10.);
  RooRealVar n2("n2","n2",20,1.,50.);   
  n2.setConstant(kTRUE);
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

  

  RooFitResult *fitres = (RooFitResult*)model.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8),Save(true));

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

//   RooArgSet selParms(sigma);
//   model.paramOn(xframe,Parameters(selParms));


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

  covQual[0] = fitres->covQual();

}

void fitSignalHighMassShapeW(int massBin,int id, int channels, int year, 
                             float lumi, bool doSfLepton,double rangeLow, double rangeHigh,
                             double bwSigma,
                             double fitValues[8], double fitErrors[8], double covQual[1], int nSigmaUp){

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


  stringstream ggFileName,vbfFileName;
  if(year==2012) ggFileName << "/cmsrm/pc21_2/emanuele/data/hzz4l/HZZ4L_53X_S1_V11_S2_V02/MC/hzzTree_id" << id << ".root"; 
  else if(year==2011) ggFileName << "/cmsrm/pc21_2/emanuele/data/hzz4l/HZZ4L_42X_S1_V11_S2_V02/MC/7TeV/yesRegrYesCalibYesMu/hzzTree_id" << id << ".root";
  else {
    cout << "Wrong year." << endl;
    return;
  }
  
  cout << "Using " << ggFileName.str() << endl;

  TFile* ggFile = new TFile(ggFileName.str().c_str()); 

  TTree* ggTree = (TTree*) ggFile->Get("zz4lTree/probe_tree");


  float channel;
  //float ntrueint;
  float nTrueInt,rho,mass,genhiggsmass,m4l;

  float pt1,eta1,id1,pt2,eta2,id2,pt3,eta3,id3,pt4,eta4,id4;
  float lsW;
  float sfLepton(1.);

  const char *base=getenv("CMSSW_BASE");
  std::string fileweights(base);
  stringstream asciimass;
  asciimass << massBin;
  if(year==2011) fileweights += "/src/WWAnalysis/AnalysisStep/data/HiggsMassReweighting/mZZ_Higgs"+asciimass.str()+"_7TeV_Lineshape+Interference.txt";
  else fileweights += "/src/WWAnalysis/AnalysisStep/data/HiggsMassReweighting/mZZ_Higgs"+asciimass.str()+"_8TeV_Lineshape+Interference.txt";
  cout << "Weighting mass with file: " << fileweights << endl;
  HiggsMassWeightProvider hmWProvider( fileweights.c_str() );


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
    if(nSigmaUp==1) localW = localW * hmWProvider.getWeightUp(genhiggsmass);
    else if(nSigmaUp==-1) localW = localW * hmWProvider.getWeightDown(genhiggsmass);
    else localW = localW * hmWProvider.getWeight(genhiggsmass);

    ntupleVarSet.setRealValue("myW",localW);
    if(x.getVal()>xMin && x.getVal()<xMax)
      dataset.add(ntupleVarSet, localW);
    
  }
  //---------

  //cout << "dataset n entries: " << dataset.sumEntries() << endl;
  //cout << "reduced data n entries: " << reducedData->sumEntries() << endl;


  //--- simple CrystalBall
  float meanR = (xInit<400) ? 5.0 : 100;
  
  // Gamma is floating. Fix the CB sigma to the extrapolation at 400 GeV
  float sigmaCBval;
  if(channels==0) sigmaCBval=4.0;
  if(channels==1) sigmaCBval=2.0;
  if(channels==2) sigmaCBval=3.5;
  RooRealVar mean("mean","mean of gaussian",0,-20.,20.) ;
  RooRealVar sigma("sigma","width of gaussian",sigmaCBval,0,200); 
  RooRealVar a1("a1","a1",1.46,0.5,5.);
  RooRealVar n1("n1","n1",1.92,0.,10.);   
  RooRealVar a2("a2","a2",1.46,1.,10.);
  RooRealVar n2("n2","n2",20,1.,50.);   
  n2.setConstant(kTRUE);

  RooDoubleCB DCBall("DCBall","Double Crystal ball",x,mean,sigma,a1,n1,a2,n2);

  // to evaluate the systematics on the BW, let's fix the resolution function to the nominal one
  RooRealVar masshiggs       ("MH", "", massBin);
  RooRealVar zerosyst        ("zerosyst", "", 0);
  
  RooFormulaVar ggh_mean_CB  ("sig_ggh_mean_CB"  , getSignalCBMeanString (massBin, channels, do7TeV, doFFT).c_str()                             , RooArgList(masshiggs,zerosyst));
  RooFormulaVar ggh_sigma_CB ("sig_ggh_sigma_CB" , getSignalCBSigmaString(massBin, channels, do7TeV).c_str()                             , RooArgList(masshiggs,zerosyst));
  RooFormulaVar ggh_alphaL   ("sig_ggh_alphaL"   , getSignalCBAlphaLString(massBin, channels, do7TeV).c_str()                            , RooArgList(masshiggs));
  RooFormulaVar ggh_alphaR   ("sig_ggh_alphaR"   , getSignalCBAlphaRString(massBin, channels, do7TeV).c_str()                            , RooArgList(masshiggs));
  RooFormulaVar ggh_nL       ("sig_ggh_nL"       , getSignalCBNLString(massBin, channels, do7TeV).c_str()                                , RooArgList(masshiggs));
  RooFormulaVar ggh_nR       ("sig_ggh_nR"       , getSignalCBNRString(massBin, channels, do7TeV).c_str()                                , RooArgList(masshiggs));

  if(nSigmaUp!=0) {
    
    float media= ggh_mean_CB.getVal();
    cout << "media = " << media << endl;
    
    mean.setVal(ggh_mean_CB.getVal());
    sigma.setVal(ggh_sigma_CB.getVal());
    a1.setVal(ggh_alphaL.getVal()); 
    a2.setVal(ggh_alphaR.getVal());
    n1.setVal(ggh_nL.getVal()); 
    n2.setVal(ggh_nR.getVal());

    mean.setConstant(kTRUE);
    sigma.setConstant(kTRUE);
    a1.setConstant(kTRUE);
    a2.setConstant(kTRUE);
    n1.setConstant(kTRUE);
    n2.setConstant(kTRUE);

  }

  //--- Breit-Wigner
  RooRealVar mean3("mean3","mean3",xInit) ;
  float bwSigmaMin = (massBin<600) ? 40. : 0.;
  RooRealVar sigma3("sigma3","width3",bwSigma,40,1000); 
  RooRealVar scale3("scale3","scale3 ",1.); 
  RooHighMassBW bw("bw","bw",x,mean3,sigma3);

  x.setBins(10000,"fft");
  RooFFTConvPdf model("model","model",x,bw,DCBall);

  // this is constrain the resolution to be lower than the Higgs width (know a priori for mH>400 GeV)
  RooFormulaVar penalty("penalty","100*(@0>@1)",RooArgSet(sigma,sigma3));
  RooAbsReal* nll = model.createNLL(dataset,Range(xMin,xMax),NumCPU(8));
  
  RooAddition nllp("nllp","nllp",RooArgSet(*nll,penalty));
  RooMinuit m(nllp) ;
  m.migrad();
  m.hesse();

  //  RooFitResult *fitres = model.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8),Save(true));

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

  //  covQual[0] = fitres->covQual();
  covQual[0] = -1;
}



void validateInterpolation(int massBin,int id, int channels, int year, 
                           float lumi, bool doSfLepton,double rangeLow, double rangeHigh,
                           double bwSigma) {
   // DUPLICATION OF CODE HERE! BUT IT IS NOT POG, SO IT IS FINE...
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

  stringstream ggFileName,vbfFileName;
  if(year==2012) ggFileName << "/cmsrm/pc21_2/emanuele/data/hzz4l/HZZ4L_53X_S1_V11_S2_V02/MC/hzzTree_id" << id << ".root"; 
  else if(year==2011) ggFileName << "/cmsrm/pc21_2/emanuele/data/hzz4l/HZZ4L_42X_S1_V11_S2_V02/MC/7TeV/yesRegrYesCalibYesMu/hzzTree_id" << id << ".root";
  else {
    cout << "Wrong year." << endl;
    return;
  }    

  cout << "Using " << ggFileName.str() << endl;

  TFile* ggFile = new TFile(ggFileName.str().c_str()); 

  TTree* ggTree = (TTree*) ggFile->Get("zz4lTree/probe_tree");


  float channel;
  //float ntrueint;
  float nTrueInt,rho,mass,m4l;

  float pt1,eta1,id1,pt2,eta2,id2,pt3,eta3,id3,pt4,eta4,id4;
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

    ntupleVarSet.setRealValue("myW",localW);
    if(x.getVal()>xMin && x.getVal()<xMax)
      dataset.add(ntupleVarSet, localW);
    
  }
  //---------

  RooRealVar masshiggs       ("MH", "", massBin);
  RooRealVar zerosyst        ("zerosyst", "", 0);

  RooFormulaVar ggh_mean_CB  ("sig_ggh_mean_CB"  , getSignalCBMeanString (massBin, channels, do7TeV, doFFT).c_str()                             , RooArgList(masshiggs,zerosyst));
  RooFormulaVar ggh_sigma_CB ("sig_ggh_sigma_CB" , getSignalCBSigmaString(massBin, channels, do7TeV).c_str()                             , RooArgList(masshiggs,zerosyst));
  RooFormulaVar ggh_alphaL   ("sig_ggh_alphaL"   , getSignalCBAlphaLString(massBin, channels, do7TeV).c_str()                            , RooArgList(masshiggs));
  RooFormulaVar ggh_alphaR   ("sig_ggh_alphaR"   , getSignalCBAlphaRString(massBin, channels, do7TeV).c_str()                            , RooArgList(masshiggs));
  RooFormulaVar ggh_nL       ("sig_ggh_nL"       , getSignalCBNLString(massBin, channels, do7TeV).c_str()                                , RooArgList(masshiggs));
  RooFormulaVar ggh_nR       ("sig_ggh_nR"       , getSignalCBNRString(massBin, channels, do7TeV).c_str()                                , RooArgList(masshiggs));

  RooDoubleCB  signalCB_ggH   ("signalCB_ggH", "", x, ggh_mean_CB,ggh_sigma_CB,ggh_alphaL,ggh_nL,ggh_alphaR,ggh_nR);

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
  RooFFTConvPdf model("model","model",x,bw,signalCB_ggH);


  cout << "ggh_mean_CB = " << ggh_mean_CB.getVal() << endl;
  cout << "sig_ggh_nR = " << ggh_nR.getVal() << endl;
  cout << "sig_ggh_alphaR = " << ggh_alphaR.getVal() << endl;

  stringstream frameTitle;
  if(channels==0){frameTitle << "4#mu, m_{H} = ";}
  if(channels==1){frameTitle << "4e, m_{H} = ";}
  if(channels==2){frameTitle << "2#mu 2e, m_{H} = ";}
  frameTitle << massBin << " GeV.         NOT A FIT: INTERPOLATED PDF";

  RooPlot* xframe = x.frame(Title(frameTitle.str().c_str() )) ;
  dataset.plotOn(xframe,DataError(RooAbsData::SumW2) );
  int col;
  if(channels==0) col=kOrange+7;
  if(channels==1) col=kAzure+2;
  if(channels==2) col=kGreen+3;
  model.plotOn(xframe,LineColor(col),LineStyle(kDashed));

  stringstream nameFile;
  nameFile << "interpolPdfM" << massBin << "_channel" << channels << ".pdf";
  xframe->Draw(); gPad->Update(); gPad->Print(nameFile.str().c_str());
 

}

void validateInterpolationHighMass(int massBin,int id, int channels, int year, 
                                   float lumi, bool doSfLepton,double rangeLow, double rangeHigh,
                                   double bwSigma, int nSigmaUp) {

  // DUPLICATION OF CODE HERE! BUT IT IS NOT POG, SO IT IS FINE...
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

  stringstream ggFileName,vbfFileName;
  if(year==2012) ggFileName << "/cmsrm/pc21_2/emanuele/data/hzz4l/HZZ4L_53X_S1_V11_S2_V02/MC/hzzTree_id" << id << ".root"; 
  else if(year==2011) ggFileName << "/cmsrm/pc21_2/emanuele/data/hzz4l/HZZ4L_42X_S1_V11_S2_V02/MC/7TeV/yesRegrYesCalibYesMu/hzzTree_id" << id << ".root";
  else {
    cout << "Wrong year." << endl;
    return;
  }

  cout << "Using " << ggFileName.str() << endl;

  TFile* ggFile = new TFile(ggFileName.str().c_str()); 

  TTree* ggTree = (TTree*) ggFile->Get("zz4lTree/probe_tree");


  float channel;
  //float ntrueint;
  float nTrueInt,rho,mass,m4l,genhiggsmass;

  float pt1,eta1,id1,pt2,eta2,id2,pt3,eta3,id3,pt4,eta4,id4;
  float sfLepton(1.);

  const char *base=getenv("CMSSW_BASE");
  std::string fileweights(base);
  stringstream asciimass;
  asciimass << massBin;
  if(year==2011) fileweights += "/src/WWAnalysis/AnalysisStep/data/HiggsMassReweighting/mZZ_Higgs"+asciimass.str()+"_7TeV_Lineshape+Interference.txt";
  else fileweights += "/src/WWAnalysis/AnalysisStep/data/HiggsMassReweighting/mZZ_Higgs"+asciimass.str()+"_8TeV_Lineshape+Interference.txt";
  cout << "Weighting mass with file: " << fileweights << endl;
  HiggsMassWeightProvider hmWProvider( fileweights.c_str() );

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

  //--- rooFit part
  double xMin,xMax,xInit;
  xInit = (double) massBin;
  xMin = rangeLow;
  xMax = rangeHigh ;
  cout << "Fit range: [" << xMin << " , " << xMax << "]. Init value = " << xInit << endl;

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
    if(nSigmaUp==1) localW = localW * hmWProvider.getWeightUp(genhiggsmass);
    else if(nSigmaUp==-1) localW = localW * hmWProvider.getWeightDown(genhiggsmass);
    else localW = localW * hmWProvider.getWeight(genhiggsmass);

    ntupleVarSet.setRealValue("myW",localW);
    if(x.getVal()>xMin && x.getVal()<xMax)
      dataset.add(ntupleVarSet, localW);
    
  }
  //---------

  RooRealVar masshiggs       ("MH", "", massBin);
  RooRealVar zerosyst        ("zerosyst", "", 0);

  RooFormulaVar ggh_mean_CB  ("sig_ggh_mean_CB"  , getSignalCBMeanString (massBin, channels, do7TeV, doFFT).c_str()                             , RooArgList(masshiggs,zerosyst));
  RooFormulaVar ggh_sigma_CB ("sig_ggh_sigma_CB" , getSignalCBSigmaString(massBin, channels, do7TeV).c_str()                             , RooArgList(masshiggs,zerosyst));
  RooFormulaVar ggh_alphaL   ("sig_ggh_alphaL"   , getSignalCBAlphaLString(massBin, channels, do7TeV).c_str()                            , RooArgList(masshiggs));
  RooFormulaVar ggh_alphaR   ("sig_ggh_alphaR"   , getSignalCBAlphaRString(massBin, channels, do7TeV).c_str()                            , RooArgList(masshiggs));
  RooFormulaVar ggh_nL       ("sig_ggh_nL"       , getSignalCBNLString(massBin, channels, do7TeV).c_str()                                , RooArgList(masshiggs));
  RooFormulaVar ggh_nR       ("sig_ggh_nR"       , getSignalCBNRString(massBin, channels, do7TeV).c_str()                                , RooArgList(masshiggs));

  RooDoubleCB  signalCB_ggH   ("signalCB_ggH", "", x, ggh_mean_CB,ggh_sigma_CB,ggh_alphaL,ggh_nL,ggh_alphaR,ggh_nR);

  //--- Breit-Wigner
  RooFormulaVar ggh_gamma_BW ("sig_ggh_gamma_BW" , getSignalBWGammaString(massBin, channels, do7TeV).c_str()                             , RooArgList(masshiggs,zerosyst));
  
  RooRealVar mean3("mean3","mean3",xInit) ;
  float bwSigmaMin = (massBin<600) ? 40. : 0.;
  RooRealVar sigma3("sigma3","width3",bwSigma,40,1000); 
  RooRealVar scale3("scale3","scale3 ",1.); 
  RooHighMassBW bw("bw","bw",x,mean3,ggh_gamma_BW);

  x.setBins(10000,"fft");
  RooFFTConvPdf model("model","model",x,bw,signalCB_ggH);


  cout << "ggh_mean_CB = " << ggh_mean_CB.getVal() << endl;
  cout << "sig_ggh_nR = " << ggh_nR.getVal() << endl;
  cout << "sig_ggh_alphaR = " << ggh_alphaR.getVal() << endl;

  stringstream frameTitle;
  if(channels==0){frameTitle << "4#mu, m_{H} = ";}
  if(channels==1){frameTitle << "4e, m_{H} = ";}
  if(channels==2){frameTitle << "2#mu 2e, m_{H} = ";}
  frameTitle << massBin << " GeV.         NOT A FIT: INTERPOLATED PDF";

  RooPlot* xframe = x.frame(Title(frameTitle.str().c_str() )) ;
  dataset.plotOn(xframe,DataError(RooAbsData::SumW2) );
  int col;
  if(channels==0) col=kOrange+7;
  if(channels==1) col=kAzure+2;
  if(channels==2) col=kGreen+3;
  model.plotOn(xframe,LineColor(col),LineStyle(kDashed));

  stringstream nameFile;
  nameFile << "interpolPdfM" << massBin << "_channel" << channels << ".pdf";
  xframe->Draw(); gPad->Update(); gPad->Print(nameFile.str().c_str());
 

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

