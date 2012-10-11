#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TF1.h"

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

using namespace RooFit;



#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

void interpolateAgain(int channel, bool highmass=false) {

  gStyle->SetOptFit(1111);

  stringstream filename;
  filename << "signalPdfs/YesRegrYesCalibV14/" << ((highmass) ? "HighMass" : "LowMass") << "/parameters_channel" << channel << ".root";
  TFile *resultfile = TFile::Open(filename.str().c_str());
  TGraphErrors *gA1 = (TGraphErrors*)resultfile->Get("gA1");
  TGraphErrors *gA2 = (TGraphErrors*)resultfile->Get("gA2");
  TGraphErrors *gN1 = (TGraphErrors*)resultfile->Get("gN1");
  TGraphErrors *gN2 = (TGraphErrors*)resultfile->Get("gN2");
  TGraphErrors *gMeanCB = (TGraphErrors*)resultfile->Get("gMeanCB");
  TGraphErrors *gSigmaCB = (TGraphErrors*)resultfile->Get("gSigmaCB");
  TGraphErrors *gSigmaBW = 0;
  if(highmass) gSigmaBW = (TGraphErrors*)resultfile->Get("gSigmaBW");

  Int_t n = gA1->GetN();

  Double_t *x = gA1->GetX();  
  Double_t *yA1 = gA1->GetY();  Double_t *yEA1 = gA1->GetEY();
  Double_t *yA2 = gA2->GetY();  Double_t *yEA2 = gA2->GetEY();
  Double_t *yN1 = gN1->GetY();  Double_t *yEN1 = gN1->GetEY();
  Double_t *yN2 = gN2->GetY();  Double_t *yEN2 = gN2->GetEY();
  Double_t *yMeanCB = gMeanCB->GetY();  Double_t *yEMeanCB = gMeanCB->GetEY();
  Double_t *ySigmaCB = gSigmaCB->GetY();  Double_t *yESigmaCB = gSigmaCB->GetEY();
  Double_t *ySigmaBW, *yESigmaBW;
  if(highmass) { ySigmaBW = gSigmaBW->GetY();  yESigmaBW = gSigmaBW->GetEY(); }

  // skip some mass points where we did not have MC                                                   
  std::vector<int> pointstoskip;
  pointstoskip.push_back(900);
  pointstoskip.push_back(1000);

  Int_t m = n-pointstoskip.size();
  TGraphErrors *cA1 = new TGraphErrors(m);
  TGraphErrors *cA2 = new TGraphErrors(m);
  TGraphErrors *cN1 = new TGraphErrors(m);
  TGraphErrors *cN2 = new TGraphErrors(m);
  TGraphErrors *cMeanCB = new TGraphErrors(m);
  TGraphErrors *cSigmaCB = new TGraphErrors(m);
  TGraphErrors *cSigmaBW = 0;
  if(highmass) cSigmaBW = new TGraphErrors(m);

  cA1->SetName("cA1");
  cA2->SetName("cA2");
  cN1->SetName("cN1");
  cN2->SetName("cN2");
  cMeanCB->SetName("cMeanCB");
  cSigmaCB->SetName("cSigmaCB");
  if(highmass) cSigmaBW->SetName("cSigmaBW");

  cA1->GetYaxis()->SetRangeUser(-1,6);
  cA2->GetYaxis()->SetRangeUser(-1,6);
  cN1->GetYaxis()->SetRangeUser(0,100);
  cN2->GetYaxis()->SetRangeUser(0,100);
  cMeanCB->GetYaxis()->SetRangeUser(-2,2);
  if(highmass) cMeanCB->GetYaxis()->SetRangeUser(-20,10);
  cSigmaCB->GetYaxis()->SetRangeUser(0,5);
  if(highmass) cSigmaCB->GetYaxis()->SetRangeUser(0,200);
  if(highmass) cSigmaBW->GetYaxis()->SetRangeUser(-10,300);

  int k=0; 
  float xMin=7000; float xMax=-100;
  for(int i=0;i<n;++i) {
    bool skip=false;
    for(int j=0;j<(int)pointstoskip.size();++j) {
      if(x[i]==int(pointstoskip[j])) { cout << "skipping point " << i << " m(H) = " << x[i] << endl; skip=true; break; }
    }
    if(skip) continue;
    if(x[i]<xMin) xMin=x[i];
    if(x[i]>xMax) xMax=x[i];
    cA1->SetPoint(k,x[i],yA1[i]);            cA1->SetPointError(k,0,yEA1[i]);
    cA2->SetPoint(k,x[i],yA2[i]);            cA2->SetPointError(k,0,yEA2[i]);
    cN1->SetPoint(k,x[i],yN1[i]);            cN1->SetPointError(k,0,yEN1[i]);
    cN2->SetPoint(k,x[i],yN2[i]);            cN2->SetPointError(k,0,yEN2[i]);
    cMeanCB->SetPoint(k,x[i],yMeanCB[i]);    cMeanCB->SetPointError(k,0,yEMeanCB[i]);
    cSigmaCB->SetPoint(k,x[i],ySigmaCB[i]);  cSigmaCB->SetPointError(k,0,yESigmaCB[i]);
    if(highmass) cSigmaBW->SetPoint(k,x[i],ySigmaBW[i]);  cSigmaBW->SetPointError(k,0,yESigmaBW[i]);
    k++;
  }

  cout << k << "   " << m << endl;
  cout << "xMin = " << xMin << "  xMax = " << xMax << endl;

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();

  stringstream channame;
  channame << "_channel" << channel << ".pdf";
  cA1->Fit("pol2","","",xMin,xMax); cA1->Draw("Ap"); gPad->Update(); gPad->Print((string("gA1")+channame.str()).c_str());
  cA2->Fit("pol2","","",xMin,xMax); cA2->Draw("Ap"); gPad->Update(); gPad->Print((string("gA2")+channame.str()).c_str());
  if(!highmass) {
    cN1->Fit("pol1","","",xMin,xMax); cN1->Draw("Ap"); gPad->Update(); gPad->Print((string("gN1")+channame.str()).c_str());
    cN2->Fit("pol1","","",xMin,xMax); cN2->Draw("Ap"); gPad->Update(); gPad->Print((string("gN2")+channame.str()).c_str());
  } else {
    cN1->Fit("pol0","","",xMin,xMax); cN1->Draw("Ap"); gPad->Update(); gPad->Print((string("gN1")+channame.str()).c_str());
    cN2->Fit("pol0","","",xMin,xMax); cN2->Draw("Ap"); gPad->Update(); gPad->Print((string("gN2")+channame.str()).c_str());
  }
  cMeanCB->Fit("pol2","","",xMin,xMax); cMeanCB->Draw("Ap"); gPad->Update(); gPad->Print((string("gMeanCB")+channame.str()).c_str());
  cSigmaCB->Fit("pol1","","",xMin,xMax); cSigmaCB->Draw("Ap"); gPad->Update(); gPad->Print((string("gSigmaCB")+channame.str()).c_str());
  if(highmass) { cSigmaBW->Fit("pol1","","",xMin,xMax); cSigmaBW->Draw("Ap"); gPad->Update(); gPad->Print((string("gSigmaBW")+channame.str()).c_str()); }

  TF1 *fA1 = (TF1*)cA1->GetFunction("pol2");
  TF1 *fA2 = (TF1*)cA2->GetFunction("pol2");
  TF1 *fN1 = (highmass) ? (TF1*)cN1->GetFunction("pol0") : (TF1*)cN1->GetFunction("pol1");
  TF1 *fN2 = (highmass) ? (TF1*)cN2->GetFunction("pol0") : (TF1*)cN2->GetFunction("pol1");
  TF1 *fMeanCB = (TF1*)cMeanCB->GetFunction("pol2");
  TF1 *fSigmaCB = (TF1*)cSigmaCB->GetFunction("pol1");
  TF1 *fSigmaBW;
  if(highmass) fSigmaBW = (TF1*)cSigmaBW->GetFunction("pol1");

  std::vector<string> names;
  std::vector<TF1*> fcns;
  fcns.push_back(fA1);       names.push_back("fA1");
  fcns.push_back(fA2);       names.push_back("fA2");
  fcns.push_back(fN1);       names.push_back("fN1");
  fcns.push_back(fN2);       names.push_back("fN2");
  fcns.push_back(fMeanCB);   names.push_back("fMeanCB");
  fcns.push_back(fSigmaCB);  names.push_back("fSigmaCB");
  if(highmass) fcns.push_back(fSigmaBW);  names.push_back("fSigmaBW");

  for(int fcn=0;fcn<(int)fcns.size();++fcn) {

    std::stringstream ss;
    
    for (int i = 0; i < fcns[fcn]->GetNumberFreeParameters(); i++) {
      if (i != 0) ss << " + (";
      else ss << "(";
      ss << fcns[fcn]->GetParameter(i);
      for (int j = 0; j < i; j++) {
        ss << "*@0";
      }
      ss << ")";
    }
    cout << "RooFormulaVar for " << names[fcn] << " = " << ss.str() << endl;
  }

}
