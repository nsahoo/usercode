#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"

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
  filename << "signalPdfs2011/YesRegrYesCalibV14/HighMass/parameters_channel" << channel << ".root";
  TFile *resultfile = TFile::Open(filename.str().c_str());
  TGraphErrors *gA1 = (TGraphErrors*)resultfile->Get("gA1");
  TGraphErrors *gA2 = (TGraphErrors*)resultfile->Get("gA2");
  TGraphErrors *gN1 = (TGraphErrors*)resultfile->Get("gN1");
  TGraphErrors *gN2 = (TGraphErrors*)resultfile->Get("gN2");
  TGraphErrors *gMeanCB = (TGraphErrors*)resultfile->Get("gMeanCB");
  TGraphErrors *gSigmaCB = (TGraphErrors*)resultfile->Get("gSigmaCB");
  TGraphErrors *gSigmaBW = 0;
  if(highmass) gSigmaBW = (TGraphErrors*)resultfile->Get("gSigmaBW");

  gA1->GetYaxis()->SetRangeUser(-1,6);
  gA2->GetYaxis()->SetRangeUser(-1,6);
  gN1->GetYaxis()->SetRangeUser(0,100);
  gN2->GetYaxis()->SetRangeUser(0,100);
  gMeanCB->GetYaxis()->SetRangeUser(-2,2);
  if(highmass) gMeanCB->GetYaxis()->SetRangeUser(-20,10);
  gSigmaCB->GetYaxis()->SetRangeUser(0,5);
  if(highmass) gSigmaCB->GetYaxis()->SetRangeUser(0,200);
  if(highmass) gSigmaBW->GetYaxis()->SetRangeUser(-10,300);

  // skip some mass points where we did not have MC                                                   

  std::vector<int> pointstoskip;
  //  pointstoskip.push_back(1);
  pointstoskip.push_back(2);
  pointstoskip.push_back(3);
  //  pointstoskip.push_back(4);
  //  pointstoskip.push_back(8);
  //pointstoskip.push_back(11);
  //pointstoskip.push_back(14);
  //pointstoskip.push_back(15);

  for(int i=0;i<(int)pointstoskip.size();++i) {
    int p=pointstoskip[i];
    gA1->SetPointError(p,0,1000);
    gA2->SetPointError(p,0,1000);
    gN1->SetPointError(p,0,1000);
    gN2->SetPointError(p,0,1000);
    gMeanCB->SetPointError(p,0,1000);
    gSigmaCB->SetPointError(p,0,1000);
    if(highmass) gSigmaBW->SetPointError(p,0,1000);
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();

  stringstream channame;
  channame << "_channel" << channel << ".pdf";
  gA1->Fit("pol2","","",400,950); gA1->Draw("Ap"); gPad->Update(); gPad->Print((string("gA1")+channame.str()).c_str());
  gA2->Fit("pol2","","",400,950); gA2->Draw("Ap"); gPad->Update(); gPad->Print((string("gA2")+channame.str()).c_str());
  gN1->Fit("pol2","","",400,950); gN1->Draw("Ap"); gPad->Update(); gPad->Print((string("gN1")+channame.str()).c_str());
  gN2->Fit("pol2","","",400,950); gN2->Draw("Ap"); gPad->Update(); gPad->Print((string("gN2")+channame.str()).c_str());
  gMeanCB->Fit("pol2","","",400,950); gMeanCB->Draw("Ap"); gPad->Update(); gPad->Print((string("gMeanCB")+channame.str()).c_str());
  gSigmaCB->Fit("pol1","","",400,950); gSigmaCB->Draw("Ap"); gPad->Update(); gPad->Print((string("gSigmaCB")+channame.str()).c_str());
  if(highmass) gSigmaBW->Fit("pol3","","",400,950); gSigmaBW->Draw("Ap"); gPad->Update(); gPad->Print((string("gSigmaBW")+channame.str()).c_str());
}
