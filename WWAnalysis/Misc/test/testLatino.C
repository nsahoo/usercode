#ifndef __CINT__
#include "LatinoPlot.C"
#include "TROOT.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"

int main(int argc, const char * argv[])
#endif
{

    gInterpreter->ExecuteMacro("scripts/LatinoStyle.C");
#ifdef __CINT__
    gROOT->LoadMacro("scripts/LatinoPlot.C");
#endif 
    TFile* file = new TFile("hists.160.root", "read");

    LatinoPlot myPlot;
    myPlot.setLumi(1300);
    myPlot.setLabel("#Delta#phi");
    myPlot.setUnits("#circ");
    myPlot.setMCHist(iWW   ,(TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/022.GluGluToWWTo4L"));
    myPlot.setMCHist(iZJets,(TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/004.DYtoMuMu"));
    myPlot.setMCHist(iTop  ,(TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/011.tTtoBLNu"));
    myPlot.setMCHist(iWZ   ,(TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/001.WZtoAny"));
    myPlot.setMCHist(iWJets,(TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/026.WJetsToLNuMad"));

    // Just as a test
    TH1F *t = (TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/022.GluGluToWWTo4L");
    myPlot.setDataHist((TH1F*)t->Clone("data"));

    // and just a test
    t = (TH1F*)file->Get("eventHists/bycut/wwmumuIPLHT/dPhill/06/101160.ggToH160toWWto2L2Nu");
    t->Scale(0.05);
    myPlot.setMCHist(iHWW  ,(TH1F*)t->Clone("hww"));

    TCanvas* canvas = new TCanvas("canvas", "canvas");
    myPlot.Draw(10);  // can pass a rebin 

    canvas->GetFrame()->DrawClone();
    canvas->SaveAs("test.png");
}


