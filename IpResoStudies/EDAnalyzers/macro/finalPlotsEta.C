
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"

#include <sstream>

void setStyle();

void cmsPrel(double intLumi);

void setLabels(int type,
	       TString& xLabel,TString& yLabel,TString& yLabelRes,TString& yLabelErr,
	       TString& text1a,TString& text1b,TString& leg1a,TString& leg1b,
	       TString& text2a,TString& text2b,TString& leg2a,TString& leg2b,
	       double& lRange,double& hRange
	       );

void finalPlotsEta(int type,bool skip3rd=false){
  //-----
  //TString inputFolderName="./setupZ3/";
  //TString outputFolderName="./finalPlotsZ3/";  

  TString xLabel;
  TString yLabel;
  TString yLabelResp;
  TString yLabelErr;
  
  double lRange,hRange;

  TString text1a,text1b,leg1a,leg1b;
  TString text2a,text2b,leg2a,leg2b;

  setLabels(type,
	    xLabel,yLabel,yLabelResp,yLabelErr,
	    text1a,text1b,leg1a,leg1b,
	    text2a,text2b,leg2a,leg2b,
	    lRange,hRange);

  //setStyle();
  //----
  stringstream stream;  stream << type;
  TString counter = stream.str();

  //TString fileName1="./finalPlotsZ1/histos"+counter+".root";
  //TString fileName2="./finalPlotsZ2/histos"+counter+".root";

  //TString fileName1="./finalPlots1GeV/histos"+counter+".root";
  //TString fileName2="./finalPlots3GeV/histos"+counter+".root";
  //TString fileName3="./finalPlots8GeV/histos"+counter+".root";

  TString fileName1="./newer1GeV/histos"+counter+".root";
  TString fileName2="./newer3GeV/histos"+counter+".root";
  TString fileName3="./newer8GeV/histos"+counter+".root";

  TString outputName1="final_type"+counter+".pdf";
  TString outputName2="final_type"+counter+".png";
  TString outputName3="final_type"+counter+".eps";

  TFile* f1 = new TFile(fileName1);
  TH1F* hData1 = (TH1F*) f1->Get("plotDataFit");
  TH1F* hSim1  = (TH1F*) f1->Get("plotSimFit");
  //TH1F* hSim1   = (TH1F*) f1->Get("plotReso");

  TFile* f2 = new TFile(fileName2);
  TH1F* hData2 = (TH1F*) f2->Get("plotDataFit");
  TH1F* hSim2  = (TH1F*) f2->Get("plotSimFit");
  //TH1F* hSim2   = (TH1F*) f2->Get("plotReso");

  TFile* f3;
  TH1F* hData3;
  TH1F* hSim3;

  if(!skip3rd){
    f3 = new TFile(fileName3);
    hData3 = (TH1F*) f3->Get("plotDataFit");
    hSim3  = (TH1F*) f3->Get("plotSimFit");
    //hSim3   = (TH1F*) f3->Get("plotReso");
  }  

  //------
  TCanvas* canvas1 = new TCanvas("canvas1","canvas1",0,0,500,500);
  canvas1->cd();
  if(type==5) canvas1->SetLogy(1);

  hData1->SetMaximum(hRange);
  hData1->SetMinimum(lRange);

  hSim1->SetMaximum(hRange);
  hSim1->SetMinimum(lRange);

  //MarkersSyle:
  //20,21,22: bullet,square,triangle
  //24,25,26: as above, but empty symbols
  
  //Colors:
  //1,2,4: black, red, blue


  hData1->SetMarkerStyle(20); hData1->SetMarkerColor(1); hData1->SetMarkerSize(1.5);
  hSim1->SetMarkerStyle(24);  hSim1->SetMarkerColor(1);  hSim1->SetMarkerSize(1.5);

  hData2->SetMarkerStyle(21); hData2->SetMarkerColor(2); hData2->SetMarkerSize(1.5);
  hSim2->SetMarkerStyle(25);  hSim2->SetMarkerColor(2);  hSim2->SetMarkerSize(1.5);

  if(!skip3rd){
    hData3->SetMarkerStyle(22); hData3->SetMarkerColor(4); hData3->SetMarkerSize(1.5);
    hSim3->SetMarkerStyle(26);  hSim3->SetMarkerColor(4);  hSim3->SetMarkerSize(1.5);
  }

  //setStyle();
  hSim1->SetYTitle(yLabel);
  hSim1->SetXTitle(xLabel);

  hSim1->Draw("E19");
  hData1->SetLineStyle(2);
  //hData1->Draw("sameHE1");
  hData1->Draw("sameE1");

  hSim2->Draw("sameE19");
  //hData2->Draw("sameHE19");
  hData2->Draw("sameE19");

  if(!skip3rd){
    hSim3->Draw("sameE19");
    hData3->Draw("sameE19");
  }

  if(skip3rd){
    TLegend* leg   = new TLegend(0.17,0.75,0.25,0.90);
    leg->AddEntry(hData1,"  Data","p");
    leg->AddEntry(hData2,"  Data","p");
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->SetShadowColor(0);
    leg->SetTextSize(0.035);
    leg->SetTextAlign(12);
    leg->Draw();

    TLegend* leg2   = new TLegend(0.33,0.75,0.35,0.90);
    leg2->AddEntry(hSim1,"","p");
    leg2->AddEntry(hSim2,"","p");
    leg2->SetFillColor(0);
    leg2->SetLineColor(0);
    leg2->SetShadowColor(0);
    leg2->SetTextSize(0.035);
    leg2->SetTextAlign(12);



    TLegend* leg3   = new TLegend(0.24,0.74,0.92,0.89);
    leg3->AddEntry(hSim1,"Simulation  p_{T} in (1.0 #pm 0.1) GeV/c, |#eta| < 0.4","");
    leg3->AddEntry(hSim2,"Simulation  p_{T} in (3.0 #pm 0.2) GeV/c, |#eta| < 0.4","");
    //leg3->AddEntry(hSim1,"MC-Truth  p_{T} in (1.0 #pm 0.1) GeV/c, |#eta| < 0.4","");
    //leg3->AddEntry(hSim2,"MC-Truth  p_{T} in (3.0 #pm 0.2) GeV/c, |#eta| < 0.4","");
    leg3->SetFillColor(0);
    leg3->SetLineColor(0);
    leg3->SetShadowColor(0);
    leg3->SetTextSize(0.035);
    leg3->SetTextAlign(22);
    leg3->Draw();
    leg2->Draw();
    leg->Draw();

  }else{
    TLegend* leg   = new TLegend(0.25,0.75,0.35,0.90);
    leg->AddEntry(hData1,"   Data","p");
    leg->AddEntry(hData2,"   Data","p");
    leg->AddEntry(hData3,"   Data","p");
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->SetShadowColor(0);
    leg->SetTextSize(0.035);    leg->SetTextAlign(12);
    leg->Draw();

    TLegend* leg2   = new TLegend(0.35,0.75,0.82,0.90);
    leg2->AddEntry(hSim1,"Simulation  p_{T} in (1.0 #pm 0.1) GeV/c","p");
    leg2->AddEntry(hSim2,"Simulation  p_{T} in (3.0 #pm 0.2) GeV/c","p");
    leg2->AddEntry(hSim3,"Simulation  p_{T} in (8.0 #pm 1.0) GeV/c","p");
    //leg2->AddEntry(hSim1,"MC-Truth  p_{T} in (1.0 #pm 0.1) GeV/c","p");
    //leg2->AddEntry(hSim2,"MC-Truth  p_{T} in (3.0 #pm 0.2) GeV/c","p");
    //leg2->AddEntry(hSim3,"MC-Truth  p_{T} in (8.0 #pm 1.0) GeV/c","p");
    leg2->SetFillColor(0);
    leg2->SetLineColor(0);
    leg2->SetShadowColor(0);
    leg2->SetTextSize(0.035);
    leg2->SetTextAlign(12);
    leg2->Draw();
  }
  cmsPrel(0);

  gPad->Update();
  gPad->Print(outputName1);
  gPad->Print(outputName2);
  gPad->Print(outputName3);
}

void cmsPrel(double intLumi){
   TLatex latex;
   latex.SetNDC();
   latex.SetTextSize(0.04);

   latex.SetTextAlign(31); // align right
   //latex.DrawLatex(0.98,0.965,"#font[22]{#sqrt{s} = 7 TeV}");
   latex.DrawLatex(0.98,0.965,"#sqrt{s} = 7 TeV");
   if (intLumi > 0.) {
     latex.SetTextAlign(31); // align right
     latex.DrawLatex(0.98,0.88,Form("#int #font[12]{L}dt = %.1f  nb^{-1}",intLumi));
   }
   latex.SetTextAlign(11); // align left
   //latex.DrawLatex(0.02,0.965,"#font[22]{CMS preliminary 2010}");
   latex.DrawLatex(0.02,0.965,"CMS preliminary 2010");
   return;
}


void setLabels(int type,
	       TString& xLabel,TString& yLabel,TString& yLabelResp,TString& yLabelErr,
	       TString& text1a,TString& text1b,TString& leg1a,TString& leg1b,
	       TString& text2a,TString& text2b,TString& leg2a,TString& leg2b,
	       double& lRange,double& hRange
	       ){

  cout << "in setLabels, type: " << type << endl;


  if(type==2){
    xLabel = "Track #eta";
    yLabel = "Transv.Impact Parameter Resolution (#mum)";
    //yLabel = "#font[22]{Transv.Impact Parameter Resolution (#mum)}";
    yLabelResp = "Smearing of the TIP due to Vertex (#mum)";
    yLabelErr  = "Transv.Impact Parameter Error (#mum)";    

    lRange = 0.;
    hRange = 260.;

    text1a = "CMS Preliminary 2010";
    text1b = "7 TeV Runs";
    leg1a  = "DATA";
    leg1b  = "MC";

    text2a = "CMS Preliminary 2010";
    text2b = "";
    leg2a  = "MC (using MC-truth)";
    leg2b  = "MC (using only reco)";
  }


  if(type==5){
    gStyle->SetOptLogy(1);
    xLabel = "Track #eta";
    yLabel = "Longit.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the LIP due to Vertex (#mum)";
    yLabelErr  = "Longit.Impact Parameter Error (#mum)";    

    lRange = 25.;
    hRange = 3000.;

    text1a = "CMS Preliminary 2010";
    text1b = "7 TeV Runs";
    leg1a  = "DATA";
    leg1b  = "MC";

    text2a = "CMS Preliminary 2010";
    text2b = "";
    leg2a  = "MC (using MC-truth)";
    leg2b  = "MC (using only reco)";
  }

}


void setStyle(){
  //=========  settings ====================
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  //gStyle->SetPadGridX(kTRUE);
  //gStyle->SetPadGridY(kTRUE);
  gStyle->SetOptStat("");
  gStyle->SetOptFit(0000000000);

  //
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.15);
  //gStyle->SetLabelOffset(1.4,"Y");
  //gStyle->SetTitleYOffset();
  
  //
  gStyle->SetTitle("");
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetTitleFontSize(0.04);

  //
  gStyle->SetStatFontSize(0.04);
  gStyle->SetTextSize(0.04);
  gStyle->SetFrameFillColor(0);

  gStyle->SetHistFillColor(0);
  gStyle->SetHistLineColor(1);
  gStyle->SetHistLineWidth(2);

  gStyle->SetFillColor(0);
  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(1);
  
  gStyle->SetMarkerSize(0.5);
  gStyle->SetErrorX(0);
  
  // For the axis: 
  gStyle->SetPadTickX(1);  //to get tick marks on the opposite side
  gStyle->SetPadTickY(1);
}
