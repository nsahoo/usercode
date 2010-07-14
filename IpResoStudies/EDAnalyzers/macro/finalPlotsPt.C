
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaveText.h"

#include <iostream>
#include <iomanip>
#include <sstream>

//void setStyle();

void cmsPrel(double intLumi);

void setLabels(int type,
	       TString& xLabel,TString& yLabel,TString& yLabelRes,TString& yLabelErr,
	       TString& text1a,TString& text1b,TString& leg1a,TString& leg1b,
	       TString& text2a,TString& text2b,TString& leg2a,TString& leg2b,
	       double& lRange,double& hRange
	       );

void finalPlotsPt(int type,bool skip3rd=false){
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


  //TString fileName1="./finalPlotsR27th/histos"+counter+".root";
  TString fileName1="./finalPlotsR27thSol6/histos"+counter+".root";
  //TString fileName1="./finalPlotsR27thSol7/histos"+counter+".root";

  TString outputName1="final_type"+counter+".pdf";
  TString outputName2="final_type"+counter+".png";
  TString outputName3="final_type"+counter+".eps";

  TFile* f1 = new TFile(fileName1);
  TH1F* hData1 = (TH1F*) f1->Get("plotDataFit");
  TH1F* hSim1  = (TH1F*) f1->Get("plotSimFit");
  //TH1F* hMC1   = (TH1F*) f1->Get("plotReso");


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
  //1,2: black, red


  hData1->SetMarkerStyle(20); hData1->SetMarkerColor(1); hData1->SetMarkerSize(1.5);
  hSim1->SetMarkerStyle(20);  hSim1->SetMarkerColor(2);  hSim1->SetMarkerSize(1.5);

  //setStyle();
  hSim1->SetYTitle(yLabel);
  hSim1->SetXTitle(xLabel);

  gStyle->SetErrorX(0.5);

  /*
  //TF1* finalFit = new TF1("finalFit","sqrt( ([0]/x)*([0]/x)+[1]*[1])",1.0,10);
  TF1* finalFit = new TF1("finalFit","[0]/x + [1]");
  finalFit->SetParameter(0,6);
  finalFit->SetParameter(1,20);
  finalFit->SetParName(0,"MultipleScattering term");
  finalFit->SetParName(1,"SinglePointResolution term");
  finalFit->SetLineColor(4);
  finalFit->SetLineWidth(2.5);
  finalFit->SetLineStyle(1);


  //hSim1->Fit("finalFit","","",1.0,10);
  //hData1->Fit("finalFit","","",1.0,10);


  double sA,sB,dA,dB;
  double sAErr,sBErr,dAErr,dBErr;

  hSim1->Fit("finalFit","N");
  sA = finalFit->GetParameter(0);
  sB = finalFit->GetParameter(1);

  hData1->Fit("finalFit");
  dA = finalFit->GetParameter(0);
  dB = finalFit->GetParameter(1);
  dAErr = finalFit->GetParError(0);
  dBErr = finalFit->GetParError(1);


  cout << "sA,sB: " << setiosflags(ios::fixed) << setprecision(1) << sA << " , " << sB << endl;
  cout << "dA,dB: " << setiosflags(ios::fixed) << setprecision(1) << dA << " , " << dB << endl;
  cout << "chi2/NDF: " << finalFit->GetChisquare()/finalFit->GetNDF() << endl;

  TString leb,den;
  if(type==1 || type==4) {
    if(type==1)
      leb="d_{0}";
    else
      leb="d_{0}";
    den="p_{T}";
  }
  if(type==7 || type==8) {
    if(type==7)
      leb="d_{0}";
    else
      leb="d_{Z}";
    den="p";  
  }

  stringstream stream2;  
  stream2 << "For Data, #sigma_{" << leb << "} = #left(#frac{" << setiosflags(ios::fixed) << setprecision(1) << dA <<" #pm " 
	  << dAErr << "}{"
	  << den << "}"<< " + (" << dB << " #pm " << dBErr << ") #right) #mum";
  cout << stream2.str() << endl;

  stringstream stream3;  
  stream3 << "Sim., #sigma(" << leb << ") = #frac{(" << setiosflags(ios::fixed) << setprecision(1) << sA << den << sB << ") #mum";
  cout << stream3.str() << endl;
  */

  TLegend* leg   = new TLegend(0.25,0.55,0.85,0.70);
  leg->AddEntry(hData1,"Data         track | #eta | < 0.4","p");
  leg->AddEntry(hSim1, "Simulation   track | #eta | < 0.4","p");

  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.04);
  leg->SetTextAlign(32);
  
  hSim1->Draw("E1");
  hData1->Draw("sameE1");
  leg->Draw();

  /*
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.03);
  latex.DrawLatex(0.35,0.60,stream2.str().c_str());
  */

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
   //latex.DrawLatex(0.20,0.20,"-0.4 < #eta < 0.4");
   return;
}


void setLabels(int type,
	       TString& xLabel,TString& yLabel,TString& yLabelResp,TString& yLabelErr,
	       TString& text1a,TString& text1b,TString& leg1a,TString& leg1b,
	       TString& text2a,TString& text2b,TString& leg2a,TString& leg2b,
	       double& lRange,double& hRange
	       ){

  cout << "in setLabels, type: " << type << endl;

  text1a = "CMS Preliminary 2010";
  text1b = "7 TeV Runs";
  leg1a  = "DATA";
  leg1b  = "MC";

  text2a = "CMS Preliminary 2010";
  text2b = "";
  leg2a  = "MC (using MC-truth)";
  leg2b  = "MC (using only reco)";
  
  if(type==1){
    xLabel = "Track p_{T} (GeV/c)";
    //yLabel = "Transv.Impact Parameter Resolution   #sigma_{d_{0}}  (#mum)";
    yLabel = "Transv.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the TIP due to Vertex (#mum)";
    yLabelErr  = "Transv.Impact Parameter Error (#mum)";    

    lRange = 10;
    hRange = 110.;
  }

  if(type==2){
    xLabel = "Track #eta";
    yLabel = "Transv.Impact Parameter Resolution (#mum)";
    //yLabel = "#font[22]{Transv.Impact Parameter Resolution (#mum)}";
    yLabelResp = "Smearing of the TIP due to Vertex (#mum)";
    yLabelErr  = "Transv.Impact Parameter Error (#mum)";    

    lRange = 0.;
    hRange = 260.;
  }

  if(type==3){
    xLabel = "Track #phi";
    yLabel = "Transv.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the TIP due to Vertex (#mum)";
    yLabelErr  = "Transv.Impact Parameter Error (#mum)";    

    lRange = 0.;
    hRange = 200.;
  }

  if(type==4){
    xLabel = "Track p_{T} (GeV/c)";
    yLabel = "Longit.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the LIP due to Vertex (#mum)";
    yLabelErr  = "Longit.Impact Parameter Error (#mum)";    

    lRange = 30.;
    hRange = 130.;
  }


  if(type==5){
    gStyle->SetOptLogy(1);
    xLabel = "Track #eta";
    yLabel = "Longit.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the LIP due to Vertex (#mum)";
    yLabelErr  = "Longit.Impact Parameter Error (#mum)";    

    lRange = 15.;
    hRange = 1500.;
  }


  if(type==6){
    gStyle->SetOptLogy(0);
    xLabel = "Track #phi";
    yLabel = "Longit.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the LIP due to Vertex (#mum)";
    yLabelErr  = "Longit.Impact Parameter Error (#mum)";    

    lRange = 0.;
    hRange = 210.;
  }

  if(type==7){
    xLabel = "Track momentum (GeV/c)";
    yLabel = "Transv.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the TIP due to Vertex (#mum)";
    yLabelErr  = "Transv.Impact Parameter Error (#mum)";    
    
    lRange = 0.;
    hRange = 140.;
  }
  
  if(type==8){
    xLabel = "Track momentum (GeV/c)";
    yLabel = "Longit.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the LIP due to Vertex (#mum)";
    yLabelErr  = "Longit.Impact Parameter Error (#mum)";    

    lRange = 0.;
    hRange = 140.;
  }


}

/*
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
*/
