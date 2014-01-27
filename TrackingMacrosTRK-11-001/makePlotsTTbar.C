#include "TPaveStats.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFitResultPtr.h" 
#include "TFitResult.h" 

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooAddPdf.h"
#include "rooDoubleCB.h"

using namespace RooFit;
#include <sstream>


void cmsLabel(double intLumi){
  // -1 ==> print only CMS Simulation
  //  0 ==> print CMS Simulation and sqrt(s)
  // >0 ==> print CMS, sqrt(s) and lumi

  //0.02 ; 0.98

   TLatex latex;
   latex.SetNDC();
   latex.SetTextSize(0.04);

   latex.SetTextAlign(31); // align right
   if(intLumi >=0) latex.DrawLatex(0.92,0.965,"#sqrt{s} = 7 TeV");
   if (intLumi > 0.) {
     latex.SetTextAlign(31); // align right
     latex.DrawLatex(0.60,0.965,Form("L = %.1f  fb^{-1}",intLumi));
   }
   latex.SetTextAlign(11); // align left
   if(intLumi <=0) latex.DrawLatex(0.14,0.965,"CMS simulation");
   else latex.DrawLatex(0.15,0.965,"CMS");
   return;
}
 

void printCanvas(TCanvas* canvas,string outputName,int type=0,int format=3){
  /*
  TPaveText *pt;
  if(cmsLabelLocation==1) pt= new TPaveText(.16,.15,.38,.25,"NDC");  //bottom-left
  if(cmsLabelLocation==2) pt = new TPaveText(.16,.84,.38,.94,"NDC"); //bottom-right
  if(cmsLabelLocation==3) pt = new TPaveText(.65,.15,.85,.25,"NDC"); //top-left
  if(cmsLabelLocation==4) pt = new TPaveText(.65,.84,.85,.94,"NDC"); //top-right


  pt->AddText("CMS Simulation");
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  pt->SetLineColor(0);
  pt->Draw();
  */

  cmsLabel(0);
  string outputFolder;
  if(type==0) outputFolder="output/gt_vs_hp/";
  else if(type==1) outputFolder="output/pu_vs_nopu/";

  if(format>=1) canvas->Print((outputFolder.c_str()+outputName+".png").c_str());
  if(format>=2) canvas->Print((outputFolder.c_str()+outputName+".pdf").c_str());
  if(format>=3) canvas->Print((outputFolder.c_str()+outputName+".eps").c_str());
}

void setLegend1(TLegend* l,
		TH1F* h1,
		TString label1){
  l->SetTextSize(0.04);
  l->SetLineColor(0);
  l->SetFillColor(0);
  l->AddEntry(h1,label1,"P");
  l->SetShadowColor(0);
}


void setLegend2(TLegend* l,
		TH1F* h1,TH1F* h2,
		TString label1, TString label2){
  l->SetTextSize(0.04);
  l->SetLineColor(0);
  l->SetFillColor(0);
  l->AddEntry(h1,label1,"P");
  l->AddEntry(h2,label2,"P");
  l->SetShadowColor(0);
}


void setLegend4(TLegend* l,
		TH1F* h1,TH1F* h2, TH1F* h3,TH1F* h4,
		TString label1, TString label2,TString label3, TString label4){
  l->SetTextSize(0.04);
  l->SetLineColor(0);
  l->SetFillColor(0);
  l->AddEntry(h1,label1,"P");
  l->AddEntry(h2,label2,"P");
  l->AddEntry(h3,label3,"P");
  l->AddEntry(h4,label4,"P");
  l->SetShadowColor(0);
}



void setHistoStyle(TH1F* newH,TH1F* refH,bool isRMS=false){
  if(isRMS){
     newH->SetMarkerStyle(20);
     refH->SetMarkerStyle(24);
  }else{
     newH->SetMarkerStyle(21);
     refH->SetMarkerStyle(25);
  }
     newH->SetMarkerColor(1);
     refH->SetMarkerColor(1);
     newH->SetMarkerSize(1.0);
     refH->SetMarkerSize(1.0);
     newH->SetLineColor(1);
     refH->SetLineColor(1);
     newH->SetLineWidth(1);
     refH->SetLineWidth(1);
}


int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
}



void makePlotsTTbar(int input=0)
{
  // +++++++++++++ BEGIN EDITME ++++++++++++++
  TString refLabel;
  TString newLabel;
  TString labelTTbar;
  TString collname2;
  TString collname1;
  TString refFile;
  TString newFile;
  TString label1Fit;
  TString label2Fit;
  TString label1RMS;
  TString label2RMS;

  if(input==0){
    // -------- for gt vs hp -----------  
    refLabel="All tracks";
    newLabel="Only high-purity tracks";
    labelTTbar = "TTbar with pileup";
    collname2="cutsRecoHp_AssociatorByHits";  
    collname1="general_AssociatorByHits";
    refFile="input/ttbar/pugt.root";
    newFile="input/ttbar/puhp.root"; 
    label1Fit = "";
    label2Fit = "";
    label1RMS = "";
    label2RMS = "";  
  }else if(input==1){
    // -------- for pu vs nopu ----------- 
    refLabel = "Without pileup";
    newLabel = "With pileup";
    labelTTbar = "TTbar with pileup";
    label1Fit = "without pileup (sigma)";
    label2Fit = "Guassian fit, sigma";
    label1RMS = "without pileup (RMS)";
    label2RMS = "RMS";
    collname2="cutsRecoHp_AssociatorByHits"; 
    collname1="cutsRecoHp_AssociatorByHits"; 
    refFile="input/ttbar/nopuhp.root";
    newFile="input/ttbar/puhp.root";
  }else{
    return;
  }

  // -------

  // +++++++++++++ END EDITME ++++++++++++++
  // +++++++++++++++++++++++++++++++++++++++


  gROOT ->Reset();
  //gROOT ->SetBatch();

 //=========  settings ====================
 //gROOT->SetStyle("Plain");
 gStyle->SetPadGridX(kTRUE);
 gStyle->SetPadGridY(kTRUE);
 gStyle->SetPadRightMargin(0.07);
 gStyle->SetPadLeftMargin(0.13);
  gStyle->SetOptStat("");
 //gStyle->SetTitleXSize(0.07); 
 //gStyle->SetTitleXOffset(0.6); 
 //tyle->SetTitleYSize(0.3);
 //gStyle->SetLabelSize(0.6) 
 //gStyle->SetTextSize(0.5);
 //char* refLabel("REF_LABEL, REF_RELEASE REFSELECTION");
 //char* newLabel("NEW_LABEL, NEW_RELEASE NEWSELECTION");



 //=============================================


 delete gROOT->GetListOfFiles()->FindObject(refFile);
 delete gROOT->GetListOfFiles()->FindObject(newFile); 

 TFile * sfile = new TFile(newFile);
 TDirectory * sdir=gDirectory;
 TFile * rfile = new TFile(refFile);
 TDirectory * rdir=gDirectory;

 if(sfile->GetDirectory("DQMData/Run 1/RecoTrackV")) sfile->cd("DQMData/Run 1/RecoTrackV/Run summary/Track");
 else if(sfile->GetDirectory("DQMData/RecoTrackV/Track"))sfile->cd("DQMData/RecoTrackV/Track");
 else if(sfile->GetDirectory("DQMData/Run 1/Tracking")) sfile->cd("DQMData/Run 1/Tracking/Run summary/Track");
 else if(sfile->GetDirectory("DQMData/Tracking/Track"))sfile->cd("DQMData/Tracking/Track");
 sdir=gDirectory;


 if(rfile->GetDirectory("DQMData/Run 1/RecoTrackV")) rfile->cd("DQMData/Run 1/RecoTrackV/Run summary/Track");
 else if(rfile->GetDirectory("DQMData/RecoTrackV/Track"))rfile->cd("DQMData/RecoTrackV/Track");
 else if(rfile->GetDirectory("DQMData/Run 1/Tracking")) rfile->cd("DQMData/Run 1/Tracking/Run summary/Track");
 else if(rfile->GetDirectory("DQMData/Tracking/Track"))rfile->cd("DQMData/Tracking/Track");
 rdir=gDirectory;

 TCanvas *canvas;

 TH1F *sh1,*rh1;
 TH1F *sh2,*rh2;
 TH1F *sh3,*rh3;
 TH1F *sh4,*rh4;
 TH1F *sh5,*rh5;
 TH1F *sh6,*rh6;





 string output;
 TLegend* l;


 //===== building
 rdir->GetObject(collname1+"/effic",rh1);
 sdir->GetObject(collname2+"/effic",sh1);
 rh1->GetYaxis()->SetRangeUser(0.5,1.0);
 sh1->GetYaxis()->SetRangeUser(0.5,1.0);
 rh1->GetXaxis()->SetTitle("#eta");
 rh1->GetXaxis()->SetLabelOffset(0.012);
 rh1->GetYaxis()->SetTitle("Efficiency");
 rh1->GetYaxis()->SetTitleSize(0.05);
 rh1->GetYaxis()->SetTitleOffset(1.2);
 rdir->GetObject(collname1+"/fakerate",rh2);
 sdir->GetObject(collname2+"/fakerate",sh2);
 rh2->GetYaxis()->SetRangeUser(0.,0.30);
 sh2->GetYaxis()->SetRangeUser(0.,0.30);
 rh2->GetXaxis()->SetTitle("#eta");
 rh2->GetXaxis()->SetLabelOffset(0.012);
 rh2->GetYaxis()->SetTitle("Fake rate");
 rh2->GetYaxis()->SetTitleSize(0.05);
 rh2->GetYaxis()->SetTitleOffset(1.2);
 

 rdir->GetObject(collname1+"/efficPt",rh3);
 sdir->GetObject(collname2+"/efficPt",sh3);
 rh3->GetYaxis()->SetRangeUser(0.,1.0);
 sh3->GetYaxis()->SetRangeUser(0.,1.0);
 rh3->GetXaxis()->SetRangeUser(0,150);
 sh3->GetXaxis()->SetRangeUser(0,150);
 rh3->GetXaxis()->SetTitle("p_{T} (GeV)");
 rh3->GetYaxis()->SetTitle("Efficiency");
 rh3->GetYaxis()->SetTitleSize(0.05);
 rh3->GetYaxis()->SetTitleOffset(1.3);
 rh3->GetXaxis()->SetTitleOffset(1.3);
 rh3->SetTitle("");
 rdir->GetObject(collname1+"/fakeratePt",rh4);
 sdir->GetObject(collname2+"/fakeratePt",sh4);
 rh4->SetTitle("");
 rh4->GetXaxis()->SetTitle("p_{T} (GeV)");
 rh4->GetYaxis()->SetTitle("Fake rate");
 rh4->GetYaxis()->SetTitleSize(0.05);
 rh4->GetYaxis()->SetTitleOffset(1.3);
 rh4->GetXaxis()->SetTitleOffset(1.3);
 rh4->GetYaxis()->SetRangeUser(0.,.30);
 sh4->GetYaxis()->SetRangeUser(0.,.30);
 rh4->GetXaxis()->SetRangeUser(0.2,150);
 sh4->GetXaxis()->SetRangeUser(0.2,150);
 
 
 rdir->GetObject(collname1+"/effic_vs_hit",rh5);
 sdir->GetObject(collname2+"/effic_vs_hit",sh5);
 rh5->GetXaxis()->SetTitle("hits");
 rh5->GetYaxis()->SetTitle("efficiency vs hits");
 rh5->GetYaxis()->SetTitleSize(0.05);
 rh5->GetYaxis()->SetTitleOffset(1.2);
 //rh3->GetXaxis()->SetRangeUser(0,30);
 //sh3->GetXaxis()->SetRangeUser(0,30);
 rdir->GetObject(collname1+"/fakerate_vs_hit",rh6);
 sdir->GetObject(collname2+"/fakerate_vs_hit",sh6);
 rh6->GetYaxis()->SetRangeUser(0.,1.0);
 rh6->GetYaxis()->SetRangeUser(0.,1.0);
 rh6->GetXaxis()->SetTitle("hits");
 rh6->GetYaxis()->SetTitle("fakerate vs hits");
 rh6->GetYaxis()->SetTitleSize(0.05);
 rh6->GetYaxis()->SetTitleOffset(1.2);
 
 TH1F * r[6]={rh1,rh2,rh3,rh4,rh5,rh6};
 TH1F * s[6]={sh1,sh2,sh3,sh4,sh5,sh6};
 
 for(int i=0; i<6;i++){
   setHistoStyle(s[i],r[i]);
 }
 
 canvas = new TCanvas("Tracks","Tracks: efficiency & fakerate",500,500);
 canvas->cd();
 
 // ----
 r[0]->Draw();
 s[0]->Draw("sames");
 
   

 l = new TLegend(0.40,0.20,0.80,0.35);
 if(input==0)  l = new TLegend(0.30,0.20,0.85,0.35);
 setLegend2(l,rh1,sh1,refLabel,newLabel);
 l->Draw();
 output="efficiencyVsEta"; printCanvas(canvas,output,input);
 delete l;
 
 
 // ----
 r[1]->Draw();
 s[1]->Draw("sames");
 
 l = new TLegend(0.37,0.77,0.77,0.92);
 if(input==0)  l = new TLegend(0.37,0.77,0.87,0.90);
 setLegend2(l,rh1,sh1,refLabel,newLabel);
 l->Draw();
 output="fakerateVsEta"; printCanvas(canvas,output,input);
 delete l;
 

 // ----
 canvas->SetLogx();
 r[2]->Draw();
 s[2]->Draw("sames");
 
 l = new TLegend(0.40,0.60,0.80,0.75);
 if(input==0)  l = new TLegend(0.35,0.60,0.90,0.75);
 setLegend2(l,rh1,sh1,refLabel,newLabel);
 l->Draw();
 output="efficiencyVsPt"; printCanvas(canvas,output,input);
 delete l;
 

 // ----
 canvas->SetLogx(1);
 r[3]->Draw();
 s[3]->Draw("sames");
 
 l = new TLegend(0.30,0.65,0.70,0.80); //all vs highpurity
 if(input==0)  l = new TLegend(0.27,0.65,0.80,0.80);
 setLegend2(l,rh1,sh1,refLabel,newLabel);
 l->Draw();
 output="fakerateVsPt"; printCanvas(canvas,output,input);
 delete l;
   

}


