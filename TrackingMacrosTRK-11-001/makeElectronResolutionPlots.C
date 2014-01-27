#include "TPaveStats.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"
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
//#include "MinimizerOptions.h"
#include <Math/MinimizerOptions.h> 

using namespace RooFit;
using namespace std;


#include <sstream>

TH1D * TH1Subset(TH1D *h,double lowEdge,double highEdge);



void extractReso(TH2F* th2,
		 TH1F*& hFit,TH1F*& hRMS,
		 TH1F*& hFitMean, TH1F*& hAverage,
		 int tailSign=-1,string name="", int type=0,bool variableBinning=false, bool wait=false);

void setHistoLabels(TH1* h,float ySize, float yOffset,float xSize, float xOffset,
		    TString yLabel,TString xLabel){
  h->SetTitle("");
  h->GetYaxis()->SetTitleSize(ySize);
  h->GetYaxis()->SetTitleOffset(yOffset);
  h->GetXaxis()->SetTitleSize(xSize);
  h->GetXaxis()->SetTitleOffset(xOffset);
  h->GetYaxis()->SetTitle(yLabel);
  h->GetXaxis()->SetTitle(xLabel);
}

void setLegend3(TLegend* l,
		TH1* h1,TH1* h2,TH1* h3,
		TString label1, TString label2, TString label3){
  l->SetTextSize(0.04);
  l->SetLineColor(0);
  l->SetFillColor(0);
  l->AddEntry(h1,label1,"P");
  l->AddEntry(h2,label2,"P");
  l->AddEntry(h3,label3,"P");
  l->SetShadowColor(0);
}

void setLegend2(TLegend* l,
		TH1* h2,TH1* h3,
		TString label2, TString label3){
  l->SetTextSize(0.04);
  l->SetLineColor(0);
  l->SetFillColor(0);
  l->AddEntry(h2,label2,"P");
  l->AddEntry(h3,label3,"P");
  l->SetShadowColor(0);
}


void setHistoStyle(TH1* h1,TH1* h2,TH1* h3,bool empty=false){
  if(empty){
     h1->SetMarkerStyle(24);
     h2->SetMarkerStyle(26);
     h3->SetMarkerStyle(25);
  }else{
     h1->SetMarkerStyle(20);
     h2->SetMarkerStyle(22);
     h3->SetMarkerStyle(21);
  }
     h1->SetMarkerColor(1);
     h2->SetMarkerColor(4);
     h3->SetMarkerColor(2);
     h1->SetMarkerSize(1.0);
     h2->SetMarkerSize(1.0);
     h3->SetMarkerSize(1.0);
     h1->SetLineColor(1);
     h2->SetLineColor(1);
     h3->SetLineColor(1);
     h1->SetLineWidth(1);
     h2->SetLineWidth(1);
     h3->SetLineWidth(1);
}


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

  cmsLabel(-1);

  string outputFolder;
  if(type==0) outputFolder="output/mu/";
  else if(type==1) outputFolder="output/pi/";
  else if(type==2) outputFolder="output/el/";

  if(format>=1) canvas->Print((outputFolder+outputName.c_str()+".png").c_str());
  if(format>=2) canvas->Print((outputFolder+outputName.c_str()+".pdf").c_str());
  if(format>=3) canvas->Print((outputFolder+outputName.c_str()+".eps").c_str());
}


int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
}


void makeElectronResolutionPlots(int input=0)
{ 
  /*
    input:
    0 --> muons
    1 --> pions
    2 --> electrons
   */

  ROOT::Math::MinimizerOptions minimizerOptions; 
  minimizerOptions.SetDefaultPrintLevel(-1);


  bool eff=0;
  bool hitEff=0;
  bool resVsEta=1;
  bool resVsPt=0;

  TString fileName1;
  TString fileName2;
  TString fileName3;
  TString fileName4;
  TString fileName5;
  TString fileName6;
  TString fileName7;
  TString fileName8;
  TString fileName9;


  TString label1;
  TString label2;
  TString label3;
  TString label4;
  TString label5;
  TString label6;

  //in microns
  double resDPhiVsEtaMax,resDPhiVsEtaMin;
  double resD0vsEtaMax,resD0vsEtaMin;
  double resDzVsEtaMax,resDzVsEtaMin;
  double resDThetaVsEtaMax,resDThetaVsEtaMin;
  double resPtVsEtaMax,resPtVsEtaMin;

  //in microns and 10-3 units
  double resDPhiVsPtMax,resDPhiVsPtMin;
  double resD0vsPtMax,resD0vsPtMin;
  double resDzVsPtMax,resDzVsPtMin;
  double resDThetaVsPtMax,resDThetaVsPtMin;
  double resPtVsPtMax,resPtVsPtMin; 
  
  int tailSignPhi;
  int tailSignTheta;
  int tailSignD0;
  int tailSignDz;
  int tailSignPt;

  int tailSignPhi_pT;
  int tailSignTheta_pT;
  int tailSignD0_pT;
  int tailSignDz_pT;
  int tailSignPt_pT;

  TString collname1;
  TString collname2;
  TString collname3;
  
  TString collname4;
  TString collname5;
  TString collname6;

  bool skipPt1(false);

  if(input==0){
    return;
  }else if(input==1){
    return;
  }else if(input==2){
    // ------ for electrons ---------    
    skipPt1 = true;

    /*
    fileName1="input/elMode/elPt10Mode.root";
    fileName2="input/elMode/elPt10Mode.root";
    fileName3="input/elMode/elPt100Mode.root";

    //for efficiency vs pt
    fileName4="input/elMode/elFlatPtBarrelForEff.root";
    fileName5="input/elMode/elFlatPtTransitionForEff.root";
    fileName6="input/elMode/elFlatPtEndcapForEff.root";
    //for fakes vs pt
    fileName7="input/elMode/elFlatPtBarrelForFake.root";
    fileName8="input/elMode/elFlatPtTransitionForFake.root";
    fileName9="input/elMode/elFlatPtEndcapForFake.root";
    */


    
    //for dz,cotTheta plots
    //fileName1="input/elKF/elPt1Tk.root";
    //fileName2="input/elKF/elPt10Tk.root";
    //fileName3="input/elKF/elPt100Tk.root";


    //for d0,phi,pt plots
    fileName1="input/elKF/elPt1Tk_largerRange.root";
    fileName2="input/elKF/elPt10Tk_largerRange.root";
    fileName3="input/elKF/elPt100Tk_largerRange.root";


    //for efficiency vs pt
    fileName4="input/elKF/elFlatPtBarrelForEff.root";
    fileName5="input/elKF/elFlatPtTransitionForEff.root";
    fileName6="input/elKF/elFlatPtEndcapForEff.root";
    //for fakes vs pt
    fileName7="input/elKF/elFlatPtBarrelForFake.root";
    fileName8="input/elKF/elFlatPtTransitionForFake.root";
    fileName9="input/elKF/elFlatPtEndcapForFake.root";



    collname1="DQMData/Tracking/Track/cutsRecoHp_AssociatorByHits";  
    collname2="DQMData/Tracking/Track/cutsRecoHp_AssociatorByHits";  
    collname3="DQMData/Tracking/Track/cutsRecoHp_AssociatorByHits";  
    collname4="DQMData/Tracking/Track/cutsRecoHp_AssociatorByHits";  
    collname5="DQMData/Tracking/Track/cutsRecoHp_AssociatorByHits";  
    collname6="DQMData/Tracking/Track/cutsRecoHp_AssociatorByHits";  
    
 

    label1="el, p_{T} = 1 GeV";
    label2="el, p_{T} = 10 GeV";
    label3="el, p_{T} = 100 GeV";
    label4="el, barrel region";
    label5="el, transition region"; 
    label6="el, endcap region"; 

    resDPhiVsEtaMax=5000;resDPhiVsEtaMin=0.05;
    resD0vsEtaMax=100000;resD0vsEtaMin=0.1;
    resDzVsEtaMax=150000;resDzVsEtaMin=2;
    resDThetaVsEtaMax=500;resDThetaVsEtaMin=0.1;
    resPtVsEtaMax=600;resPtVsEtaMin=2; //for el
    //resPtVsEtaMax=400;resPtVsEtaMin=0.5; //for el
  
    //in microns and 10-3 units
    resDPhiVsPtMax=80;resDPhiVsPtMin=0.05;
    resD0vsPtMax=1500;resD0vsPtMin=5;
    resDzVsPtMax=1500;resDzVsPtMin=20;
    resDThetaVsPtMax=100;resDThetaVsPtMin=0.2;
    resPtVsPtMax=50;resPtVsPtMin=0.5; //for mu and pi

    tailSignPhi =-1;
    tailSignTheta =-1;
    tailSignD0 =+1;
    tailSignDz =-1;
    tailSignPt =-1;
 
    tailSignPhi_pT =-1;
    tailSignTheta_pT =-1;
    tailSignD0_pT =+1;
    tailSignDz_pT =-1;
    tailSignPt_pT =-1;

  }else{
    return;
  }



  //gROOT->ProcessLine(".x HistoCompare_Tracks.C");
  //gROOT ->Reset();
 //gROOT ->SetBatch();

 //=========  settings ====================
 //gROOT->SetStyle("Plain");
 gStyle->SetPadGridX(kTRUE);
 gStyle->SetPadGridY(kTRUE);
 gStyle->SetPadRightMargin(0.07);
 gStyle->SetPadLeftMargin(0.13);
 //gStyle->SetOptStat("eour");
 gStyle->SetOptStat("");
 //gStyle->SetTitleXSize(0.07); 
 //gStyle->SetTitleXOffset(0.6); 
 //tyle->SetTitleYSize(0.3);
 //gStyle->SetLabelSize(0.6) 
 //gStyle->SetTextSize(0.5);



 //=============================================


 TCanvas *canvas;
 

 //==============

 string output;
 if (eff){
   TFile * file1 = new TFile(fileName1); TDirectory * dir1=file1->GetDirectory("");
   TFile * file2 = new TFile(fileName2); TDirectory * dir2=file2->GetDirectory("");
   TFile * file3 = new TFile(fileName3); TDirectory * dir3=file3->GetDirectory("");
   TFile * file4 = new TFile(fileName4); TDirectory * dir4=file4->GetDirectory("");
   TFile * file5 = new TFile(fileName5); TDirectory * dir5=file5->GetDirectory("");
   TFile * file6 = new TFile(fileName6); TDirectory * dir6=file6->GetDirectory("");
   TFile * file7 = new TFile(fileName7); TDirectory * dir7=file7->GetDirectory("");
   TFile * file8 = new TFile(fileName8); TDirectory * dir8=file8->GetDirectory("");
   TFile * file9 = new TFile(fileName9); TDirectory * dir9=file9->GetDirectory("");

   TH1F *f1h1,*f2h1,*f3h1;
   TH1F *f1h2,*f2h2,*f3h2;
   TH1F *f1h3,*f2h3,*f3h3;
   TH1F *f1h4,*f2h4,*f3h4;
   TH1F *f1h5,*f2h5,*f3h5;
   TH1F *f1h6,*f2h6,*f3h6;

   dir1->GetObject(collname1+"/effic",f1h1);
   dir2->GetObject(collname2+"/effic",f2h1);
   dir3->GetObject(collname3+"/effic",f3h1);
   f1h1->GetYaxis()->SetRangeUser(0.5,1.05);
   f2h1->GetYaxis()->SetRangeUser(0.5,1.05);
   //f1h1->GetYaxis()->SetRangeUser(0.,1.05);
   //f2h1->GetYaxis()->SetRangeUser(0.,1.05);
   setHistoLabels(f1h1,0.05,1.2,0.07,0.6,"Efficiency","#eta");
   setHistoLabels(f2h1,0.05,1.2,0.07,0.6,"Efficiency","#eta");


   // ---
   TH1F  *xFakeNum1,*xFakeNum2,*xFakeNum3;
   TH1F  *xFakeDen1,*xFakeDen2,*xFakeDen3;
   dir1->GetObject(collname1+"/num_assoc(recoToSim)_eta",xFakeNum1);
   dir2->GetObject(collname2+"/num_assoc(recoToSim)_eta",xFakeNum2);
   dir3->GetObject(collname3+"/num_assoc(recoToSim)_eta",xFakeNum3);
   dir1->GetObject(collname1+"/num_reco_eta",xFakeDen1);
   dir2->GetObject(collname2+"/num_reco_eta",xFakeDen2);
   dir3->GetObject(collname3+"/num_reco_eta",xFakeDen3);

   xFakeNum1->Add(xFakeDen1,xFakeNum1,1.,-1.);
   xFakeNum2->Add(xFakeDen2,xFakeNum2,1.,-1.);
   xFakeNum3->Add(xFakeDen3,xFakeNum3,1.,-1.);

   dir1->GetObject(collname1+"/fakerate",f1h2);
   dir2->GetObject(collname2+"/fakerate",f2h2);
   dir3->GetObject(collname3+"/fakerate",f3h2);
   
   // replace fake rate histograms
   f1h2->Divide(xFakeNum1,xFakeDen1,1,1,"B");
   f2h2->Divide(xFakeNum2,xFakeDen2,1,1,"B");
   f3h2->Divide(xFakeNum3,xFakeDen3,1,1,"B");


   //

   f1h2->GetYaxis()->SetRangeUser(0.,0.30);
   f2h2->GetYaxis()->SetRangeUser(0.,0.30);
   setHistoLabels(f1h2,0.05,1.2,0.07,0.6,"Fake Rate","#eta");
   setHistoLabels(f2h2,0.05,1.2,0.07,0.6,"Fake Rate","#eta");

   // --
   dir4->GetObject(collname4+"/efficPt",f1h3);
   dir5->GetObject(collname5+"/efficPt",f2h3);
   dir6->GetObject(collname6+"/efficPt",f3h3);
   if(!skipPt1){
     f1h3->GetXaxis()->SetRangeUser(0,150);
     f2h3->GetXaxis()->SetRangeUser(0,150);
     f3h3->GetXaxis()->SetRangeUser(0,150);
   }else{
     f1h3->GetXaxis()->SetRangeUser(1.,150);
     f2h3->GetXaxis()->SetRangeUser(1.,150);
     f3h3->GetXaxis()->SetRangeUser(1.,150);
   }
   f1h3->GetYaxis()->SetRangeUser(0.,1.05);
   setHistoLabels(f1h3,0.05,1.3,0.055,1.0,"Efficiency","p_{T} (GeV)");

   //f2h3->GetXaxis()->SetRangeUser(0,150);
   f2h3->GetYaxis()->SetRangeUser(0.,1.05);
   setHistoLabels(f2h3,0.05,1.3,0.055,1.0,"Efficiency","p_{T} (GeV)");

   // --- 
   TH1F  *xFakePtNum1,*xFakePtNum2,*xFakePtNum3;
   TH1F  *xFakePtDen1,*xFakePtDen2,*xFakePtDen3;
   dir7->GetObject(collname1+"/num_assoc(recoToSim)_eta",xFakePtNum1);
   dir8->GetObject(collname2+"/num_assoc(recoToSim)_eta",xFakePtNum2);
   dir9->GetObject(collname3+"/num_assoc(recoToSim)_eta",xFakePtNum3);
   dir7->GetObject(collname1+"/num_reco_eta",xFakePtDen1);
   dir8->GetObject(collname2+"/num_reco_eta",xFakePtDen2);
   dir9->GetObject(collname3+"/num_reco_eta",xFakePtDen3);

   xFakePtNum1->Add(xFakePtDen1,xFakePtNum1,1.,-1.);
   xFakePtNum2->Add(xFakePtDen2,xFakePtNum2,1.,-1.);
   xFakePtNum3->Add(xFakePtDen3,xFakePtNum3,1.,-1.);

   dir7->GetObject(collname1+"/fakeratePt",f1h4);
   dir8->GetObject(collname2+"/fakeratePt",f2h4);
   dir9->GetObject(collname3+"/fakeratePt",f3h4);
   /*
   f1h4->Divide(xFakePtNum1,xFakePtDen1,1,1,"B");
   f2h4->Divide(xFakePtNum2,xFakePtDen2,1,1,"B");
   f3h4->Divide(xFakePtNum3,xFakePtDen3,1,1,"B");
   */

   f1h4->GetYaxis()->SetRangeUser(0.,0.30);
   if(!skipPt1){
     f1h4->GetXaxis()->SetRangeUser(0.2,150);
     f2h4->GetXaxis()->SetRangeUser(0.2,150);
     f3h4->GetXaxis()->SetRangeUser(0.2,150);
   }else{
     f1h4->GetYaxis()->SetRangeUser(0.,0.50);
     f1h4->GetXaxis()->SetRangeUser(1,150);
     f2h4->GetXaxis()->SetRangeUser(1,150);
     f3h4->GetXaxis()->SetRangeUser(1,150);
   }


   setHistoLabels(f1h4,0.05,1.3,0.055,1.0,"Fake Rate","p_{T} (GeV)");

   //f2h4->GetYaxis()->SetRangeUser(0.,0.30);
   //f2h4->GetXaxis()->SetRangeUser(0.2,150);
   setHistoLabels(f2h4,0.05,1.3,0.055,1.0,"Fake Rate","p_{T} (GeV)");


   dir1->GetObject(collname1+"/effic_vs_hit",f1h5);
   dir2->GetObject(collname2+"/effic_vs_hit",f2h5);
   dir3->GetObject(collname3+"/effic_vs_hit",f3h5);
   setHistoLabels(f1h5,0.05,1.2,0.07,1.2,"Efficiency","hits");
   
   dir1->GetObject(collname1+"/fakerate_vs_hit",f1h6);
   dir2->GetObject(collname2+"/fakerate_vs_hit",f2h6);
   dir3->GetObject(collname3+"/fakerate_vs_hit",f3h6);

   f1h6->GetYaxis()->SetRangeUser(0.,1.0);
   setHistoLabels(f1h6,0.05,1.2,0.07,1.2,"Fake Rate","hits");
   f2h6->GetYaxis()->SetRangeUser(0.,1.0);
   setHistoLabels(f2h6,0.05,1.2,0.07,1.2,"Fake Rate","hits");




   TH1 * f1histos[6]={f1h1,f1h2,f1h3,f1h4,f1h5,f1h6};
   TH1 * f2histos[6]={f2h1,f2h2,f2h3,f2h4,f2h5,f2h6};
   TH1 * f3histos[6]={f3h1,f3h2,f3h3,f3h4,f3h5,f3h6};

   for(int i=0; i<6;i++){
     setHistoStyle(f1histos[i],f2histos[i],f3histos[i]);
   }

   canvas = new TCanvas("Tracks","Tracks: efficiency & fakerate",500,500);
   canvas->cd();

   // ----
   TLegend* l = new TLegend(0.40,0.20,0.80,0.35);
   if(!skipPt1){
     f1histos[0]->Draw();
     f2histos[0]->Draw("same");
     f3histos[0]->Draw("same");
     setLegend3(l,f1h1,f2h1,f3h3,label1,label2,label3);}
   else{
     f2histos[0]->Draw();
     f3histos[0]->Draw("same");
     setLegend2(l,f2h1,f3h3,label2,label3);
   }
   
   l->Draw();
   output="efficiencyVsEta";  printCanvas(canvas,output,input);
   delete l;


   // ----
   l = new TLegend(0.40,0.75,0.80,0.90);
   if(!skipPt1){
     f1histos[1]->Draw();
     f2histos[1]->Draw("sames");
     f3histos[1]->Draw("sames");
     setLegend3(l,f1h1,f2h1,f3h3,label1,label2,label3);}   
   else{
     f2histos[1]->Draw();
     f3histos[1]->Draw("sames");
     setLegend2(l,f2h1,f3h3,label2,label3);
   }

   l->Draw();
   output="fakerateVsEta"; printCanvas(canvas,output,input);
   delete l;


   // ----
   canvas->SetLogx();
   l = new TLegend(0.40,0.40,0.80,0.55);
   f1histos[2]->Draw();
   f2histos[2]->Draw("same");
   f3histos[2]->Draw("same");


   setLegend3(l,f1h1,f2h1,f3h3,label4,label5,label6);
   l->Draw();
   output="efficiencyVsPt"; printCanvas(canvas,output,input);
   delete l;


   // ----
   canvas->SetLogx(1);
   f1histos[3]->Draw();
   f2histos[3]->Draw("same");
   f3histos[3]->Draw("same");

   l = new TLegend(0.30,0.65,0.70,0.80); //all vs highpurity
   setLegend3(l,f1h1,f2h1,f3h3,label4,label5,label6);
   l->Draw();
   output="fakerateVsPt"; printCanvas(canvas,output,input);
   delete l;



 }
 
 if(hitEff){
   TFile * file1 = new TFile(fileName1); TDirectory * dir1=file1->GetDirectory("");
   TFile * file2 = new TFile(fileName2); TDirectory * dir2=file2->GetDirectory("");
   TFile * file3 = new TFile(fileName3); TDirectory * dir3=file3->GetDirectory("");

   TH1 *f1h7,*f2h7,*f3h7;
   TH1 *f1h8,*f2h8,*f3h8;


   //--- hit-finding plots
   TH2F* tmpH2;
   dir1->GetObject(collname1+"/hitFindingEff_vs_eta",tmpH2); f1h7 = (TH1*) tmpH2->ProfileX()->DrawCopy()->Clone();
   dir2->GetObject(collname1+"/hitFindingEff_vs_eta",tmpH2); f2h7 = (TH1*) tmpH2->ProfileX()->DrawCopy()->Clone();
   dir3->GetObject(collname1+"/hitFindingEff_vs_eta",tmpH2); f3h7 = (TH1*) tmpH2->ProfileX()->DrawCopy()->Clone();
   f1h7->GetYaxis()->SetRangeUser(0.5,1.01);
   f2h7->GetYaxis()->SetRangeUser(0.5,1.01);
   f3h7->GetYaxis()->SetRangeUser(0.5,1.01);
   setHistoLabels(f1h7,0.05,1.2,0.07,0.6,"Hit-finding Efficiency","#eta");
   setHistoLabels(f2h7,0.05,1.2,0.07,0.6,"Hit-finding Efficiency","#eta");

   dir1->GetObject(collname1+"/hitPurity_vs_eta",tmpH2); f1h8 = (TH1*) tmpH2->ProfileX()->DrawCopy()->Clone();
   dir2->GetObject(collname1+"/hitPurity_vs_eta",tmpH2); f2h8 = (TH1*) tmpH2->ProfileX()->DrawCopy()->Clone();
   dir3->GetObject(collname1+"/hitPurity_vs_eta",tmpH2); f3h8 = (TH1*) tmpH2->ProfileX()->DrawCopy()->Clone();
   f1h8->GetYaxis()->SetRangeUser(0.5,1.01);
   f2h8->GetYaxis()->SetRangeUser(0.5,1.01);
   f3h8->GetYaxis()->SetRangeUser(0.5,1.01);
   setHistoLabels(f1h8,0.05,1.2,0.07,0.6,"Hits-on-Track Purity","#eta");
   setHistoLabels(f2h8,0.05,1.2,0.07,0.6,"Hits-on-Track Purity","#eta");


   TH1 * f1histos[2]={f1h7,f1h8};
   TH1 * f2histos[2]={f2h7,f2h8};
   TH1 * f3histos[2]={f3h7,f3h8};

   for(int i=0; i<2;i++){
     setHistoStyle(f1histos[i],f2histos[i],f3histos[i]);
   }

   canvas = new TCanvas("Tracks","hit-finding efficiency",500,500);
   canvas->cd();


   // ----
   TLegend* l = new TLegend(0.50,0.15,0.90,0.35); 
   canvas->SetLogx(0);
   if(!skipPt1){
     f1histos[0]->Draw();
     f2histos[0]->Draw("sames");
     f3histos[0]->Draw("sames");
     setLegend3(l,f1h7,f2h7,f3h7,label1,label2,label3);}
   else{
     f2histos[0]->Draw();
     f3histos[0]->Draw("sames");
     setLegend2(l,f2h7,f3h7,label2,label3);
   }

   l->Draw();
   output="hitFindingEffVsEta"; printCanvas(canvas,output,input);
   delete l;


   // ----
   l = new TLegend(0.50,0.15,0.90,0.35); 
   canvas->SetLogx(0);
   if(!skipPt1){
     f1histos[1]->Draw();
     f2histos[1]->Draw("sames");
     f3histos[1]->Draw("sames");
     setLegend3(l,f1h7,f2h7,f3h7,label1,label2,label3);}
   else{
     f2histos[1]->Draw();
     f3histos[1]->Draw("sames");
     setLegend2(l,f2h7,f3h7,label2,label3);
   }

   l->Draw();
   output="hitPurityVsEta"; printCanvas(canvas,output,input);
   delete l;

 }


 if(resVsEta){
   //===== resolutions vs eta
   TFile * file1 = new TFile(fileName1); TDirectory * dir1=file1->GetDirectory("");
   TFile * file2 = new TFile(fileName2); TDirectory * dir2=file2->GetDirectory("");
   TFile * file3 = new TFile(fileName3); TDirectory * dir3=file3->GetDirectory("");
      


   TH1F *f1h1(0),*f2h1,*f3h1;
   TH1F *f1h2,*f2h2,*f3h2;
   TH1F *f1h3,*f2h3,*f3h3;
   TH1F *f1h4,*f2h4,*f3h4;
   TH1F *f1h5,*f2h5,*f3h5;
    
   TH1F *f1h1RMS,*f2h1RMS,*f3h1RMS;
   TH1F *f1h2RMS,*f2h2RMS,*f3h2RMS;
   TH1F *f1h3RMS,*f2h3RMS,*f3h3RMS;
   TH1F *f1h4RMS,*f2h4RMS,*f3h4RMS;
   TH1F *f1h5RMS,*f2h5RMS,*f3h5RMS;

   TH1F *f1h1Mean,*f2h1Mean,*f3h1Mean;
   TH1F *f1h2Mean,*f2h2Mean,*f3h2Mean;
   TH1F *f1h3Mean,*f2h3Mean,*f3h3Mean;
   TH1F *f1h4Mean,*f2h4Mean,*f3h4Mean;
   TH1F *f1h5Mean,*f2h5Mean,*f3h5Mean;

   TH1F *f1h1Avg,*f2h1Avg,*f3h1Avg;
   TH1F *f1h2Avg,*f2h2Avg,*f3h2Avg;
   TH1F *f1h3Avg,*f2h3Avg,*f3h3Avg;
   TH1F *f1h4Avg,*f2h4Avg,*f3h4Avg;
   TH1F *f1h5Avg,*f2h5Avg,*f3h5Avg;

   f1h1 = new TH1F();    f2h1 = new TH1F();    f3h1 = new TH1F();
   f1h2 = new TH1F();    f2h2 = new TH1F();    f3h2 = new TH1F();
   f1h3 = new TH1F();    f2h3 = new TH1F();    f3h3 = new TH1F();
   f1h4 = new TH1F();    f2h4 = new TH1F();    f3h4 = new TH1F();
   f1h5 = new TH1F();    f2h5 = new TH1F();    f3h5 = new TH1F();

   f1h1RMS = new TH1F();    f2h1RMS = new TH1F();    f3h1RMS = new TH1F();
   f1h2RMS = new TH1F();    f2h2RMS = new TH1F();    f3h2RMS = new TH1F();
   f1h3RMS = new TH1F();    f2h3RMS = new TH1F();    f3h3RMS = new TH1F();
   f1h4RMS = new TH1F();    f2h4RMS = new TH1F();    f3h4RMS = new TH1F();
   f1h5RMS = new TH1F();    f2h5RMS = new TH1F();    f3h5RMS = new TH1F();

   f1h1Mean = new TH1F();    f2h1Mean = new TH1F();    f3h1Mean = new TH1F();
   f1h2Mean = new TH1F();    f2h2Mean = new TH1F();    f3h2Mean = new TH1F();
   f1h3Mean = new TH1F();    f2h3Mean = new TH1F();    f3h3Mean = new TH1F();
   f1h4Mean = new TH1F();    f2h4Mean = new TH1F();    f3h4Mean = new TH1F();
   f1h5Mean = new TH1F();    f2h5Mean = new TH1F();    f3h5Mean = new TH1F();

   f1h1Avg = new TH1F();    f2h1Avg = new TH1F();    f3h1Avg = new TH1F();
   f1h2Avg = new TH1F();    f2h2Avg = new TH1F();    f3h2Avg = new TH1F();
   f1h3Avg = new TH1F();    f2h3Avg = new TH1F();    f3h3Avg = new TH1F();
   f1h4Avg = new TH1F();    f2h4Avg = new TH1F();    f3h4Avg = new TH1F();
   f1h5Avg = new TH1F();    f2h5Avg = new TH1F();    f3h5Avg = new TH1F();


   // --- 
   TH2F* tmpH2;
   string name;


   
   if(!skipPt1){
     name="phires_vs_eta_pt1";
     dir1->GetObject(collname1+"/phires_vs_eta",tmpH2);  extractReso(tmpH2,f1h1,f1h1RMS,f1h1Mean,f1h1Avg,
								     tailSignPhi,name,input,false,false);
   }
   name="phires_vs_eta_pt10";
   dir2->GetObject(collname2+"/phires_vs_eta",tmpH2);  extractReso(tmpH2,f2h1,f2h1RMS,f2h1Mean,f2h1Avg,
								     tailSignPhi,name,input,false,false); 
   name="phires_vs_eta_pt100";
   dir3->GetObject(collname3+"/phires_vs_eta",tmpH2);  extractReso(tmpH2,f3h1,f3h1RMS,f3h1Mean,f3h1Avg,
   								     tailSignPhi,name,input,false,false); 


   if(!skipPt1){
     name="dxyres_vs_eta_pt1";
     dir1->GetObject(collname1+"/dxyres_vs_eta",tmpH2);   extractReso(tmpH2,f1h3,f1h3RMS,f1h3Mean,f1h3Avg,
								      tailSignD0,name,input,false,false); 
   }
   name="dxyres_vs_eta_pt10";
   dir2->GetObject(collname2+"/dxyres_vs_eta",tmpH2);   extractReso(tmpH2,f2h3,f2h3RMS,f2h3Mean,f2h3Avg,
								      tailSignD0,name,input,false,false); 
   name="dxyres_vs_eta_pt100";
   dir3->GetObject(collname3+"/dxyres_vs_eta",tmpH2);   extractReso(tmpH2,f3h3,f3h3RMS,f3h3Mean,f3h3Avg,
   tailSignD0,name,input,false,false); 



   if(!skipPt1){
   name="ptres_vs_eta_pt1";
   dir1->GetObject(collname1+"/ptres_vs_eta",tmpH2);   extractReso(tmpH2,f1h5,f1h5RMS,f1h5Mean,f1h5Avg,
                                                                   tailSignPt,name,input,false,false); 
   }

   name="ptres_vs_eta_pt10";
   dir2->GetObject(collname2+"/ptres_vs_eta",tmpH2);   extractReso(tmpH2,f2h5,f2h5RMS,f2h5Mean,f2h5Avg,
   								     tailSignPt,name,input,false,false); 
   name="ptres_vs_eta_pt100";
   dir3->GetObject(collname3+"/ptres_vs_eta",tmpH2);   extractReso(tmpH2,f3h5,f3h5RMS,f3h5Mean,f3h5Avg,
   								   tailSignPt,name,input,false,false); 

   							   

   /*
   if(!skipPt1){
     name="cotThetares_vs_eta_pt1";
     dir1->GetObject(collname1+"/cotThetares_vs_eta",tmpH2);  extractReso(tmpH2,f1h2,f1h2RMS,f1h2Mean,f1h2Avg,
									  tailSignTheta,name,input,false,false); 
   }
   name="cotThetares_vs_eta_pt10";
   dir2->GetObject(collname2+"/cotThetares_vs_eta",tmpH2);  extractReso(tmpH2,f2h2,f2h2RMS,f2h2Mean,f2h2Avg,
									  tailSignTheta,name,input,false,false);  
   name="cotThetares_vs_eta_pt100";
   dir3->GetObject(collname3+"/cotThetares_vs_eta",tmpH2);  extractReso(tmpH2,f3h2,f3h2RMS,f3h2Mean,f3h2Avg,
									  tailSignTheta,name,input,false,false);  

 
   if(!skipPt1){
     name="dzres_vs_eta_pt1";
     dir1->GetObject(collname1+"/dzres_vs_eta",tmpH2);   extractReso(tmpH2,f1h4,f1h4RMS,f1h4Mean,f1h4Avg,
								     tailSignDz,name,input,false,false); 
   }
   name="dzres_vs_eta_pt10";
   dir2->GetObject(collname2+"/dzres_vs_eta",tmpH2);   extractReso(tmpH2,f2h4,f2h4RMS,f2h4Mean,f2h4Avg,
								     tailSignDz,name,input,false,false); 
   name="dzres_vs_eta_pt100";
   dir3->GetObject(collname3+"/dzres_vs_eta",tmpH2);   extractReso(tmpH2,f3h4,f3h4RMS,f3h4Mean,f3h4Avg,
   								     tailSignDz,name,input,false,false); 
   */

   // ---

   
   //

   f1h1->Scale(1000.);    f2h1->Scale(1000.);    f3h1->Scale(1000.); 
   f1h1RMS->Scale(1000.); f2h1RMS->Scale(1000.); f3h1RMS->Scale(1000.); 
   f1h1Mean->Scale(1000.);  f2h1Mean->Scale(1000.);  f3h1Mean->Scale(1000.); 
   f1h1Avg->Scale(1000.);   f2h1Avg->Scale(1000.);   f3h1Avg->Scale(1000.);
   f1h1->GetYaxis()->SetRangeUser(resDPhiVsEtaMin,resDPhiVsEtaMax);
   f2h1->GetYaxis()->SetRangeUser(resDPhiVsEtaMin,resDPhiVsEtaMax);
   f1h1Mean->GetYaxis()->SetRangeUser(-20,10);
   f2h1Mean->GetYaxis()->SetRangeUser(-20,10);
   setHistoLabels(f1h1,0.05,1.2,0.07,0.7,"Resolution in #phi (10^{-3}radians)","#eta");
   setHistoLabels(f2h1,0.05,1.2,0.07,0.7,"Resolution in #phi (10^{-3}radians)","#eta");
   setHistoLabels(f1h1Mean,0.05,1.2,0.07,0.7,"Bias in #phi (10^{-3}radians)","#eta");
   setHistoLabels(f2h1Mean,0.05,1.2,0.07,0.7,"Bias in #phi (10^{-3}radians)","#eta");

   
   f1h2->Scale(1000.);    f2h2->Scale(1000.);      f3h2->Scale(1000.); 
   f1h2RMS->Scale(1000.); f2h2RMS->Scale(1000.);   f3h2RMS->Scale(1000.); 
   f1h2Mean->Scale(1000.);  f2h2Mean->Scale(1000.);  f3h2Mean->Scale(1000.); 
   f1h2Avg->Scale(1000.);   f2h2Avg->Scale(1000.);   f3h2Avg->Scale(1000.);
   f1h2->GetYaxis()->SetRangeUser(resDThetaVsEtaMin,resDThetaVsEtaMax);
   f2h2->GetYaxis()->SetRangeUser(resDThetaVsEtaMin,resDThetaVsEtaMax);
   f1h2Mean->GetYaxis()->SetRangeUser(-4,8);
   f2h2Mean->GetYaxis()->SetRangeUser(-4,8);
   setHistoLabels(f1h2,0.05,1.2,0.07,0.7,"Resolution in cot(#theta) (10^{-3})","#eta");
   setHistoLabels(f2h2,0.05,1.2,0.07,0.7,"Resolution in cot(#theta) (10^{-3})","#eta");
   setHistoLabels(f1h2Mean,0.05,1.2,0.07,0.7,"Bias in cot(#theta) (10^{-3})","#eta");
   setHistoLabels(f2h2Mean,0.05,1.2,0.07,0.7,"Bias in cot(#theta) (10^{-3})","#eta");


   f1h3->Scale(10000.);    f2h3->Scale(10000.);    f3h3->Scale(10000.); 
   f1h3RMS->Scale(10000.); f2h3RMS->Scale(10000.); f3h3RMS->Scale(10000.); 
   f1h3Mean->Scale(10000.);  f2h3Mean->Scale(10000.);  f3h3Mean->Scale(10000.); 
   f1h3Avg->Scale(10000.);   f2h3Avg->Scale(10000.);   f3h3Avg->Scale(10000.);
   f1h3->GetYaxis()->SetRangeUser(resD0vsEtaMin,resD0vsEtaMax);  
   f2h3->GetYaxis()->SetRangeUser(resD0vsEtaMin,resD0vsEtaMax);  
   f1h3Mean->GetYaxis()->SetRangeUser(-100,600);
   f2h3Mean->GetYaxis()->SetRangeUser(-100,600);
   setHistoLabels(f1h3,0.05,1.2,0.07,0.7,"Resolution in d_{0} (#mum)","#eta");
   setHistoLabels(f2h3,0.05,1.2,0.07,0.7,"Resolution in d_{0} (#mum)","#eta");
   setHistoLabels(f1h3Mean,0.05,1.2,0.07,0.7,"Bias in d_{0} (#mum)","#eta");
   setHistoLabels(f2h3Mean,0.05,1.2,0.07,0.7,"Bias in d_{0} (#mum)","#eta");



   f1h4->Scale(10000.);    f2h4->Scale(10000.);    f3h4->Scale(10000.); 
   f1h4RMS->Scale(10000.); f2h4RMS->Scale(10000.); f3h4RMS->Scale(10000.); 
   f1h4Mean->Scale(10000.);  f2h4Mean->Scale(10000.);  f3h4Mean->Scale(10000.); 
   f1h4Avg->Scale(10000.);   f2h4Avg->Scale(10000.);   f3h4Avg->Scale(10000.);
   f1h4->GetYaxis()->SetRangeUser(resDzVsEtaMin,resDzVsEtaMax);
   f2h4->GetYaxis()->SetRangeUser(resDzVsEtaMin,resDzVsEtaMax);
   f1h4Mean->GetYaxis()->SetRangeUser(-200,200);
   f2h4Mean->GetYaxis()->SetRangeUser(-200,200);
   setHistoLabels(f1h4,0.05,1.2,0.07,0.7,"Resolution in z_{0} (#mum)","#eta");
   setHistoLabels(f2h4,0.05,1.2,0.07,0.7,"Resolution in z_{0} (#mum)","#eta");
   setHistoLabels(f1h4Mean,0.05,1.2,0.07,0.7,"Bias in z_{0} (#mum)","#eta");
   setHistoLabels(f2h4Mean,0.05,1.2,0.07,0.7,"Bias in z_{0} (#mum)","#eta");


   f1h5->Scale(100.);    f2h5->Scale(100.);    f3h5->Scale(100.); 
   f1h5RMS->Scale(100.); f2h5RMS->Scale(100.); f3h5RMS->Scale(100.); 
   f1h5Mean->Scale(100.);  f2h5Mean->Scale(100.);  f3h5Mean->Scale(100.); 
   f1h5Avg->Scale(100.);   f2h5Avg->Scale(100.);   f3h5Avg->Scale(100.);
   f1h5->GetYaxis()->SetRangeUser(resPtVsEtaMin,resPtVsEtaMax);
   f2h5->GetYaxis()->SetRangeUser(resPtVsEtaMin,resPtVsEtaMax);
   f1h5Mean->GetYaxis()->SetRangeUser(-55,30);
   f2h5Mean->GetYaxis()->SetRangeUser(-55,30);
   setHistoLabels(f1h5,0.05,1.2,0.07,0.7,"(Resolution in p_{T})/p_{T} (%)","#eta");
   setHistoLabels(f2h5,0.05,1.2,0.07,0.7,"(Resolution in p_{T})/p_{T} (%)","#eta");
   setHistoLabels(f1h5Mean,0.05,1.2,0.07,0.7,"(Bias in p_{T})/p_{T} (%)","#eta");
   setHistoLabels(f2h5Mean,0.05,1.2,0.07,0.7,"(Bias in p_{T})/p_{T} (%)","#eta");



   TH1F* f1histos2[5] = {f1h1,f1h2,f1h3,f1h4,f1h5};
   TH1F* f2histos2[5] = {f2h1,f2h2,f2h3,f2h4,f2h5};
   TH1F* f3histos2[5] = {f3h1,f3h2,f3h3,f3h4,f3h5};
   TH1F* f1histos2RMS[5] = {f1h1RMS,f1h2RMS,f1h3RMS,f1h4RMS,f1h5RMS};
   TH1F* f2histos2RMS[5] = {f2h1RMS,f2h2RMS,f2h3RMS,f2h4RMS,f2h5RMS};
   TH1F* f3histos2RMS[5] = {f3h1RMS,f3h2RMS,f3h3RMS,f3h4RMS,f3h5RMS};

   TH1F* f1hMeans[5] = {f1h1Mean,f1h2Mean,f1h3Mean,f1h4Mean,f1h5Mean};
   TH1F* f2hMeans[5] = {f2h1Mean,f2h2Mean,f2h3Mean,f2h4Mean,f2h5Mean};
   TH1F* f3hMeans[5] = {f3h1Mean,f3h2Mean,f3h3Mean,f3h4Mean,f3h5Mean};
   TH1F* f1hAvgs[5] =  {f1h1Avg,f1h2Avg,f1h3Avg,f1h4Avg,f1h5Avg,};
   TH1F* f2hAvgs[5] =  {f2h1Avg,f2h2Avg,f2h3Avg,f2h4Avg,f2h5Avg,};
   TH1F* f3hAvgs[5] =  {f3h1Avg,f3h2Avg,f3h3Avg,f3h4Avg,f3h5Avg,};



   for(int i=0; i<5;i++){
     setHistoStyle(f1histos2[i],f2histos2[i],f3histos2[i]);
     setHistoStyle(f1histos2RMS[i],f2histos2RMS[i],f3histos2RMS[i],true);
     
     setHistoStyle(f1hMeans[i],f2hMeans[i],f3hMeans[i]);
     setHistoStyle(f1hAvgs[i],f2hAvgs[i],f3hAvgs[i],true);
   }

   /* //what is this about ?? 
   setHistoStyle(f1h5,f1h5,f1h5);
   setHistoStyle(f1h5RMS,f1h5RMS,f1h5RMS,true);
   setHistoStyle(f1h5Mean,f1h5Mean,f1h5Mean);
   setHistoStyle(f1h5Avg,f1h5Avg,f1h5Avg,true);
   */
   
   canvas = new TCanvas("Tracks","Resolution vs Eta",500,500);
   canvas->cd();


   // -----
   canvas->SetLogx(0);
   canvas->SetLogy(1);
   TLegend* l = new TLegend(0.35,0.77,0.75,0.92); //single mu

   if(!skipPt1){
     f1histos2[0]->Draw();         f1histos2RMS[0]->Draw("sames");
     f2histos2[0]->Draw("sames");  f2histos2RMS[0]->Draw("sames");
     f3histos2[0]->Draw("sames");  f3histos2RMS[0]->Draw("sames");
     setLegend3(l,f1histos2[0],f2histos2[0],f3histos2[0],label1,label2,label3);}
   else{
     f2histos2[0]->Draw();  f2histos2RMS[0]->Draw("same");
     f3histos2[0]->Draw("sames");  f3histos2RMS[0]->Draw("sames");
     setLegend2(l,f2histos2[0],f3histos2[0],label2,label3);
   }
 
   l->Draw();
   output="resolutionPhiVsEta"; printCanvas(canvas,output,input);
   delete l;



   // -----
   canvas->SetLogx(0);
   canvas->SetLogy(1);
   l = new TLegend(0.35,0.75,0.75,0.90); //single mu

   if(!skipPt1){
     f1histos2[1]->Draw();        f1histos2RMS[1]->Draw("sames");
     f2histos2[1]->Draw("sames"); f2histos2RMS[1]->Draw("sames");
     f3histos2[1]->Draw("sames"); f3histos2RMS[1]->Draw("sames");
     setLegend3(l,f1histos2[1],f2histos2[1],f3histos2[1],label1,label2,label3);}
   else{
     f2histos2[1]->Draw(); f2histos2RMS[1]->Draw("same");
     f3histos2[1]->Draw("sames"); f3histos2RMS[1]->Draw("sames");
     setLegend2(l,f2histos2[1],f3histos2[1],label2,label3);
   }

   l->Draw();
   output="resolutionThetaVsEta"; printCanvas(canvas,output,input);
   delete l;

   // -----
   canvas->SetLogx(0);
   canvas->SetLogy(1);

   //l = new TLegend(0.35,0.75,0.75,0.90); //single mu
   l = new TLegend(0.35,0.15,0.75,0.30); //single mu
   if(!skipPt1){
     f1histos2[2]->Draw();        f1histos2RMS[2]->Draw("sames");
     f2histos2[2]->Draw("sames"); f2histos2RMS[2]->Draw("sames");
     f3histos2[2]->Draw("sames"); f3histos2RMS[2]->Draw("sames");
     setLegend3(l,f1histos2[2],f2histos2[2],f3histos2[2],label1,label2,label3);
   }else{
     f2histos2[2]->Draw(); f2histos2RMS[2]->Draw("same");
     f3histos2[2]->Draw("sames"); f3histos2RMS[2]->Draw("sames");
     setLegend2(l,f2histos2[2],f3histos2[2],label2,label3);
   }

   l->Draw();
   output="resolutionD0VsEta"; printCanvas(canvas,output,input);
   delete l;

   // -----
   canvas->SetLogx(0);
   canvas->SetLogy(1);
   l = new TLegend(0.50,0.15,0.90,0.30); //single mu

   if(!skipPt1){
     f1histos2[3]->Draw();        f1histos2RMS[3]->Draw("sames");
     f2histos2[3]->Draw("sames"); f2histos2RMS[3]->Draw("sames");
     f3histos2[3]->Draw("sames"); f3histos2RMS[3]->Draw("sames");
     setLegend3(l,f1histos2[3],f2histos2[3],f3histos2[3],label1,label2,label3);
   }else{
     f2histos2[3]->Draw(); f2histos2RMS[3]->Draw("same");
     f3histos2[3]->Draw("sames"); f3histos2RMS[3]->Draw("sames");
     setLegend2(l,f2histos2[3],f3histos2[3],label2,label3);
   }

   l->Draw();
   output="resolutionDzVsEta"; printCanvas(canvas,output,input);
   delete l;

   
   // -----
   canvas->SetLogx(0);
   canvas->SetLogy(1);
   l = new TLegend(0.40,0.77,0.80,0.92); //single mu   
   if(!skipPt1){
     //f1h5->Draw(); f1h5RMS->Draw("same");
     f1histos2[4]->Draw();        f1histos2RMS[4]->Draw("sames");
     f2histos2[4]->Draw("sames"); f2histos2RMS[4]->Draw("sames");
     f3histos2[4]->Draw("sames"); f3histos2RMS[4]->Draw("sames");
     setLegend3(l,f1histos2[4],f2histos2[4],f3histos2[4],label1,label2,label3);
   }else{
     f1h5->Draw(); f1h5RMS->Draw("same");
     f2histos2[4]->Draw(); f2histos2RMS[4]->Draw("same");
     f3histos2[4]->Draw("sames"); f3histos2RMS[4]->Draw("sames");
     setLegend2(l,f2histos2[4],f3histos2[4],label2,label3);
   }

   l->Draw();
   output="resolutionPtVsEta"; printCanvas(canvas,output,input);
   delete l;



   // ======= here print the bias plots =======
   canvas->SetLogx(0);
   canvas->SetLogy(0);
   l = new TLegend(0.35,0.75,0.75,0.90); 
   if(!skipPt1){
     f1hMeans[0]->Draw();         f1hAvgs[0]->Draw("sames");
     f2hMeans[0]->Draw("sames");  f2hAvgs[0]->Draw("sames");
     f3hMeans[0]->Draw("sames");  f3hAvgs[0]->Draw("sames");
     setLegend3(l,f1hMeans[0],f2hMeans[0],f3hMeans[0],label1,label2,label3);}
   else{
     f2hMeans[0]->Draw();         f2hAvgs[0]->Draw("same");
     f3hMeans[0]->Draw("sames");  f3hAvgs[0]->Draw("sames");
     setLegend2(l,f2hMeans[0],f3hMeans[0],label2,label3);
   }
   l->Draw();
   output="biasPhiVsEta"; printCanvas(canvas,output,input);
   delete l;

   l = new TLegend(0.35,0.77,0.75,0.92); 
   if(!skipPt1){
     f1hMeans[1]->Draw();         f1hAvgs[1]->Draw("sames");
     f2hMeans[1]->Draw("sames");  f2hAvgs[1]->Draw("sames");
     f3hMeans[1]->Draw("sames");  f3hAvgs[1]->Draw("sames");
     setLegend3(l,f1hMeans[1],f2hMeans[1],f3hMeans[1],label1,label2,label3);}
   else{
     f2hMeans[1]->Draw();         f2hAvgs[1]->Draw("same");
     f3hMeans[1]->Draw("sames");  f3hAvgs[1]->Draw("sames");
     setLegend2(l,f2hMeans[1],f3hMeans[1],label2,label3);
   }
   l->Draw();
   output="biasThetaVsEta"; printCanvas(canvas,output,input);
   delete l;

   l = new TLegend(0.35,0.75,0.75,0.90); 
   if(!skipPt1){
     f1hMeans[2]->Draw();         f1hAvgs[2]->Draw("sames");
     f2hMeans[2]->Draw("sames");  f2hAvgs[2]->Draw("sames");
     f3hMeans[2]->Draw("sames");  f3hAvgs[2]->Draw("sames");
     setLegend3(l,f1hMeans[2],f2hMeans[2],f3hMeans[0],label1,label2,label3);}
   else{
     f2hMeans[2]->Draw();         f2hAvgs[2]->Draw("same");
     f3hMeans[2]->Draw("sames");  f3hAvgs[2]->Draw("sames");
     setLegend2(l,f2hMeans[2],f3hMeans[2],label2,label3);
   }
   l->Draw();
   output="biasD0VsEta"; printCanvas(canvas,output,input);
   delete l;

   l = new TLegend(0.35,0.75,0.75,0.90); 
   if(!skipPt1){
     f1hMeans[3]->Draw();         f1hAvgs[3]->Draw("sames");
     f2hMeans[3]->Draw("sames");  f2hAvgs[3]->Draw("sames");
     f3hMeans[3]->Draw("sames");  f3hAvgs[3]->Draw("sames");
     setLegend3(l,f1hMeans[3],f2hMeans[3],f3hMeans[0],label1,label2,label3);}
   else{
     f2hMeans[3]->Draw();         f2hAvgs[3]->Draw("same");
     f3hMeans[3]->Draw("sames");  f3hAvgs[3]->Draw("sames");
     setLegend2(l,f2hMeans[3],f3hMeans[3],label2,label3);
   }
   l->Draw();
   output="biasDzVsEta"; printCanvas(canvas,output,input);
   delete l;

   l = new TLegend(0.35,0.75,0.75,0.90); 
   if(!skipPt1){
     f1hMeans[4]->Draw();         f1hAvgs[4]->Draw("sames");
     f2hMeans[4]->Draw("sames");  f2hAvgs[4]->Draw("sames");
     f3hMeans[4]->Draw("sames");  f3hAvgs[4]->Draw("sames");
     setLegend3(l,f1hMeans[4],f2hMeans[4],f3hMeans[0],label1,label2,label3);}
   else{
     f2hMeans[4]->Draw();         f2hAvgs[4]->Draw("same");
     f3hMeans[4]->Draw("sames");  f3hAvgs[4]->Draw("sames");
     setLegend2(l,f2hMeans[4],f3hMeans[4],label2,label3);
   }
   l->Draw();
   output="biasPtVsEta"; printCanvas(canvas,output,input);
   delete l;


 }

 if(resVsPt){
   //===== resolutions vs pt
   TFile * file4 = new TFile(fileName4); TDirectory * dir1=file4->GetDirectory("");
   TFile * file5 = new TFile(fileName5); TDirectory * dir2=file5->GetDirectory("");
   TFile * file6 = new TFile(fileName6); TDirectory * dir3=file6->GetDirectory("");
      
   
   TH1F *f1h1,*f2h1,*f3h1;
   TH1F *f1h2,*f2h2,*f3h2;
   TH1F *f1h3,*f2h3,*f3h3;
   TH1F *f1h4,*f2h4,*f3h4;
   TH1F *f1h5,*f2h5,*f3h5;

   
   TH1F *f1h1RMS,*f2h1RMS,*f3h1RMS;
   TH1F *f1h2RMS,*f2h2RMS,*f3h2RMS;
   TH1F *f1h3RMS,*f2h3RMS,*f3h3RMS;
   TH1F *f1h4RMS,*f2h4RMS,*f3h4RMS;
   TH1F *f1h5RMS,*f2h5RMS,*f3h5RMS;

   TH1F *f1h1Mean,*f2h1Mean,*f3h1Mean;
   TH1F *f1h2Mean,*f2h2Mean,*f3h2Mean;
   TH1F *f1h3Mean,*f2h3Mean,*f3h3Mean;
   TH1F *f1h4Mean,*f2h4Mean,*f3h4Mean;
   TH1F *f1h5Mean,*f2h5Mean,*f3h5Mean;

   TH1F *f1h1Avg,*f2h1Avg,*f3h1Avg;
   TH1F *f1h2Avg,*f2h2Avg,*f3h2Avg;
   TH1F *f1h3Avg,*f2h3Avg,*f3h3Avg;
   TH1F *f1h4Avg,*f2h4Avg,*f3h4Avg;
   TH1F *f1h5Avg,*f2h5Avg,*f3h5Avg;

   
   // --- 
   TH2F* tmpH2;

   string name;

   name="dxyres_vs_pt_barrel";
   dir1->GetObject(collname4+"/dxyres_vs_pt",tmpH2);   extractReso(tmpH2,f1h3,f1h3RMS,f1h3Mean,f1h3Avg,
								     tailSignD0_pT,name,input,true); 
   name="dxyres_vs_pt_transition";
   dir2->GetObject(collname5+"/dxyres_vs_pt",tmpH2);   extractReso(tmpH2,f2h3,f2h3RMS,f2h3Mean,f2h3Avg,
								     tailSignD0_pT,name,input,true); 
   name="dxyres_vs_pt_endcap";
   dir3->GetObject(collname6+"/dxyres_vs_pt",tmpH2);   extractReso(tmpH2,f3h3,f3h3RMS,f3h3Mean,f3h3Avg,
								     tailSignD0_pT,name,input,true); 
   name="phires_vs_pt_barrel";
   dir1->GetObject(collname4+"/phires_vs_pt",tmpH2);  extractReso(tmpH2,f1h1,f1h1RMS,f1h1Mean,f1h1Avg,
								    tailSignPhi_pT,name,input,true);
   name="phires_vs_pt_transition";
   dir2->GetObject(collname5+"/phires_vs_pt",tmpH2);  extractReso(tmpH2,f2h1,f2h1RMS,f2h1Mean,f2h1Avg,
								    tailSignPhi_pT,name,input,true);
   name="phires_vs_pt_endcap";
   dir3->GetObject(collname6+"/phires_vs_pt",tmpH2);  extractReso(tmpH2,f3h1,f3h1RMS,f3h1Mean,f3h1Avg,
								    tailSignPhi_pT,name,input,true);

   name="ptres_vs_pt_barrel";
   dir1->GetObject(collname4+"/ptres_vs_pt",tmpH2);   extractReso(tmpH2,f1h5,f1h5RMS,f1h5Mean,f1h5Avg,
								    tailSignPt_pT,name,input,true); 
   name="ptres_vs_pt_transition";
   dir2->GetObject(collname5+"/ptres_vs_pt",tmpH2);   extractReso(tmpH2,f2h5,f2h5RMS,f2h5Mean,f2h5Avg,
								    tailSignPt_pT,name,input,true); 
   name="ptres_vs_pt_endcap";
   dir3->GetObject(collname6+"/ptres_vs_pt",tmpH2);   extractReso(tmpH2,f3h5,f3h5RMS,f3h5Mean,f3h5Avg,
								    tailSignPt_pT,name,input,true); 



   name="cotThetares_vs_pt_barrel";
   dir1->GetObject(collname4+"/cotThetares_vs_pt",tmpH2);  extractReso(tmpH2,f1h2,f1h2RMS,f1h2Mean,f1h2Avg,
									 tailSignTheta_pT,name,input,true); 
   name="cotThetares_vs_pt_transition";
   dir2->GetObject(collname5+"/cotThetares_vs_pt",tmpH2);  extractReso(tmpH2,f2h2,f2h2RMS,f2h2Mean,f2h2Avg,
									 tailSignTheta_pT,name,input,true);  
   name="cotThetares_vs_pt_endcap";
   dir3->GetObject(collname6+"/cotThetares_vs_pt",tmpH2);  extractReso(tmpH2,f3h2,f3h2RMS,f3h2Mean,f3h2Avg,
									 tailSignTheta_pT,name,input,true);  

   name="dzres_vs_pt_barrel";
   dir1->GetObject(collname4+"/dzres_vs_pt",tmpH2);   extractReso(tmpH2,f1h4,f1h4RMS,f1h4Mean,f1h4Avg,
								    tailSignDz_pT,name,input,true); 
   name="dzres_vs_pt_transition";
   dir2->GetObject(collname5+"/dzres_vs_pt",tmpH2);   extractReso(tmpH2,f2h4,f2h4RMS,f2h4Mean,f2h4Avg,
								    tailSignDz_pT,name,input,true); 
   name="dzres_vs_pt_endcap";
   dir3->GetObject(collname6+"/dzres_vs_pt",tmpH2);   extractReso(tmpH2,f3h4,f3h4RMS,f3h4Mean,f3h4Avg,
								    tailSignDz_pT,name,input,true); 

   // ---

   f1h1->Scale(1000.); f2h1->Scale(1000.); f3h1->Scale(1000.); 
   f1h1RMS->Scale(1000.); f2h1RMS->Scale(1000.); f3h1RMS->Scale(1000.); 
   f1h1->GetYaxis()->SetRangeUser(resDPhiVsPtMin,resDPhiVsPtMax);
   f1h1->GetXaxis()->SetRangeUser(0.2,150);
   f2h1->GetXaxis()->SetRangeUser(0.2,150);
   f3h1->GetXaxis()->SetRangeUser(0.2,150);
   setHistoLabels(f1h1,0.05,1.2,0.055,1.1,"resolution in #phi [10^{-3}]","p_{T} (GeV)");


   f1h2->Scale(1000.); f2h2->Scale(1000.); f3h2->Scale(1000.); 
   f1h2RMS->Scale(1000.); f2h2RMS->Scale(1000.); f3h2RMS->Scale(1000.); 
   f1h2->GetYaxis()->SetRangeUser(resDThetaVsPtMin,resDThetaVsPtMax);
   f1h2->GetXaxis()->SetRangeUser(0.2,150);
   f2h2->GetXaxis()->SetRangeUser(0.2,150);
   f3h2->GetXaxis()->SetRangeUser(0.2,150);
   setHistoLabels(f1h2,0.05,1.2,0.055,1.1,"resolution in cot(#theta) [10^{-3}]","p_{T} (GeV)");


   f1h3->Scale(10000.); f2h3->Scale(10000.); f3h3->Scale(10000.); 
   f1h3RMS->Scale(10000.); f2h3RMS->Scale(10000.); f3h3RMS->Scale(10000.); 
   f1h3->GetYaxis()->SetRangeUser(resD0vsPtMin,resD0vsPtMax);  
   f1h3->GetXaxis()->SetRangeUser(0.2,150);
   f2h3->GetXaxis()->SetRangeUser(0.2,150);
   f3h3->GetXaxis()->SetRangeUser(0.2,150);
   setHistoLabels(f1h3,0.05,1.2,0.055,1.1,"resolution in d_{0} [#mum]","p_{T} (GeV)");
 

   f1h4->Scale(10000.); f2h4->Scale(10000.); f3h4->Scale(10000.); 
   f1h4RMS->Scale(10000.); f2h4RMS->Scale(10000.); f3h4RMS->Scale(10000.); 
   f1h4->GetYaxis()->SetRangeUser(resDzVsPtMin,resDzVsPtMax);
   f1h4->GetXaxis()->SetRangeUser(0.2,150);
   f2h4->GetXaxis()->SetRangeUser(0.2,150);
   f3h4->GetXaxis()->SetRangeUser(0.2,150);
   setHistoLabels(f1h4,0.05,1.2,0.055,1.1,"resolution in z_{0} [#mum]","p_{T} (GeV)");


   f1h5->Scale(100.); f2h5->Scale(100.); f3h5->Scale(100.); 
   f1h5RMS->Scale(100.); f2h5RMS->Scale(100.); f3h5RMS->Scale(100.); 
   f1h5->GetYaxis()->SetRangeUser(resPtVsPtMin,resPtVsPtMax);
   f1h5->GetXaxis()->SetRangeUser(0.2,150);
   f2h5->GetXaxis()->SetRangeUser(0.2,150);
   f3h5->GetXaxis()->SetRangeUser(0.2,150);
   setHistoLabels(f1h5,0.05,1.2,0.055,1.1,"(resolution in p_{T})/p_{T} [%]","p_{T} (GeV)");



   TH1F* f1histos3[5] = {f1h1,f1h2,f1h3,f1h4,f1h5};
   TH1F* f2histos3[5] = {f2h1,f2h2,f2h3,f2h4,f2h5};
   TH1F* f3histos3[5] = {f3h1,f3h2,f3h3,f3h4,f3h5};

   TH1F* f1histos3RMS[5] = {f1h1RMS,f1h2RMS,f1h3RMS,f1h4RMS,f1h5RMS};
   TH1F* f2histos3RMS[5] = {f2h1RMS,f2h2RMS,f2h3RMS,f2h4RMS,f2h5RMS};
   TH1F* f3histos3RMS[5] = {f3h1RMS,f3h2RMS,f3h3RMS,f3h4RMS,f3h5RMS};

   for(int i=0; i<5;i++){
     setHistoStyle(f1histos3[i],f2histos3[i],f3histos3[i]);
     setHistoStyle(f1histos3RMS[i],f2histos3RMS[i],f3histos3RMS[i],true);
   }

   canvas = new TCanvas("Tracks","Resolution vs Pt",500,500);
   canvas->cd();

   // -----
   canvas->SetLogx(1);
   canvas->SetLogy(1);
   f1histos3[0]->Draw();         f1histos3RMS[0]->Draw("sames");
   f2histos3[0]->Draw("sames");  f2histos3RMS[0]->Draw("sames");
   f3histos3[0]->Draw("sames");  f3histos3RMS[0]->Draw("sames");
   
   TLegend* l = new TLegend(0.40,0.65,0.83,0.90);
   setLegend3(l,f1histos3[0],f2histos3[0],f3histos3[0],label4,label5,label6);
   l->Draw();
   output="resolutionPhiVsPt"; printCanvas(canvas,output,input);
   delete l;


   // -----
   canvas->SetLogx(1);
   canvas->SetLogy(1);
   f1histos3[1]->Draw();         f1histos3RMS[1]->Draw("sames");
   f2histos3[1]->Draw("sames");  f2histos3RMS[1]->Draw("sames");
   f3histos3[1]->Draw("sames");  f3histos3RMS[1]->Draw("sames");

   l = new TLegend(0.44,0.65,0.87,0.90);
   setLegend3(l,f1histos3[1],f2histos3[1],f3histos3[1],label4,label5,label6);
   l->Draw();
   output="resolutionThetaVsPt"; printCanvas(canvas,output,input);
   delete l;


   // -----
   canvas->SetLogx(1);
   canvas->SetLogy(1);

   f1histos3[2]->Draw();         f1histos3RMS[2]->Draw("sames");
   f2histos3[2]->Draw("sames");  f2histos3RMS[2]->Draw("sames");
   f3histos3[2]->Draw("sames");  f3histos3RMS[2]->Draw("sames");

   l = new TLegend(0.45,0.65,0.88,0.90);
   setLegend3(l,f1histos3[2],f2histos3[2],f3histos3[2],label4,label5,label6);
   l->Draw();
   output="resolutionD0VsPt"; printCanvas(canvas,output,input);
   delete l;


   // -----
   canvas->SetLogx(1);
   canvas->SetLogy(1);

   f1histos3[3]->Draw();         f1histos3RMS[3]->Draw("sames");
   f2histos3[3]->Draw("sames");  f2histos3RMS[3]->Draw("sames");
   f3histos3[3]->Draw("sames");  f3histos3RMS[3]->Draw("sames");
   
   l = new TLegend(0.44,0.65,0.86,0.90);
   setLegend3(l,f1histos3[3],f2histos3[3],f3histos3[3],label4,label5,label6);
   l->Draw();
   output="resolutionDzVsPt"; printCanvas(canvas,output,input);
   delete l;

   // -----
   canvas->SetLogx(1);
   canvas->SetLogy(1);
   f1histos3[4]->Draw();         f1histos3RMS[4]->Draw("sames");
   f2histos3[4]->Draw("sames");  f2histos3RMS[4]->Draw("sames");
   f3histos3[4]->Draw("sames");  f3histos3RMS[4]->Draw("sames");
   
   l = new TLegend(0.28,0.65,0.70,0.90);
   setLegend3(l,f1histos3[4],f2histos3[4],f3histos3[4],label4,label5,label6);
   l->Draw();
   output="resolutionPtVsPt"; printCanvas(canvas,output,input);
   delete l;
 }


 //delete canvas;
}



// ----------------------------------




void extractReso(TH2F* th2, 
		 TH1F*& hFit, TH1F*& hRMS,
		 TH1F*& hFitMean, TH1F*& hAverage,
		 int tailSign,string name, int type,bool variableBinning,bool wait)
{
  cout << "--- calling extracReso for input " << type << " ---" << endl;
  //int rangeType(0);  //inclusive
  //int rangeType(+1);    //run on side which is gaussian-like
  int rangeType(-1);  //run on side which has non-gaussian tail
  

  

  TCanvas* canv = new TCanvas("tmp","tmp",500,500); canv->cd();
  //static double sigma;


  int xBins = th2->GetNbinsX();
  //cout << "xBins: " << xBins << endl; 



  //for resolution vs pt
  if(variableBinning){
    hFit = new TH1F("tmp","tmp",xBins,th2->GetXaxis()->GetXbins()->GetArray());
    hRMS = new TH1F("tmp2","tmp2",xBins,th2->GetXaxis()->GetXbins()->GetArray());
    hFitMean = new TH1F("tmp3","tmp3",xBins,th2->GetXaxis()->GetXbins()->GetArray());
    hAverage = new TH1F("tmp4","tmp4",xBins,th2->GetXaxis()->GetXbins()->GetArray());
  }else{
    //for resolution vs eta
    hFit = new TH1F("tmp","tmp",xBins,th2->GetBinLowEdge(1),th2->GetBinLowEdge(xBins+1));
    hRMS = new TH1F("tmp2","tmp2",xBins,th2->GetBinLowEdge(1),th2->GetBinLowEdge(xBins+1));
    hFitMean = new TH1F("tmp3","tmp4",xBins,th2->GetBinLowEdge(1),th2->GetBinLowEdge(xBins+1));
    hAverage = new TH1F("tmp3","tmp4",xBins,th2->GetBinLowEdge(1),th2->GetBinLowEdge(xBins+1));
  }


  double previousMean,previousSigma,previousA,previousN,previousADx,previousNDx;

  for(int i=1;i<=xBins; i++){
    cout << "i,xLeft: " << i << " , " << th2->GetBinLowEdge(i) << endl;
    TH1D* proj = th2->ProjectionY("prova",i,i+1);
    double rms = proj->GetRMS();

    double tmpMean;
    double tmpSigma;


    // --------- new ---
    double firstRangeLeft;    
    double firstRangeRight;

    bool isPtPlot(false),isDzTheta(false),isD0(false);
    if(name.find("phires_vs_eta") != string::npos){
      firstRangeLeft = -0.005;    
      firstRangeRight = 0.005;
    }
    if(name.find("cotThetares_vs_eta") != string::npos){
      firstRangeLeft = -0.01;    
      firstRangeRight = 0.01;
      isDzTheta = true;
    }
    if(name.find("dxyres_vs_eta") != string::npos){
      firstRangeLeft = -0.05;    
      firstRangeRight = 0.05;
      isD0 = true;
    }
    if(name.find("dzres_vs_eta") != string::npos){
      firstRangeLeft = -0.1;    
      firstRangeRight = 0.1;
      isDzTheta = true;
    }
    if(name.find("ptres_vs_eta") != string::npos){
      firstRangeLeft = -1;    
      firstRangeRight = 2.0;
      isPtPlot = true;
    }
    TF1* f1 = new TF1("f1","gaus",firstRangeLeft,firstRangeRight);
    // ------------

    TFitResultPtr r = proj->Fit(f1,"SMRL Q");   
    proj->Draw(); gPad->Update(); 
    stringstream outputName;
    if(i<10)      outputName << name << "_raw_bin_" << "0" << i ; else      outputName << name << "_bin" << i ;
    //printCanvas(canv,outputName.str(),type,1);    
    
    if(wait) Wait();
    tmpMean = r->Parameter(1);
    tmpSigma = r->Parameter(2);    

    double xMin,xMax;
    xMin = tmpMean - tmpSigma*5.;
    xMax = tmpMean + tmpSigma*5.;



    RooRealVar x("x","x",xMin,xMax) ;

    double meanRangeMin;
    double meanRangeMax;
    
    cout << "=== tmpMean: " << tmpMean << endl;

    if(tmpMean<0){
      meanRangeMin = 5.0*tmpMean; //new
      meanRangeMax = -2*tmpMean;
      meanRangeMax = -5*tmpMean;
    }else{
      meanRangeMin = -2*tmpMean;
      meanRangeMin = -5*tmpMean;
      meanRangeMax = 5.0*tmpMean; //new
    }

    RooRealVar meanRoo("mean","mean of gaussian",tmpMean,meanRangeMin,meanRangeMax) ;
    if(name.find("ptres_vs_eta") != string::npos){ //new
      meanRoo.setRange(firstRangeLeft,firstRangeRight) ;
      x.setRange(firstRangeLeft,firstRangeRight) ;
    }

    RooRealVar sigmaRoo("sigma","width of gaussian",tmpSigma,tmpSigma*0.5,tmpSigma*1.5); 
    if(isPtPlot){
      sigmaRoo.setMin(tmpSigma*0.1);
      sigmaRoo.setMax(tmpSigma*5.0);
    }


    RooRealVar a("a","a",3.,0.3,10.); //new
    RooRealVar aDx("aDx","aDx",3.,1.,10.);
    RooRealVar n("n","n",5.,0.,10.);   
    RooRealVar nDx("nDx","nDx",5.,0.,10.);   

    if(tailSign!=0){
      if(tailSign<0){
	a.setVal(1); a.setMin(0.1); a.setMax(3.);
	x.setRange(firstRangeLeft,firstRangeRight/3.) ; //31/08/13
      }else{
	aDx.setVal(1); aDx.setMin(0.1); aDx.setMax(3.);
	x.setRange(firstRangeLeft/3.,firstRangeRight) ; //31/08/13
      }
    }

    if(i>1){
      meanRoo.setVal(previousMean);
      sigmaRoo.setVal(previousSigma);
      a.setVal(previousA);
      n.setVal(previousN);
      aDx.setVal(previousADx);
      nDx.setVal(previousNDx);
    }

    RooDoubleCB func1("cb","cb PDF",x,meanRoo,sigmaRoo,a,n,aDx,nDx) ;
    RooPlot* xframe = x.frame(Title("Gaussian p.d.f.")) ;
    
    //RooDataHist  dh("dh","dh",x,Import(*proj));

    // ---- new ----
    TH1D* subProj =  TH1Subset(proj,xMin,xMax);
    if(isPtPlot) subProj->Rebin(2);
    RooDataHist  dh("dh","dh",x,Import(*subProj));
    delete subProj;
    // ---------

    func1.fitTo(dh,NumCPU(2),RooFit::PrintLevel(0));
    dh.plotOn(xframe);
    func1.plotOn(xframe) ;
    //

    xframe->Draw(); gPad->Update();
    stringstream outputName2; 
    if(i<10)      outputName2 << name << "_bin" << "0" << i ; else      outputName2 << name << "_bin" << i ;
    //printCanvas(canv,outputName2.str(),type,1);    
    if(wait) Wait();
    
    tmpMean = meanRoo.getVal();
    tmpSigma = sigmaRoo.getVal();
 
    double sigma,sigmaErr;
    double mean,meanErr;
    sigma= sigmaRoo.getVal();    sigmaErr = sigmaRoo.getError();     
    mean= meanRoo.getVal(); meanErr = meanRoo.getError();


    //save values for next fit
    previousMean = meanRoo.getVal();
    previousSigma = sigmaRoo.getVal();
    previousA = a.getVal();
    previousN = n.getVal(); 
    previousADx = aDx.getVal();
    previousNDx = nDx.getVal();

    /*
    proj->SetAxisRange(-5*sigma,5*sigma);
    rms = proj->GetRMS();
    double rmsErr = proj->GetRMSError();
    hFit->SetBinContent(i,sigma);
    hFit->SetBinError(i,sigmaErr);
    hRMS->SetBinContent(i,rms);
    hRMS->SetBinError(i,rmsErr);
    delete proj;
    */

    // --- HERE I get the 68% and 95% range

    // --------- new ----------
    double fullAverage = proj->GetMean();
    double step = proj->GetBinWidth(1);
    int peakBin = proj->FindBin(tmpMean);
    double fullIntegral,fullIntegralSx,fullIntegralDx;
    fullIntegral = proj->Integral(0,proj->GetNbinsX()+1);
    fullIntegralSx = proj->Integral(0,peakBin);
    fullIntegralDx = proj->Integral(peakBin,proj->GetNbinsX()+1);


    //
    bool found68(false),found68Sx(false),found68Dx(false);
    bool found95(false),found95Sx(false),found95Dx(false);
    double range68(0),range68Sx(0),range68Dx(0), uncert68(0.);
    double range95(0),range95Sx(0),range95Dx(0), uncert95(0.);
    bool touchUnder(false),touchOver(false);
    double fraction(0),fractionSx(0),fractionDx(0);
    for(int j=0; j<proj->GetNbinsX()/2 ; j++){
      if((peakBin-j)<1 && !found95Sx) {
	if(!touchUnder) {
	  cout << "WARNING: get to under-flow" << endl;
	  touchUnder=true;
	  if(!found68Sx) {range68Sx = step*(2*j+1)*0.5; found68Sx=true;}
	  if(!found95Sx) {range95Sx = step*(2*j+1)*0.5; found95Sx=true;}
	}
	if(touchUnder && touchOver) break;
      }
      if( (peakBin+j) > proj->GetNbinsX() && !found95Dx) {
	if(!touchOver){
	  cout << "WARNING: get to over-flow" << endl; 
	  touchOver = true;
	  if(!found68Dx) {range68Dx = step*(2*j+1)*0.5; found68Dx=true;}
	  if(!found95Dx) {range95Dx = step*(2*j+1)*0.5; found95Dx=true;}
	}
	if(touchUnder && touchOver) break;
      }
      if((peakBin-j)<1 && (peakBin+j) > proj->GetNbinsX()) {
	if(!found68) range68 = step*(2*j+1)*0.5;
	if(!found95) range95 = step*(2*j+1)*0.5;
	break;
      }


      fraction= proj->Integral(peakBin-j,peakBin+j)/fullIntegral;
      fractionSx= proj->Integral(peakBin-j,peakBin)/fullIntegralSx;
      fractionDx= proj->Integral(peakBin,peakBin+j)/fullIntegralDx;    

      //double range2nd(0.954);
      //double range2nd(0.80);
      //double range2nd(0.85);
      double range2nd(0.90);

      if(fraction>0.682 && !found68){ //1sigma range
	found68=true;
	range68 = step*(2*j+1)*0.5;
	double averageBinContent = (proj->GetBinContent(peakBin-j) + proj->GetBinContent(peakBin+j))/2.0 ; 
	uncert68 = sqrt(0.682*(1-0.682)/fullIntegral)/(averageBinContent/step/fullIntegral);
      }
      if(fractionSx>0.682 && !found68Sx){ //1sigma range
	found68Sx=true;
	range68Sx = step*(2*j+1)*0.5;
	double averageBinContent = proj->GetBinContent(peakBin-j) ; 
	uncert68 = sqrt(0.682*(1-0.682)/fullIntegralSx)/(averageBinContent/step/fullIntegralSx);
      }
      if(fractionDx>0.682 && !found68Dx){ //1sigma range
	found68Dx=true;
	range68Dx = step*(2*j+1)*0.5;
	double averageBinContent = proj->GetBinContent(peakBin+j) ; 
	uncert68 = sqrt(0.682*(1-0.682)/fullIntegralSx)/(averageBinContent/step/fullIntegralSx);
      }

      if(fraction>range2nd && !found95){ //2sigma range
	found95=true;
	range95 = step*(2*j+1)*0.5;
	double averageBinContent = (proj->GetBinContent(peakBin-j) + proj->GetBinContent(peakBin+j))/2.0 ; 
	uncert95 = sqrt(range2nd*(1-range2nd)/fullIntegral)/(averageBinContent/step/fullIntegral);	
      }
      if(fractionSx>range2nd && !found95Sx){ //2sigma range
	found95Sx=true;
	range95Sx = step*(2*j+1)*0.5;
	double averageBinContent = proj->GetBinContent(peakBin-j) ; 
	uncert95 = sqrt(range2nd*(1-range2nd)/fullIntegralSx)/(averageBinContent/step/fullIntegralSx);	
      }
      if(fractionDx>range2nd && !found95Dx){ //2sigma range
	found95Dx=true;
	range95Dx = step*(2*j+1)*0.5;
	double averageBinContent = proj->GetBinContent(peakBin+j) ; 
	uncert95 = sqrt(range2nd*(1-range2nd)/fullIntegralDx)/(averageBinContent/step/fullIntegralDx);	
      }

      if(found68 && found68Sx && found68Dx && found95 && found95Sx && found95Dx) break;     

    }
    // --------- end new ----------



    // ------- new --------
    cout << "range68,range68sx,range68dx: " 
	 << range68 << " , " 
	 << range68Sx << " , " 
	 << range68Dx << endl;

    cout << "range95,range95sx,range95dx: " 
	 << range95 << " , " 
	 << range95Sx << " , " 
	 << range95Dx << endl;


    //68 plus 95
    if(rangeType == +1){
      //---
      if(isDzTheta){ 
	if(i<=25){
	  hFit->SetBinContent(i,range68Dx);
	  hRMS->SetBinContent(i,range95Dx);
	}else{
	  hFit->SetBinContent(i,range68Sx);
	  hRMS->SetBinContent(i,range95Sx);
	}
      }else      if(isD0){
	hFit->SetBinContent(i,range68Sx);
	hRMS->SetBinContent(i,range95Sx);
      }else{
	hFit->SetBinContent(i,range68Dx);
	hRMS->SetBinContent(i,range95Dx);
      }
      //---

      hFit->SetBinError(i,1/1000000.); 
      hRMS->SetBinError(i,1/1000000.);
    }else if(rangeType == -1){
      //---
      if(isDzTheta){ 
	if(i<=25){
	  hFit->SetBinContent(i,range68Sx);
	  hRMS->SetBinContent(i,range95Sx);
	}else{
	  hFit->SetBinContent(i,range68Dx);
	  hRMS->SetBinContent(i,range95Dx);
	}
      }else      if(isD0){
	hFit->SetBinContent(i,range68Dx);
	hRMS->SetBinContent(i,range95Dx);
      }else{
	hFit->SetBinContent(i,range68Sx);
	hRMS->SetBinContent(i,range95Sx);
      }
      //---
      hFit->SetBinError(i,1/1000000.);
      hRMS->SetBinError(i,1/1000000.);
    }else{
      hFit->SetBinContent(i,range68);
      hFit->SetBinError(i,1/1000000.);
      hRMS->SetBinContent(i,range95);
      hRMS->SetBinError(i,1/1000000.);
    }
    

    //Set real bin uncertainty
    hFit->SetBinError(i,uncert68); 
    hRMS->SetBinError(i,uncert95);

    hFitMean->SetBinContent(i,mean);
    hFitMean->SetBinError(i,meanErr);
    hAverage->SetBinContent(i,fullAverage);
    hAverage->SetBinError(i,0.00001);
    delete proj;
    // ------- end new --------



  }
  hFit->SetDirectory(gROOT);
  hRMS->SetDirectory(gROOT);
  delete canv;
  //hFit->Draw();gPad->Update(); Wait();
}


TH1D * TH1Subset(TH1D *h,double lowEdge,double highEdge)
{
  int initialNBins = h->GetNbinsX();
  //cout << "initial nbins: " << initialNBins << endl;
  double hmin = h->GetBinLowEdge(1);
  double hmax = h->GetBinLowEdge(initialNBins);
  int binLow  = h->FindBin(max(lowEdge,hmin));
  int binHigh = h->FindBin(min(highEdge,hmax));
  int nBins = binHigh - binLow +1;
  //cout << "lowEdge,highEdge: " << lowEdge << " , " << highEdge << endl;
  //cout << "binLow,binHigh,nBins: " << binLow << " , " << binHigh << " , " << nBins << endl;
  double newLow = h->GetBinLowEdge(binLow);
  double newHigh = h->GetBinLowEdge(binHigh+1);

  TH1D *newH = new TH1D("subSetH","subSetH",nBins, newLow,newHigh);
       
  Double_t content;
  
  for (Int_t i=0; i<nBins; i++) {
      content = h->GetBinContent(binLow+i);
      newH->SetBinContent(i+1, content);
  }
  
  return newH;
}
