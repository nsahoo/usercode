#include "TPaveStats.h"
#include "TF1.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
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


using namespace std;
using namespace RooFit;
#include <sstream>


void projectAndFit(TH2F* th2,TH1F*& hFit,TH1F*& hRMS,int tailSign=-1,string name="", bool variableBinning=false, bool wait=false);


void setHistoLabels(TH1F* h,float ySize, float yOffset,float xSize, float xOffset,
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
		TH1F* h1,TH1F* h2,TH1F* h3,
		TString label1, TString label2, TString label3){
  l->SetTextSize(0.04);
  l->SetLineColor(0);
  l->SetFillColor(0);
  l->AddEntry(h1,label1,"P");
  l->AddEntry(h2,label2,"P");
  l->AddEntry(h3,label3,"P");
  l->SetShadowColor(0);
}

void setHistoStyle(TH1F* h1,TH1F* h2,TH1F* h3,bool empty=false){
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

void printCanvas(TCanvas* canvas,string outputName,int format=3){

  cmsLabel(0);

  string outputFolder;
  outputFolder="output/ttbarSplit/";

  if(format>=1) canvas->Print((outputFolder.c_str()+outputName+".png").c_str());
  if(format>=2) canvas->Print((outputFolder.c_str()+outputName+".pdf").c_str());
  if(format>=3) canvas->Print((outputFolder.c_str()+outputName+".eps").c_str());
}


int Wait() {
     cout << " Continue [<RET>|q]?  "; 
     char x;
     x = getchar();
     if ((x == 'q') || (x == 'Q')) return 1;
     return 0;
}


void makePlotsTTbarSplit()
{ 
  /*
    input:
    0 --> muons
    1 --> pions
    2 --> electrons
   */

  //TString fileName1;
  //TString fileName2;
  //TString fileName3;
  TString fileName4;
  TString fileName5;
  TString fileName6;
  //TString fileName7;
  //TString fileName8;
  //TString fileName9;


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


  // ------ for pions and TTbar ---------
  fileName4="input/ttbar/puhp_barrel.root";
  fileName5="input//ttbar/puhp_trans.root";
  fileName6="input//ttbar/puhp_endcap.root";



  label1="#pi, p_{T} = 1 GeV";
  label2="#pi, p_{T} = 10 GeV";
  label3="#pi, p_{T} = 100 GeV";
  label4="#pi, Barrel region";
  label5="#pi, Transition region"; 
  label6="#pi, Endcap region"; 

  label1=" p_{T} = 1 GeV";
  label2="p_{T} = 10 GeV";
  label3="p_{T} = 100 GeV";
  label4="Barrel region";
  label5="Transition region"; 
  label6="Endcap region"; 


  resDPhiVsEtaMax=30;resDPhiVsEtaMin=0.02;
  resD0vsEtaMax=1000;resD0vsEtaMin=2;
  resDzVsEtaMax=2000;resDzVsEtaMin=5;
  resDThetaVsEtaMax=50;resDThetaVsEtaMin=0.1;
  resPtVsEtaMax=30;resPtVsEtaMin=0.5; //for mu and pi
  
  //in microns and 10-3 units
  resDPhiVsPtMax=200;resDPhiVsPtMin=0.05;
  resD0vsPtMax=5000;resD0vsPtMin=5;
  resDzVsPtMax=8000;resDzVsPtMin=20;
  resDThetaVsPtMax=500;resDThetaVsPtMin=0.2;
  resPtVsPtMax=100;resPtVsPtMin=0.5; //for mu and pi
    

  tailSignPhi =0;
  tailSignTheta =0;
  tailSignD0 =0;
  tailSignDz =0;
  tailSignPt =0;
  
  tailSignPhi_pT =0;
  tailSignTheta_pT =0;
  tailSignD0_pT =0;
  tailSignDz_pT =0;
  tailSignPt_pT =0;

  //bool eff=0;
  //bool resVsEta=0;
  bool resVsPt=1;


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

 if(resVsPt){
   //===== resolutions vs pt
   TFile * file4 = new TFile(fileName4); TDirectory * dir1=file4->GetDirectory("");
   TFile * file5 = new TFile(fileName5); TDirectory * dir2=file5->GetDirectory("");
   TFile * file6 = new TFile(fileName6); TDirectory * dir3=file6->GetDirectory("");
   
   TString collname4="DQMData/Tracking/Track/cutsRecoHp_AssociatorByHits";  
   TString collname5="DQMData/Tracking/Track/cutsRecoHp_AssociatorByHits";  
   TString collname6="DQMData/Tracking/Track/cutsRecoHp_AssociatorByHits";  
   
   
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
   
   // --- 
   TH2F* tmpH2;

   string name;
   name="phires_vs_pt_barrel";
   dir1->GetObject(collname4+"/phires_vs_pt",tmpH2);  projectAndFit(tmpH2,f1h1,f1h1RMS,tailSignPhi_pT,name,true);
   name="phires_vs_pt_transition";
   dir2->GetObject(collname5+"/phires_vs_pt",tmpH2);  projectAndFit(tmpH2,f2h1,f2h1RMS,tailSignPhi_pT,name,true);
   name="phires_vs_pt_endcap";
   dir3->GetObject(collname6+"/phires_vs_pt",tmpH2);  projectAndFit(tmpH2,f3h1,f3h1RMS,tailSignPhi_pT,name,true);

   name="cotThetares_vs_pt_barrel";
   dir1->GetObject(collname4+"/cotThetares_vs_pt",tmpH2);  projectAndFit(tmpH2,f1h2,f1h2RMS,tailSignTheta_pT,name,true); 
   name="cotThetares_vs_pt_transition";
   dir2->GetObject(collname5+"/cotThetares_vs_pt",tmpH2);  projectAndFit(tmpH2,f2h2,f2h2RMS,tailSignTheta_pT,name,true);  
   name="cotThetares_vs_pt_endcap";
   dir3->GetObject(collname6+"/cotThetares_vs_pt",tmpH2);  projectAndFit(tmpH2,f3h2,f3h2RMS,tailSignTheta_pT,name,true);  

   name="dxyres_vs_pt_barrel";
   dir1->GetObject(collname4+"/dxyres_vs_pt",tmpH2);   projectAndFit(tmpH2,f1h3,f1h3RMS,tailSignD0_pT,name,true); 
   name="dxyres_vs_pt_transition";
   dir2->GetObject(collname5+"/dxyres_vs_pt",tmpH2);   projectAndFit(tmpH2,f2h3,f2h3RMS,tailSignD0_pT,name,true); 
   name="dxyres_vs_pt_endcap";
   dir3->GetObject(collname6+"/dxyres_vs_pt",tmpH2);   projectAndFit(tmpH2,f3h3,f3h3RMS,tailSignD0_pT,name,true); 

   name="dzres_vs_pt_barrel";
   dir1->GetObject(collname4+"/dzres_vs_pt",tmpH2);   projectAndFit(tmpH2,f1h4,f1h4RMS,tailSignDz_pT,name,true); 
   name="dzres_vs_pt_transition";
   dir2->GetObject(collname5+"/dzres_vs_pt",tmpH2);   projectAndFit(tmpH2,f2h4,f2h4RMS,tailSignDz_pT,name,true); 
   name="dzres_vs_pt_endcap";
   dir3->GetObject(collname6+"/dzres_vs_pt",tmpH2);   projectAndFit(tmpH2,f3h4,f3h4RMS,tailSignDz_pT,name,true); 

   name="ptres_vs_pt_barrel";
   dir1->GetObject(collname4+"/ptres_vs_pt",tmpH2);   projectAndFit(tmpH2,f1h5,f1h5RMS,tailSignPt_pT,name,true); 
   name="ptres_vs_pt_transition";
   dir2->GetObject(collname5+"/ptres_vs_pt",tmpH2);   projectAndFit(tmpH2,f2h5,f2h5RMS,tailSignPt_pT,name,true); 
   name="ptres_vs_pt_endcap";
   dir3->GetObject(collname6+"/ptres_vs_pt",tmpH2);   projectAndFit(tmpH2,f3h5,f3h5RMS,tailSignPt_pT,name,true); 
   // ---


   f1h1->Scale(1000.); f2h1->Scale(1000.); f3h1->Scale(1000.); 
   f1h1RMS->Scale(1000.); f2h1RMS->Scale(1000.); f3h1RMS->Scale(1000.); 
   f1h1->GetYaxis()->SetRangeUser(resDPhiVsPtMin,resDPhiVsPtMax);
   f1h1->GetXaxis()->SetRangeUser(0.2,100);
   setHistoLabels(f1h1,0.05,1.2,0.055,1.1,"Resolution in #phi (10^{-3}radians)","p_{T} (GeV)");


   f1h2->Scale(1000.); f2h2->Scale(1000.); f3h2->Scale(1000.); 
   f1h2RMS->Scale(1000.); f2h2RMS->Scale(1000.); f3h2RMS->Scale(1000.); 
   f1h2->GetYaxis()->SetRangeUser(resDThetaVsPtMin,resDThetaVsPtMax);
   f1h2->GetXaxis()->SetRangeUser(0.2,100);
   setHistoLabels(f1h2,0.05,1.2,0.055,1.1,"Resolution in cot(#theta) (10^{-3})","p_{T} (GeV)");


   f1h3->Scale(10000.); f2h3->Scale(10000.); f3h3->Scale(10000.); 
   f1h3RMS->Scale(10000.); f2h3RMS->Scale(10000.); f3h3RMS->Scale(10000.); 
   f1h3->GetYaxis()->SetRangeUser(resD0vsPtMin,resD0vsPtMax);  
   f1h3->GetXaxis()->SetRangeUser(0.2,100);
   setHistoLabels(f1h3,0.05,1.2,0.055,1.1,"Resolution in d_{0} (#mum)","p_{T} (GeV)");


   f1h4->Scale(10000.); f2h4->Scale(10000.); f3h4->Scale(10000.); 
   f1h4RMS->Scale(10000.); f2h4RMS->Scale(10000.); f3h4RMS->Scale(10000.); 
   f1h4->GetYaxis()->SetRangeUser(resDzVsPtMin,resDzVsPtMax);
   f1h4->GetXaxis()->SetRangeUser(0.2,100);
   setHistoLabels(f1h4,0.05,1.2,0.055,1.1,"Resolution in z_{0} (#mum)","p_{T} (GeV)");


   f1h5->Scale(100.); f2h5->Scale(100.); f3h5->Scale(100.); 
   f1h5RMS->Scale(100.); f2h5RMS->Scale(100.); f3h5RMS->Scale(100.); 
   f1h5->GetYaxis()->SetRangeUser(resPtVsPtMin,resPtVsPtMax);
   f1h5->GetXaxis()->SetRangeUser(0.2,100);
   setHistoLabels(f1h5,0.05,1.2,0.055,1.1,"(Resolution in p_{T})/p_{T} (%)","p_{T} (GeV)");



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
   
   TLegend* l = new TLegend(0.40,0.67,0.83,0.92);
   setLegend3(l,f1histos3[0],f2histos3[0],f3histos3[0],label4,label5,label6);
   l->Draw();
   output="resolutionPhiVsPt"; printCanvas(canvas,output);
   delete l;


   // -----
   canvas->SetLogx(1);
   canvas->SetLogy(1);
   f1histos3[1]->Draw();         f1histos3RMS[1]->Draw("sames");
   f2histos3[1]->Draw("sames");  f2histos3RMS[1]->Draw("sames");
   f3histos3[1]->Draw("sames");  f3histos3RMS[1]->Draw("sames");

   l = new TLegend(0.42,0.67,0.85,0.92);
   setLegend3(l,f1histos3[1],f2histos3[1],f3histos3[1],label4,label5,label6);
   l->Draw();
   output="resolutionThetaVsPt"; printCanvas(canvas,output);
   delete l;


   // -----
   canvas->SetLogx(1);
   canvas->SetLogy(1);

   f1histos3[2]->Draw();         f1histos3RMS[2]->Draw("sames");
   f2histos3[2]->Draw("sames");  f2histos3RMS[2]->Draw("sames");
   f3histos3[2]->Draw("sames");  f3histos3RMS[2]->Draw("sames");

   l = new TLegend(0.40,0.67,0.83,0.92);
   setLegend3(l,f1histos3[2],f2histos3[2],f3histos3[2],label4,label5,label6);
   l->Draw();
   output="resolutionD0VsPt"; printCanvas(canvas,output);
   delete l;


   // -----
   canvas->SetLogx(1);
   canvas->SetLogy(1);

   f1histos3[3]->Draw();         f1histos3RMS[3]->Draw("sames");
   f2histos3[3]->Draw("sames");  f2histos3RMS[3]->Draw("sames");
   f3histos3[3]->Draw("sames");  f3histos3RMS[3]->Draw("sames");
   
   l = new TLegend(0.44,0.67,0.86,0.92);
   setLegend3(l,f1histos3[3],f2histos3[3],f3histos3[3],label4,label5,label6);
   l->Draw();
   output="resolutionDzVsPt"; printCanvas(canvas,output);
   delete l;

   // -----
   canvas->SetLogx(1);
   canvas->SetLogy(1);
   f1histos3[4]->Draw();         f1histos3RMS[4]->Draw("sames");
   f2histos3[4]->Draw("sames");  f2histos3RMS[4]->Draw("sames");
   f3histos3[4]->Draw("sames");  f3histos3RMS[4]->Draw("sames");
   
   l = new TLegend(0.20,0.67,0.62,0.92);
   setLegend3(l,f1histos3[4],f2histos3[4],f3histos3[4],label4,label5,label6);
   l->Draw();
   output="resolutionPtVsPt"; printCanvas(canvas,output);
   delete l;
 }


 //delete canvas;
}



// ----------------------------------
//======= Doing real fit and producing Fit/RMS resolution histograms =======
void projectAndFit(TH2F* th2, TH1F*& hFit, TH1F*& hRMS,int tailSign,string name, bool variableBinning,bool wait){
  TCanvas* canv = new TCanvas("tmp","tmp",500,500); canv->cd();
  //static double sigma;


  int xBins = th2->GetNbinsX();
  //cout << "xBins: " << xBins << endl; 



  //for resolution vs pt
  if(variableBinning){
    hFit = new TH1F("tmp","tmp",xBins,th2->GetXaxis()->GetXbins()->GetArray());
    hRMS = new TH1F("tmp2","tmp2",xBins,th2->GetXaxis()->GetXbins()->GetArray());
  }else{
    //for resolution vs eta
    hFit = new TH1F("tmp","tmp",xBins,th2->GetBinLowEdge(1),th2->GetBinLowEdge(xBins+1));
    hRMS = new TH1F("tmp2","tmp2",xBins,th2->GetBinLowEdge(1),th2->GetBinLowEdge(xBins+1));
  }


  for(int i=1;i<=xBins; i++){
    cout << "i,xLeft: " << i << " , " << th2->GetBinLowEdge(i) << endl;
    TH1D* proj = th2->ProjectionY("prova",i,i+1);
    if(proj->GetEntries()<20){
      //cout << "N entries < 20. Should I continue? " << endl; Wait();
      hRMS->SetBinContent(i,0);
      hRMS->SetBinError(i,0);
      hFit->SetBinContent(i,0);
      hFit->SetBinError(i,0);
      delete proj;
      continue;
    }
    double rms = proj->GetRMS();

    double tmpMean;
    double tmpSigma;


    // quick fit with standard gauss
    double leftBound,rightBound;
    if(tailSign<0){
      leftBound = -1.5*rms;
      rightBound = +1.0*rms;
      cout << "negative tailSign gives: " << leftBound << " , " << rightBound << endl;
    } else if(tailSign>0){
      leftBound = -1.0*rms;
      rightBound = +1.5*rms;
      cout << "positive tailSign gives: " << leftBound << " , " << rightBound << endl;
    }else{
      leftBound = -1.0*rms;
      rightBound = +1.0*rms;
      cout << "null tailSign gives: " << leftBound << " , " << rightBound << endl;
    }
    
    TF1* f1 = new TF1("f1","gaus",leftBound,rightBound);
    TFitResultPtr r = proj->Fit(f1,"S M R L");   
    proj->Draw(); gPad->Update(); 
    stringstream outputName;
    if(i<10)      outputName << name << "_raw_bin_" << "0" << i ; else      outputName << name << "_bin" << i ;
    //printCanvas(canv,outputName.str(),1);    
    
    if(wait) Wait();
    tmpMean = r->Parameter(1);
    tmpSigma = r->Parameter(2);    

    double xMin,xMax;
    xMin = tmpMean - tmpSigma*5.;
    xMax = tmpMean + tmpSigma*5.;



    RooRealVar x("x","x",xMin,xMax) ;

    double meanRangeMin;
    double meanRangeMax;
    if(tmpMean<0){
      meanRangeMin = 1.6*tmpMean;
      meanRangeMax = 0.1*tmpMean;
    }else{
      meanRangeMin = 0.1*tmpMean;
      meanRangeMax = 1.6*tmpMean;
    }

    RooRealVar meanRoo("mean","mean of gaussian",tmpMean,meanRangeMin,meanRangeMax) ;
    RooRealVar sigmaRoo("sigma","width of gaussian",tmpSigma,tmpSigma*0.5,tmpSigma*1.5); 

    RooRealVar a("a","a",3.,2.,10.);
    RooRealVar aDx("aDx","aDx",3.,2.,10.);
    RooRealVar n("n","n",5.,0.,10.);   
    RooRealVar nDx("nDx","nDx",5.,0.,10.);   

    if(tailSign!=0){
      if(tailSign<0){
	a.setVal(1); a.setMin(0.5); a.setMax(3.);
      }else{
	aDx.setVal(1); aDx.setMin(0.5); aDx.setMax(3.);
      }
    }


    RooDoubleCB func1("cb","cb PDF",x,meanRoo,sigmaRoo,a,n,aDx,nDx) ;
    RooPlot* xframe = x.frame(Title("Gaussian p.d.f.")) ;
    
    RooDataHist  dh("dh","dh",x,Import(*proj));
    func1.fitTo(dh);
    dh.plotOn(xframe);
    func1.plotOn(xframe) ;
    //

    xframe->Draw(); gPad->Update();
    stringstream outputName2; 
    if(i<10)      outputName2 << name << "_bin" << "0" << i ; else      outputName2 << name << "_bin" << i ;
    //printCanvas(canv,outputName2.str(),1);    
    if(wait) Wait();
    
    tmpMean = meanRoo.getVal();
    tmpSigma = sigmaRoo.getVal();
 
    double sigma;
    double sigmaErr;
    sigma= sigmaRoo.getVal();
    sigmaErr = sigmaRoo.getError();     

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


    // ---- NEW IMPLEMENTATION
    double fullIntegral = proj->Integral(0,proj->GetNbinsX()+1);
    //double fullIntegral = proj->Integral();
    //double fullAverage = proj->GetMean();
    double step = proj->GetBinWidth(1);;
    int peakBin = proj->FindBin(tmpMean);

    //double range2nd(0.954);
    double range2nd(0.90);


    //
    bool found68(false);
    double range68(0), uncert68(0.);
    double range95(0), uncert95(0.);
    for(int j=0; j<proj->GetNbinsX()/2 ; j++){
      if((peakBin-j)<1 || (peakBin+j) > proj->GetNbinsX()) break;
      double fraction = proj->Integral(peakBin-j,peakBin+j)/fullIntegral;
      if(fraction>0.682 && !found68){ //2sigma range
	found68=true;
	range68 = step*(2*j+1)*0.5;
	double averageBinContent = (proj->GetBinContent(peakBin-j) + proj->GetBinContent(peakBin+j))/2.0 ; 
	uncert68 = sqrt(0.682*(1-0.682)/fullIntegral)/(averageBinContent/step/fullIntegral);
      }
      if(fraction>range2nd){ //3sigma range
	range95 = step*(2*j+1)*0.5;
	double averageBinContent = (proj->GetBinContent(peakBin-j) + proj->GetBinContent(peakBin+j))/2.0 ; 
	uncert95 = sqrt(range2nd*(1-range2nd)/fullIntegral)/(averageBinContent/step/fullIntegral);
	break;
      }
    }


    // --- OLD IMPLEMENTATION
    /*
    double fullIntegral = proj->Integral();
    double step = sigma/20.;
    //
    double range68(0);
    for(int j=2; j<20000;j++){
      xMin = tmpMean - step*j; 
      xMax = tmpMean + step*j;  
      proj->SetAxisRange(xMin,xMax);
      double fraction = proj->Integral()/fullIntegral;
      if(fraction>0.682){ //2sigma range
	range68 = step*j;
	break;
      }
    }
    //
    double range95(0);
    for(int j=2; j<20000;j++){
      xMin = tmpMean - step*j; 
      xMax = tmpMean + step*j;  
      proj->SetAxisRange(xMin,xMax);
      double fraction = proj->Integral()/fullIntegral;
      if(fraction>range2nd){ //2sigma range
	range95 = step*j;
	break;
      }
    }
    uncert68 = step;
    uncert95 = step;
    */




    //--------
    //hFit->SetBinContent(i,sigma);
    //hFit->SetBinError(i,sigmaErr);
    hFit->SetBinContent(i,range68);
    hFit->SetBinError(i,uncert68);
    hRMS->SetBinContent(i,range95);
    hRMS->SetBinError(i,uncert95);
    delete proj;


    
  }
  hFit->SetDirectory(gROOT);
  hRMS->SetDirectory(gROOT);
  delete canv;
  //hFit->Draw();gPad->Update(); Wait();
}

