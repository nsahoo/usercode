#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"

#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"

#include <sstream>


using namespace RooFit ;

TString xLabel;
TString yLabel;
TString yLabelResp;


bool wantMore2() {
  // ask if user wants more
  fprintf(stderr,"Type <CR> to continue or q to quit ==> ");
  // read first char
  int readch = getchar(), answer = readch;
  // poll out remaining chars from buffer
  while (readch != '\n' && readch != EOF) readch = getchar();
  // check first char
  return !(answer == 'q' || answer == 'Q');  
}


double getFitRange(TH1* h);
void getResponse(TH1* h_reso,TH1* h_resp,TH1* h_bef,TH1* h_befData,
		 double& reso, double& resp1, double& resp2, double& frac,
		 double& fit,double& dataFit,
		 double& resoErr, double& fitErr, double& dataFitErr, 
		 double input1,double input2,double input3,double rangeMax,
		 TCanvas* canvas1,int typePlot=0);

void getResponse2(double resp1, //estimated resp1 from MC (previous step)
		  double resp2, //estimated resp1 from MC (previous step)
		  double& frac, //do I need this??
		  TH1* h_resp,TH1* h_bef,TH1* h_befData,     	  
		  double& fit,double& dataFit,
		  double& fitErr, double& dataFitErr, 
		  double input1,
		  double input2,
		  double input3,
		  double rangeMax,
		  TCanvas* canvas1,
		  int typePlot=0);


void plotHistos(TH1F* plotReso, TH1F* plotFit, TH1F* plotDataFit,TH1F* plotResp1,TH1F* plotResp2,
		double* resolutions, double* fits, double* dataFits,
		double* responses1,double* responses2,
		double* resolutionsErr, double* fitsErr, double* dataFitsErr,
		TString yLabel, TString xLabel, double hRange, double lRange,
		TCanvas* canvas2,int type=0);


void ipUnfold(int type=1,int typePlot=0,int typePlot2=1,int typePlot3=0)
{
  //type(1-6): d0vsPt,d0vsEta,d0vsPhi, dzVsPt,dzVsEta,dzVsPhi
  //typePlot1(1-4): MC-resolutions(1), responses(2), measured resolutions sim(3), meas. res data (4) 
  //typePlot2(0-3): nothing, MC-resolutions, responses, measured resolutions 


  //=========  settings ====================
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  //gStyle->SetPadGridX(kTRUE);
  //gStyle->SetPadGridY(kTRUE);
  gStyle->SetOptStat("");
  gStyle->SetOptFit(0000000000);

  TString fnResp;
  TString fnReso;
  TString fnRawSim;
  TString fnRawData;

  double lRange,hRange;

  TString text1a,text1b,leg1a,leg1b;
  TString text2a,text2b,leg2a,leg2b;
  
  TString outputName1,outputName2;

  const unsigned int nbins_const = 60;
  //const unsigned int nbins_const = 50;
  unsigned int nbins;
  double xlow,xhigh;
  double fitRangeMax(0),initialResp1(0),initialResp2(0),initialFrac(0);


  TString hname1, hname2, hname3,hname4;

  TString rootFolderNameRaw="./raw/";
  TString rootFolderNameResp="./resp/";
  
  if(type==1 || type==4){//d0 and dz vs pt
    initialResp1 = 20.;
    initialResp2 = 20.;
    initialFrac = 0.6;
  }

  if(type==2){//d0 vs eta
    initialResp1 = 20.;
    initialResp2 = 20.;
    initialFrac = 0.6;
  }

  if(type==3 || type==6){//d0/dz vs phi
    //STILL TO BE DEFINED
  }

  if(type==5){//dz vs eta
    initialResp1 = 400.;
    initialResp2 = 200.;
    initialFrac = 0.65;
  }

  switch ( type )
    {
    case 1: //d0 vs pt
      hname1 ="sim_reso_dxyReso_vsPt_n"; //resolution
      hname2 ="sim_resp_dxyResp_vsPt_n"; //response 
      hname3 ="sim_raw_d0_vsPt_n"; //raw MC
      hname4 ="data_raw_d0_vsPt_n"; //raw DATA
      fitRangeMax = 400;
      break;

    case 2: //d0 vs eta
      hname1 ="sim_reso_dxyReso_vsEta_n"; //resolution
      hname2 ="sim_resp_dxyResp_vsEta_n"; //response 
      hname3 ="sim_raw_d0_vsEta_n"; //raw MC
      hname4 ="data_raw_d0_vsEta_n"; //raw DATA
      fitRangeMax = 600;
      break;

    case 3: //d0 vs phi
      hname1 ="sim_reso_dxyReso_vsPhi_n"; //resolution
      hname2 ="sim_resp_dxyResp_vsPhi_n"; //response 
      hname3 ="sim_raw_d0_vsPhi_n"; //raw MC
      hname4 ="data_raw_d0_vsPhi_n"; //raw DATA
      fitRangeMax = 400;
      break;

    case 4: //dz vs pt
      hname1 ="sim_reso_dzReso_vsPt_n"; //resolution
      hname2 ="sim_resp_dzResp_vsPt_n"; //response 
      hname3 ="sim_raw_dz_vsPt_n"; //raw MC
      hname4 ="data_raw_dz_vsPt_n"; //raw DATA
      fitRangeMax = 1500;
      break;

    case 5: //dz vs eta
      hname1 ="sim_reso_dzReso_vsEta_n"; //resolution
      hname2 ="sim_resp_dzResp_vsEta_n"; //response 
      hname3 ="sim_raw_dz_vsEta_n"; //raw MC
      hname4 ="data_raw_dz_vsEta_n"; //raw DATA
      fitRangeMax = 1500;
      break;

    case 6: //dz vs phi
      hname1 ="sim_reso_dzReso_vsPhi_n"; //resolution
      hname2 ="sim_resp_dzResp_vsPhi_n"; //response 
      hname3 ="sim_raw_dz_vsPhi_n"; //raw MC
      hname4 ="data_raw_dz_vsPhi_n"; //raw DATA
      fitRangeMax = 1500;
      break;

      
    default:
      break;
    }
  

  if(type==1){
    fnRawSim  = rootFolderNameRaw+"sim.raw.d0.vsPt.60bins.root";
    fnRawData = rootFolderNameRaw+"data.raw.d0.vsPt.60bins.root";
    fnResp    = rootFolderNameResp+"sim.resp.dxyResp.vsPt.60bins.root";
    fnReso    = rootFolderNameResp+"sim.reso.dxyReso.vsPt.60bins.root";

    xLabel = "Track p_{T} (GeV/c)";
    yLabel = "Transv.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the TIP due to Vertex (#mum)";
    
    lRange = 0.;
    hRange = 210.;

    text1a = "CMS Preliminary 2009";
    text1b = "900 GeV Runs";
    leg1a  = "DATA";
    leg1b  = "MC";

    text2a = "CMS Preliminary 2009";
    text2b = "";
    leg2a  = "MC (using MC-truth)";
    leg2b  = "MC (using only reco)";

    outputName1 = "resoD0_vs_pt.DATA.png";
    outputName2 = "resoD0_vs_pt.MC.png";

    nbins = 60;
    xlow=0.7;
    xhigh=2.2;
  }

  if(type==2){
    fnRawSim = rootFolderNameRaw+"sim.raw.d0.vsEta.50bins.root";
    fnRawData = rootFolderNameRaw+"data.raw.d0.vsEta.50bins.root";
    fnResp = rootFolderNameResp+"sim.resp.dxyResp.vsEta.50bins.root";
    fnReso    = rootFolderNameResp+"sim.reso.dxyReso.vsEta.50bins.root";

    xLabel = "Track #eta";
    yLabel = "Transv.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the TIP due to Vertex (#mum)";

    lRange = 0.;
    hRange = 250.;
    text1a = "CMS Preliminary 2009";
    text1b = "900 GeV Runs";
    leg1a  = "DATA";
    leg1b  = "MC";

    text2a = "CMS Preliminary 2009";
    text2b = "";
    leg2a  = "MC (using MC-truth)";
    leg2b  = "MC (using only reco)";

    outputName1 = "resoD0_vs_eta.Pt08.DATA.png";
    outputName2 = "resoD0_vs_eta.Pt08.MC.png";

    nbins = 50;
    xlow= -2.5;
    xhigh= 2.5;
  }

  if(type==3){
    //STILL TO BE DEFINED
  }

  if(type==4){
    fnRawSim  = rootFolderNameRaw+"sim.raw.dz.vsPt.60bins.root";
    fnRawData = rootFolderNameRaw+"data.raw.dz.vsPt.60bins.root";
    fnResp    = rootFolderNameResp+"sim.resp.dzResp.vsPt.60bins.root";
    fnReso    = rootFolderNameResp+"sim.reso.dzReso.vsPt.60bins.root";

    xLabel = "Track p_{T} (GeV/c)";
    yLabel = "Longit.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the LIP due to Vertex (#mum)";

    lRange = 0.;
    hRange = 210.;
    text1a = "CMS Preliminary 2009";
    text1b = "900 GeV Runs";
    leg1a  = "DATA";
    leg1b  = "MC";

    text2a = "CMS Preliminary 2009";
    text2b = "";
    leg2a  = "MC (using MC-truth)";
    leg2b  = "MC (using only reco)";

    outputName1 = "resoDz_vs_pt.DATA.png";
    outputName2 = "resoDz_vs_pt.MC.png";

    nbins = 60;
    xlow=0.7;
    xhigh=2.2;
  }


  if(type==5){
    gStyle->SetOptLogy(1);
    fnRawSim  = rootFolderNameRaw+"sim.raw.dz.vsEta.50bins.root";
    fnRawData = rootFolderNameRaw+"data.raw.dz.vsEta.50bins.root";
    fnResp    = rootFolderNameResp+"sim.resp.dzResp.vsEta.50bins.root";
    fnReso    = rootFolderNameResp+"sim.reso.dzReso.vsEta.50bins.root";

    xLabel = "Track #eta";
    yLabel = "Longit.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the LIP due to Vertex (#mum)";

    lRange = 15.;
    hRange = 3000.;
    text1a = "CMS Preliminary 2009";
    text1b = "900 GeV Runs";
    leg1a  = "DATA";
    leg1b  = "MC";

    text2a = "CMS Preliminary 2009";
    text2b = "";
    leg2a  = "MC (using MC-truth)";
    leg2b  = "MC (using only reco)";

    outputName1 = "resoDz_vs_eta.Pt08.DATA.png";
    outputName2 = "resoDz_vs_eta.Pt08.MC.png";

    nbins = 50;
    xlow= -2.5;
    xhigh= 2.5;
  }


  if(type==6){
    //STILL TO BE DEFINED
  }




  TFile* fileRawSim  = new TFile(fnRawSim);
  TFile* fileRawData = new TFile(fnRawData);
  TFile* fileResp    = new TFile(fnResp);
  TFile* fileReso    = new TFile(fnReso);
  
  TCanvas* canvas1 = new TCanvas("canvas1","canvas1",0,0,500,500);
  TCanvas* canvas2 = new TCanvas("canvas2","canvas2",500,0,500,500);
  canvas2->SetRightMargin(0.05);
  canvas2->SetLeftMargin(0.15);


  TH1F* plotReso    = new TH1F("plotReso","plotReso",nbins,xlow,xhigh);
  TH1F* plotSimFit  = new TH1F("plotSimFit","plotSimFit",nbins,xlow,xhigh);
  TH1F* plotDataFit = new TH1F("plotDataFit","plotDataFit",nbins,xlow,xhigh);

  TH1F* plotResp1  = new TH1F("plotResp1","plotResp2",nbins,xlow,xhigh);
  TH1F* plotResp2  = new TH1F("plotResp1","plotResp2",nbins,xlow,xhigh);
  TH1F* plotResp1B = new TH1F("plotResp1B","plotResp2B",nbins,xlow,xhigh);
  TH1F* plotResp2B = new TH1F("plotResp1B","plotResp2B",nbins,xlow,xhigh);

  double resolutions[nbins_const+1]; resolutions[0]=0;
  double responses1[nbins_const+1];  responses1[0]=0;
  double responses2[nbins_const+1];  responses1[0]=0;
  double fracs[nbins_const+1];  fracs[0]=0;

  double fits[nbins_const+1]; fits[0]=0;
  double dataFits[nbins_const+1]; dataFits[0]=0;

  double resolutionsErr[nbins_const+1]; resolutionsErr[0]=0;
  double fitsErr[nbins_const+1]; fitsErr[0]=0;
  double dataFitsErr[nbins_const+1]; dataFitsErr[0]=0;

  for(unsigned int i=0; i<(nbins+1); ++i){
    resolutions[i]=0;  resolutionsErr[i]=0;
    fits[i]=0;         fitsErr[i]=0;
    dataFits[i]=0;     dataFitsErr[i]=0;

    responses1[i]=0;
    responses2[i]=0;
    fracs[i]=0;
  }

  for(unsigned int i=1; i<(nbins+1); ++i){ // loop on [1,nbins]
    stringstream stream;  stream << i;
    TString counter = stream.str();

    TString hn1,hn2,hn3,hn4;
    hn1 = hname1+counter;
    hn2 = hname2+counter;
    hn3 = hname3+counter;
    hn4 = hname4+counter;
  
    cout << "hn1: " << hn1 << endl;
    cout << "hn2: " << hn2 << endl;
    cout << "hn3: " << hn3 << endl;
    cout << "hn4: " << hn4 << endl;
    //wantMore2();


    TH1* h_bef      = (TH1*) fileRawSim->Get(hn3);
    TH1* h_befData  = (TH1*) fileRawData->Get(hn4);
    TH1* h_resp     = (TH1*) fileResp->Get(hn2);
    TH1* h_reso     = (TH1*) fileReso->Get(hn1);

  
    double reso,resp1,resp2,frac,fit,dataFit;
    double resoErr,fitErr,dataFitErr;

    double inputResp1(initialResp1),inputResp2(initialResp2),inputFrac(initialFrac);
    
    if(i>1){
      inputResp1 = responses1[i-1];
      inputResp2 = responses2[i-1];
      inputFrac = fracs[i-1];
    }

    //cout << "counter: " << counter << endl; wantMore2();
    getResponse(h_reso,h_resp,h_bef,h_befData, //input projection histos 
		reso,resp1,resp2,frac,         //MC-reso and resp to be fitted
		fit,dataFit,                   //reso to be measured
		resoErr,fitErr,dataFitErr,
		inputResp1,inputResp2,inputFrac,fitRangeMax,
		canvas1,typePlot);

    resolutions[i] = reso;
    responses1[i]  = resp1;
    responses2[i]  = resp2;
    fracs[i]       = frac;

    fits[i]        = fit;
    dataFits[i]    = dataFit;

    resolutionsErr[i] = resoErr;
    fitsErr[i]        = fitErr;
    dataFitsErr[i]    = dataFitErr;

    
    plotHistos(plotReso,plotSimFit,plotDataFit,plotResp1,plotResp2,
	       resolutions, fits, dataFits,responses1,responses2,
	       resolutionsErr,fitsErr,dataFitsErr,
	       yLabel,xLabel,hRange,lRange,
	       canvas2,typePlot2);    
  }
  

  cout << endl;
  for(unsigned int i=1; i<nbins+1; ++i){
    //cout << "sigmaResolution(" << i << "): " << resolutions[i] << endl;
    //cout << "sigmaResponse1(" << i << "): "<< responses1[i] << endl;
    //cout << "sigmaResponse2(" << i << "): "<< responses2[i] << endl;
    cout << "fracs(" << i << "): " << fracs[i] << endl;
  }
  
  return;
  //==================== I"m currently skipping the rest of the class =================
  //===================================================================================

  /*
  plotResp1->SetMaximum(120);
  plotResp1->SetMinimum(0);
  
  TF1* myResp1Fit;
  TF1* myResp2Fit;

  if(type==5){//dz vs eta
    myResp1Fit = new TF1("myResp1Fit","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",0.5,2.0);
    myResp2Fit = new TF1("myResp2Fit","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",0.5,2.0);
  }else{
    myResp1Fit = new TF1("myResp1Fit","[0]+[1]*x+[2]*x*x",0.5,2.0);
    myResp2Fit = new TF1("myResp2Fit","[0]+[1]*x+[2]*x*x",0.5,2.0);
  }

  //TF1* myResp1Fit = new TF1("myResp1Fit","[0]+[1]*x",0.5,2.0);
  //TF1* myResp2Fit = new TF1("myResp2Fit","[0]+[1]*x",0.5,2.0);

  myResp1Fit->SetParName(0,"q");
  myResp1Fit->SetParName(1,"p");
  myResp1Fit->SetParName(2,"p2");
  myResp1Fit->SetLineWidth(2);
  myResp1Fit->SetLineColor(2);

  myResp2Fit->SetParName(0,"q");
  myResp2Fit->SetParName(1,"p");
  myResp2Fit->SetParName(2,"p2");
  myResp2Fit->SetLineWidth(2);
  myResp2Fit->SetLineColor(2);

  plotResp1->Fit("myResp1Fit");   
  plotResp2->Fit("myResp2Fit");   
  
  plotResp1->Draw("E1");   gPad->Update();     wantMore2();
  plotResp2->Draw("E1same");  gPad->Update();  wantMore2();

  //TH1* fit1Histo = myResp1Fit->GetHistogram();
  //TH1* fit2Histo = myResp2Fit->GetHistogram();
  
  for(unsigned int i=1; i!=nbins+1; ++i){
    double xvalue = plotResp1->GetBinLowEdge(i)+0.5*plotResp1->GetBinWidth(i); 
    double tmp = myResp1Fit->Eval(xvalue);
    cout << "x,xvalue: " << xvalue << " , " << tmp << endl;
    plotResp1B->SetBinContent(i,tmp);

    xvalue = plotResp2->GetBinLowEdge(i)+0.5*plotResp2->GetBinWidth(i);
    tmp = myResp2Fit->Eval(xvalue);
    cout << "x,xvalue: " << xvalue << " , " << tmp << endl;
    plotResp2B->SetBinContent(i,tmp);

  }


  canvas1->cd(); canvas1->Clear();gPad->Update();
  plotResp1->SetMaximum(hRange);
  plotResp1->SetMinimum(lRange);


  plotResp1->Draw("E1");   gPad->Update();     wantMore2();
  plotResp2->Draw("E1same");  gPad->Update();  wantMore2();

  canvas2->cd(); canvas2->Clear();gPad->Update();

  //myResp1Fit->Draw("E1");   gPad->Update();     wantMore2();
  //myResp2Fit->Draw("E1same");  gPad->Update();  wantMore2();


  //canv->cd(); canv->Clear(); gPad->Update();

  
  canvas1->cd();
  fit1Histo->SetMaximum(120);
  fit1Histo->SetMinimum(0);
  
  fit1Histo->Draw("E1");   gPad->Update();     wantMore2();
  fit2Histo->Draw("E1same");  gPad->Update();  wantMore2();


  for(unsigned int i=1; i<nbins+1; ++i){
    //resolutions[i]=0;
    responses1[i]=0;
    responses2[i]=0;
    fracs[i]=0;
    fits[i]=0;
    dataFits[i]=0;
    //resolutionsErr[i]=0;
    fitsErr[i]=0;
    dataFitsErr[i]=0;
  }



  // fit again resolutions after stabilization of response functions
  for(unsigned int i=1; i<nbins+1; ++i){
    stringstream stream;  stream << i;
    TString counter = stream.str();

    cout << "counter: " << counter << endl;

    TString hn1,hn2,hn3,hn4;
    hn1 = hname1+counter;
    hn2 = hname2+counter;
    hn3 = hname3+counter;
    hn4 = hname4+counter;
  

    TH1* h_bef     = (TH1*) fileRawSim->Get(hn3);
    TH1* h_befData = (TH1*) fileRawData->Get(hn4);
    TH1* h_resp    = (TH1*) fileResp->Get(hn2);
    TH1* h_reso    = (TH1*) fileReso->Get(hn1);
  
    double frac;
    double fit,dataFit;
    double fitErr,dataFitErr;

    double inputResp1(initialResp1),inputResp2(initialResp2),inputFrac(initialFrac);
    
    if(i>1){
      //inputResp1 = fit1Histo->GetBinContent(i-1);
      //inputResp2 = fit2Histo->GetBinContent(i-1);
      inputResp1 = plotResp1B->GetBinContent(i-1);
      inputResp2 = plotResp2B->GetBinContent(i-1);
      inputFrac = fracs[i-1];
    }

    //double prev_resp1 = fit1Histo->GetBinContent(i);
    //double prev_resp2 = fit2Histo->GetBinContent(i);
    double prev_resp1 = plotResp1B->GetBinContent(i);
    double prev_resp2 = plotResp2B->GetBinContent(i);

    
    getResponse2(prev_resp1,prev_resp2,
		 frac,
		 h_resp,h_bef,h_befData,
		 fit,dataFit,
		 fitErr,dataFitErr,
		 inputResp1,inputResp2,inputFrac,
		 fitRangeMax,
		 canvas1,
		 typePlot3);

    responses1[i]= prev_resp1;
    responses2[i]= prev_resp2;

    fracs[i]       = frac;
    fits[i]        = fit;
    dataFits[i]    = dataFit;

    fitsErr[i]        = fitErr;
    dataFitsErr[i]    = dataFitErr;
    
   
    plotHistos(plotReso,plotSimFit,plotDataFit,plotResp1,plotResp2,
	       resolutions, fits, dataFits,responses1,responses2,
	       resolutionsErr,fitsErr,dataFitsErr,
	       yLabel,xLabel,hRange,lRange,
	       canvas2,typePlot);
   
  }

  


  //plotDataFit->Draw("E1");
  //plotSimFit->Draw("sameE1");
  //gPad->Update();

  TPaveText* txt =  new TPaveText(0.45,0.80,0.75,0.90,"NDC");
  TLegend* leg   = new TLegend(0.50,0.65,0.80,0.79);

  TPaveText* txt2 =  new TPaveText(0.45,0.80,0.75,0.90,"NDC");
  TLegend* leg2   = new TLegend(0.50,0.65,0.80,0.79);
			     
  leg->SetShadowColor(0);
  leg->SetFillColor(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.035);
  
  txt->SetShadowColor(0);
  txt->SetFillColor(0);
  txt->SetLineColor(0);
  txt->SetTextSize(0.035);

  leg2->SetShadowColor(0);
  leg2->SetFillColor(0);
  leg2->SetLineColor(0);
  leg2->SetTextSize(0.035);
  
  txt2->SetShadowColor(0);
  txt2->SetFillColor(0);
  txt2->SetLineColor(0);
  txt2->SetTextSize(0.035);

  leg->AddEntry(plotDataFit,leg1a,"lp");
  leg->AddEntry(plotSimFit,leg1b,"lp");
  txt->AddText(text1a);
  txt->AddText(text1b);

  leg2->AddEntry(plotReso,leg2a,"lp");
  leg2->AddEntry(plotSimFit,leg2b,"lp");
  txt2->AddText(text2a);


  leg->Draw();
  txt->Draw();

  gPad->Print(outputName1);
  wantMore2();

  gPad->Clear();
  //plotSimFit->Draw("E1");
  //plotReso->Draw("sameE1");

  leg2->Draw();
  txt2->Draw();
  gPad->Update();

  gPad->Print(outputName2);



  */
}


void getResponse(TH1* h_reso,TH1* h_resp,TH1* h_bef,TH1* h_befData, //input projection histos 
		 double& reso, double& resp1, double& resp2,double& frac,
		 double& fit,double& dataFit,
		 double& resoErr, double& fitErr, double& dataFitErr,
		 double inputResp1,double inputResp2,double inputFrac,double rangeMax,
		 TCanvas* canvas1,int typePlot)
{ 
  canvas1->cd();
  //============ RESOLUTION
  double range = getFitRange(h_reso);
  RooRealVar x1("x1","MC-based track-TP residual",-range*2.,+range*2.);
  RooDataHist datah("datah","dataset from h",x1,h_reso) ;

  RooPlot* frame1 = x1.frame(); 
  datah.plotOn(frame1);


  RooRealVar mean1("mean1","mean of gaussian1",0,-10,10) ;
  RooRealVar sigma1("sigma1","width of gaussian1",inputResp1,10,rangeMax) ;
  RooGaussian gauss1("gauss1","gaussian PDF 1",x1,mean1,sigma1) ;  

  gauss1.fitTo(datah);
  gauss1.plotOn(frame1);
  reso = sigma1.getVal();
  resoErr = sigma1.getError();

  // --- printouts
  cout << "+++ fitting resolution+++" << endl; 
  cout << "range: " << -2.*range << " , " << 2.*range << endl;
  cout << "MC-reso, MC-resoErr: " << reso << " , " << resoErr << endl; 
  cout << "MC-mean, MC-meanErr: " 
       << mean1.getVal() << " , "
       << mean1.getError() << endl;

  // --- draw histograms
  if(typePlot==1){
    frame1->Draw(); gPad->Update(); 
    wantMore2();
  }



  

  //============ RESPONSE
  range = getFitRange(h_resp);
  RooRealVar x2("x2","Vertex smearing",-range*2.,+range*2.);
  x2.setRange(-range*2.,+range*2.);
  RooDataHist datah2("datah2","dataset from h",x2,h_resp) ;

  RooPlot* frame2 = x2.frame(); 
  datah2.plotOn(frame2);
 
  RooRealVar mean2("mean2","mean of gaussian2",0,-10,10) ;
  RooRealVar sigma2("sigma2","width of gaussian2",inputResp1,0,rangeMax) ;
  RooGaussian gauss2("gauss2","gaussian PDF 2",x2,mean2,sigma2) ;  


  gauss2.fitTo(datah2);

  // --- printouts
  resp1 = sigma2.getVal();
  cout << "=======" << endl;
  cout << "resp1: " << resp1 << endl;


  // --- draw histograms
  if(typePlot==2){
    gauss2.plotOn(frame2,LineColor(kBlue));
    frame2->Draw(); gPad->Update(); 
    wantMore2();
  }



  //============ CONVOLUTION
  //range = getFitRange(h_bef);
  range = getFitRange(h_reso);
  RooRealVar x3("x3","raw residual sim",-range*2.,+range*2.);
  x3.setRange(-range*2.,+range*2.);
  RooDataHist datah3("datah3","dataset from h",x3,h_bef) ;
  RooPlot* frame3 = x3.frame(); 


  RooRealVar mean2b("mean2b","mean of gaussian2b",0,-20,20) ;
  RooRealVar sigma2b("sigma2b","width of gaussian2b",sigma2.getVal(),0,rangeMax) ;
  RooGaussian gauss2b("gauss2b","gaussian PDF 2b",x3,mean2b,sigma2b) ;  

  RooRealVar mean3("mean3","mean of gaussian3",0,-20,20) ;
  RooRealVar sigma3("sigma3","width of gaussian3",reso,0,rangeMax) ;
  RooGaussian gauss3("gauss3","gaussian PDF 3",x3,mean3,sigma3) ;  
  RooFFTConvPdf conv("conv"," gauss2 (X) gauss3",x3,gauss2b,gauss3) ;

  //fix values of the vertex smearing pdf
  mean2b.setConstant(kTRUE);
  sigma2b.setConstant(kTRUE);

  //avoid to fit the mean
  mean3.setConstant(kTRUE);

  conv.fitTo(datah3) ;


  // --- printouts
  fit = sigma3.getVal();
  fitErr = sigma3.getError();
  cout << "+++ fitting conv1 +++" << endl; 
  cout << "range: " << -range*2. << " , " << range*2. << endl;
  cout << "reso: " << reso << endl;
  cout << "    mean2,mean3: " << mean2b.getVal() << " , " << mean3.getVal() << endl;
  cout << "    sigma2,sigma3: " << sigma2b.getVal() << " , " << sigma3.getVal() << endl;



  // --- draw histograms
  if(typePlot==3){
    datah3.plotOn(frame3);
    conv.plotOn(frame3,LineColor(kBlue)) ;
    frame3->Draw(); gPad->Update(); 
    wantMore2();
  }


  

  //============ CONVOLUTION 2  
  //range = getFitRange(h_befData);
  //cout << "++++++++++++++++++++ range befData: " << range << endl;
  range = getFitRange(h_reso);
  RooRealVar x4("x4","raw residual DATA",-range*2.,+range*2.);
  x4.setRange(-range*2.,+range*2.);

  RooDataHist datah4("datah4","dataset from h",x4,h_befData) ;
  RooPlot* frame4 = x4.frame(); 

  RooRealVar mean2c("mean2c","mean of gaussian2c",0,-20,20) ;
  RooRealVar sigma2c("sigma2c","width of gaussian2c",sigma2.getVal(),0,rangeMax) ;
  RooGaussian gauss2c("gauss2c","gaussian PDF 2c",x4,mean2c,sigma2c) ;  

  RooRealVar mean4("mean4","mean of gaussian4",0,-20.,20.) ;
  RooRealVar sigma4("sigma4","width of gaussian4",reso,0,rangeMax) ;
  RooGaussian gauss4("gauss4","gaussian PDF 4",x4,mean4,sigma4) ;  
  RooFFTConvPdf conv2("conv2"," gauss2c (X) gauss4",x4,gauss2c,gauss4) ;

  //fix values of the vertex smearing pdf
  mean2c.setConstant(kTRUE);
  sigma2c.setConstant(kTRUE);

  //avoid to fit the mean
  mean4.setConstant(kTRUE);

  conv2.fitTo(datah4) ;

  
  // --- printouts
  cout << "+++ fitting conv2+++" << endl;
  dataFit = sigma4.getVal();
  dataFitErr = sigma4.getError();


  // --- draw histograms
  if(typePlot==4){
    datah4.plotOn(frame4);
    conv2.plotOn(frame4,LineColor(kBlue)) ;
    frame4->Draw(); gPad->Update(); 
    wantMore2();
  }
}




void getResponse2(double resp1, //estimated resp1 from MC (previous step)
		  double resp2, //estimated resp1 from MC (previous step)
		  double& frac, //do I need this??
		  TH1* h_resp,TH1* h_bef,TH1* h_befData,     	  
		  double& fit,double& dataFit,
		  double& fitErr, double& dataFitErr, 
		  double inputResp1, //fitted resolution (MC) from previous step
		  double inputResp2, //fitted resolution (DATA) from previous step
		  double inputFrac, //fitted fraction  from previous step ??
		  double rangeMax,
		  TCanvas* canvas1,
		  int typePlot)
{ 
  canvas1->cd();
  //=========== again response
  double range = getFitRange(h_resp);
  RooRealVar x("x","x",-range*2.,+range*2.);
  x.setRange(-range*2.,+range*2.);
  RooDataHist datah2("datah2","dataset from h",x,h_resp) ;

  RooPlot* frame2 = x.frame(); 
  datah2.plotOn(frame2);

  



  //=========== preparing 2-gaus sum function
  RooRealVar mean1("mean1","mean of gaussian1",0,-10,10) ;
  RooRealVar sigma1("sigma1","width of gaussian1",resp1,10,rangeMax) ;
  RooGaussian gauss1("gauss1","gaussian PDF 1",x,mean1,sigma1) ;  

  RooRealVar mean2("mean2","mean of gaussian2",0,-10,10) ;
  RooRealVar sigma2("sigma2","width of gaussian2",resp2,5,rangeMax) ;
  RooGaussian gauss2("gauss2","gaussian PDF 2",x,mean2,sigma2) ;  

  RooRealVar fraction("fraction","fraction second gaussian",inputFrac,0.55,0.80);//for pt plot
  RooAddPdf twoGausSum("twoGausSum","twoGausSum pdf",RooArgList(gauss1,gauss2),fraction);

  mean1.setConstant(kTRUE);
  sigma1.setConstant(kTRUE);
  mean2.setConstant(kTRUE);
  sigma2.setConstant(kTRUE);

  twoGausSum.fitTo(datah2);
  if(typePlot==1){
    twoGausSum.plotOn(frame2,LineColor(kBlack),LineStyle(kDashed));
    gauss1.plotOn(frame2,LineColor(kGreen));
    gauss2.plotOn(frame2,LineColor(kBlue));
    frame2->Draw(); gPad->Update(); 
    cout << "=======fraction: " << fraction.getVal() << endl; //wantMore2();
  }

  frac  = fraction.getVal();

  //============ CONVOLUTION
  range = getFitRange(h_bef);
  x.setRange(-range*2.,+range*2.);
  RooDataHist datah3("datah3","dataset from h",x,h_bef) ;

  RooPlot* frame3 = x.frame(); 
  datah3.plotOn(frame3);

  

  RooRealVar mean3("mean3","mean of gaussian3",0,-20,20) ;
  //RooRealVar sigma3("sigma3","width of gaussian3",20,0,200) ;
  RooRealVar sigma3("sigma3","width of gaussian3",inputResp1,0,rangeMax) ;
  RooGaussian gauss3("gauss3","gaussian PDF 3",x,mean3,sigma3) ;  
  RooFFTConvPdf conv("conv"," twoGausSum (X) gauss3",x,twoGausSum,gauss3) ;

  fraction.setConstant(kTRUE);
  mean1.setConstant(kTRUE);
  sigma1.setConstant(kTRUE);
  mean2.setConstant(kTRUE);
  sigma2.setConstant(kTRUE);
  mean3.setConstant(kTRUE);

  conv.fitTo(datah3) ;
  conv.plotOn(frame3,LineColor(kBlack),LineStyle(kDashed)) ;
  gauss3.plotOn(frame3,LineColor(kRed));
  //twoGausSum.plotOn(frame3,LineColor(kBlack));
  gauss1.plotOn(frame3,LineColor(kGreen));
  gauss2.plotOn(frame3,LineColor(kBlue));

  
  fit = sigma3.getVal();
  fitErr = sigma3.getError();

  //frame3->Draw(); gPad->Update();
  cout << "+++ fitting conv1+++" << endl; //wantMore2();

  //============ CONVOLUTION 2
  
  //range = getFitRange(h_befData);
  x.setRange(-range*2.,+range*2.);
  RooDataHist datah4("datah4","dataset from h",x,h_befData) ;

  RooPlot* frame4 = x.frame(); 
  datah4.plotOn(frame4);


  RooRealVar mean4("mean4","mean of gaussian4",0,-5.,5.) ;
  RooRealVar sigma4("sigma4","width of gaussian4",inputResp2,0,rangeMax) ;
  RooGaussian gauss4("gauss4","gaussian PDF 4",x,mean4,sigma4) ;  
  RooFFTConvPdf conv2("conv2"," twoGausSum (X) gauss4",x,twoGausSum,gauss4) ;


  fraction.setConstant(kTRUE);
  mean1.setConstant(kTRUE);
  sigma1.setConstant(kTRUE);
  mean2.setConstant(kTRUE);
  sigma2.setConstant(kTRUE);
  //mean4.setConstant(kTRUE);

  conv2.fitTo(datah4) ;
  if(typePlot==2){
    conv2.plotOn(frame4,LineColor(kBlack),LineStyle(kDashed)) ;    
    gauss4.plotOn(frame4,LineColor(kRed));
    //gauss2.plotOn(frame2,LineColor(kBlue));
    frame4->Draw(); gPad->Update(); wantMore2();
  }




  //frame4->Draw(); gPad->Update(); 
  cout << "+++ fitting conv2+++" << endl;

  dataFit = sigma4.getVal();
  dataFitErr = sigma4.getError();


  //dataFit = 10;
  //dataFitErr = 5; 

}





double getFitRange(TH1* h){
  h->GetListOfFunctions()->Delete();
  h->Fit("gaus","Q");gPad->Clear();
  TF1* firstFit = (TF1*) h->GetListOfFunctions()->Last();
  //firstFit->Print();

  double firstSigma = firstFit->GetParameter(2);
  //cout << "firstSigma: " << firstSigma << endl;
  //h->Draw(); gPad->Update();

  double nSigma=1.;
  return nSigma*firstSigma;

}


void plotHistos(TH1F* plotReso, TH1F* plotSimFit, TH1F* plotDataFit,TH1F* plotResp1,TH1F* plotResp2,
		double* resolutions, double* fits, double* dataFits, double* responses1,double* responses2,
		double* resolutionsErr, double* fitsErr, double* dataFitsErr,
		TString yLabel, TString xLabel, double hRange, double lRange,
		TCanvas* canvas2, int typePlot2){

  // --- set bin contents
  plotReso->SetContent(resolutions);   plotReso->SetError(resolutionsErr);
  plotSimFit->SetContent(fits);           plotSimFit->SetError(fitsErr);
  plotDataFit->SetContent(dataFits);   plotDataFit->SetError(dataFitsErr);
  plotResp1->SetContent(responses1);    
  plotResp2->SetContent(responses2);           


  // --- set labels
  plotReso->SetTitle("");
  plotReso->GetYaxis()->SetTitleSize(0.04);
  plotReso->GetXaxis()->SetTitleSize(0.04);
  plotReso->GetYaxis()->SetLabelSize(0.04);
  plotReso->GetXaxis()->SetLabelSize(0.04);
  plotReso->GetYaxis()->SetTitleOffset(1.4);
  plotReso->SetYTitle(yLabel);
  plotReso->SetXTitle(xLabel);

  plotResp1->SetTitle("");
  plotResp1->GetYaxis()->SetTitleSize(0.04);
  plotResp1->GetXaxis()->SetTitleSize(0.04);
  plotResp1->GetYaxis()->SetLabelSize(0.04);
  plotResp1->GetXaxis()->SetLabelSize(0.04);
  plotResp1->GetYaxis()->SetTitleOffset(1.4);
  plotResp1->SetYTitle(yLabelResp);
  plotResp1->SetXTitle(xLabel);

  plotDataFit->SetTitle("");
  plotDataFit->GetYaxis()->SetTitleSize(0.04);
  plotDataFit->GetXaxis()->SetTitleSize(0.04);
  plotDataFit->GetYaxis()->SetLabelSize(0.04);
  plotDataFit->GetXaxis()->SetLabelSize(0.04);
  plotDataFit->GetYaxis()->SetTitleOffset(1.4);
  plotDataFit->SetYTitle(yLabel);
  plotDataFit->SetXTitle(xLabel);

  plotSimFit->SetTitle("");
  plotSimFit->GetYaxis()->SetTitleSize(0.04);
  plotSimFit->GetXaxis()->SetTitleSize(0.04);
  plotSimFit->GetYaxis()->SetLabelSize(0.04);
  plotSimFit->GetXaxis()->SetLabelSize(0.04);
  plotSimFit->GetYaxis()->SetTitleOffset(1.4);
  plotSimFit->SetYTitle(yLabel);
  plotSimFit->SetXTitle(xLabel);


  // --- set ranges
  plotReso->SetMaximum(hRange);
  plotReso->SetMinimum(lRange);

  plotResp1->SetMaximum(hRange);
  plotResp1->SetMinimum(lRange);

  plotDataFit->SetMaximum(hRange);
  plotDataFit->SetMinimum(lRange);

  plotSimFit->SetMaximum(hRange);
  plotSimFit->SetMinimum(lRange);
  


  plotReso->SetMarkerStyle(20);
  plotReso->SetLineColor(1);
  plotReso->SetMarkerColor(1);
  plotReso->SetLineWidth(1);

  plotResp1->SetMarkerStyle(20);
  plotResp1->SetLineColor(1);
  plotResp1->SetMarkerColor(3);
  plotResp1->SetLineWidth(1);

  plotSimFit->SetMarkerStyle(21);
  plotSimFit->SetLineColor(1);
  plotSimFit->SetMarkerColor(2);
  plotSimFit->SetLineWidth(1);

  plotDataFit->SetMarkerStyle(20);
  plotDataFit->SetLineColor(1);
  plotDataFit->SetMarkerColor(1);
  plotDataFit->SetLineWidth(1);

  canvas2->cd();

  switch ( typePlot2 ) {
  case 0 :
    //draw nothing
    break;
    
  case 1 :
    plotReso->Draw("E1"); gPad->Update();
    break;
      
  case 2 :
    plotResp1->Draw("E1"); gPad->Update();
    plotResp2->Draw("sameE1"); gPad->Update();
    break;
      
  case 3 :
    plotReso->Draw("E1"); gPad->Update();
    plotSimFit->Draw("sameE1"); gPad->Update();
    break;

  case 4 :
    //plotSimFit->Draw("E1"); gPad->Update();
    plotDataFit->Draw("E1"); gPad->Update();
    break;
      
  default :
    break;
  }
}
