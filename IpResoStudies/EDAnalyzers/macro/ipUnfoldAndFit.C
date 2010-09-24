#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
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

void fitRawErrors(TH1* h_rErrSim,TH1* h_rErrData,
		  double& rErrSim, double& rErrData,
		  double& rErrSimUncert, double& rErrDataUncert,
		  TCanvas* canv);

void fitResponse(TH1* h_resp, int type, double xlow, double xhigh);

void plotHistos(TH1F* plotReso, TH1F* plotFit, TH1F* plotDataFit,TH1F* plotResp1,TH1F* plotResp2,
		double* resolutions, double* fits, double* dataFits,
		double* responses1,double* responses2,
		double* resolutionsErr, double* fitsErr, double* dataFitsErr,
		TString yLabel, TString xLabel, double hRange, double lRange,
		TCanvas* canvas2,int type=0);

void plotRawErrHistos(TH1* plotErrSim, TH1* plotErrData,
		      double* errSim, double* errData,
		      double* errSimUncert, double* errDataUncert,
		      TString xLabel,
		      double hRange, double lRange,
		      TCanvas* canv);

void finalFit(TH1F* histoToFit1, TH1F* histoToFit2,  TString xLabel, TString yLabel, TPaveText *txt, TString fitOutputfile);

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

  gStyle->SetFillColor(0);
  gStyle->SetLineWidth(1);
  gStyle->SetLineColor(0);
  


  TString fnResp;
  TString fnReso;
  TString fnRawSim;
  TString fnRawErrSim;
  TString fnRawData;
  TString fnRawErrData;

  double lRange,hRange;

  TString text1a,text1b,leg1a,leg1b;
  TString text2a,text2b,leg2a,leg2b;
  
  TString outputName1,outputName2,outputName3,outputName4;

  const unsigned int nbins_const = 60;
  //const unsigned int nbins_const = 50;
  unsigned int nbins;
  double xlow,xhigh;
  double fitRangeMax(0),initialResp1(0),initialResp2(0),initialFrac(0);


  TString hname1, hname2, hname3,hname4,hname5,hname6;

  //TString inputFolderName="./lastVersionOutputProjection/";
  TString inputFolderName="./input/";
  TString outputFolderName="./finalPlots/";  

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
      hname5 ="sim_raw_d0FitErr_vsPt_n"; //raw error MC
      hname6 ="data_raw_d0FitErr_vsPt_n"; //raw error DATA
      fitRangeMax = 400;
      break;

    case 2: //d0 vs eta
      hname1 ="sim_reso_dxyReso_vsEta_n"; //resolution
      hname2 ="sim_resp_dxyResp_vsEta_n"; //response 
      hname3 ="sim_raw_d0_vsEta_n"; //raw MC
      hname4 ="data_raw_d0_vsEta_n"; //raw DATA
      hname5 ="sim_raw_d0FitErr_vsEta_n"; //raw error MC
      hname6 ="data_raw_d0FitErr_vsEta_n"; //raw error DATA
      fitRangeMax = 600;
      break;

    case 3: //d0 vs phi
      hname1 ="sim_reso_dxyReso_vsPhi_n"; //resolution
      hname2 ="sim_resp_dxyResp_vsPhi_n"; //response 
      hname3 ="sim_raw_d0_vsPhi_n"; //raw MC
      hname4 ="data_raw_d0_vsPhi_n"; //raw DATA
      hname5 ="sim_raw_d0FitErr_vsPhi_n"; //raw error MC
      hname6 ="data_raw_d0FitErr_vsPhi_n"; //raw error DATA
      fitRangeMax = 400;
      break;

    case 4: //dz vs pt
      hname1 ="sim_reso_dzReso_vsPt_n"; //resolution
      hname2 ="sim_resp_dzResp_vsPt_n"; //response 
      hname3 ="sim_raw_dz_vsPt_n"; //raw MC
      hname4 ="data_raw_dz_vsPt_n"; //raw DATA
      hname5 ="sim_raw_dzFitErr_vsPt_n"; //raw error MC
      hname6 ="data_raw_dzFitErr_vsPt_n"; //raw error DATA
      fitRangeMax = 1500;
      break;

    case 5: //dz vs eta
      hname1 ="sim_reso_dzReso_vsEta_n"; //resolution
      hname2 ="sim_resp_dzResp_vsEta_n"; //response 
      hname3 ="sim_raw_dz_vsEta_n"; //raw MC
      hname4 ="data_raw_dz_vsEta_n"; //raw DATA
      hname5 ="sim_raw_dzFitErr_vsEta_n"; //raw error MC
      hname6 ="data_raw_dzFitErr_vsEta_n"; //raw error DATA
      fitRangeMax = 1500;
      break;

    case 6: //dz vs phi
      hname1 ="sim_reso_dzReso_vsPhi_n"; //resolution
      hname2 ="sim_resp_dzResp_vsPhi_n"; //response 
      hname3 ="sim_raw_dz_vsPhi_n"; //raw MC
      hname4 ="data_raw_dz_vsPhi_n"; //raw DATA
      hname5 ="sim_raw_dzFitErr_vsPhi_n"; //raw error MC
      hname6 ="data_raw_dzFitErr_vsPhi_n"; //raw error DATA
      fitRangeMax = 1500;
      break;

      
    default:
      break;
    }
  

  if(type==1){
    fnRawSim     = inputFolderName+"sim.raw.d0.vsPt.60bins.root";
    fnRawErrSim  = inputFolderName+"sim.raw.d0FitErr.vsPt.60bins.root";
    fnRawData    = inputFolderName+"data.raw.d0.vsPt.60bins.root";
    fnRawErrData = inputFolderName+"data.raw.d0FitErr.vsPt.60bins.root";
    fnResp    = inputFolderName+"sim.resp.dxyResp.vsPt.60bins.root";
    fnReso    = inputFolderName+"sim.reso.dxyReso.vsPt.60bins.root";

    xLabel = "Track p_{T} (GeV/c)";
    yLabel = "Transv.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the TIP due to Vertex (#mum)";
    
    lRange = 0.;
    hRange = 210.;

    text1a = "CMS Preliminary 2010";
    text1b = "7 TeV Runs";
    leg1a  = "DATA";
    leg1b  = "MC";

    text2a = "CMS Preliminary 2010";
    text2b = "";
    leg2a  = "MC (using MC-truth)";
    leg2b  = "MC (using only reco)";

    outputName1 = outputFolderName+"resoD0_vs_pt.DATA.png";
    outputName2 = outputFolderName+"resoD0_vs_pt.MC.png";
    outputName3 = outputFolderName+"errorD0_vs_pt.png";
    outputName4 = outputFolderName+"resoErrorRatio_D0_vs_pt.png";

    nbins = 60;
    xlow=0.7;
    xhigh=2.2;
  }

  if(type==2){
    fnRawSim     = inputFolderName+"sim.raw.d0.vsEta.50bins.root";
    fnRawErrSim  = inputFolderName+"sim.raw.d0FitErr.vsEta.50bins.root";
    fnRawData    = inputFolderName+"data.raw.d0.vsEta.50bins.root";
    fnRawErrData = inputFolderName+"data.raw.d0FitErr.vsEta.50bins.root";
    fnResp    = inputFolderName+"sim.resp.dxyResp.vsEta.50bins.root";
    fnReso    = inputFolderName+"sim.reso.dxyReso.vsEta.50bins.root";

    xLabel = "Track #eta";
    yLabel = "Transv.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the TIP due to Vertex (#mum)";

    lRange = 0.;
    //hRange = 250.;
    hRange = 350.;
    text1a = "CMS Preliminary 2010";
    text1b = "7 TeV Runs";
    leg1a  = "DATA";
    leg1b  = "MC";

    text2a = "CMS Preliminary 2010";
    text2b = "";
    leg2a  = "MC (using MC-truth)";
    leg2b  = "MC (using only reco)";

    outputName1 = outputFolderName+"resoD0_vs_eta.Pt08.DATA.png";
    outputName2 = outputFolderName+"resoD0_vs_eta.Pt08.MC.png";
    outputName3 = outputFolderName+"errorD0_vs_eta.png";
    outputName4 = outputFolderName+"resoErrorRatio_D0_vs_eta.png";

    nbins = 50;
    xlow= -2.5;
    xhigh= 2.5;
  }

  if(type==3){
    //STILL TO BE DEFINED
  }

  if(type==4){
    fnRawSim     = inputFolderName+"sim.raw.dz.vsPt.60bins.root";
    fnRawErrSim  = inputFolderName+"sim.raw.dzFitErr.vsPt.60bins.root";
    fnRawData    = inputFolderName+"data.raw.dz.vsPt.60bins.root";
    fnRawErrData = inputFolderName+"data.raw.dzFitErr.vsPt.60bins.root";
    fnResp    = inputFolderName+"sim.resp.dzResp.vsPt.60bins.root";
    fnReso    = inputFolderName+"sim.reso.dzReso.vsPt.60bins.root";

    xLabel = "Track p_{T} (GeV/c)";
    yLabel = "Longit.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the LIP due to Vertex (#mum)";

    lRange = 0.;
    hRange = 210.;
    text1a = "CMS Preliminary 2010";
    text1b = "7 TeV Runs";
    leg1a  = "DATA";
    leg1b  = "MC";

    text2a = "CMS Preliminary 2010";
    text2b = "";
    leg2a  = "MC (using MC-truth)";
    leg2b  = "MC (using only reco)";

    outputName1 = outputFolderName+"resoDz_vs_pt.DATA.png";
    outputName2 = outputFolderName+"resoDz_vs_pt.MC.png";
    outputName3 = outputFolderName+"errorDz_vs_pt.png";
    outputName4 = outputFolderName+"resoErrorRatio_Dz_vs_pt.png";


    nbins = 60;
    xlow=0.7;
    xhigh=2.2;
  }


  if(type==5){
    gStyle->SetOptLogy(1);
    fnRawSim     = inputFolderName+"sim.raw.dz.vsEta.50bins.root";
    fnRawErrSim  = inputFolderName+"sim.raw.dzFitErr.vsEta.50bins.root";
    fnRawData    = inputFolderName+"data.raw.dz.vsEta.50bins.root";
    fnRawErrData = inputFolderName+"data.raw.dzFitErr.vsEta.50bins.root";
    fnResp    = inputFolderName+"sim.resp.dzResp.vsEta.50bins.root";
    fnReso    = inputFolderName+"sim.reso.dzReso.vsEta.50bins.root";

    xLabel = "Track #eta";
    yLabel = "Longit.Impact Parameter Resolution (#mum)";
    yLabelResp = "Smearing of the LIP due to Vertex (#mum)";

    lRange = 15.;
    hRange = 3000.;
    text1a = "CMS Preliminary 2010";
    text1b = "7 TeV Runs";
    leg1a  = "DATA";
    leg1b  = "MC";

    text2a = "CMS Preliminary 2010";
    text2b = "";
    leg2a  = "MC (using MC-truth)";
    leg2b  = "MC (using only reco)";

    outputName1 = outputFolderName+"resoDz_vs_eta.Pt08.DATA.png";
    outputName2 = outputFolderName+"resoDz_vs_eta.Pt08.MC.png";
    outputName3 = outputFolderName+"errorDz_vs_eta.png";
    outputName4 = outputFolderName+"resoErrorRatio_Dz_vs_eta.png";

    nbins = 50;
    xlow= -2.5;
    xhigh= 2.5;
  }


  if(type==6){
    //STILL TO BE DEFINED
  }




  TFile* fileRawSim      = new TFile(fnRawSim);
  TFile* fileRawErrSim   = new TFile(fnRawErrSim);
  TFile* fileRawData     = new TFile(fnRawData);
  TFile* fileRawErrData  = new TFile(fnRawErrData);
  TFile* fileResp    = new TFile(fnResp);
  TFile* fileReso    = new TFile(fnReso);
  
  TCanvas* canvas1 = new TCanvas("canvas1","canvas1",0,0,500,500);
  TCanvas* canvas2 = new TCanvas("canvas2","canvas2",500,0,500,500);


  TH1F* plotReso    = new TH1F("plotReso","",nbins,xlow,xhigh);
  TH1F* plotSimFit      = new TH1F("plotSimFit","",nbins,xlow,xhigh);
  TH1F* plotSimErrFit   = new TH1F("plotSimErrFit","",nbins,xlow,xhigh);
  TH1F* plotDataFit     = new TH1F("plotDataFit","",nbins,xlow,xhigh);
  TH1F* plotDataErrFit  = new TH1F("plotDataErrFit","",nbins,xlow,xhigh);

  TH1F* plotResp1  = new TH1F("plotResp1","",nbins,xlow,xhigh);
  TH1F* plotResp2  = new TH1F("plotResp1","",nbins,xlow,xhigh);
  TH1F* plotResp1B = new TH1F("plotResp1B","",nbins,xlow,xhigh);
  TH1F* plotResp2B = new TH1F("plotResp1B","",nbins,xlow,xhigh);

  double resolutions[nbins_const+1]; 
  double resolutionsUncert[nbins_const+1]; 

  double responses1[nbins_const+1];  
  double responses2[nbins_const+1];  
  double fracs[nbins_const+1];  

  double fits[nbins_const+1]; 
  double fitsUncert[nbins_const+1]; 
  double dataFits[nbins_const+1]; 
  double dataFitsUncert[nbins_const+1]; 

  double simErrFits[nbins_const+1]; 
  double simErrFitsUncert[nbins_const+1]; 
  double dataErrFits[nbins_const+1]; 
  double dataErrFitsUncert[nbins_const+1]; 

  for(unsigned int i=0; i<(nbins+1); ++i){
    resolutions[i]=0;  resolutionsUncert[i]=0;

    responses1[i]=0;
    responses2[i]=0;
    fracs[i]=0;

    fits[i]=0;         fitsUncert[i]=0;
    dataFits[i]=0;     dataFitsUncert[i]=0;    

    simErrFits[i]=0;      simErrFitsUncert[i]=0;
    dataErrFits[i]=0;     dataErrFitsUncert[i]=0;


  }

  for(unsigned int i=1; i<(nbins+1); ++i){ // loop on [1,nbins]
    stringstream stream;  stream << i;
    TString counter = stream.str();

    TString hn1,hn2,hn3,hn4,hn5,hn6;
    hn1 = hname1+counter;
    hn2 = hname2+counter;
    hn3 = hname3+counter;
    hn4 = hname4+counter;
    hn5 = hname5+counter;
    hn6 = hname6+counter;

  
    TH1* h_bef       = (TH1*) fileRawSim->Get(hn3);
    TH1* h_befData   = (TH1*) fileRawData->Get(hn4);
    TH1* h_rErrSim   = (TH1*) fileRawErrSim->Get(hn5);
    TH1* h_rErrData  = (TH1*) fileRawErrData->Get(hn6);

    TH1* h_resp     = (TH1*) fileResp->Get(hn2);
    TH1* h_reso     = (TH1*) fileReso->Get(hn1);

  
    double reso,resp1,resp2,frac,fit,dataFit;
    double resoUncert,fitUncert,dataFitUncert;
    
    double simErrFit,simErrFitUncert,dataErrFit,dataErrFitUncert;
    
    double inputResp1(initialResp1),inputResp2(initialResp2),inputFrac(initialFrac);
    
    if(i>1){
      inputResp1 = responses1[i-1];
      inputResp2 = responses2[i-1];
      inputFrac = fracs[i-1];
    }

    gStyle->SetLineColor(1);
    getResponse(h_reso,h_resp,h_bef,h_befData, //input projection histos 
		reso,resp1,resp2,frac,         //MC-reso and resp to be fitted
		fit,dataFit,                   //reso to be measured
		resoUncert,fitUncert,dataFitUncert,
		inputResp1,inputResp2,inputFrac,fitRangeMax,
		canvas1,typePlot);

    fitRawErrors(h_rErrSim,h_rErrData,
		 simErrFit, dataErrFit,
		 simErrFitUncert,dataErrFitUncert,
		 canvas1);

    resolutions[i] = reso;
    resolutionsUncert[i] = resoUncert;

    responses1[i]  = resp1;
    responses2[i]  = resp2;
    fracs[i]       = frac;

    fits[i]              = fit;
    dataFits[i]          = dataFit;
    fitsUncert[i]        = fitUncert;
    dataFitsUncert[i]    = dataFitUncert;

    simErrFits[i]         = simErrFit;
    simErrFitsUncert[i]   = simErrFitUncert;
    dataErrFits[i]        = dataErrFit;
    dataErrFitsUncert[i]  = dataErrFitUncert;
    


    // Plot the final plots as their bin content is set.
    // NB: this function also *set* the bin content of the final plots. 
    plotHistos(plotReso,plotSimFit,plotDataFit,plotResp1,plotResp2,
	       resolutions, fits, dataFits,responses1,responses2,
	       resolutionsUncert,fitsUncert,dataFitsUncert,
	       yLabel,xLabel,hRange,lRange,
	       canvas2,typePlot2);    

    
    plotRawErrHistos(plotSimErrFit, plotDataErrFit, 
		     simErrFits, dataErrFits,
		     simErrFitsUncert, dataErrFitsUncert,
		     xLabel,
		     hRange,lRange,
		     canvas2);		           
  }
  
  gStyle->SetLineColor(0);
  // --- Here are the final VertexSmearing functions
  cout << "Are you ready for the final plots?" << endl; wantMore2();
  
  canvas1->cd();  gPad->Clear();
  plotResp1->Draw("E1"); gPad->Update(); wantMore2();

  fitResponse(plotResp1,type,xlow,xhigh);

  canvas2->cd();  gPad->Clear();
  plotResp1->Draw("E1"); gPad->Update(); wantMore2();

  // --- Here are the final plot
  canvas1->cd();  gPad->Clear();
			     
  TPaveText* txt =  new TPaveText(0.50,0.75,0.75,0.85,"NDC");
  TLegend* leg   = new TLegend(0.50,0.60,0.75,0.74);
  leg->AddEntry(plotDataFit,leg1a,"lp");
  leg->AddEntry(plotSimFit,leg1b,"lp");
  txt->AddText(text1a);
  txt->AddText(text1b);

  TPaveText* txt2 =  new TPaveText(0.50,0.75,0.75,0.85,"NDC");
  TLegend* leg2   = new TLegend(0.50,0.60,0.75,0.74);
  leg2->AddEntry(plotReso,leg2a,"lp");
  leg2->AddEntry(plotSimFit,leg2b,"lp");
  txt2->AddText(text2a);

  leg->SetTextSize(0.035);  
  txt->SetTextSize(0.035);
  leg2->SetTextSize(0.035);  
  txt2->SetTextSize(0.035);

  // --- plots DATA vs MC
  plotReso->SetMarkerStyle(21);      plotReso->SetMarkerColor(2);
  plotDataFit->SetMarkerStyle(20);   plotDataFit->SetMarkerColor(1);

  plotReso->Draw("E1");
  plotDataFit->Draw("sameE1");
  gPad->Update();

  leg->Draw();
  txt->Draw();

  finalFit(plotReso,plotDataFit,xLabel,yLabel,txt,outputFolderName+"resoD0_vs_pt.DATA.fit.png");

  gPad->Print(outputName1);
  wantMore2();


  // --- plots SIM vs MC-Truth-based
  canvas2->cd();  gPad->Clear();

  plotReso->SetMarkerStyle(21);      plotReso->SetMarkerColor(2);
  plotSimFit->SetMarkerStyle(20);    plotSimFit->SetMarkerColor(1);

  plotReso->Draw("E1");
  plotSimFit->Draw("sameE1");
  gPad->Update();
 
  leg2->Draw();
  txt2->Draw();
  gPad->Update();

  // finalFit(plotSimFit,plotReso);

  gPad->Print(outputName2);
  wantMore2();


  // --- plots DATA vs Sim, <error>
  canvas1->cd();  gPad->Clear();
  plotSimErrFit->SetMarkerStyle(21);      plotSimErrFit->SetMarkerColor(2);
  plotDataErrFit->SetMarkerStyle(20);     plotDataErrFit->SetMarkerColor(1);

  TPaveText* txt3 =  new TPaveText(0.50,0.75,0.75,0.85,"NDC");
  TLegend* leg3   = new TLegend(0.50,0.60,0.75,0.74);
  leg3->AddEntry(plotSimErrFit,leg1b,"lp");
  leg3->AddEntry(plotDataErrFit,leg1a,"lp");
  txt3->AddText(text1b);
  txt3->AddText(text1a);
  leg3->SetTextSize(0.035);  
  txt3->SetTextSize(0.035);

  plotSimErrFit->Draw("E1");
  plotDataErrFit->Draw("sameE1");
  gPad->Update();

  leg3->Draw();
  txt3->Draw();

  gPad->Print(outputName3);
  wantMore2();


  // --- plots DATA vs Sim, <resol>/<error>
  canvas2->cd();  gPad->Clear();

  canvas2->SetLogy(0);
  plotSimFit->Divide(plotSimErrFit);
  plotDataFit->Divide(plotDataErrFit);

  plotSimFit->SetMaximum(2.);
  plotSimFit->SetMinimum(0.);
  plotSimFit->SetYTitle(" resolution / < error > ");

  plotSimFit->SetMarkerStyle(21);      plotSimFit->SetMarkerColor(2);
  plotDataFit->SetMarkerStyle(20);     plotDataFit->SetMarkerColor(1);
  
  TPaveText* txt4 =  new TPaveText(0.50,0.75,0.75,0.85,"NDC");
  TLegend* leg4   = new TLegend(0.50,0.60,0.75,0.74);
  leg4->AddEntry(plotSimFit,leg1b,"lp");
  leg4->AddEntry(plotDataFit,leg1a,"lp");
  txt4->AddText(text1b);
  txt4->AddText(text1a);
  leg4->SetTextSize(0.035);  
  txt4->SetTextSize(0.035);

  plotSimFit->Draw("E1");
  plotDataFit->Draw("sameE1");
  gPad->Update();

  leg4->Draw();
  txt4->Draw();
  gPad->Update();

  gPad->Print(outputName4);


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


double getFitRange(TH1* h){
  h->GetListOfFunctions()->Delete();
  h->Fit("gaus","Q");gPad->Clear(); 
  TF1* firstFit = (TF1*) h->GetListOfFunctions()->Last();
  //firstFit->Print();

  double firstSigma = firstFit->GetParameter(2);
  //cout << "firstSigma: " << firstSigma << endl;
  //h->Draw(); gPad->Update();
  //wantMore2();

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
  plotReso->SetYTitle(yLabel);
  plotReso->SetXTitle(xLabel);

  plotResp1->SetYTitle(yLabelResp);
  plotResp1->SetXTitle(xLabel);

  plotDataFit->SetYTitle(yLabel);
  plotDataFit->SetXTitle(xLabel);

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
  

  // --- set style
  plotReso->SetMarkerStyle(20);
  plotReso->SetMarkerColor(1);

  plotResp1->SetMarkerStyle(20);
  plotResp1->SetMarkerColor(3);

  plotSimFit->SetMarkerStyle(21);
  plotSimFit->SetMarkerColor(2);

  plotDataFit->SetMarkerStyle(20);
  plotDataFit->SetMarkerColor(1);

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


void fitResponse(TH1* h_resp, int type, double xlow, double xhigh){
  TF1* myResp1Fit;

  if(type==1 || type==2 || type==4)
    myResp1Fit = new TF1("myResp1Fit","[0]+[1]*x",xlow,xhigh);
  
  if(type==5)
    myResp1Fit = new TF1("myResp1Fit","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",xlow,xhigh);
  
  myResp1Fit->SetParName(0,"q");
  myResp1Fit->SetParName(1,"p");
  myResp1Fit->SetParName(2,"p2");
  myResp1Fit->SetParName(3,"p3");
  myResp1Fit->SetParName(4,"p4");
  myResp1Fit->SetLineWidth(2);
  myResp1Fit->SetLineColor(2);

  h_resp->Fit("myResp1Fit");   
  
  int nbins = h_resp->GetNbinsX();
  for(unsigned int i=1; i!=nbins+1; ++i){
    double xvalue = h_resp->GetBinLowEdge(i)+0.5*h_resp->GetBinWidth(i); 
    double tmp = myResp1Fit->Eval(xvalue);
    //cout << "x,xvalue: " << xvalue << " , " << tmp << endl;
    h_resp->SetBinContent(i,tmp);
  }
  
}


void fitRawErrors(TH1* h_rErrSim,TH1* h_rErrData,
		  double& rErrSim, double& rErrData,
		  double& rErrSimUncert, double& rErrDataUncert,
		  TCanvas* canv){

  /*
  rErrSim  = h_rErrSim->GetMean();   rErrSimUncert  = h_rErrSim->GetMeanError();
  rErrData = h_rErrData->GetMean();  rErrDataUncert = h_rErrData->GetMeanError();
  */

  double range;
  // --- section about SIM
  h_rErrSim->Fit("gaus","Q"); gPad->Clear();
  TF1* firstFit = (TF1*) h_rErrSim->GetListOfFunctions()->Last();
  double firstMean = firstFit->GetParameter(1);
  double firstSigma = firstFit->GetParameter(2);
  cout << "firstMean,firstSigma: " << firstMean << " , " << firstSigma << endl;
  

  RooRealVar x1("x1","rawError Sim",firstMean-firstSigma*2.0,firstMean+firstSigma);
  RooDataHist datah("datah","dataset from h",x1,h_rErrSim) ;

  RooPlot* frame1 = x1.frame(); 
  datah.plotOn(frame1);

  RooRealVar mean1("mean1","mean of gaussian1",firstMean,0,2000) ;
  RooRealVar sigma1("sigma1","width of gaussian1",firstSigma,0,300) ;
  RooGaussian gauss1("gauss1","gaussian PDF 1",x1,mean1,sigma1) ;  

  gauss1.fitTo(datah);
  gauss1.plotOn(frame1);
  rErrSim = mean1.getVal();
  rErrSimUncert = mean1.getError();

  //frame1->Draw(); gPad->Update(); 
  //wantMore2();

  // --- section about DATA
  RooRealVar x2("x2","rawError DATA",firstMean-firstSigma*2.0,firstMean+firstSigma);
  RooDataHist datah2("datah2","dataset from h2",x2,h_rErrData) ;

  RooPlot* frame2 = x2.frame(); 
  datah2.plotOn(frame2);

  RooRealVar mean2("mean2","mean of gaussian2",firstMean,0,2000) ;
  RooRealVar sigma2("sigma2","width of gaussian2",firstSigma,0,300) ;
  RooGaussian gauss2("gauss2","gaussian PDF 2",x2,mean2,sigma2) ;  

  gauss2.fitTo(datah2);
  gauss2.plotOn(frame2);
  rErrData = mean2.getVal();
  rErrDataUncert = mean2.getError();

  frame2->Draw(); gPad->Update(); 
  wantMore2();
}


void plotRawErrHistos(TH1* plotErrSim, TH1* plotErrData,
		      double* errSim, double* errData,
		      double* errSimUncert, double* errDataUncert,
		      TString xLabel,
		      double hRange, double lRange,
		      TCanvas* canv){
  // --- set bin contents
  plotErrSim->SetContent(errSim); plotErrSim->SetError(errSimUncert);
  plotErrData->SetContent(errData); plotErrData->SetError(errDataUncert);


  // --- set labels
  //plotErrSim->SetTitle("");
  plotErrSim->SetYTitle("Impact Parameter Error (#mum)");
  plotErrSim->SetXTitle(xLabel);

  
  // --- set ranges
  plotErrSim->SetMaximum(hRange);
  plotErrSim->SetMinimum(lRange);

  
  // --- set style
  plotErrSim->SetMarkerStyle(21);
  plotErrSim->SetMarkerColor(2);
  
  plotErrData->SetMarkerStyle(20);
  plotErrData->SetMarkerColor(1);


  canv->cd();
  plotErrSim->Draw("E1"); gPad->Update();
  plotErrData->Draw("sameE1"); gPad->Update();
  //cout << "before exiting plotRawErrHistos()" << endl; 
  //cout << "hRange,lRange: " << hRange << " , " << lRange << endl;
  //wantMore2();
}

void finalFit(TH1F* histoToFit1, TH1F* histoToFit2, TString xLabel, TString yLabel, TPaveText *txt, TString fitOutputfile){
  RooRealVar pt("pt","pt",0.7,2.2);
  pt.setBins(60);

  RooDataHist _histoToFit1("_histoToFit1","_histoToFit1",RooArgList(pt),histoToFit1);
  RooDataHist _histoToFit2("_histoToFit2","_histoToFit2",RooArgList(pt),histoToFit2);

  RooRealVar a("a","constant",0.,1.);
  RooRealVar b("b","inverse",0.,1.);

  RooFormulaVar t("t","1/@0",RooArgList(pt));
  RooPolynomial resolFit("resolFit","resolFit",t,RooArgList(a,b));
  resolFit.fitTo(_histoToFit1);

  TString fitRes_a1, fitRes_b1, fitResErr_a1, fitResErr_b1;
 
  stringstream stream;  
  stream << setprecision(3) << a.getVal();
  fitRes_a1 = stream.str();
  stream.str("");

  stream << a.getError();
  fitResErr_a1 = stream.str();
  stream.str("");

  stream << b.getVal();
  fitRes_b1 = stream.str();
  stream.str("");

  stream << b.getError();
  fitResErr_b1 = stream.str();
  stream.str("");

  TString fitResult1 = fitRes_a1+" #pm "+fitResErr_a1+" + #frac{"+fitRes_b1+" #pm "+fitResErr_b1+"}{p_{T}}";

  RooPlot *xplot = pt.frame();
  _histoToFit1.plotOn(xplot,MarkerColor(2),MarkerStyle(21));
  resolFit.plotOn(xplot);
  xplot->getAttLine()->SetLineColor(2);
  xplot->SetMaximum(200.);
  xplot->SetTitle("");

  resolFit.fitTo(_histoToFit2);
  _histoToFit2.plotOn(xplot,MarkerColor(1),MarkerStyle(20));
  resolFit.plotOn(xplot);
  xplot->getAttLine()->SetLineColor(1);

  TCanvas* fitCanvas = new TCanvas("fitCanvas","fitCanvas",0,0,500,500);
  fitCanvas->cd();

  xplot->GetXaxis()->SetTitle(xLabel);
  xplot->GetYaxis()->SetTitle(yLabel);
  xplot->Draw();
  txt->Draw();

  TString fitRes_a2, fitRes_b2, fitResErr_a2, fitResErr_b2;
 
  stream  << a.getVal();
  fitRes_a2 = stream.str();
  stream.str("");
  
  stream << b.getVal();
  fitRes_b2 = stream.str();
  stream.str("");
  
  stream << a.getError();
  fitResErr_a2 = stream.str();
  stream.str("");

  stream << b.getError();
  fitResErr_b2 = stream.str();
  stream.str("");

  TString fitResult2 = fitRes_a2+" #pm "+fitResErr_a2+" + #frac{"+fitRes_b2+" #pm "+fitResErr_b2+"}{p_{T}}";

  TLegend *leg2 = new TLegend(0.2056452,0.1652542,0.6028226,0.3580508,NULL,"brNDC");
  leg2->AddEntry(histoToFit1,"MC  "+fitResult1,"lp");
  leg2->AddEntry(histoToFit2,"DATA  "+fitResult2,"lp");
  leg2->Draw();
 
  fitCanvas->SaveAs(fitOutputfile);

}
