#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TF1.h"
#include "TH2D.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "RooHZZStyle.C"
#include "contours.cxx"

#include <sstream>
#include <iostream>

using namespace std;

TLatex *CMSPreliminary(float lumi7TeV=5.1, float lumi8TeV=19.8) {

  stringstream line;
  line << "CMS Preliminary  #sqrt{s}=7 TeV, L=" << lumi7TeV << " fb^{-1}  #sqrt{s}=8 TeV, L=" << lumi8TeV << " fb^{-1}";
  TLatex* CP = new TLatex(0.15,0.96, line.str().c_str());
  CP->SetNDC(kTRUE);
  CP->SetTextSize(0.032);

  return CP;

}

TLatex *CMS(float lumi7TeV=5.1, float lumi8TeV=19.8) {

  stringstream line;
  line << "CMS                   #sqrt{s}=7 TeV, L=" << lumi7TeV << " fb^{-1}  #sqrt{s}=8 TeV, L=" << lumi8TeV << " fb^{-1}";
  TLatex* CP = new TLatex(0.15,0.96, line.str().c_str());
  CP->SetNDC(kTRUE);
  CP->SetTextSize(0.032);

  return CP;

}

TPaveText *text(const char *txt, float x1, float y1, float x2, float y2) {
  TPaveText *text = new TPaveText(x1,y1,x2,y2,"brNDC");
  text->AddText(txt);
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(42);
  text->SetTextSize(0.05);
  return text;
}

float getMode(TH1F *h) {
  float contMax=0.;
  int bin=0;
  for(int b=1;b<h->GetNbinsX();++b) {
    float cont=h->GetBinContent(b);
    if(cont > contMax) {
      contMax=cont;
      bin=b;
    }
  }
  return h->GetBinCenter(bin);
}

void analyzeMassToys(int fitdimensions, int channel) {

  string chstr;
  if(channel==0) chstr = "4mu";
  if(channel==1) chstr = "4e";
  if(channel==2) chstr = "2e2mu";
  if(channel==3) chstr = "comb";

  // fit limits
  float massMin=120;
  float massMax=130;

  stringstream oss;
  if(fitdimensions==1) oss << "histos1D_" << chstr << ".root";
  else if(fitdimensions==2) oss << "histos2D_" << chstr << ".root";
  else if(fitdimensions==3) oss << "histos3D_" << chstr << ".root";
  
  TFile *histos = TFile::Open(oss.str().c_str(),"recreate");
        
  stringstream fss;
  if(fitdimensions==1) fss << "results/c1DNoMassErr/toys_" << chstr << "_new.root";
  else if(fitdimensions==2) fss << "results/c1DMassErr/toys_" << chstr << "_new.root";
  else if(fitdimensions==3) fss << "results/c2DMassErr/toys_" << chstr << "_new.root";
  else {
    cout << "fitdimensions should be 1,2,3" << endl;
    return;
  }

  TFile *tfile = TFile::Open(fss.str().c_str());
  TTree *toys = (TTree*)tfile->Get("limit");

  float errMax;
  switch (channel) {
  case 0: errMax=1.5; break;
  case 1: errMax=4.0; break;
  case 2: errMax=2.0; break;
  case 3: errMax=1.5; break;
  default: errMax=10; break;
  }
  TH1F *massh = new TH1F("massh","",25,123.,129.);
  TH1F *pullh = new TH1F("pullh","",25,-5.,5.);
  TH1F *errh = new TH1F("errh","",30,0.2,errMax);
  TH1F *errl = new TH1F("errl","",30,0.2,errMax);

  errh->GetXaxis()->SetTitle("mass error (GeV)");
  errl->GetXaxis()->SetTitle("mass error (GeV)");
  pullh->GetXaxis()->SetTitle("pull");

  errh->GetYaxis()->SetTitle("n. toys");
  errl->GetYaxis()->SetTitle("n. toys");
  pullh->GetYaxis()->SetTitle("n. toys");

  double genmh;
  float MH;
  float quantileExpected;
  toys->SetBranchAddress("mh", &genmh);
  toys->SetBranchAddress("MH", &MH);
  toys->SetBranchAddress("quantileExpected", &quantileExpected);

  float mass, lowere, uppere;
  mass=lowere=uppere=0;
  for(int i=0; i<(int)toys->GetEntries(); ++i) {
    toys->GetEntry(i);
    // this to distinguish -1 sigma / +1 sigma
    if(quantileExpected==1) {
      massh->Fill(MH);
      mass = MH;
    } else {
      if((i-1)%4==0) {
        lowere = mass-MH;
	if(mass!=massMin && mass!=massMax && lowere>0.05) errl->Fill(lowere);
      }
      if((i-2)%4==0) {
        uppere = MH-mass;
	if(mass!=massMin && mass!=massMax && uppere>0.05) errh->Fill(uppere);
      }
      if((i-3)%4==0) {
        float avgerr = (lowere + uppere)/2.;
	if(mass!=massMin && mass!=massMax && lowere>0.05 && uppere>0.05) pullh->Fill((mass-genmh)/avgerr);
      }
    }
  }

  histos->cd();
  massh->Write();
  pullh->Write();
  errh->Write();
  errl->Write();
  histos->Close();

}

void finalizeMassToys(int channel) {

  // data fits (uppere)
  //  float dm1h=0.572;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  string chstr;
  if(channel==0) chstr = "4mu";
  if(channel==1) chstr = "4e";
  if(channel==2) chstr = "2e2mu";
  if(channel==3) chstr = "comb";

  stringstream oss1D, oss2D, oss3D;
  oss1D << "histos1D_" << chstr << ".root";
  oss2D << "histos2D_" << chstr << ".root";
  oss3D << "histos3D_" << chstr << ".root";

  TFile *fit1D = TFile::Open(oss1D.str().c_str());
  TH1F *pull1D = (TH1F*)fit1D->Get("pullh");
  TH1F *errh1D = (TH1F*)fit1D->Get("errh");

  TFile *fit2D = TFile::Open(oss2D.str().c_str());
  TH1F *pull2D = (TH1F*)fit2D->Get("pullh");
  TH1F *errh2D = (TH1F*)fit2D->Get("errh");

  TFile *fit3D = TFile::Open(oss3D.str().c_str());
  TH1F *pull3D = (TH1F*)fit3D->Get("pullh");
  TH1F *errh3D = (TH1F*)fit3D->Get("errh");

  float m1=getMode(errh1D);
  float m2=getMode(errh2D);
  float m3=getMode(errh3D);

  stringstream m1s,m2s,m3s;
  m1s.precision(2); m2s.precision(2); m3s.precision(2);
  m1s << std::fixed << "m_{4l}. Mode = " << m1 << " GeV";
  m2s << std::fixed << "m_{4l}, #sigma_{m}. Mode = " << m2 << " GeV";
  m3s << std::fixed << "m_{4l}, #sigma_{m}, KD. Mode = " << m3 << " GeV";

  std::vector<TH1F*> pulls,errs;
  pulls.push_back(pull1D);
  pulls.push_back(pull2D);
  pulls.push_back(pull3D);
  errs.push_back(errh1D);
  errs.push_back(errh2D);
  errs.push_back(errh3D);

  
  TPaveText *text = new TPaveText(0.50,0.50,0.85,0.60,"brNDC");
  text->AddText(m1s.str().c_str());
  text->AddText(m2s.str().c_str());
  text->AddText(m3s.str().c_str());
  text->SetBorderSize(0);
  text->SetFillStyle(0);
  text->SetTextAlign(12);
  text->SetTextFont(42);
  text->SetTextSize(0.03);


  // draw the legend
  TLegend *legend = new TLegend(0.50,0.65,0.85,0.85,NULL,"brNDC");
  legend->SetBorderSize(     0);
  legend->SetFillColor (     0);
  legend->SetTextAlign (    12);
  legend->SetTextFont  (    42);
  legend->SetTextSize  (0.05);
  legend->AddEntry(errh1D, "m_{4l}","f");
  legend->AddEntry(errh2D, "m_{4l},#sigma_{m}","f");
  legend->AddEntry(errh3D, "m_{4l},#sigma_{m},KD","f");

  errs[0]->SetLineColor(kAzure+5);
  errs[0]->SetFillColor(kAzure+5);
  errs[1]->SetLineColor(kOrange-3);
  errs[1]->SetFillColor(kOrange-3);
  errs[2]->SetLineColor(kGreen+3);
  errs[2]->SetFillColor(kGreen+3);

  errs[1]->SetFillStyle(3001);
  errs[2]->SetFillStyle(3001);

  TCanvas *c1 = new TCanvas("c1","",750,750);

  float yMax;
  switch (channel) {
  case 0: yMax=400; break;
  case 1: yMax=500; break;
  case 2: yMax=100; break;
  case 3: yMax=400; break;
  default: yMax=100; break;
  }

  errs[0]->SetMaximum(yMax);
  for(int i=0;i<3;++i) {
    if(i==0) errs[i]->Draw();
    else {
      errs[i]->Scale(errs[0]->Integral()/errs[i]->Integral());
      errs[i]->Draw("same");
    }
  }
  legend->Draw();
  text->Draw();
  stringstream errc;
  errc << "errors_ch" << channel << ".pdf";
  c1->SaveAs(errc.str().c_str());

  for(int i=0;i<3;++i) {
    pulls[i]->Draw();
    pulls[i]->Fit("gaus");
    stringstream pullpdf,pullpng;
    pullpdf << "pull_" << i+1 << "D_ch" << channel << ".pdf";
    pullpng << "pull_" << i+1 << "D_ch" << channel << ".png";
    c1->SaveAs(pullpdf.str().c_str());
    c1->SaveAs(pullpng.str().c_str());
  }  
 
}

void drawAllToys() {

  for(int ch=0;ch<4;++ch) {
    for(int d=1;d<4;++d) {
      analyzeMassToys(d,ch);
    }
    finalizeMassToys(ch);
  }
}

void plotScanByChannel(int ndim, bool fullerror) {

  float xmin=122;
  float xmax=132;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  
  TStyle *mystyle = RooHZZStyle("ZZ");
  mystyle->cd();

  stringstream fss;
  if(ndim==1) fss << "results/c1DNoMassErr/scanMass-MH_";
  else if(ndim==2) fss << "results/c1DMassErr/scanMass-MH_";
  else if(ndim==3) fss << "results/c2DMassErr/scanMass-MH_";
  else {
    cout << "ndim should be 1,2,3" << endl;
    return;
  }

  stringstream file4mu,file4e,file2e2mu,filecomb;
  if(fullerror) {
    file4mu    << fss.str() << "4mu_new.root";
    file4e     << fss.str() << "4e_new.root";
    file2e2mu  << fss.str() << "2e2mu_new.root";
    filecomb   << fss.str() << "comb_new.root";
  } else {
    file4mu    << fss.str() << "4mu_nosyst.root";
    file4e     << fss.str() << "4e_nosyst.root";
    file2e2mu  << fss.str() << "2e2mu_nosyst.root";
    filecomb   << fss.str() << "comb_nosyst.root";
  }

  TFile *fit4mu = TFile::Open(file4mu.str().c_str());
  TTree *tree4mu = (TTree*)fit4mu->Get("limit");

  TFile *fit4e = TFile::Open(file4e.str().c_str());
  TTree *tree4e = (TTree*)fit4e->Get("limit");

  TFile *fit2e2mu = TFile::Open(file2e2mu.str().c_str());
  TTree *tree2e2mu = (TTree*)fit2e2mu->Get("limit");

  TFile *fitcomb = TFile::Open(filecomb.str().c_str());
  TTree *treecomb = (TTree*)fitcomb->Get("limit");

  TGraph *g4mu = new TGraph();
  TGraph *g4e = new TGraph();
  TGraph *g2e2mu = new TGraph();
  TGraph *gcomb = new TGraph();

  vector<TTree*> trees;
  trees.push_back(tree4mu);
  trees.push_back(tree4e);
  trees.push_back(tree2e2mu);
  trees.push_back(treecomb);

  vector<TGraph*> graphs;
  graphs.push_back(g4mu);
  graphs.push_back(g4e);
  graphs.push_back(g2e2mu);
  graphs.push_back(gcomb);

  for(int cha=0; cha<(int)trees.size(); ++cha) {
    
    cout << "Analyzing scan for channel = " << cha << endl;

    float MH;
    float deltaNLL;
    trees[cha]->SetBranchAddress("MH", &MH);
    trees[cha]->SetBranchAddress("deltaNLL", &deltaNLL);

    if(cha<3) graphs[cha]->SetLineWidth(2);
    else graphs[cha]->SetLineWidth(4);

    for(int i=1; i<(int)trees[cha]->GetEntries();++i) {
       trees[cha]->GetEntry(i);
       graphs[cha]->SetPoint(i-1,MH,2*deltaNLL);
    }
    graphs[cha]->Sort();
  }

  graphs[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  graphs[0]->GetYaxis()->SetRangeUser(0,10);
  graphs[0]->GetXaxis()->SetTitle("m_{X} (GeV)");
  graphs[0]->GetYaxis()->SetTitle("-2 #Delta lnL");

  graphs[0]->SetMarkerColor(kRed);
  graphs[1]->SetMarkerColor(kGreen+2);
  graphs[2]->SetMarkerColor(kBlue);
  graphs[3]->SetMarkerColor(kBlack);

  graphs[0]->SetLineColor(kRed);
  graphs[1]->SetLineColor(kGreen+2);
  graphs[2]->SetLineColor(kBlue);
  graphs[3]->SetLineColor(kBlack);

  if(!fullerror) 
    for(int i=0;i<4;++i) graphs[i]->SetLineStyle(kDashed);

  if(fullerror) graphs[0]->Draw("l");
  else  graphs[0]->Draw("al");
  for(int i=1;i<4;++i) graphs[i]->Draw("l");

  if(fullerror) {
    // draw the legend
    TLegend *legend = new TLegend(0.55,0.65,0.85,0.90,NULL,"brNDC");
    legend->SetBorderSize(     0);
    legend->SetFillColor (     0);
    legend->SetTextAlign (    12);
    legend->SetTextFont  (    42);
    legend->SetTextSize  (0.05);
    legend->AddEntry(graphs[3], "Combined","l");
    legend->AddEntry(graphs[1], "H #rightarrow ZZ #rightarrow 4e ","l");
    legend->AddEntry(graphs[0], "H #rightarrow ZZ #rightarrow 4#mu ","l");
    legend->AddEntry(graphs[2], "H #rightarrow ZZ #rightarrow 2e2#mu ","l");
    legend->Draw();
  }

}

void plotScanByChannel(int ndim) {

  TCanvas *c1 = new TCanvas("c1","",750,750);

  plotScanByChannel(ndim,false);
  plotScanByChannel(ndim,true);

  TLine *line1 = new TLine(122,1,132,1);
  line1->SetLineColor(kRed);
  line1->SetLineWidth(3.0);
  line1->Draw();

  TLine *line2 = new TLine(122,4,132,4);
  line2->SetLineColor(kRed);
  line2->SetLineWidth(1.5);
  line2->Draw();

  TLatex *CP = CMS();
  CP->Draw();

  TPaveText *comment;
  if(ndim==1) comment = text("H #rightarrow ZZ, 1D",0.20,0.90,0.40,0.90);
  if(ndim==2) comment = text("H #rightarrow ZZ, 2D",0.20,0.90,0.40,0.90);
  if(ndim==3) comment = text("H #rightarrow ZZ, 3D",0.20,0.90,0.40,0.90);
  comment->Draw();

  stringstream outnamepdf;
  outnamepdf << "scanmass_bychannel_" << ndim << "D.pdf";
  c1->SaveAs(outnamepdf.str().c_str());
  stringstream outnamepng;
  outnamepng << "scanmass_bychannel_" << ndim << "D.png";
  c1->SaveAs(outnamepng.str().c_str());


}

void plotScanByDim(int channel, bool fullerror) {

  float xmin=122;
  float xmax=132;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  
  TStyle *mystyle = RooHZZStyle("ZZ");
  mystyle->cd();

  string channelstr;
  if (channel==0) channelstr = "4mu";
  if (channel==1) channelstr = "4e";
  if (channel==2) channelstr = "2e2mu";
  if (channel==3) channelstr = "comb";

  stringstream file1D,file2D,file3D;
  if(fullerror) {
    file1D << "results/c1DNoMassErr/scanMass-MH_" << channelstr << "_new.root";
    file2D << "results/c1DMassErr/scanMass-MH_" << channelstr << "_new.root";
    file3D << "results/c2DMassErr/scanMass-MH_" << channelstr << "_new.root";
  } else {
    file1D << "results/c1DNoMassErr/scanMass-MH_" << channelstr << "_nosyst.root";
    file2D << "results/c1DMassErr/scanMass-MH_" << channelstr << "_nosyst.root";
    file3D << "results/c2DMassErr/scanMass-MH_" << channelstr << "_nosyst.root";
  }

  TFile *fit1D = TFile::Open(file1D.str().c_str());
  TTree *tree1D = (TTree*)fit1D->Get("limit");

  TFile *fit2D = TFile::Open(file2D.str().c_str());
  TTree *tree2D = (TTree*)fit2D->Get("limit");

  TFile *fit3D = TFile::Open(file3D.str().c_str());
  TTree *tree3D = (TTree*)fit3D->Get("limit");

  TGraph *g1D = new TGraph();
  TGraph *g2D = new TGraph();
  TGraph *g3D = new TGraph();

  vector<TTree*> trees;
  trees.push_back(tree1D);
  trees.push_back(tree2D);
  trees.push_back(tree3D);

  vector<TGraph*> graphs;
  graphs.push_back(g1D);
  graphs.push_back(g2D);
  graphs.push_back(g3D);

  for(int dim=0; dim<(int)trees.size(); ++dim) {
    
    cout << "Analyzing scan for n. dim. fit = " << dim << endl;

    float MH;
    float deltaNLL;
    trees[dim]->SetBranchAddress("MH", &MH);
    trees[dim]->SetBranchAddress("deltaNLL", &deltaNLL);

    graphs[dim]->SetLineWidth(2);

    for(int i=1; i<(int)trees[dim]->GetEntries();++i) {
       trees[dim]->GetEntry(i);
       graphs[dim]->SetPoint(i-1,MH,2*deltaNLL);
    }
    graphs[dim]->Sort();
  }

  graphs[0]->SetMarkerColor(kRed);
  graphs[1]->SetMarkerColor(kGreen+2);
  graphs[2]->SetMarkerColor(kBlue);

  graphs[0]->SetLineColor(kRed);
  graphs[1]->SetLineColor(kGreen+2);
  graphs[2]->SetLineColor(kBlue);

  if(!fullerror) 
    for(int i=0;i<3;++i) graphs[i]->SetLineStyle(kDashed);

  graphs[0]->GetXaxis()->SetRangeUser(xmin,xmax);
  graphs[0]->GetYaxis()->SetRangeUser(0,10);
  graphs[0]->GetXaxis()->SetTitle("m_{X} (GeV)");
  graphs[0]->GetYaxis()->SetTitle("-2 #Delta lnL");

  if(fullerror) graphs[0]->Draw("l");
  else graphs[0]->Draw("al");
  for(int i=1;i<3;++i) graphs[i]->Draw("l");

  if(fullerror) {
    // draw the legend
    TLegend *legend = new TLegend(0.75,0.65,0.90,0.90,NULL,"brNDC");
    legend->SetBorderSize(     0);
    legend->SetFillColor (     0);
    legend->SetTextAlign (    12);
    legend->SetTextFont  (    42);
    legend->SetTextSize  (0.05);
    legend->AddEntry(graphs[0], "1D fit","l");
    legend->AddEntry(graphs[1], "2D fit","l");
    legend->AddEntry(graphs[2], "3D fit","l");
    legend->Draw();
  }

}

void plotScanByDim(int channel) {

  TCanvas *c1 = new TCanvas("c1","",750,750);

  plotScanByDim(channel,false);  
  plotScanByDim(channel,true);  

  TLine *line1 = new TLine(122,1,132,1);
  line1->SetLineColor(kRed);
  line1->SetLineWidth(3.0);
  line1->Draw();

  TLine *line2 = new TLine(122,4,132,4);
  line2->SetLineColor(kRed);
  line2->SetLineWidth(1.5);
  line2->Draw();

  TLatex *CP = CMS();
  CP->Draw();

  TPaveText *comment;
  if(channel==0) comment = text("H #rightarrow ZZ #rightarrow 4#mu",0.20,0.90,0.40,0.90);
  if(channel==1) comment = text("H #rightarrow ZZ #rightarrow 4e",0.20,0.90,0.40,0.90);
  if(channel==2) comment = text("H #rightarrow ZZ #rightarrow 2e2#mu",0.20,0.90,0.40,0.90);
  if(channel==3) comment = text("H #rightarrow ZZ #rightarrow 4l",0.20,0.90,0.40,0.90);
  comment->Draw();

  string channelstr;
  if (channel==0) channelstr = "4mu";
  if (channel==1) channelstr = "4e";
  if (channel==2) channelstr = "2e2mu";
  if (channel==3) channelstr = "comb";

  stringstream outnamepdf;
  outnamepdf << "scanmass_bydim_" << channelstr << ".pdf";
  c1->SaveAs(outnamepdf.str().c_str());

  stringstream outnamepng;
  outnamepng << "scanmass_bydim_" << channelstr << ".png";
  c1->SaveAs(outnamepng.str().c_str());
  
}
    
void cccPlot(int ndim) {

  float fitval[8], fiterrl[8], fiterrh[8];
  float fitstaterrl[8], fitstaterrh[8], fitsysterrl[8], fitsysterrh[8];
  for(int cha=0; cha<8; ++cha) {
    fiterrl[cha]=0;
    fiterrh[cha]=0;
    fitstaterrl[cha]=0;
    fitstaterrh[cha]=0;
    fitsysterrl[cha]=0;
    fitsysterrh[cha]=0;
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  
  TStyle *mystyle = RooHZZStyle("ZZ");
  mystyle->cd();

  stringstream fss;
  if(ndim==1) fss << "results/c1DNoMassErr/scanMass-MH_";
  else if(ndim==2) fss << "results/c1DMassErr/scanMass-MH_";
  else if(ndim==3) fss << "results/c2DMassErr/scanMass-MH_";
  else {
    cout << "ndim should be 1,2,3" << endl;
    return;
  }

  stringstream file4mu,file4e,file2e2mu,filecomb;
  file4mu    << fss.str() << "4mu_new.root";
  file4e     << fss.str() << "4e_new.root";
  file2e2mu  << fss.str() << "2e2mu_new.root";
  filecomb   << fss.str() << "comb_new.root";

  stringstream file4mu_nosyst,file4e_nosyst,file2e2mu_nosyst,filecomb_nosyst;
  file4mu_nosyst    << fss.str() << "4mu_nosyst.root";
  file4e_nosyst     << fss.str() << "4e_nosyst.root";
  file2e2mu_nosyst  << fss.str() << "2e2mu_nosyst.root";
  filecomb_nosyst   << fss.str() << "comb_nosyst.root";

  TFile *fit4mu = TFile::Open(file4mu.str().c_str());
  TTree *tree4mu = (TTree*)fit4mu->Get("limit");

  TFile *fit4e = TFile::Open(file4e.str().c_str());
  TTree *tree4e = (TTree*)fit4e->Get("limit");

  TFile *fit2e2mu = TFile::Open(file2e2mu.str().c_str());
  TTree *tree2e2mu = (TTree*)fit2e2mu->Get("limit");

  TFile *fitcomb = TFile::Open(filecomb.str().c_str());
  TTree *treecomb = (TTree*)fitcomb->Get("limit");


  TFile *fit4mu_nosyst = TFile::Open(file4mu_nosyst.str().c_str());
  TTree *tree4mu_nosyst = (TTree*)fit4mu_nosyst->Get("limit");

  TFile *fit4e_nosyst = TFile::Open(file4e_nosyst.str().c_str());
  TTree *tree4e_nosyst = (TTree*)fit4e_nosyst->Get("limit");

  TFile *fit2e2mu_nosyst = TFile::Open(file2e2mu_nosyst.str().c_str());
  TTree *tree2e2mu_nosyst = (TTree*)fit2e2mu_nosyst->Get("limit");

  TFile *fitcomb_nosyst = TFile::Open(filecomb_nosyst.str().c_str());
  TTree *treecomb_nosyst = (TTree*)fitcomb_nosyst->Get("limit");

  vector<TTree*> trees;
  trees.push_back(tree4mu);
  trees.push_back(tree4e);
  trees.push_back(tree2e2mu);
  trees.push_back(treecomb);
  trees.push_back(tree4mu_nosyst);
  trees.push_back(tree4e_nosyst);
  trees.push_back(tree2e2mu_nosyst);
  trees.push_back(treecomb_nosyst);

  for(int cha=0; cha<(int)trees.size(); ++cha) {
    
    cout << "Analyzing scan for channel = " << cha << endl;

    float MH;
    float deltaNLL;
    trees[cha]->SetBranchAddress("MH", &MH);
    trees[cha]->SetBranchAddress("deltaNLL", &deltaNLL);
    
    float dNLLMinus=1000; float dNLLPlus=1000;
    for(int i=0; i<(int)trees[cha]->GetEntries();++i) {
      trees[cha]->GetEntry(i);
      if(i==0) fitval[cha]=MH;
      else {
	if(fabs(2*deltaNLL-1)<dNLLMinus && MH<fitval[cha]) {
	  if(cha<4) fiterrl[cha]=MH;
	  else fitstaterrl[cha]=MH;
	  dNLLMinus=fabs(2*deltaNLL-1);
	}
	if(fabs(2*deltaNLL-1)<dNLLPlus && MH>fitval[cha]) {
	  if(cha<4) fiterrh[cha]=MH;
	  else fitstaterrh[cha]=MH;
	  dNLLPlus=fabs(2*deltaNLL-1);
	}
      }
    }
  }

  for(int cha=0; cha<4; ++cha) {
    fiterrl[cha]=fitval[cha]-fiterrl[cha];
    fitstaterrl[cha]=fitval[cha]-fitstaterrl[cha+4];
    fitsysterrl[cha]=sqrt(fabs(pow(fiterrl[cha],2)-pow(fitstaterrl[cha],2)));
    fiterrh[cha]=fiterrh[cha]-fitval[cha];
    fitstaterrh[cha]=fitstaterrh[cha+4]-fitval[cha];
    fitsysterrh[cha]=sqrt(fabs(pow(fiterrh[cha],2)-pow(fitstaterrh[cha],2)));
    // patch if the scan arrested too early
    if(fiterrh[cha]==(-fitval[cha])) fiterrh[cha]=fiterrl[cha];
    
    std::cout << std::fixed;
    cout << std::setprecision(2) << "Mass for channel " << cha << " = " << fitval[cha] 
	 << " -" << fiterrl[cha] << " +" << fiterrh[cha] << " GeV (tot)" << endl;
    cout << std::setprecision(2) << "Mass for channel " << cha << " = " << fitval[cha] 
	 << " -" << fitstaterrl[cha] << " +" << fitstaterrh[cha] << " GeV (stat)" 
	 << " -" << fitsysterrl[cha] << " +" << fitsysterrh[cha] << " GeV (syst)" 
	 << endl;
  }


    TLatex l; l.SetTextFont(43); l.SetNDC(); l.SetTextSize(25);

    TCanvas *c1 = new TCanvas("c1","",0,44,540,750);
    c1->SetLeftMargin(0.16);
    c1->SetGridx(1);

    int nChann = 3;
    TH2F frame("frame",";best fit m_{X} (GeV);",1,123,129,nChann,0,nChann);

    TGraphAsymmErrors points(nChann);
    TGraphAsymmErrors points_nosyst(nChann);
    for (int cha=0; cha<3; ++cha) {
      TString channame("");
      if (cha==0) channame+=" 4#mu";
      if (cha==1) channame+=" 4e";
      if (cha==2) channame+=" 2e2#mu";
      points.SetPoint(cha,       fitval[cha],  cha+0.5);
      points.SetPointError(cha,  max(fiterrl[cha],fitstaterrl[cha]), max(fiterrh[cha],fitstaterrh[cha]), 0, 0);
      points_nosyst.SetPoint(cha,       fitval[cha],  cha+0.5);
      points_nosyst.SetPointError(cha,  fitstaterrl[cha], fitstaterrh[cha], 0, 0);
      frame.GetYaxis()->SetBinLabel(cha+1, channame);
    }
    points.SetLineColor(kRed);
    points.SetLineWidth(3);
    points.SetMarkerStyle(21);

    points_nosyst.SetLineColor(kBlack);
    points_nosyst.SetLineWidth(3);
    points.SetMarkerStyle(21);

    frame.GetXaxis()->SetNdivisions(6,kFALSE);
    frame.GetXaxis()->SetTitleSize(0.05);
    frame.GetXaxis()->SetLabelSize(0.04);
    frame.GetYaxis()->SetLabelSize(0.06);
    frame.Draw(); gStyle->SetOptStat(0);
    TBox globalFitBand(fitval[3]-fiterrl[3], 0, fitval[3]+fiterrh[3], nChann);
    globalFitBand.SetFillStyle(3013);
    globalFitBand.SetFillColor(65);
    globalFitBand.SetLineStyle(0);
    globalFitBand.DrawClone();
    TLine globalFitLine(fitval[3], 0, fitval[3], nChann);
    globalFitLine.SetLineWidth(4);
    globalFitLine.SetLineColor(214);
    globalFitLine.DrawClone();
    points.Draw("P SAME");
    points.Draw("[] SAME");
    points_nosyst.Draw("[] SAME");
    if(ndim==1) l.DrawLatex(0.20, 0.90, Form("H #rightarrow ZZ, 1D"));
    if(ndim==2) l.DrawLatex(0.20, 0.90, Form("H #rightarrow ZZ, 2D"));
    if(ndim==3) l.DrawLatex(0.20, 0.90, Form("H #rightarrow ZZ, 3D"));

    stringstream outnamepdf;
    outnamepdf << "bestfit_bychannel_" << ndim << "D.pdf";
    c1->SaveAs(outnamepdf.str().c_str());
    stringstream outnamepng;
    outnamepng << "bestfit_bychannel_" << ndim << "D.png";
    c1->SaveAs(outnamepng.str().c_str());
    stringstream outnamec;
    outnamec << "bestfit_bychannel_" << ndim << "D.C";
    c1->SaveAs(outnamec.str().c_str());

}


void plot2DScan() {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);

  TStyle *mystyle = RooHZZStyle("ZZ");
  mystyle->SetPalette(1);
  mystyle->cd();

  stringstream file78TeV;
  TFile *fit78TeV = TFile::Open("results/c2DMassErr/scanMass-MHMu_comb_new7TeV.root");

  TCanvas *c1 = new TCanvas("c1", "",0,22,763,622);
  gStyle->SetOptStat(0);
  c1->Range(106.8241,-0.3755458,161.5748,2.122271);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.1493289);
  c1->SetRightMargin(0.2114094);
  c1->SetTopMargin(0.04895105);
  c1->SetBottomMargin(0.1503496);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);

  contour("MH",122,132,"r",0,3);
  TH2 *h2d = (TH2*) gROOT->FindObject("h2d");
  h2d->SetContour(1000);
  h2d->SetTitle("");
  h2d->GetXaxis()->SetTitle("M_{H} (GeV)");
  h2d->GetYaxis()->SetTitle("#sigma/#sigma_{SM}");
  h2d->GetZaxis()->SetTitle("-2 #Delta lnL");
  h2d->GetZaxis()->SetLabelFont(42);
  h2d->GetZaxis()->SetLabelOffset(0.007);
  h2d->GetZaxis()->SetLabelSize(0.045);
  h2d->GetZaxis()->SetTitleSize(0.05);
  h2d->GetZaxis()->SetTitleFont(42);

  h2d->GetZaxis()->SetRangeUser(0,20);
  h2d->GetXaxis()->SetRangeUser(123,128.5);

  TLatex *CP = CMSPreliminary();
  CP->Draw();

  TPaveText *comment;
  comment = text("H #rightarrow ZZ",0.20,0.90,0.40,0.90);

  c1->SaveAs("mh_scan2D_full_hzz4l.pdf");
  c1->SaveAs("mh_scan2D_full_hzz4l.png");

}
