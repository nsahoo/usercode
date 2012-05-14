#include <iostream>
#include <vector>

#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TString.h>
#include <TLegend.h>
#include <TF1.h>
#include <TPaveText.h>

using namespace std;


void prep(TH1* hist, TF1* f, TPaveText* t, TCanvas* canvas, int pad, const char* leptoninfo, const char* isoinfo) {
    hist->SetTitle("");
    hist->GetXaxis()->SetLabelSize(0.035);
    hist->GetYaxis()->SetLabelSize(0.035);
    hist->GetXaxis()->SetTitleSize(0.035);
    hist->GetYaxis()->SetTitleSize(0.035);
    hist->GetXaxis()->SetTitleOffset(1.2);
    hist->GetYaxis()->SetTitleOffset(0.8);
    hist->GetXaxis()->SetTitle("P_{T} [GeV/c]");
    hist->GetYaxis()->SetTitle("Fake Rate (Reco#rightarrow ID+Isolation)");
    hist->SetMarkerStyle(20);
    hist->SetMarkerColor(4);
    hist->SetMarkerSize(0.7);
    hist->SetLineColor(4);
    hist->SetLineWidth(1);
    hist->SetStats(0);
    hist->GetYaxis()->SetRangeUser(0.0,0.4);

    f->SetLineWidth(1);
    f->SetLineStyle(2);
    f->SetLineColor(kBlue);

    t->AddText("#sqrt{s} = 7 TeV 2011, L = 4.9 fb^{-1}");
    t->AddText(leptoninfo);
    t->AddText(isoinfo);
    t->SetBorderSize(0);
    t->SetFillStyle(0);
    t->SetTextAlign(12);
    t->SetTextFont(42);
    t->SetTextSize(0.06);

    canvas->cd(pad);
    for (int i = 0; i < f->GetNumberFreeParameters(); ++i) f->SetParameter(i,0.1);
    hist->Fit(f, "", "", 5.0, 80.0);
    hist->Draw();
    t->Draw();

    for (int i = 0; i < f->GetNumberFreeParameters(); ++i) std::cerr << "\tfloat par" << i << " = " << f->GetParameter(i) << ";" << std::endl;
    std::cerr << std::endl;

}

void drawOneComparison(vector<TH1F*> histos, vector<TString> descr, TString xaxislabel, const char *filename) {

  if(histos.size()>3) {
    cout << "more than 3 histos not implemented." << endl;
    return;
  }

  if(histos.size()!=descr.size()) {
    cout << "description not complete!" << endl;
    return;
  }

  vector<int> colors;
  colors.push_back(kRed+1);
  colors.push_back(kAzure-6);
  colors.push_back(kTeal+3);

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  c1->SetGridx();
  c1->SetGridy();  
  TLegend* legend = new TLegend(0.20, 0.70, 0.43, 0.86);
  legend->SetBorderSize(   0);
  legend->SetFillColor (   0);
  legend->SetTextAlign (  12);
  legend->SetTextFont  (  42);
  legend->SetTextSize  (0.05);

  for(int i=0;i<(int)histos.size();++i) {
    
    histos[i]->SetMinimum(0);
    histos[i]->SetMaximum(0.4);
    histos[i]->SetMarkerSize(2);
    histos[i]->SetMarkerStyle(20);
    histos[i]->SetMarkerColor(colors[i]);
    histos[i]->SetLineColor(colors[i]);
    histos[i]->SetTitle("");
    if(TString(histos[i]->GetName()).Contains("PU")) {
      histos[i]->Fit("pol1","","same",4,35);
      histos[i]->GetFunction("pol1")->SetLineColor(colors[i]);
    }
    histos[i]->GetXaxis()->SetTitle(xaxislabel);
    histos[i]->GetYaxis()->SetTitle("efficiency");

    legend->AddEntry(histos[i],descr[i]);

    if(i==0) { 
      histos[i]->Draw("pe1");
      for(int bin=1;bin<=histos[i]->GetNbinsX();bin++) {
	//cout << "bin i = " << bin << ": " <<  histos[i]->GetBinContent(bin) << " +/- " << histos[i]->GetBinError(bin) << endl;
	cout << "m35_fakeRate[" << bin << "] = " << histos[i]->GetBinContent(bin) << ";" << endl;
      }
      for(int bin=1;bin<=histos[i]->GetNbinsX();bin++) {
	cout << "m35_fakeRate_err[" << bin << "] = " << histos[i]->GetBinError(bin) << ";" << endl;
      }
    }
    else histos[i]->Draw("same pe1");
  }
  legend->Draw();

  TString basename(filename);
  basename.ReplaceAll("Eff","FR");
  c1->SaveAs(basename+TString(".png"));
  c1->SaveAs(basename+TString(".pdf"));

}

void drawOneToOne(vector<TH1F*> set1, vector<TH1F*> set2, const char* desc1, const char* desc2, const char* xaxislabel) {
  if(set1.size()!=set2.size()) {
    cout << "first set and second set of histos have different sizes! ERROR! " << endl;
    return;
  }
  for(int i=0;i<(int)set1.size();++i) {
    vector<TH1F*> histos;
    vector<TString> desc;
    histos.push_back(set1[i]);
    histos.push_back(set2[i]);
    desc.push_back(TString(desc1));
    desc.push_back(TString(desc2));
    if(set1[i]==0 || set2[i]==0) {
      cout << "histogram not found!" << endl;
      continue;
    }
    drawOneComparison(histos,desc,TString(xaxislabel),set1[i]->GetName());
  }

}

void drawOneToTwo(vector<TH1F*> set1, vector<TH1F*> set2, vector<TH1F*> set3, const char* desc1, const char* desc2, const char* desc3, const char* xaxislabel) {
  if(set1.size()!=set2.size() || set1.size()!=set3.size()) {
    cout << "first set and second/third set of histos have different sizes! ERROR! " << endl;
    return;
  }
  for(int i=0;i<(int)set1.size();++i) {
    vector<TH1F*> histos;
    vector<TString> desc;
    histos.push_back(set1[i]);
    histos.push_back(set2[i]);
    histos.push_back(set3[i]);
    desc.push_back(TString(desc1));
    desc.push_back(TString(desc2));
    desc.push_back(TString(desc3));
    if(set1[i]==0 || set2[i]==0 || set3[i]==0) {
      cout << "histogram not found!" << endl;
      continue;
    }
    drawOneComparison(histos,desc,TString(xaxislabel),set1[i]->GetName());
  }

}

void compareElectronIdsUnbiased() {

  TFile *file = TFile::Open("fakerates_zll1e.root");

  // eta
  // std iso
  TH1F *NoTrgElectronhzzPfIsoEtaHighPt = (TH1F*)file->Get("NoTrgElectronhzzPfIsoEtaHighPt_Eff");
  TH1F *NoTrgElectronhzzPfIsoEtaLowPt = (TH1F*)file->Get("NoTrgElectronhzzPfIsoEtaLowPt_Eff");
  // MVA iso
  TH1F *NoTrgElectronhzzMvaLooseEtaHighPt = (TH1F*)file->Get("NoTrgElectronhzzMvaLooseEtaHighPt_Eff");
  TH1F *NoTrgElectronhzzMvaLooseEtaLowPt = (TH1F*)file->Get("NoTrgElectronhzzMvaLooseEtaLowPt_Eff");
  TH1F *NoTrgElectronhzzMvaTightEtaHighPt = (TH1F*)file->Get("NoTrgElectronhzzMvaTightEtaHighPt_Eff");
  TH1F *NoTrgElectronhzzMvaTightEtaLowPt = (TH1F*)file->Get("NoTrgElectronhzzMvaTightEtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2, etaSet3;
  etaSet1.push_back(NoTrgElectronhzzPfIsoEtaHighPt);
  etaSet1.push_back(NoTrgElectronhzzPfIsoEtaLowPt);
  etaSet2.push_back(NoTrgElectronhzzMvaLooseEtaHighPt);
  etaSet2.push_back(NoTrgElectronhzzMvaLooseEtaLowPt);
  etaSet3.push_back(NoTrgElectronhzzMvaTightEtaHighPt);
  etaSet3.push_back(NoTrgElectronhzzMvaTightEtaLowPt);
 
  drawOneToTwo(etaSet1,etaSet2,etaSet3,"H#rightarrow ZZ MVA ref.", "H#rightarrow ZZ MVA Loose","H#rightarrow ZZ MVA Tight","#eta");

  // pt
  // std iso
  TH1F *NoTrgElectronhzzPfIsoPtBarrel = (TH1F*)file->Get("NoTrgElectronhzzPfIsoPtBarrel_Eff");
  TH1F *NoTrgElectronhzzPfIsoPtEndcap = (TH1F*)file->Get("NoTrgElectronhzzPfIsoPtEndcap_Eff");
  // MVA iso
  TH1F *NoTrgElectronhzzMvaLoosePtBarrel = (TH1F*)file->Get("NoTrgElectronhzzMvaLoosePtBarrel_Eff");
  TH1F *NoTrgElectronhzzMvaLoosePtEndcap = (TH1F*)file->Get("NoTrgElectronhzzMvaLoosePtEndcap_Eff");
  TH1F *NoTrgElectronhzzMvaTightPtBarrel = (TH1F*)file->Get("NoTrgElectronhzzMvaTightPtBarrel_Eff");
  TH1F *NoTrgElectronhzzMvaTightPtEndcap = (TH1F*)file->Get("NoTrgElectronhzzMvaTightPtEndcap_Eff");

  vector<TH1F*> ptSet1, ptSet2, ptSet3;
  ptSet1.push_back(NoTrgElectronhzzPfIsoPtBarrel);
  ptSet1.push_back(NoTrgElectronhzzPfIsoPtEndcap);
  ptSet2.push_back(NoTrgElectronhzzMvaLoosePtBarrel);
  ptSet2.push_back(NoTrgElectronhzzMvaLoosePtEndcap);
  ptSet3.push_back(NoTrgElectronhzzMvaTightPtBarrel);
  ptSet3.push_back(NoTrgElectronhzzMvaTightPtEndcap);

  drawOneToTwo(ptSet1,ptSet2,ptSet3,"H#rightarrow ZZ MVA ref","H#rightarrow ZZ MVA Loose","H#rightarrow ZZ MVA Tight","p_{T} [GeV]");

  // PU
  // TH1F *NoTrgElectronhzzPfIsoPUBarrel = (TH1F*)file->Get("NoTrgElectronhzzPfIsoPUBarrel_Eff");
  // TH1F *NoTrgElectronhzzPfIsoPUEndcap = (TH1F*)file->Get("NoTrgElectronhzzPfIsoPUEndcap_Eff");
  // // MVA iso
  // TH1F *NoTrgElectronhzzMvaLoosePUBarrel = (TH1F*)file->Get("NoTrgElectronhzzMvaLoosePUBarrel_Eff");
  // TH1F *NoTrgElectronhzzMvaLoosePUEndcap = (TH1F*)file->Get("NoTrgElectronhzzMvaLoosePUEndcap_Eff");
  // TH1F *NoTrgElectronhzzMvaTightPUBarrel = (TH1F*)file->Get("NoTrgElectronhzzMvaTightPUBarrel_Eff");
  // TH1F *NoTrgElectronhzzMvaTightPUEndcap = (TH1F*)file->Get("NoTrgElectronhzzMvaTightPUEndcap_Eff");

  // vector<TH1F*> puSet1, puSet2, puSet3;
  // puSet1.push_back(NoTrgElectronhzzPfIsoPUBarrel);
  // puSet1.push_back(NoTrgElectronhzzPfIsoPUEndcap);
  // puSet2.push_back(NoTrgElectronhzzMvaLoosePUBarrel);
  // puSet2.push_back(NoTrgElectronhzzMvaLoosePUEndcap);
  // puSet3.push_back(NoTrgElectronhzzMvaTightPUBarrel);
  // puSet3.push_back(NoTrgElectronhzzMvaTightPUEndcap);

  // drawOneToTwo(puSet1,puSet2,puSet3,"H#rightarrow ZZ MVA ref","H#rightarrow ZZ MVA Loose","H#rightarrow ZZ MVA Tight","# vertices");


}


void compareMuonIdsUnbiased() {

  TFile *file = TFile::Open("fakerates_zll1m.root");

  // eta
  TH1F *NoTrgMuonhzzPfIsoEtaHighPt = (TH1F*)file->Get("NoTrgMuonhzzPfIsoEtaHighPt_Eff");
  TH1F *NoTrgMuonhzzPfIsoEtaLowPt = (TH1F*)file->Get("NoTrgMuonhzzPfIsoEtaLowPt_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *NoTrgMuonhzzMvaPfIsoEtaHighPt = (TH1F*)file->Get("NoTrgMuonhzzMvaPfIsoEtaHighPt_Eff");
  TH1F *NoTrgMuonhzzMvaPfIsoEtaLowPt = (TH1F*)file->Get("NoTrgMuonhzzMvaPfIsoEtaLowPt_Eff");

  vector<TH1F*> etaSet1, etaSet2;
  etaSet1.push_back(NoTrgMuonhzzPfIsoEtaHighPt);
  etaSet1.push_back(NoTrgMuonhzzPfIsoEtaLowPt);
  etaSet2.push_back(NoTrgMuonhzzMvaPfIsoEtaHighPt);
  etaSet2.push_back(NoTrgMuonhzzMvaPfIsoEtaLowPt);
  
  drawOneToOne(etaSet1,etaSet2,"H#rightarrow ZZ PF iso","H#rightarrow ZZ MVA PF Iso","#eta");

  // pt
  TH1F *NoTrgMuonhzzPfIsoPtBarrel = (TH1F*)file->Get("NoTrgMuonhzzPfIsoPtBarrel_Eff");
  TH1F *NoTrgMuonhzzPfIsoPtEndcap = (TH1F*)file->Get("NoTrgMuonhzzPfIsoPtEndcap_Eff");
  // ---> HZZ id + EA corrected isolation
  TH1F *NoTrgMuonhzzMvaPfIsoPtBarrel = (TH1F*)file->Get("NoTrgMuonhzzMvaPfIsoPtBarrel_Eff");
  TH1F *NoTrgMuonhzzMvaPfIsoPtEndcap = (TH1F*)file->Get("NoTrgMuonhzzMvaPfIsoPtEndcap_Eff");

  vector<TH1F*> ptSet1, ptSet2;
  ptSet1.push_back(NoTrgMuonhzzPfIsoPtBarrel);
  ptSet1.push_back(NoTrgMuonhzzPfIsoPtEndcap);
  ptSet2.push_back(NoTrgMuonhzzMvaPfIsoPtBarrel);
  ptSet2.push_back(NoTrgMuonhzzMvaPfIsoPtEndcap);

  drawOneToOne(ptSet1,ptSet2,"H#rightarrow ZZ PF iso","H#rightarrow ZZ MVA PF iso","p_{T} [GeV]");

  // PU
  // TH1F *NoTrgMuonhzzPfIsoPUBarrel = (TH1F*)file->Get("NoTrgMuonhzzPfIsoPUBarrel_Eff");
  // TH1F *NoTrgMuonhzzPfIsoPUEndcap = (TH1F*)file->Get("NoTrgMuonhzzPfIsoPUEndcap_Eff");
  // // ---> HZZ id + EA corrected isolation
  // TH1F *NoTrgMuonhzzMvaPfIsoPUBarrel = (TH1F*)file->Get("NoTrgMuonhzzMvaPfIsoPUBarrel_Eff");
  // TH1F *NoTrgMuonhzzMvaPfIsoPUEndcap = (TH1F*)file->Get("NoTrgMuonhzzMvaPfIsoPUEndcap_Eff");

  // vector<TH1F*> puSet1, puSet2;
  // puSet1.push_back(NoTrgMuonhzzPfIsoPUBarrel);
  // puSet1.push_back(NoTrgMuonhzzPfIsoPUEndcap);
  // puSet2.push_back(NoTrgMuonhzzMvaPfIsoPUBarrel);
  // puSet2.push_back(NoTrgMuonhzzMvaPfIsoPUEndcap);

  // drawOneToOne(puSet1,puSet2,"H#rightarrow ZZ PF iso","H#rightarrow ZZ MVA PF iso","# vertices");


}

enum {
  kMuPfIso = 0,
  kMuMvaPfIso
};

void fitMuonIdsUnbiased(int wp) {

  if(wp>1) {
    cout << "wp must be kPfIso / kMvaPfIso" << endl;
    return;
  }
  TString wpstr = (wp==kMuPfIso) ? TString("PfIso") : TString("MvaPfIso");

  TFile *file = TFile::Open("fakerates_zll1m.root");

  TH1F *barrel = (TH1F*)file->Get(TString("NoTrgMuonhzz")+wpstr+TString("PtBarrel_Eff"));
  TH1F *endcap = (TH1F*)file->Get(TString("NoTrgMuonhzz")+wpstr+TString("PtEndcap_Eff"));

  std::vector<TH1F*> fakerates1D;
  fakerates1D.push_back(barrel);
  fakerates1D.push_back(endcap);

  TCanvas* canvas = new TCanvas("fakerate","Fake Rates",1000,1000);
  canvas->Divide(1,2);

  TPaveText* text1  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");
  TPaveText* text2  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");

  TF1 * f1 = new TF1("f1","[0]+[1]*(exp(-[2]*x) - exp(-[3]*x))");
  for (int i = 0; i <=f1->GetNumberFreeParameters(); ++i) f1->SetParameter(i,1);
  f1->SetParLimits(1,0,5);

  prep(barrel, f1, text1, canvas, 1, "Barrel", "");
  prep(endcap, f1, text2, canvas, 2, "Endcap", "");
  
  TString pdffile = (wp==kMuPfIso) ? TString("fakerates_zll1m_PfIso.pdf") : TString("fakerates_zll1m_MvaPfIso.pdf");
  canvas->SaveAs(pdffile);

  TString pngfile = (wp==kMuPfIso) ? TString("fakerates_zll1m_PfIso.png") : TString("fakerates_zll1m_MvaPfIso.png");
  canvas->SaveAs(pngfile);

}

enum {
  kElePfIso = 0,
  kEleMvaLoose,
  kEleMvaTight
};

void fitElectronIdsUnbiased(int wp) {

  if(wp>2) {
    cout << "wp must be kElePfIso / kMvaLoose / kMvaTight" << endl;
    return;
  }
  TString wpstr;
  if(wp==kElePfIso) wpstr = TString("PfIso");
  else if(wp==kEleMvaLoose) wpstr = TString("MvaLoose");
  else  wpstr = TString("MvaTight");

  TFile *file = TFile::Open("fakerates_zll1e.root");

  TH1F *barrel = (TH1F*)file->Get(TString("NoTrgElectronhzz")+wpstr+TString("PtBarrel_Eff"));
  TH1F *endcap = (TH1F*)file->Get(TString("NoTrgElectronhzz")+wpstr+TString("PtEndcap_Eff"));

  std::vector<TH1F*> fakerates1D;
  fakerates1D.push_back(barrel);
  fakerates1D.push_back(endcap);

  TCanvas* canvas = new TCanvas("fakerate","Fake Rates",1000,1000);
  canvas->Divide(1,2);

  TPaveText* text1  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");
  TPaveText* text2  = new TPaveText(0.15, 0.9, 0.8, 0.7, "ndc");

  TF1 * f1 = new TF1("f1","[0]+[1]*(exp(-[2]*x) - exp(-[3]*x))");
  for (int i = 0; i <=f1->GetNumberFreeParameters(); ++i) f1->SetParameter(i,1);

  prep(barrel, f1, text1, canvas, 1, "Barrel", "");
  prep(endcap, f1, text2, canvas, 2, "Endcap", "");
  
  TString pdffile;
  if(wp==kElePfIso) pdffile=TString("fakerates_zll1e_PfIso.pdf");
  else if(wp==kEleMvaLoose) pdffile=TString("fakerates_zll1e_MvaLoosePfIso.pdf");
  else pdffile=TString("fakerates_zll1e_MvaTightPfIso.pdf");

  canvas->SaveAs(pdffile);

  TString pngfile;
  if(wp==kElePfIso) pngfile=TString("fakerates_zll1e_PfIso.png");
  else if(wp==kEleMvaLoose) pngfile=TString("fakerates_zll1e_MvaLoosePfIso.png");
  else pngfile=TString("fakerates_zll1e_MvaTightPfIso.png");

  canvas->SaveAs(pngfile);

}
