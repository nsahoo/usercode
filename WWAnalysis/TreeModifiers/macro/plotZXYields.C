#include <iostream>
#include <sstream>

#include "TMath.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"

using namespace std;

double wavg(double x, double y, double ex, double ey) {
  return (x/ex/ex + y/ey/ey)/(1/ex/ex + 1/ey/ey);
}

pair<double,double> werror(double x, double y, double ex, double ey) {
  double max = TMath::Max(x+ex,y+ey);
  double min = TMath::Min(x-ex,y-ey);

  double val = wavg(x,y,ex,ey);
  double eu=fabs(max-val);
  double ed=fabs(val-min);
  
  return std::make_pair<double,double>(ed,eu);
}

void plotZXYields(bool do7TeV) {

  // SS method
  double ssyield_7TeV[3] = {0.92, 2.52, 3.10};
  double ssyielde_7TeV[3] = {0.20, 0.30, 0.48};
  double ssyield_8TeV[3] = {3.75, 9.52, 12.8};
  double ssyielde_8TeV[3] = {0.52, 0.64, 0.91};

  // OS method
  double osyield_7TeV[3] = {1.18, 0.47, 0.89};
  double osyielde_7TeV[3] = {0.42, 0.21, 0.58};
  double osyield_8TeV[3] = {1.52, 3.98, 5.18};
  double osyielde_8TeV[3] = {0.71, 0.61, 1.02};

  // average
  double  yield_7TeV[3], yieldeup_7TeV[3], yieldedn_7TeV[3], yield_8TeV[3], yieldeup_8TeV[3], yieldedn_8TeV[3];
  
  for(int i=0;i<3;++i) {
    yield_7TeV[i] = wavg(ssyield_7TeV[i],osyield_7TeV[i],ssyielde_7TeV[i],osyielde_7TeV[i]);
    pair<double,double> err=werror(ssyield_7TeV[i],osyield_7TeV[i],ssyielde_7TeV[i],osyielde_7TeV[i]);
    yieldedn_7TeV[i] = err.first;
    yieldeup_7TeV[i] = err.second;

    yield_8TeV[i] = wavg(ssyield_8TeV[i],osyield_8TeV[i],ssyielde_8TeV[i],osyielde_8TeV[i]);
    pair<double,double> error=werror(ssyield_8TeV[i],osyield_8TeV[i],ssyielde_8TeV[i],osyielde_8TeV[i]);
    yieldedn_8TeV[i] = error.first;
    yieldeup_8TeV[i] = error.second;

    cout << "channel = " << i << "\tyield(7 TeV) = " << yield_7TeV[i] << " + " << yieldeup_7TeV[i] << " - " << yieldedn_7TeV[i] << endl;
    cout << "channel = " << i << "\tyield(8 TeV) = " << yield_8TeV[i] << " + " << yieldeup_8TeV[i] << " - " << yieldedn_8TeV[i] << endl;
  }

  gStyle->SetLineWidth(1.);
  gStyle->SetMarkerStyle(8);

  TGraphAsymmErrors *grSS = new TGraphAsymmErrors(3);
  TGraphAsymmErrors *grOS = new TGraphAsymmErrors(3);
  TGraphAsymmErrors *grAvg = new TGraphAsymmErrors(3);

  TGraph *grDummy = new TGraph(2);
  grDummy->SetPoint(0,-1.,6);
  grDummy->SetPoint(1,3,6);
  grDummy->SetMarkerColor(0);
  grDummy->GetXaxis()->SetBinLabel(0,"");
  grDummy->SetTitle("");

  for(int ch=0;ch<3;++ch) {
    if(do7TeV) {
      grSS->SetPoint(ch,ch-0.1,ssyield_7TeV[ch]);
      grSS->SetPointError(ch,0,0,ssyielde_7TeV[ch],ssyielde_7TeV[ch]);

      grOS->SetPoint(ch,ch+0.1,osyield_7TeV[ch]);
      grOS->SetPointError(ch,0,0,osyielde_7TeV[ch],osyielde_7TeV[ch]);

      grAvg->SetPoint(ch,ch,yield_7TeV[ch]);
      grAvg->SetPointError(ch,0,0,yieldedn_7TeV[ch],yieldeup_7TeV[ch]);
    } else {
      grSS->SetPoint(ch,ch-0.1,ssyield_8TeV[ch]);
      grSS->SetPointError(ch,0,0,ssyielde_8TeV[ch],ssyielde_8TeV[ch]);

      grOS->SetPoint(ch,ch+0.1,osyield_8TeV[ch]);
      grOS->SetPointError(ch,0,0,osyielde_8TeV[ch],osyielde_8TeV[ch]);

      grAvg->SetPoint(ch,ch,yield_8TeV[ch]);
      grAvg->SetPointError(ch,0,0,yieldedn_8TeV[ch],yieldeup_8TeV[ch]);
    }
  }

  grSS->SetMarkerColor(kBlue+2);
  grOS->SetMarkerColor(kCyan+2);
  grAvg->SetMarkerColor(kViolet-2);

  grSS->SetLineColor(kBlue+2);
  grOS->SetLineColor(kCyan+2);
  grAvg->SetLineColor(kViolet-2);

  TLegend* leg = new TLegend(0.2, 0.7, 0.5, 0.9);
  leg->SetFillColor(0);
  leg->AddEntry(grOS, "OS method","pl");
  leg->AddEntry(grSS, "SS method","pl");
  leg->AddEntry(grAvg, "Average","pl");

  TLatex* CP = new TLatex(-1.,do7TeV ? 10.5 : 20.5, do7TeV ? "CMS Preliminary                                                  #sqrt{s} = 7 TeV, L = 5.05 fb^{-1}" : "CMS Preliminary                                              #sqrt{s} = 8 TeV, L = 12.2 fb^{-1}");
  CP->SetTextSize(0.035);

  TLatex* xax0 = new TLatex(0-.1,-1, "4#mu");
  TLatex* xax1 = new TLatex(1-.1,-1, "4e");
  TLatex* xax2 = new TLatex(2-.1,-1, "2e2#mu");

  TCanvas *c1 = new TCanvas("c1","c1");
  grDummy->GetYaxis()->SetRangeUser(0,((do7TeV) ? 10 : 20));
  grDummy->Draw("ap");
  grSS->Draw("PE");
  grOS->Draw("PE");
  grAvg->Draw("PE");
  CP->Draw();
  xax0->Draw();
  xax1->Draw();
  xax2->Draw();
  leg->Draw();

  stringstream fss;
  if(do7TeV) fss << "ZXYields-7TeV.pdf";
  else fss << "ZXYields-8TeV.pdf";
  c1->SaveAs(fss.str().c_str());

}
