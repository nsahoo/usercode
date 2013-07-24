#include <TString.h>
#include <TStyle.h>
#include <TChain.h>
#include <TPaveText.h>
#include <RooWorkspace.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooHistPdf.h>
#include <RooProdPdf.h>
#include <RooPlot.h>
#include <RooBinning.h>
#include <RooMsgService.h>
#include <RooParametricStepFunction.h>
#include <TArrayD.h>
#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooFitResult.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TF2.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TFile.h>
#include <iostream>
#include <cstdio>
#include <cmath>

#include "MuScleFit/Calibration/src/MuScleFitCorrector.cc"
#include "MuScleFit/Calibration/src/Functions.cc"

MuScleFitCorrector musclecor_data("MuScleFit/Calibration/data/MuScleFit_2011_DATA_42X.txt");
MuScleFitCorrector musclecor_mc("MuScleFit/Calibration/data/MuScleFit_2011_MC_42X.txt");
MuScleFitCorrector musclecor2012_dataABC("MuScleFit/Calibration/data/MuScleFit_2012ABC_DATA_ReReco_53X.txt");
MuScleFitCorrector musclecor2012_dataD("MuScleFit/Calibration/data/MuScleFit_2012D_DATA_ReReco_53X.txt");
MuScleFitCorrector musclecor2012_mc("MuScleFit/Calibration/data/MuScleFit_2012_MC_53X_smearReReco.txt");
#define MUSCLE2012_SMEARING(X,Y,Z)  X.applyPtSmearing(Y, Z, false)

int gYear = 2011;
TString label = "";
TString labelLeft = "";
TString labelRight = "";


RooWorkspace *w = new RooWorkspace("w","w");
TCanvas *c1 = new TCanvas("c1","c1");
TString prefix = "";
TF2 *ebeCorr  = 0;
TH2 *ebeCorrH = 0;
TF2 *accRej  = 0;
void initCanvas() {
    c1 = new TCanvas("c1","c1");
    c1->cd();
}

void initBaseNoParam() {
    // assumes that zmll, logdm and MREF are already defined
    w->factory("expr::sigmaEff(\"exp(@0)*@1*@2\", logdm, MREF, sigma[2,.5,4])");
}

void initBaseParam(bool doPtBinning) {
    // assumes that zmll and MREF are already defined
    w->factory("pt1[5,200]");  w->factory("eta1[1,0,2.5]"); w->factory("dm1[0,1]"); 
    w->factory("pt2[5,200]");  w->factory("eta2[1,0,2.5]"); w->factory("dm2[0,1]"); 

    const int ptbins = 6;
    const double ptvals[ptbins+1] =  { 5, 9, 15, 20, 30, 50,  100 };
#if 1
    const int etabins = 7;
    const double etavals[etabins+1] =  { 0.0, 0.5, 0.8, 1.2, 1.5, 1.8, 2.1, 2.5 };
#else
    const int etabins = 26;
    double etavals[etabins+1]; for (int i = 0; i <= etabins; ++i) etavals[i] = 0.1*i;
#endif


    w->var("pt1")->setBinning(RooBinning(ptbins,ptvals));
    w->var("eta1")->setBinning(RooBinning(etabins,etavals));

    TArrayD etavalsAD(etabins+1, etavals); RooArgList etavarshi, etavarslo;
    // implementation notice: a RooParametricStepFunction with N bins has has (N-1) parameters, and its integral is 1.0
    //                        so, for a uniform scale factor of 1.0 we expect bin values of 1.0/2.5 = 0.4.
    //                        and then we have to multiply the RooParametricStepFunction by a factor 2.5 to get 1.0
    for (int i = 0; i < etabins-1; ++i) {
        etavarshi.add( *w->factory(Form("dmScaleEta_hiPt_eta%02d_%02d[0.4,.1,2.0]", int(10*etavals[i]), int(10*etavals[i+1]))) );
        etavarslo.add( *w->factory(Form("dmScaleEta_loPt_eta%02d_%02d[0.4,.1,2.0]", int(10*etavals[i]), int(10*etavals[i+1]))) );
    }
    RooParametricStepFunction *rpsf1_hi = new RooParametricStepFunction("dm1ScaleEta_hiPt","",*w->var("eta1"),etavarshi,etavalsAD,etabins);
    RooParametricStepFunction *rpsf1_lo = new RooParametricStepFunction("dm1ScaleEta_loPt","",*w->var("eta1"),etavarslo,etavalsAD,etabins);
    RooParametricStepFunction *rpsf2_hi = new RooParametricStepFunction("dm2ScaleEta_hiPt","",*w->var("eta2"),etavarshi,etavalsAD,etabins);
    RooParametricStepFunction *rpsf2_lo = new RooParametricStepFunction("dm2ScaleEta_loPt","",*w->var("eta2"),etavarslo,etavalsAD,etabins);
    w->import(*rpsf1_hi, RooFit::Silence(1)); w->import(*rpsf2_hi, RooFit::Silence(1));
    w->import(*rpsf1_lo, RooFit::Silence(1)); w->import(*rpsf2_lo, RooFit::Silence(1));
    w->factory("dmScaleNorm_hiPt[2.5,0.5,10.0]");
    w->factory("dmScaleNorm_loPt[2.5,0.5,10.0]");
    if (doPtBinning) {
        w->factory("expr::dm1Scale(\"@0*@1*(10./(@4+5)) + @2*@3*(1-10./(@4+5))\", dm1ScaleEta_loPt, dmScaleNorm_loPt, dm1ScaleEta_hiPt, dmScaleNorm_hiPt, pt1)"); 
        w->factory("expr::dm2Scale(\"@0*@1*(10./(@4+5)) + @2*@3*(1-10./(@4+5))\", dm2ScaleEta_loPt, dmScaleNorm_loPt, dm2ScaleEta_hiPt, dmScaleNorm_hiPt, pt2)"); 
    } else {
        w->factory("expr::dm1Scale(\"@0*@1\", dm1ScaleEta_loPt, dmScaleNorm_loPt)"); 
        w->factory("expr::dm2Scale(\"@0*@1\", dm2ScaleEta_loPt, dmScaleNorm_loPt)"); 
    }

    //w->factory("expr::sigmaEff(\"@0\", sigma[2,.5,4])");
    //w->factory("expr::sigmaEff(\"sqrt(@0*@0+@1*@1)*@2*@3\", dm1,dm2,MREF, sigma[1,0.4,2.5])");
    w->factory("expr::sigmaEff(\"sqrt(@0*@0+@1*@1)*@2\", prod(dm1,dm1Scale),prod(dm2,dm2Scale),MREF)");
}

void initModelZ(int param, bool doPtBin, bool doDoubleCB=false) {
    w->factory("MREF[91.1876]");
    w->factory("zmll[75,105]"); w->var("zmll")->setBins(60);
    //w->factory("zmll[81,105]"); w->var("zmll")->setBins(30);
    switch (param) {
        case 1:
            initBaseParam(doPtBin);
            break;
        case 0:
            w->factory("logdm[-5.5,-2.5]"); w->var("logdm")->setBins(30);
            initBaseNoParam();
            break;
        case -1:
            w->factory("logdm[-5.5,-2.5]"); w->var("logdm")->setBins(30);
            w->factory("expr::sigmaEff(\"0.01*@0*@1\", MREF, sigma[2,.5,4])");
            break;
        
    }

    if (doDoubleCB) {
        w->factory("DoubleCB::resol(zmll, mean[0,-5,5], sigmaEff, alpha[3., 0.05, 50.], n[1, 0.1, 100.], alpha2[3., 0.05, 50.],   n2[1, 0.1, 100.])");
    } else {
        w->factory("CBShape::resol(zmll, mean[0,-5,5], sigmaEff, alpha[3., 0.05, 50.], n[1, 0.1, 100.])");
    }
    w->factory("BreitWigner::true(zmll, MREF, width[2.495])");
    w->factory("FCONV::zed(zmll,true,resol)");
    //w->factory("FCONV::signal(zmll,true,resol)");

    w->factory("Exponential::background(zmll, lp[-0.1,-1,0])");
    w->factory("SUM::signal(fSig[0.9,.7,1.0] * zed, background)");
    //w->Print("V");
    w->saveSnapshot("clean", w->allVars());
    w->saveSnapshot("veryClean", w->allVars());
}


void initModelJ(int param, bool doPtBin) {
    w->factory("MREF[3.09692]");
    w->factory("zmll[2.91,3.3]"); w->var("zmll")->setBins(30);
    if (param == 1) {
        initBaseParam(doPtBin);
    } else if (param == 0) {
        w->factory("logdm[-5.5,-2.5]"); w->var("logdm")->setBins(30);
        initBaseNoParam();
    } else if (param == -1) {
        w->factory("logdm[-5.5,-2.5]"); w->var("logdm")->setBins(30);
        w->factory("expr::sigmaEff(\"0.01*@0*@1\", MREF, sigma[2,.5,4])");
    }

    //w->factory("CBShape::jpsi(zmll,   sum(MREF,   mean[0,-.1,.1]), sigmaEff,   alpha[3., 0.05, 50.],   n[1, 0.1, 100.])");
    w->factory("DoubleCB::jpsi(zmll,   sum(MREF,   mean[0,-.1,.1]), sigmaEff,   alpha[3., 0.05, 50.],   n[2, 0.8, 100.], alpha2[3., 0.05, 50.],   n2[1, 0.8, 100.])");
    w->factory("Bernstein::bkg3(zmll, { c0[0,1], c1[0,1], c2[0,1] })");
    w->factory("SUM::signal(fSig[.8,0.5,1]*jpsi, bkg3)");
    w->saveSnapshot("clean", w->allVars());
    w->saveSnapshot("veryClean", w->allVars());
}

void initModelU(int mc, int param, bool doPtBin) {
    w->factory("MREF[9.460]");
    w->factory("MREF_2S[10.02]"); w->factory("SCAL_2S[1.059]");
    w->factory("MREF_3S[10.36]"); w->factory("SCAL_3S[1.095]");
    w->factory("zmll[8.7,11]"); w->var("zmll")->setBins(60);
    if (param == 1) {
        initBaseParam(doPtBin);
    } else if (param == 0) {
        w->factory("logdm[-5.5,-2.5]"); w->var("logdm")->setBins(15);
        initBaseNoParam();
    } else if (param == -1) {
        w->factory("logdm[-5.5,-2.5]"); w->var("logdm")->setBins(30);
        w->factory("expr::sigmaEff(\"0.01*@0*@1\", MREF, sigma[2,.5,4])");
    }

    w->factory("CBShape::ups1s(zmll, sum(MREF, mean[0,-.1,.1]), sigmaEff,  alpha[3., 0.05, 50.],   n[1, 0.1, 100.])");
    if (mc) {
        w->factory("Bernstein::bkg3(zmll, { c0[0,1], c1[0,1], c2[0,1] })");
        w->factory("SUM::signal(fSig[.8,0,1]*ups1s, bkg3)");
    } else {
        w->factory("Bernstein::bkg3(zmll, { c0[0,1], c1[0,1], c2[0,1], c3[0,1], c4[0,1] })");
        w->factory("CBShape::ups2s(zmll,   sum(MREF_2S, mean), prod(SCAL_2S,  sigmaEff),   alpha,   n)");
        w->factory("CBShape::ups3s(zmll,   sum(MREF_3S, mean), prod(SCAL_3S,  sigmaEff),   alpha,   n)");
        w->factory("SUM::signal(  fSig1S[.6,0,1]*ups1s,    fSig2S[.2,0,1]*ups2s,    fSig3S[.1,0,1]*ups3s,    bkg3)");
    }
    w->saveSnapshot("clean", w->allVars());
    w->saveSnapshot("veryClean", w->allVars());
}


class GioProjWData {
    public:
        GioProjWData(const RooAbsPdf &pdf, const RooRealVar &xvar, const RooAbsData &data, int rebinFactor=20) ;
        GioProjWData(const RooAbsPdf &pdf, const RooRealVar &xvar,       RooDataSet &data, int rebinFactor=20) ;
        ~GioProjWData() { delete hist_; delete pdf_; }
        const RooAbsPdf & pdf() const { return *pdf_; }
        RooAbsPdf & pdf() { return *pdf_; }
        RooPlot* plotOn(RooPlot* frame,
                          const RooCmdArg& arg1=RooCmdArg::none(), const RooCmdArg& arg2=RooCmdArg::none(),
                          const RooCmdArg& arg3=RooCmdArg::none(), const RooCmdArg& arg4=RooCmdArg::none(),
                          const RooCmdArg& arg5=RooCmdArg::none(), const RooCmdArg& arg6=RooCmdArg::none(),
                          const RooCmdArg& arg7=RooCmdArg::none(), const RooCmdArg& arg8=RooCmdArg::none(),
                          const RooCmdArg& arg9=RooCmdArg::none(), const RooCmdArg& arg10=RooCmdArg::none()
              ) const { return pdf().plotOn(frame,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10); }
    private:
        RooDataHist *hist_; 
        RooHistPdf  *pdf_;
};
GioProjWData::GioProjWData(const RooAbsPdf &pdf, const RooRealVar &xvar, const RooAbsData &data, int rebinFactor) {
    RooArgSet vars(*data.get());  vars.remove(xvar,false,true);
    RooArgSet *pars = pdf.getParameters(RooArgSet(xvar));
    TH1 *hist = 0;
    RooDataHist others("","dm",vars,data);
    for (int i = 0, n = others.numEntries(); i < n; ++i) {
        *pars = *others.get(i);
        TH1 *h = pdf.createHistogram("htemp",xvar,RooFit::Binning(rebinFactor*xvar.getBins()));
        h->SetDirectory(0);
        h->Scale(others.weight()/h->Integral());
        if (hist == 0) hist = h; else { hist->Add(h); delete h; }
    }
    hist_ = new RooDataHist("","",RooArgList(xvar),hist);
    pdf_  = new RooHistPdf("","",RooArgSet(xvar),*hist_,2);
    delete pars;
    delete hist;
}
GioProjWData::GioProjWData(const RooAbsPdf &pdf, const RooRealVar &xvar, RooDataSet &data, int rebinFactor) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    RooArgSet vars(*data.get());  vars.remove(xvar,false,true);
    RooArgSet *pars = pdf.getParameters(RooArgSet(xvar));
    TH1 *hist = 0;
    RooDataSet others("","dm",&data,vars);
    std::cout << "Will have to average over " << others.numEntries() << " entries." << std::endl;
    for (int i = 0, n = others.numEntries(); i < n; ++i) {
        *pars = *others.get(i);
        TH1 *h = pdf.createHistogram("htemp",xvar,RooFit::Binning(rebinFactor*xvar.getBins()));
        h->SetDirectory(0);
        h->Scale(others.weight()/h->Integral());
        if (hist == 0) hist = h; else { hist->Add(h); delete h; }
    }
    hist_ = new RooDataHist("","",RooArgList(xvar),hist);
    pdf_  = new RooHistPdf("","",RooArgSet(xvar),*hist_,2);
    delete pars;
    delete hist;
}

void dmScalePlotPrint(TH2 *hist, TString name, TString postfix, TString drawOpt="COLZ TEXT") {
    hist->GetXaxis()->SetTitle("muon p_{T} (GeV)");
    hist->GetXaxis()->SetMoreLogLabels(1);
    hist->GetYaxis()->SetTitle("muon |#eta|");
    hist->GetYaxis()->SetDecimals(1);
    hist->GetXaxis()->SetNoExponent(1);
    hist->GetYaxis()->SetNdivisions(505);
    hist->GetXaxis()->SetLabelOffset(0.00);
    hist->GetZaxis()->SetTitle("#sigma(obs.)/#sigma(pred.)");
    hist->GetZaxis()->SetRangeUser(0.501,1.499);
    hist->GetZaxis()->SetDecimals(1);
    hist->SetContour(100);
    hist->SetMarkerSize(2.0);
    gStyle->SetPaintTextFormat(".2f");
    gStyle->SetTextFont(42);
    c1->SetLogx(1); float rm = c1->GetRightMargin(); c1->SetRightMargin(0.2);
    hist->Draw(drawOpt);
    c1->Print(Form("plots/%d/fit2D_ptEta_",gYear)+prefix+name+"_"+postfix+".png");
    c1->SetLogx(0); c1->SetRightMargin(rm);
}

TH2* dmScalePlot(TString name, bool plot=true) {
    using namespace RooFit;
    RooAbsReal *dm1Scale = w->function("dm1Scale");
    if (dm1Scale == 0) { std::cerr << "Missing dm1Scale: cannot do dmScalePlot." << std::endl; return 0; }
    TH2 *hist = dynamic_cast<TH2*>(dm1Scale->createHistogram("dmScale", *w->var("pt1"), YVar(*w->var("eta1")), Scaling(0)));
    if (plot) dmScalePlotPrint(hist,name,"map2D");
    return hist;
}
TProfile2D * dmScalePlotErr(RooFitResult *res, TString name, int ntoys=5000) {
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    TH2 *nominal = dmScalePlot(name,true); nominal->SetDirectory(0);
    int nx = nominal->GetNbinsX(), ny = nominal->GetNbinsY();
    TAxis *ax = nominal->GetXaxis(), *ay = nominal->GetYaxis();
    double xbins[999], ybins[999];
    for (int i = 1; i <= nx; ++i) xbins[i-1] = ax->GetBinLowEdge(i);
    for (int i = 1; i <= ny; ++i) ybins[i-1] = ay->GetBinLowEdge(i);
    xbins[nx] = nominal->GetXaxis()->GetBinUpEdge(nx);
    ybins[ny] = nominal->GetYaxis()->GetBinUpEdge(ny);
    TProfile2D *toys = new TProfile2D("toys","toys", nx, xbins, ny, ybins, "S");
    for (int i = 0; i < ntoys; ++i) {
        w->allVars() = res->randomizePars();
        TH2 *toy = dmScalePlot(name,false);
        for (int ix = 1; ix <= nx; ++ix) { for (int iy = 1; iy <= ny; ++iy) {
            //if (ix == 1 && iy == 1) printf("toy %3d: %lf\n", i, toy->GetBinContent(ix,iy));
            toys->Fill( ax->GetBinCenter(ix), ay->GetBinCenter(iy), toy->GetBinContent(ix,iy) );
        } }
        delete toy;
    }
    w->allVars() = res->floatParsFinal();
    dmScalePlotPrint(toys, name, "map2DErr", "COLZ TEXTE");
    return toys;
}


void cmsprelim(double x1=0.75, double y1=0.40, double x2=0.95, double y2=0.48, const int align=12, const char *text="CMS Preliminary", double textSize=0.033, int color=1) { 
    TPaveText *cmsprel = new TPaveText(x1,y1,x2,y2,"NDC");
    cmsprel->SetTextSize(textSize);
    cmsprel->SetFillColor(0);
    cmsprel->SetFillStyle(0);
    cmsprel->SetLineStyle(2);
    cmsprel->SetLineColor(0);
    cmsprel->SetTextAlign(align);
    cmsprel->SetTextFont(42);
    cmsprel->SetTextColor(color);
    cmsprel->AddText(text);
    cmsprel->Draw("same");
}

void tinyPrelim(TString what) {
    cmsprelim(.26, .955, .40, .999, 12, "CMS Preliminary"); 
    if (what.Contains("reco53x") || what.Contains("reco42x")) {
        cmsprelim(.48, .955, .99, .999, 32, Form("#sqrt{s} = %d TeV, L = %.1f fb^{-1}", (gYear+7-2011), (gYear==2011?5.1:19.5)));
    } else if (what.Contains("mc42x") || what.Contains("mc53x")) {
        cmsprelim(.48, .955, .99, .999, 32, Form("Simulation, #sqrt{s} = %d TeV", (gYear+7-2011)));
    } 
}

void fitMass(RooAbsData *data, TString name, TString binName) { 
    w->loadSnapshot("clean");
    RooArgSet vars(*data->get());
    vars.remove(*w->var("zmll"),false,true);
   
    const RooCmdArg & cpu =  (prefix+name).Contains("unbin") ? RooFit::NumCPU(6) :  RooFit::NumCPU(4);
    RooFitResult *res = w->pdf("signal")->fitTo(*data,RooFit::ConditionalObservables(vars), RooFit::Save(1), cpu);
    RooFitResult *resmc = 0;

    std::cout << "PARAMETERS FOR " << name << std::endl;
    w->var("mean")->Print("");
    if (w->var("sigma")) w->var("sigma")->Print("");

    c1->cd(); c1->Clear();

    w->var("zmll")->SetTitle("dimuon mass (GeV)");
    RooPlot *frame = w->var("zmll")->frame();
    data->plotOn(frame);
    if ((prefix+name).Contains("unbin")) {
        GioProjWData(*w->pdf("signal"), *w->var("zmll"), dynamic_cast<RooDataSet &>(*data)).plotOn(frame,RooFit::LineColor(kBlue));
    } else {
        GioProjWData(*w->pdf("signal"), *w->var("zmll"), *data).plotOn(frame,RooFit::LineColor(kBlue));
    }
    if (TString(binName+"_"+prefix+name).Contains("noEBE")) {
        if (TString(binName+"_"+prefix+name).Contains("jpsi")) {
            w->pdf("signal")->plotOn(frame,RooFit::LineColor(4), RooFit::Components("jpsi"));
        } else if (TString(binName+"_"+prefix+name).Contains("ups")) {
            w->pdf("signal")->plotOn(frame,RooFit::LineColor(4), RooFit::Components("ups1s"));
        } else if (TString(binName+"_"+prefix+name).Contains("zMuMu")) {
            w->pdf("signal")->plotOn(frame,RooFit::LineColor(4), RooFit::Components("zed"));
            TString outRes = TString::Format("plots/%d_v7_all/fit2D_ptEta_%s_%s%s.root",gYear,binName.Data(),prefix.Data(),name.Data());
            outRes.ReplaceAll("rereco","reco"); outRes.ReplaceAll("reco","mc");
            TFile *mc = TFile::Open(outRes);
            if (mc && ( mc->Get("res") != 0)) {
                resmc = (RooFitResult *) mc->Get("res");
                RooArgSet params(*w->pdf("zed")->getParameters(data));
                params = resmc->floatParsFinal();
                w->pdf("signal")->plotOn(frame,RooFit::LineColor(2), RooFit::LineWidth(2));
                params = res->floatParsFinal();
            }
        }
    }
    frame->Draw();
    std::cout << "Max: " << frame->GetMaximum() << std::endl;
    double xoffs = 0;
    if (frame->GetMaximum() >= 10000) {
        c1->SetLeftMargin(0.21); xoffs = 0.05;
        frame->Draw();
        frame->GetYaxis()->SetTitleOffset(1.7);
    } else if (frame->GetMaximum() >= 1000) {
        c1->SetLeftMargin(0.18); xoffs = 0.02;
        frame->Draw();
        frame->GetYaxis()->SetTitleOffset(1.45);
    } else {
        c1->SetLeftMargin(0.16); // default
        frame->Draw();
        frame->GetYaxis()->SetTitleOffset(1.25);
    }
    c1->SetRightMargin(0.05); 
    frame->GetXaxis()->SetNdivisions(505);   
    frame->GetXaxis()->SetDecimals(1);   

    cmsprelim(.15+xoffs, .955, .40, .998, 12, (name.Contains("mc") ? "CMS Simulation": "CMS Preliminary")); 
    cmsprelim(.48+xoffs, .955, .97, .998, 32, Form("#sqrt{s} = %d TeV, L = %.1f fb^{-1}", (gYear+7-2011), (gYear==2011?5.0:19.5)));
    //if (label != "") cmsprelim(0.28, 0.84, 0.45, 0.91, 22, label.Data(), 0.04);
    if (name.Contains("zMuMu")) {
       if (labelLeft != "")  cmsprelim(0.18+xoffs, 0.86, 0.45, 0.91, 11, labelLeft.Data(),  0.042);
       if (labelRight != "") cmsprelim(0.18+xoffs, 0.81, 0.45, 0.85, 11, labelRight.Data(), 0.042);
    } else {
       if (labelLeft != "")  cmsprelim(0.64+xoffs, 0.86, 0.92, 0.91, 31, labelLeft.Data(),  0.042);
       if (labelRight != "") cmsprelim(0.64+xoffs, 0.81, 0.92, 0.85, 31, labelRight.Data(), 0.042);
    }

#if 0 // OLD
    if (TString(binName+"_"+prefix+name).Contains("noEBE")) {
        frame->GetYaxis()->SetNoExponent(1);
        frame->GetYaxis()->SetTitle("Events");
        frame->GetYaxis()->SetTitleOffset(1.45);
        if (resmc) {
            cmsprelim(.67,.87,.92,.92,12,Form("Data:"), 0.045, 4);
            cmsprelim(.67,.82,.92,.87,12,Form("#Delta/M = %+0.2f%%", w->var("mean")->getVal()/w->var("MREF")->getVal()*100.), 0.045, 4);
            cmsprelim(.67,.77,.92,.82,12,Form("#sigma/M =   %4.2f%%", w->var("sigma")->getVal()), 0.045, 4);
            cmsprelim(.20,.87,.42,.92,12,Form("Sim.:"), 0.045, 2);
            cmsprelim(.20,.82,.42,.87,12,Form("#Delta/M = %+0.2f%%", RooArgSet(resmc->floatParsFinal()).getRealValue("mean")/w->var("MREF")->getVal()*100.), 0.045, 2);
            cmsprelim(.20,.77,.42,.82,12,Form("#sigma/M =   %4.2f%%", RooArgSet(resmc->floatParsFinal()).getRealValue("sigma")), 0.045, 2);
        } else {
            cmsprelim(.67,.87,.92,.92,12,Form("#Delta/M = %+0.2f%%", w->var("mean")->getVal()/w->var("MREF")->getVal()*100.), 0.045);
            cmsprelim(.67,.82,.92,.87,12,Form("#sigma/M =   %4.2f%%", w->var("sigma")->getVal()), 0.045);
        }
    }
#endif
    TString fullName = Form("plots/%d_v7_all/fit2D_ptEta_%s_",gYear,binName.Data())+prefix+name;
    c1->Print(fullName+".pdf");
    c1->Print(fullName+".eps");
    TString convOpt = "-q  -dBATCH -dSAFER  -dNOPAUSE  -dAlignToPixels=0 -dEPSCrop  -dPrinted -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -sDEVICE=png16m";
    TString convCmd = Form("gs %s -sOutputFile=%s.png -q \"%s.eps\" -c showpage -c quit", convOpt.Data(), fullName.Data(), fullName.Data());
    gSystem->Exec(convCmd);

#if 0
    if ((prefix+name).Contains("_noEBE") && !(prefix+name).Contains("_unbin")) {
        RooPlot *frame = w->var("logdm")->frame();
        data->plotOn(frame);
    }
#endif

    if (!(prefix+name).Contains("unbin")) {
        FILE *mylog = fopen(Form("plots/%d_v7_all/fit2D_ptEta_%s_%s%s.txt",gYear,binName.Data(),prefix.Data(),name.Data()), "w");
        fprintf(mylog,"%s\t%s\t%g\t%g\t%g\t%g\n",name.Data(),binName.Data(),w->var("mean"  )->getVal(), w->var("mean"  )->getError(), w->var("sigma"  )->getVal(), w->var("sigma"  )->getError());
        fclose(mylog); 
        TFile *thisOut = TFile::Open(Form("plots/%d_v7_all/fit2D_ptEta_%s_%s%s.root",gYear,binName.Data(),prefix.Data(),name.Data()), "RECREATE");
        std::cout << Form("plots/%d_v7_all/fit2D_ptEta_%s_%s%s.root",gYear,binName.Data(),prefix.Data(),name.Data()) << std::endl;
        thisOut->WriteTObject(res,"res");
        thisOut->Close();
    } else {
        TFile *thisOut = TFile::Open(Form("plots/%d_v7_all/fit2D_ptEta_%s_%s%s.root",gYear,binName.Data(),prefix.Data(),name.Data()), "RECREATE");
        std::cout << "tried to save on " << thisOut->GetName() << std::endl;
        std::cout << Form("plots/%d_v7_all/fit2D_ptEta_%s_%s%s.root",gYear,binName.Data(),prefix.Data(),name.Data()) << std::endl;
        TProfile2D *prof = dmScalePlotErr(res,name);
        thisOut->WriteTObject(prof);
        thisOut->WriteTObject(res);
        thisOut->Close();
        FILE *mylog = fopen(Form("plots/%d_v7_all/fit2D_ptEta_%s_%s%s.txt",gYear,binName.Data(),prefix.Data(),name.Data()), "w");
        RooLinkedListIter it = res->floatParsFinal().iterator();
        for (RooRealVar *v = (RooRealVar *) it.Next(); v != 0; v = (RooRealVar *) it.Next()) {
            fprintf(mylog,"%s\t%s\t%s\t%g\t%g\n",name.Data(),binName.Data(),v->GetName(), v->getVal(), v->getError());
        }
        fclose(mylog);
    }
}

void fitMass2D(TH2* mass, TString name, TString binName) { 
    using namespace RooFit;
    RooDataHist *data = new RooDataHist("hzmll","zmll",RooArgList(*w->var("zmll"), *w->var("logdm")), mass);
    fitMass(data, name, binName);
    delete data;
}


void fitMassRochesterUnbin(TTree *tree, TString res, TString name, double ptMin, double ptMax, double absEtaMin, double absEtaMax, Long64_t nmax=-1, TString binName="all2D") {
#if 0
    RooRealVar *vzmll =  w->var("zmll");
    double mMin = w->var("zmll")->getMin(), mMax = w->var("zmll")->getMax();
    RooRealVar *pt1 = w->var("pt1"), *eta1 = w->var("eta1"), *dm1 = w->var("dm1");
    RooRealVar *pt2 = w->var("pt2"), *eta2 = w->var("eta2"), *dm2 = w->var("dm2");
    RooArgSet obs(*vzmll, *pt1,*eta1,*dm1, *pt2,*eta2,*dm2);
    RooDataSet *draw  = new RooDataSet("draw","draw",obs);
    RooDataSet *dcorr = new RooDataSet("dcorr","dcorr",obs);
    double mREF = w->var("MREF")->getVal();
    Float_t zmll, l1pt, l1eta, l1phi, l1charge, l2pt, l2eta, l2phi, l2charge, mErr1, mErr2; UInt_t run; Int_t pf = 1, tagPf = 1;

    tree->SetBranchAddress("mass", &zmll);
    tree->SetBranchAddress("tag_pt", &l1pt); tree->SetBranchAddress("tag_eta", &l1eta); tree->SetBranchAddress("tag_phi", &l1phi); 
    tree->SetBranchAddress("pt", &l2pt); tree->SetBranchAddress("eta", &l2eta); tree->SetBranchAddress("phi", &l2phi); 
    tree->SetBranchAddress("pair_massErr_tag",   &mErr1);
    tree->SetBranchAddress("pair_massErr_probe", &mErr2);
    tree->SetBranchAddress("tag_charge", &l1charge);
    tree->SetBranchAddress("charge", &l2charge);
    tree->SetBranchAddress("tag_PF", &tagPf); tree->SetBranchAddress("PF", &pf);
    tree->SetBranchAddress("run", &run);

    TLorentzVector l1, l2;   // leptons
    for (Long64_t i = 0, n = tree->GetEntries(); i < n && i != nmax; ++i) {
        if ((i+1) % 100000 == 0) std::cerr << "  entry " << (i+1) << "/" << (n+1) << " of " << name << std::endl;
        tree->GetEntry(i);
        if (l1pt < 5  || l2pt < 5)  continue;
        if (!pf || !tagPf) continue;
        if (l1pt < ptMin || l2pt < ptMin || l1pt > ptMax || l2pt > ptMax) continue;
        if (fabs(l1eta) < absEtaMin || fabs(l2eta) < absEtaMin || fabs(l1eta) > absEtaMax || fabs(l2eta) > absEtaMax) continue;
        if (zmll < mMin || zmll > mMax) continue;
        if (accRej != 0 && gRandom->Rndm() < accRej->Eval(0.5*(l1pt+l2pt), 0.5*(fabs(l1eta)+fabs(l2eta)))) continue;
        pt1->setVal(l1pt); eta1->setVal(fabs(l1eta)); dm1->setVal(mErr1/mREF);
        pt2->setVal(l2pt); eta2->setVal(fabs(l2eta)); dm2->setVal(mErr2/mREF);

        l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, 0.105);
        l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, 0.105);
        if (gYear == 2011) {
            if (run == 1) {
                corrector.momcor_mc(l1, l1charge < 0 ? -1 : +1, 0.,  0);
                corrector.momcor_mc(l2, l2charge < 0 ? -1 : +1, 0.,  0);
            } else {
                corrector.momcor_data(l1, l1charge < 0 ? -1 : +1, 0., run <= 173692 ? 0 : 1 );
                corrector.momcor_data(l2, l2charge < 0 ? -1 : +1, 0., run <= 173692 ? 0 : 1 );
            }
        } else {
            if (run == 1) {
                corrector2012.momcor_mc(l1, l1charge < 0 ? -1 : +1, 0.,  0);
                corrector2012.momcor_mc(l2, l2charge < 0 ? -1 : +1, 0.,  0);
            } else {
                corrector2012.momcor_data(l1, l1charge < 0 ? -1 : +1, 0., 0);
                corrector2012.momcor_data(l2, l2charge < 0 ? -1 : +1, 0., 0);
            }
        }
        vzmll->setVal(zmll);        draw->add(obs);
        vzmll->setVal((l1+l2).M()); dcorr->add(obs);
    }

    std::cout << "Will perform an unbinned fit with " << dcorr->numEntries() << " events." << std::endl;
    std::cout << "  .... go grab a cup of tea, it will take some time." << std::endl;
    //fitMass( draw,  Form("%s_ptAvg_%.1f_%.1f_etaAvg_%.1f_%.1f_raw_%s",  res.Data(), ptMin, ptMax, absEtaMin, absEtaMax, name.Data()), binName);
    fitMass( dcorr, Form("%s_ptAvg_%05.1f_%05.1f_etaAvg_%.1f_%.1f_corr_%s", res.Data(), ptMin, ptMax, absEtaMin, absEtaMax, name.Data()), binName);
    delete draw;
    delete dcorr;
#endif
}

void fitMassRochester(TTree *tree, TString res, TString name, double ptMin, double ptMax, double absEtaMin, double absEtaMax, Long64_t nmax=-1, TString binName="all2D") {
    labelRight = TString::Format("%.1f < |#kern[0.1]{#eta}(#mu)| < %.1f",absEtaMin,absEtaMax);
    labelLeft  = TString::Format("%.0f < p_{T}(#mu) < %.0f GeV",ptMin,ptMax);

    int nbins = w->var("zmll")->getBins(), nbinsy = w->var("logdm")->getBins();
    double mMin = w->var("zmll")->getMin(), mMax = w->var("zmll")->getMax();
    double dmMin = w->var("logdm")->getMin(), dmMax = w->var("logdm")->getMax();
    double mREF = w->var("MREF")->getVal();
    TH2F *hzm2d = new TH2F("hzm2d","hzm2d",nbins,mMin,mMax,nbinsy,dmMin,dmMax);
    TH2F *hzum2d = new TH2F("hzum2d","hzm2d",nbins,mMin,mMax,nbinsy,dmMin,dmMax);

    Float_t zmll, l1pt, l1eta, l1phi, l1charge, l2pt, l2eta, l2phi, l2charge, mErr1, mErr2; UInt_t run; Int_t pf = 1, tagPf = 1;
    tree->SetBranchAddress("mass", &zmll);
    tree->SetBranchAddress("tag_pt", &l1pt); tree->SetBranchAddress("tag_eta", &l1eta); tree->SetBranchAddress("tag_phi", &l1phi); 
    tree->SetBranchAddress("pt", &l2pt); tree->SetBranchAddress("eta", &l2eta); tree->SetBranchAddress("phi", &l2phi); 
    tree->SetBranchAddress("pair_massErr_tag",   &mErr1);
    tree->SetBranchAddress("pair_massErr_probe", &mErr2);
    tree->SetBranchAddress("tag_charge", &l1charge);
    tree->SetBranchAddress("charge", &l2charge);
    tree->SetBranchAddress("tag_PF", &tagPf); tree->SetBranchAddress("PF", &pf);

    tree->SetBranchAddress("run", &run);
    TLorentzVector l1, l2;   // leptons
    bool doCorr = !res.Contains("zEE");
    bool isPtAvg = binName.Contains("ptAvg");
    for (Long64_t i = 0, n = tree->GetEntries(); i < n && i != nmax; ++i) {
        if ((i+1) % 100000 == 0) std::cerr << "  entry " << (i+1) << "/" << (n+1) << " of " << name << std::endl;
        tree->GetEntry(i);
        if (isPtAvg) { // cut on average pT
            if (0.5*(l1pt+l2pt) < ptMin || 0.5*(l1pt+l2pt) > ptMax) continue;
        } else {
            if (l1pt < ptMin || l2pt < ptMin || l1pt > ptMax || l2pt > ptMax) continue;
        }
        if (fabs(l1eta) < absEtaMin || fabs(l2eta) < absEtaMin || fabs(l1eta) > absEtaMax || fabs(l2eta) > absEtaMax) continue;
        if (l1pt < 5  || l2pt < 5)  continue;
        if (!pf || !tagPf) continue;
        double weight = 1.0;
        l1.SetPtEtaPhiM(l1pt, l1eta, l1phi, 0.105);
        l2.SetPtEtaPhiM(l2pt, l2eta, l2phi, 0.105);
        if (doCorr) {
            if (gYear == 2011) {
                if (run == 1) {
                    musclecor_mc.applyPtCorrection(l1, l1charge < 0 ? -1 : +1);
                    musclecor_mc.applyPtCorrection(l2, l2charge < 0 ? -1 : +1);
                } else {
                    musclecor_data.applyPtCorrection(l1, l1charge < 0 ? -1 : +1);
                    musclecor_data.applyPtCorrection(l2, l2charge < 0 ? -1 : +1);
                }
            } else {
                if (run == 1) {
                    musclecor2012_mc.applyPtCorrection(l1, l1charge < 0 ? -1 : +1);
                    musclecor2012_mc.applyPtCorrection(l2, l2charge < 0 ? -1 : +1);
                    MUSCLE2012_SMEARING(musclecor2012_mc, l1, l1charge < 0 ? -1 : +1);
                    MUSCLE2012_SMEARING(musclecor2012_mc, l2, l2charge < 0 ? -1 : +1);
                } else {
                    if (run >= 203773) {
                        musclecor2012_dataD.applyPtCorrection(l1, l1charge < 0 ? -1 : +1);
                        musclecor2012_dataD.applyPtCorrection(l2, l2charge < 0 ? -1 : +1);
                    } else {
                        musclecor2012_dataABC.applyPtCorrection(l1, l1charge < 0 ? -1 : +1);
                        musclecor2012_dataABC.applyPtCorrection(l2, l2charge < 0 ? -1 : +1);
                    }
                }
            }
        }
        if (ebeCorr != 0) {
            mErr1 *= ebeCorr->Eval(l1pt, fabs(l1eta));
            mErr2 *= ebeCorr->Eval(l2pt, fabs(l2eta));
        }
        if (ebeCorrH != 0) {
            mErr1 *= ebeCorrH->GetBinContent(ebeCorrH->FindBin(l1pt, fabs(l1eta)));
            mErr2 *= ebeCorrH->GetBinContent(ebeCorrH->FindBin(l2pt, fabs(l2eta)));
        }
        double mErr = hypot(mErr1,mErr2);
        double logdm = std::max(dmMin+1e-5,std::min(dmMax-1e-5, log(mErr/mREF)));
        hzum2d->Fill( zmll, logdm, weight );
        hzm2d->Fill( (l1+l2).M(), logdm, weight );
    }

    //fitMass2D( hzum2d,  Form("%s_ptAvg_%.1f_%.1f_etaAvg_%.1f_%.1f_raw_%s",  res.Data(), ptMin, ptMax, absEtaMin, absEtaMax, name.Data()), binName);
    fitMass2D( hzm2d,   Form("%s_ptAvg_%05.1f_%05.1f_etaAvg_%.1f_%.1f_corr_%s", res.Data(), ptMin, ptMax, absEtaMin, absEtaMax, name.Data()), binName);
}


void fitMassCaliforniaPtEtaBinsEBEv2(int data=1, int res=0, int energy=7, int ijob = 0, int maxev=-1, const char *local=0) {
    prefix = "";
    initCanvas(); 
    gStyle->SetOptTitle(0);
    TString sRes = ((res % 100) == 0 ? "zMuMu" : "jPsi");
    if ((res % 100) == 2) sRes = "ups"; 
    if ((res % 100) == 4) sRes = "zEE"; 
    if (res / 100 == 1) sRes += "_unbin"; 
    if (res / 100 == 2) sRes += "_unbinPt"; 
    if (res / 100 == 3) sRes += "_noEBE"; 
    gYear = (energy == 7 ? 2011 : 2012);

    if ((res % 100) == 0) {
        if (res / 100 == 0) initModelZ(0,0);
        else if (res / 100 == 1) initModelZ(1,0);
        else if (res / 100 == 2) initModelZ(1,1);
        else if (res / 100 == 3) initModelZ(-1,0,true);
        TChain *reco = new TChain("tpTree/fitter_tree");
        TString name = "";
        switch(10*energy+data) {
            case 70:
                name = "mc42x";
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2", "zMuMu_lineshape_MC_42X_DYJets.root"));
                break;
            case 71:
                name = "reco42x";
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2", "zMuMu_lineshape_Data_2011A_16Jan2012.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2", "zMuMu_lineshape_Data_2011B_16Jan2012.root"));
                break;
            case 80:
                name = "mc53x";
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V3", "zMuMu_lineshape_MC_powheg.root"));
                break;
            case 81:
                name = "reco53x";
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_DoubleMu_Run2012A.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_DoubleMuParked_Run2012B.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_DoubleMuParked_Run2012C.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_DoubleMuParked_Run2012D.root"));
                break;
        }
        if (res/100 == 0 || res/100 == 3) {
            TFile *fIn = TFile::Open(Form("mcPull_%d_ptEta_ALLFine.root",gYear));
            TString bin = "all";
            if (ijob / 30 == 1) { ebeCorrH = (TH2*) fIn->Get("ALLFine_rms"); bin = "all_mcPull"; }
            if (ijob / 30 == 2) { ebeCorrH = (TH2*) fIn->Get("ALLFine_sigma"); bin = "all_mcPullSig"; }
            switch(ijob % 30) {
                case 0:
                case 1:
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 0.4, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 0.5, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 0.6, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 0.7, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 0.8, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 1.0, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 1.2, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 1.3, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 1.4, maxev, bin);
                    if (ijob % 30) break;
                case 2:
                    fitMassRochester(reco, sRes, name, 30., 80., 0.5, 1.2, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.6, 1.2, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.7, 1.2, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.8, 1.2, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.9, 1.2, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.0, 1.2, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.1, 1.5, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.2, 1.5, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.2, 1.6, maxev, bin);
                    if (ijob % 30) break;
                 case 3:
                    fitMassRochester(reco, sRes, name, 30., 80., 1.2, 2.4, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 2.4, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.6, 2.4, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.7, 2.4, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.8, 2.4, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.9, 2.4, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 2.0, 2.4, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 2.1, 2.4, maxev, bin);
                    if (ijob % 30) break;
                case 4:
                    fitMassRochester(reco, sRes, name, 30., 80., 0.8, 1.4, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.8, 1.5, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.8, 1.6, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.8, 1.7, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 1.8, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 1.9, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 2.0, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 2.1, maxev, bin);
                    if (ijob % 30) break;
                case 5:
                    fitMassRochester(reco, sRes, name, 5., 99., 0.0, 1.2, maxev, bin);
                    fitMassRochester(reco, sRes, name, 5., 99., 1.2, 2.4, maxev, bin);
                    if (ijob % 30) break;
                case 11:
                    fitMassRochester(reco, sRes, name, 5., 99., 0.0, 0.8, maxev, bin);
                    fitMassRochester(reco, sRes, name, 5., 99., 0.8, 1.6, maxev, bin);
                    fitMassRochester(reco, sRes, name, 5., 99., 1.6, 2.4, maxev, bin);
                    if (ijob % 30) break;
                case 12:
                    fitMassRochester(reco, sRes, name, 20., 99., 0.0, 0.8, maxev, bin);
                    fitMassRochester(reco, sRes, name, 20., 99., 0.8, 1.6, maxev, bin);
                    fitMassRochester(reco, sRes, name, 20., 99., 1.6, 2.4, maxev, bin);
                    if (ijob % 30) break;
                case 13:
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 0.8, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 0.8, 1.6, maxev, bin);
                    fitMassRochester(reco, sRes, name, 30., 80., 1.6, 2.4, maxev, bin);
                    if (ijob % 30) break;
            }
            fIn->Close(); ebeCorrH = 0;
        } else {
            switch(ijob) {
                case 0:
                case 1:
                    fitMassRochesterUnbin(reco, sRes, name, 30., 80., 0.0, 0.8, maxev, "all");
                    if (ijob) break;
                case 2:
                    fitMassRochesterUnbin(reco, sRes, name, 30., 80., 0.0, 2.4, maxev, "all");
                    if (ijob) break;
                case 3:
                    fitMassRochesterUnbin(reco, sRes, name, 5.0, 100., 0.0, 2.4, maxev, "all");
                    if (ijob) break;
                case 4:
                    fitMassRochesterUnbin(reco, sRes, name, 5.0, 100., 0.0, 0.8, maxev, "all");
                    if (ijob) break;
            }
        }
    } else if ((res % 100) == 1) {
        if (res / 100 == 0) initModelJ(0,0);
        else if (res / 100 == 1) initModelJ(1,0);
        else if (res / 100 == 2) initModelJ(1,1);
        else if (res / 100 == 3) initModelJ(-1,0);
        TChain *reco = new TChain("tpTree/fitter_tree");
        TString name = "";
        switch(10*energy+data) {
            case 70:
                name = "mc42x";
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2/jpsi_lineshape_MC_42X.root");
                break;
            case 71:
                name = "reco42x";
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2/oniaMuMu_lineshape_Data_2011A_05Aug2011.root");
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2/oniaMuMu_lineshape_Data_2011A_10May2011.root");
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2/oniaMuMu_lineshape_Data_2011Av4.root");
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2/oniaMuMu_lineshape_Data_2011Av6.root");
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2/oniaMuMu_lineshape_Data_2011B.root");
                break;
            case 80:
                name = "mc53x";
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2", "jPsi_lineshape_MC_53X.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2", "jPsi_lineshape_MC_53X_Bs.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2", "jPsi_lineshape_MC_53X_Bu.root"));
                break;
            case 81:
                name = "reco53x";
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_MuOnia_Run2012A.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_MuOnia_Run2012B.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_MuOnia_Run2012C.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_MuOnia_Run2012D.root"));
                break;
        }
        if (res / 100 == 0 || res/100 == 3) {
            TString mypostfix="";
            if (ijob / 30 == 1) {
                TFile *fIn = TFile::Open(Form("mcPull_%d_ptEta_ALLFine.root",gYear));
                ebeCorrH = (TH2*) fIn->Get("ALLFine_rms");
                mypostfix = "_mcPull";
            }
            switch(ijob % 30) {
                case 0:
                case 1:
                    fitMassRochester(reco, sRes, name, 5., 10., 0.0, 0.2, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 0.2, 0.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 0.5, 0.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 0.8, 1.2, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 1.2, 1.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 1.5, 1.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 1.8, 2.1, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 2.1, 2.5, maxev, "all"+mypostfix);
                    if (ijob % 30) break;
                case 2:
                    fitMassRochester(reco, sRes, name, 10., 20., 0.0, 0.2, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 0.2, 0.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 0.5, 0.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 0.8, 1.2, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 1.2, 1.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 1.5, 1.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 1.8, 2.1, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 2.1, 2.5, maxev, "all"+mypostfix);
                    if (ijob % 30) break;
                case 3:
                    fitMassRochester(reco, sRes, name, 20., 100., 0.0, 0.2, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 0.2, 0.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 0.5, 0.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 0.8, 1.2, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 1.2, 1.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 1.5, 1.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 1.8, 2.1, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 2.1, 2.5, maxev, "all"+mypostfix);
                    if (ijob % 30) break;
                case 7:
                    fitMassRochester(reco, sRes, name, 5.,  20.,  0.0, 1.2, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5.,  20.,  1.2, 2.4, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 0.0, 1.2, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 1.2, 2.4, maxev, "all"+mypostfix);
                    if (ijob % 30) break;
                case 8:
                    fitMassRochester(reco, sRes, name, 5.,  20.,  0.0, 0.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 0.0, 0.5, maxev, "all"+mypostfix);
                    if (ijob % 30) break;
                case 9:
                    fitMassRochester(reco, sRes, name, 5.,  20.,  0.5, 1.2, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 0.5, 1.2, maxev, "all"+mypostfix);
                    if (ijob % 30) break;
                case 10:
                    fitMassRochester(reco, sRes, name, 5.,  20.,  0.0, 0.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5.,  20.,  0.8, 1.6, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5.,  20.,  1.6, 2.4, maxev, "all"+mypostfix);
                    if (ijob % 30) break;
                case 11:
                    fitMassRochester(reco, sRes, name, 20., 100., 0.0, 0.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 0.8, 1.6, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 1.6, 2.4, maxev, "all"+mypostfix);
                    if (ijob % 30) break;

            }
        } else {
           switch(ijob) {
                case 0:
                case 1:
                    fitMassRochesterUnbin(reco, sRes, name, 5., 10.,  0.0, 2.4, maxev, "all");
                    if (ijob) break;
                case 2:
                    fitMassRochesterUnbin(reco, sRes, name, 10., 20., 0.0, 2.4, maxev, "all");
                    if (ijob) break;
                case 3:
                    fitMassRochesterUnbin(reco, sRes, name, 20., 80.0, 0.0, 2.4, maxev, "all");
                    if (ijob) break;
                case 20:
                    accRej = new TF2("accRej","pow(1-10/(x+5),2)+0*y",0,200,0,2.5);
                    fitMassRochesterUnbin(reco, sRes, name, 5., 80.0, 0.0, 2.4, maxev, "accRej");
                    if (ijob) break;
           }
        }
    } else if ((res % 100) == 2) {
        if (res / 100 == 0) initModelU(!data,0,0);
        else if (res / 100 == 1) initModelU(!data,1,0);
        else if (res / 100 == 2) initModelU(!data,1,1);
        else if (res / 100 == 3) initModelU(!data,-1,0);
        TChain *reco = new TChain("tpTree/fitter_tree");
        TString name = "";
        switch(10*energy+data) {
            case 70:
                name = "mc42x";
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2/upMuMu_lineshape_MC_42X.root");
                break;
            case 71:
                name = "reco42x";
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2/oniaMuMu_lineshape_Data_2011A_05Aug2011.root");
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2/oniaMuMu_lineshape_Data_2011A_10May2011.root");
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2/oniaMuMu_lineshape_Data_2011Av4.root");
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2/oniaMuMu_lineshape_Data_2011Av6.root");
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2/oniaMuMu_lineshape_Data_2011B.root");
                break;
            case 80:
                name = "mc53x";
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V2", "upMuMu_lineshape_MC_53X.root"));
                break;
            case 81:
                name = "reco53x";
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_MuOnia_Run2012A.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_MuOnia_Run2012B.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_MuOnia_Run2012C.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_MuOnia_Run2012D.root"));
                break;
        }
        if (res / 100 == 0 || res/100 == 3) {
            TString mypostfix="";
            if (ijob / 30 == 1) {
                TFile *fIn = TFile::Open(Form("mcPull_%d_ptEta_ALLFine.root",gYear));
                ebeCorrH = (TH2*) fIn->Get("ALLFine_rms");
                mypostfix = "_mcPull";
            }
            switch(ijob % 30) {
                case 0:
                case 1:
                    fitMassRochester(reco, sRes, name, 5., 10., 0.0, 0.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 0.5, 0.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 0.8, 1.2, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 1.2, 1.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 1.5, 1.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 1.8, 2.1, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5., 10., 2.1, 2.5, maxev, "all"+mypostfix);
                    if (ijob % 30) break;
                case 2:
                    fitMassRochester(reco, sRes, name, 10., 20., 0.0, 0.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 0.5, 0.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 0.8, 1.2, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 1.2, 1.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 1.5, 1.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 1.8, 2.1, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 10., 20., 2.1, 2.5, maxev, "all"+mypostfix);
                    if (ijob % 30) break;
                case 3:
                    fitMassRochester(reco, sRes, name, 20., 100., 0.0, 0.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 0.5, 0.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 0.8, 1.2, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 1.2, 1.5, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 1.5, 1.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 1.8, 2.1, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 2.1, 2.5, maxev, "all"+mypostfix);
                    if (ijob % 30) break;
                case 10:
                    fitMassRochester(reco, sRes, name, 5.,  20.,  0.0, 0.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5.,  20.,  0.8, 1.6, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 5.,  20.,  1.6, 2.4, maxev, "all"+mypostfix);
                    if (ijob % 30) break;
                case 11:
                    fitMassRochester(reco, sRes, name, 20., 100., 0.0, 0.8, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 0.8, 1.6, maxev, "all"+mypostfix);
                    fitMassRochester(reco, sRes, name, 20., 100., 1.6, 2.4, maxev, "all"+mypostfix);
                    if (ijob % 30) break;

            }
        } else {
           switch(ijob) {
                case 0:
                case 1:
                    fitMassRochesterUnbin(reco, sRes, name, 5., 10.,  0.0, 2.4, maxev, "all");
                    if (ijob) break;
                case 2:
                    fitMassRochesterUnbin(reco, sRes, name, 10., 20., 0.0, 2.4, maxev, "all");
                    if (ijob) break;
                case 3:
                    fitMassRochesterUnbin(reco, sRes, name, 20., 80.0, 0.0, 2.4, maxev, "all");
                    if (ijob) break;
                case 20:
                    accRej = new TF2("accRej","pow(1-10/(x+5),2)+0*y",0,200,0,2.5);
                    fitMassRochesterUnbin(reco, sRes, name, 5., 80.0, 0.0, 2.4, maxev, "accRej");
                    if (ijob) break;
           }
        }

    } else if ((res % 100) == 4) {
        if (res / 100 == 0) initModelZ(0,0);
        else if (res / 100 == 1) initModelZ(1,0);
        else if (res / 100 == 2) initModelZ(1,1);
        else if (res / 100 == 3) initModelZ(-1,0,true);
        TChain *reco = new TChain("tpTree/fitter_tree");
        TString name = "";
        switch(10*energy+data) {
            case 70:
                name = "mc42x";
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4/lineshape_DYJets_EE_42X_16Jan2012_LegacyPaper.root");
                break;
            case 71:
                name = "reco42x";
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4/lineshape_DoubleElectron_Run2011A.root");
                reco->AddFile("root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4/lineshape_DoubleElectron_Run2011B.root");
                break;
            case 80:
                name = "mc53x";
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_DYJets_EE_Summer2012_LegacyPaper.root"));
                break;
            case 81:
                name = "reco53x";
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_DoubleElectron_Run2012A.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_DoubleElectron_Run2012B.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_DoubleElectron_Run2012C.root"));
                reco->AddFile(Form("%s/%s", local ? local : "root://eoscms//eos/cms/store/caf/user/gpetrucc/muonScale/V4", "lineshape_DoubleElectron_Run2012D.root"));
                break;
        }
        if (res/100 == 0 || res/100 == 3) {
            switch(ijob) {
                case 0:
                case 1:
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 1.5, maxev, "all");
                    if (ijob) break;
                case 2:
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 2.5, maxev, "all");
                    if (ijob) break;
                case 3:
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 0.9, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 0.9, 1.5, maxev, "all");
                    if (ijob) break;
                case 4:
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 2.0, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 2.0, 2.5, maxev, "all");
                    if (ijob) break;
                case 5:
                    fitMassRochester(reco, sRes, name, 0., 200., 0.0, 1.5, maxev, "all");
                    if (ijob) break;
                case 6:
                    fitMassRochester(reco, sRes, name, 0., 200., 1.5, 2.5, maxev, "all");
                    if (ijob) break;
                case 7:
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 0.4, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 0.6, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 0.8, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 1.0, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 1.2, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 1.4, maxev, "all");
                    if (ijob) break;
                case 8:
                    fitMassRochester(reco, sRes, name, 30., 80., 0.8, 1.5, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 0.9, 1.5, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 1.0, 1.5, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 1.1, 1.5, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 1.2, 1.5, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 1.3, 1.5, maxev, "all");
                    if (ijob) break;
                case 9:
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 1.7, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 1.8, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 1.9, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 2.0, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 2.1, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 2.2, maxev, "all");
                    if (ijob) break;
                case 10:
                    fitMassRochester(reco, sRes, name, 30., 80., 1.7, 2.5, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 1.8, 2.5, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 1.9, 2.5, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 2.0, 2.5, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 2.1, 2.5, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 2.2, 2.5, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 2.3, 2.5, maxev, "all");
                    if (ijob) break;
                case 12:
                    fitMassRochester(reco, sRes, name, 30., 80., 0.0, 0.8, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 0.8, 1.5, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 1.5, 2.0, maxev, "all");
                    fitMassRochester(reco, sRes, name, 30., 80., 2.0, 2.5, maxev, "all");
            }
        } else {
        }
    } 
}
