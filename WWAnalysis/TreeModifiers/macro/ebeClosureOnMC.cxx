#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TChain.h>
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
#include <TH1.h>
#include <TH2.h>
#include <TF2.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TFile.h>
#include <iostream>
#include <cstdio>
#include <cmath>

int gYear = 2011;

RooWorkspace *w = new RooWorkspace("w","w");
TCanvas *c1 = new TCanvas("c1","c1");
TString prefix = "";
TF2 *ebeCorr  = 0, *ebeCorr4E  = 0;
TH2 *ebeCorrH = 0, *ebeCorrH4E = 0;
TF2 *accRej  = 0;
void initCanvas() {
    gStyle->SetOptStat(0);
    gStyle->SetCanvasDefW(850); //Width of canvas
    gStyle->SetPaperSize(850./600.*20.,20.);
    c1 = new TCanvas("c1","c1");
    c1->cd();
    c1->SetWindowSize(850 + (850 - c1->GetWw()), 600 + (600 - c1->GetWh()));
    c1->SetRightMargin(0.04);
    c1->SetGridy(0); c1->SetGridx(0);
}

void initModel(int mref=125) {
    // for pulls
    w->factory("pull[-8,8]"); w->var("pull")->setBins(160);

    w->factory("Gaussian::pullFit_g(pull, pull_mean[-2,2], pull_sigma[1,0.1,5])");
    w->factory("Voigtian::pullFit_v(pull, pull_mean, pull_bww[0.0041], pull_sigma)");
    w->factory("CBShape::pullFit_cb(pull, pull_mean, pull_sigma, alpha[3.,0.05, 50.], n[2, 1, 100.])");
    //w->factory("SUM::pullFit_dcb(fLeft[0.5,0,1]*pullFit_cb, CBShape::pullFit_cbFlip(pull, pull_mean, pull_sigma, alphaFlip[-3., -50, -0.005], nFlip[1, 0.5, 200.]) )");
    w->factory("DoubleCB::pullFit_dcb2(pull, pull_mean, pull_sigma, alpha, n, alpha2[3., 0.5, 50.], nFlip[2, 0.5, 100.])");

    w->factory(Form("MREF[%d]",mref));
    w->factory("mass[115,135]"); w->var("mass")->setBins(40);
    w->factory("logdm[-6.0,-1.0]"); w->var("logdm")->setBins(80);

    w->factory("CBShape::massFit_cb(mass, sum(MREF,pull_mean), expr::mass_sigma(\"@0*@1*exp(@2)\",MREF,pull_sigma,logdm), alpha, n)");
    w->factory("DoubleCB::massFit_dcb2(mass, sum(MREF,pull_mean), mass_sigma, alpha, n, alpha2, nFlip)");

    //w->var("pull_sigma")->setVal(1.254); w->var("pull_sigma")->setConstant(true);
    w->saveSnapshot("clean", w->allVars());
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

void fitPull(RooAbsData *data, TString pdf, TString name, TString binName) { 
    using namespace RooFit;
    w->loadSnapshot("clean");
    RooFitResult *res = w->pdf(pdf)->fitTo(*data, Save(1));

    RooPlot *frame = w->var("pull")->frame();
    data->plotOn(frame);
    w->pdf(pdf)->plotOn(frame);
    w->pdf(pdf)->paramOn(frame, 
                Format("N", FixedPrecision(4)),
                Layout(.70,.97,.97),
                Parameters(RooArgSet(*w->var("pull_mean"), *w->var("pull_sigma"))));
    frame->Draw();

    c1->Print(Form("plots/%d_v7_all/fitPull_ptEta_%s_%s_",gYear,binName.Data(),pdf.Data())+prefix+name+".png");
    c1->Print(Form("plots/%d_v7_all/fitPull_ptEta_%s_%s_",gYear,binName.Data(),pdf.Data())+prefix+name+".pdf");
    FILE *mylog = fopen(Form("plots/%d_v7_all/fitPull_ptEta_%s_%s_%s%s.txt",gYear,binName.Data(),pdf.Data(),prefix.Data(),name.Data()), "w");
    fprintf(mylog,"%s\t%s\t%g\t%g\t%g\t%g\n",name.Data(),binName.Data(),w->var("pull_mean"  )->getVal(), w->var("pull_mean"  )->getError(), w->var("pull_sigma"  )->getVal(), w->var("pull_sigma"  )->getError());
    fclose(mylog); 
}

void fitMass2D(RooAbsData *data, TString pdf, TString name, TString binName) { 
    using namespace RooFit;
    w->loadSnapshot("clean");
    RooArgSet vars(*data->get());
    vars.remove(*w->var("mass"),false,true);
   
    const RooCmdArg & cpu =  name.Contains("unbin") ? NumCPU(6) :  RooCmdArg::none();
    RooFitResult *res = w->pdf(pdf)->fitTo(*data,ConditionalObservables(vars), Save(1), cpu);

    RooPlot *frame = w->var("mass")->frame();
    data->plotOn(frame);
    GioProjWData(*w->pdf(pdf), *w->var("mass"), *data).plotOn(frame,LineColor(kBlue));
    w->pdf(pdf)->paramOn(frame, 
                Format("N", FixedPrecision(4)),
                Layout(.70,.97,.97),
                Parameters(RooArgSet(*w->var("pull_mean"), *w->var("pull_sigma"))));
    frame->Draw();

    c1->Print(Form("plots/%d_v7_all/fitMass2D_ptEta_%s_%s_",gYear,binName.Data(),pdf.Data())+prefix+name+".png");
    c1->Print(Form("plots/%d_v7_all/fitMass2D_ptEta_%s_%s_",gYear,binName.Data(),pdf.Data())+prefix+name+".pdf");
    FILE *mylog = fopen(Form("plots/%d_v7_all/fitMass2D_ptEta_%s_%s_%s%s.txt",gYear,binName.Data(),pdf.Data(),prefix.Data(),name.Data()), "w");
    fprintf(mylog,"%s\t%s\t%g\t%g\t%g\t%g\n",name.Data(),binName.Data(),w->var("pull_mean"  )->getVal(), w->var("pull_mean"  )->getError(), w->var("pull_sigma"  )->getVal(), w->var("pull_sigma"  )->getError());
    fclose(mylog); 
}



void fitPull(TH1* mass, TString pdf, TString name, TString binName) { 
    using namespace RooFit;
    RooDataHist *data = new RooDataHist("hzmll","zmll",RooArgList(*w->var("pull")), mass);
    fitPull(data, pdf, name, binName);
    delete data;
}
void fitMass2D(TH2* mass, TString pdf, TString name, TString binName) { 
    using namespace RooFit;
    RooDataHist *data = new RooDataHist("hzmll","zmll",RooArgList(*w->var("mass"), *w->var("logdm")), mass);
    fitMass2D(data, pdf, name, binName);
    delete data;
}



void fitPulls(TTree *tree, TString res, TString name, Long64_t nmax=-1, TString binName="all2D") {
    int nbinsp = w->var("pull")->getBins(), nbins = w->var("mass")->getBins(), nbinsy = w->var("logdm")->getBins();
    double pMin = w->var("pull")->getMin(), pMax = w->var("pull")->getMax();
    double mMin = w->var("mass")->getMin(), mMax = w->var("mass")->getMax();
    double dmMin = w->var("logdm")->getMin(), dmMax = w->var("logdm")->getMax();
    TH1D *hzp1d_4l  = new TH1D("hzm_4l","hzm", nbinsp,pMin,pMax);
    TH1D *hzp1d_gen = new TH1D("hzm_gen","hzm",nbinsp,pMin,pMax);
    TH1D *hzp1d_ref = new TH1D("hzm_ref","hzm",nbinsp,pMin,pMax);
    TH2F *hzm2d_4l  = new TH2F("hzm2_4l", "hzm",nbins,mMin,mMax,nbinsy,dmMin,dmMax);
    TH2F *hzm2d_gen = new TH2F("hzm2_gen","hzm",nbins,mMin,mMax,nbinsy,dmMin,dmMax);
    TH2F *hzm2d_ref = new TH2F("hzm2_ref","hzm",nbins,mMin,mMax,nbinsy,dmMin,dmMax);

    Float_t mass, massErr, m4l, gen_mass, gen_m4l; Int_t gen_match;
    Float_t l1pt, l1eta, l1phi, l1pdgId, l1massErr, l2pt, l2eta, l2phi, l2pdgId, l2massErr;
    Float_t l3pt, l3eta, l3phi, l3pdgId, l3massErr, l4pt, l4eta, l4phi, l4pdgId, l4massErr;
    Float_t pho1pt, pho2pt;
    tree->SetBranchAddress("mass", &mass);
    tree->SetBranchAddress("massErr", &massErr);
    tree->SetBranchAddress("m4l", &m4l);
    tree->SetBranchAddress("genhiggsmass", &gen_mass);
    tree->SetBranchAddress("gen_m4l", &gen_m4l);
    tree->SetBranchAddress("mc_4lany", &gen_match);
    tree->SetBranchAddress("l1pt", &l1pt); tree->SetBranchAddress("l1eta", &l1eta); tree->SetBranchAddress("l1phi", &l1phi); tree->SetBranchAddress("l1pdgId", &l1pdgId); tree->SetBranchAddress("l1massErr", &l1massErr);
    tree->SetBranchAddress("l2pt", &l2pt); tree->SetBranchAddress("l2eta", &l2eta); tree->SetBranchAddress("l2phi", &l2phi); tree->SetBranchAddress("l2pdgId", &l2pdgId); tree->SetBranchAddress("l2massErr", &l2massErr);
    tree->SetBranchAddress("l3pt", &l3pt); tree->SetBranchAddress("l3eta", &l3eta); tree->SetBranchAddress("l3phi", &l3phi); tree->SetBranchAddress("l3pdgId", &l3pdgId); tree->SetBranchAddress("l3massErr", &l3massErr);
    tree->SetBranchAddress("l4pt", &l4pt); tree->SetBranchAddress("l4eta", &l4eta); tree->SetBranchAddress("l4phi", &l4phi); tree->SetBranchAddress("l4pdgId", &l4pdgId); tree->SetBranchAddress("l4massErr", &l4massErr);
    tree->SetBranchAddress("pho1pt", &pho1pt); tree->SetBranchAddress("pho2pt", &pho2pt); 
    int reqmu = res.Contains("4mu") ? 4 : (res.Contains("4e") ? 0 : 2);
    int onlyfsr =  res.Contains("FSR");
    int withfsr =  res.Contains("fsr") || onlyfsr;
    int barrel = res.Contains("EB") - res.Contains("EE");
    for (Long64_t i = 0, n = tree->GetEntries(); i < n && i != nmax; ++i) {
        tree->GetEntry(i);
        int nmu = 2*(abs(l1pdgId) == 13) + 2*(abs(l3pdgId) == 13);
        if (nmu != reqmu) continue;
        bool hasfsr = (pho1pt + pho2pt > 0);
        if ((hasfsr && !withfsr) || (!hasfsr && onlyfsr)) continue; 
        if (barrel == 1) {
            if (abs(l1pdgId) == 11 && fabs(l1eta) > 1.5) continue;
            if (abs(l2pdgId) == 11 && fabs(l2eta) > 1.5) continue;
            if (abs(l3pdgId) == 11 && fabs(l3eta) > 1.5) continue;
            if (abs(l4pdgId) == 11 && fabs(l4eta) > 1.5) continue;
        } else if (barrel == -1) {
            if (abs(l1pdgId) == 11 && fabs(l1eta) < 1.5) continue;
            if (abs(l2pdgId) == 11 && fabs(l2eta) < 1.5) continue;
            if (abs(l3pdgId) == 11 && fabs(l3eta) < 1.5) continue;
            if (abs(l4pdgId) == 11 && fabs(l4eta) < 1.5) continue;
        }
        double mErrFSR2 = massErr*massErr - (l1massErr*l1massErr + l2massErr*l2massErr + l3massErr*l3massErr + l4massErr*l4massErr);
        if (ebeCorr != 0) {
            if (abs(l1pdgId) == 13) {
                l1massErr *= ebeCorr->Eval(l1pt, fabs(l1eta));
                l2massErr *= ebeCorr->Eval(l2pt, fabs(l2eta));
            }
            if (abs(l3pdgId) == 13) {
                l3massErr *= ebeCorr->Eval(l3pt, fabs(l3eta));
                l4massErr *= ebeCorr->Eval(l4pt, fabs(l4eta));
            }
        }
        if (ebeCorrH != 0) {
            if (abs(l1pdgId) == 13) {
                l1massErr *= ebeCorrH->GetBinContent(ebeCorrH->FindBin(l1pt, fabs(l1eta)));
                l2massErr *= ebeCorrH->GetBinContent(ebeCorrH->FindBin(l2pt, fabs(l2eta)));
            }
            if (abs(l3pdgId) == 13) {
                l3massErr *= ebeCorrH->GetBinContent(ebeCorrH->FindBin(l3pt, fabs(l3eta)));
                l4massErr *= ebeCorrH->GetBinContent(ebeCorrH->FindBin(l4pt, fabs(l4eta)));
            }
        }
        if (ebeCorr4E != 0) {
            if (abs(l1pdgId) == 11) {
                l1massErr *= ebeCorr4E->Eval(l1pt, fabs(l1eta));
                l2massErr *= ebeCorr4E->Eval(l2pt, fabs(l2eta));
            }
            if (abs(l3pdgId) == 11) {
                l3massErr *= ebeCorr4E->Eval(l3pt, fabs(l3eta));
                l4massErr *= ebeCorr4E->Eval(l4pt, fabs(l4eta));
            }
        }
        if (ebeCorrH4E != 0) {
            if (abs(l1pdgId) == 11) {
                l1massErr *= ebeCorrH4E->GetBinContent(ebeCorrH4E->FindBin(l1pt, fabs(l1eta)));
                l2massErr *= ebeCorrH4E->GetBinContent(ebeCorrH4E->FindBin(l2pt, fabs(l2eta)));
            }
            if (abs(l3pdgId) == 11) {
                l3massErr *= ebeCorrH4E->GetBinContent(ebeCorrH4E->FindBin(l3pt, fabs(l3eta)));
                l4massErr *= ebeCorrH4E->GetBinContent(ebeCorrH4E->FindBin(l4pt, fabs(l4eta)));
            }
        }

        double mErr4l = sqrt(l1massErr*l1massErr + l2massErr*l2massErr + l3massErr*l3massErr + l4massErr*l4massErr);
        double mErr   = sqrt(l1massErr*l1massErr + l2massErr*l2massErr + l3massErr*l3massErr + l4massErr*l4massErr + mErrFSR2);
        double mRef = 125.0; //floor(gen_mass+0.5);
        double pull_4l = (m4l - gen_m4l)/mErr4l, pull_g = (m4l - gen_mass)/mErr, pull_r = (m4l - mRef)/mErr;
        double logdm4l = std::max(dmMin+1e-5,std::min(dmMax-1e-5, log(mErr4l/mRef)));
        double logdm   = std::max(dmMin+1e-5,std::min(dmMax-1e-5, log(mErr  /mRef)));
        hzp1d_4l->Fill(pull_4l);
        hzp1d_gen->Fill(pull_g);
        hzp1d_ref->Fill(pull_r);
        hzm2d_4l->Fill(m4l-(gen_m4l-mRef), logdm4l);
        hzm2d_gen->Fill(mass-(gen_mass-mRef), logdm);
        hzm2d_ref->Fill(mass, logdm);
    }

#if 0
    //fitPull( hzp1d_4l,  "pullFit_g",    Form("4l_%s_%s", res.Data(), name.Data()), binName);
    fitPull( hzp1d_4l,  "pullFit_cb",   Form("4l_%s_%s", res.Data(), name.Data()), binName);
    //fitPull( hzp1d_4l,  "pullFit_dcb",  Form("4l_%s_%s", res.Data(), name.Data()), binName);
    fitPull( hzp1d_4l,  "pullFit_dcb2", Form("4l_%s_%s", res.Data(), name.Data()), binName);
    //fitPull( hzp1d_gen, "pullFit_g",    Form("gen_%s_%s", res.Data(), name.Data()), binName);
    fitPull( hzp1d_gen, "pullFit_cb",   Form("gen_%s_%s", res.Data(), name.Data()), binName);
    //fitPull( hzp1d_gen, "pullFit_dcb",  Form("gen_%s_%s", res.Data(), name.Data()), binName);
    fitPull( hzp1d_gen, "pullFit_dcb2", Form("gen_%s_%s", res.Data(), name.Data()), binName);
    //fitPull( hzp1d_ref, "pullFit_v",    Form("ref_%s_%s", res.Data(), name.Data()), binName);
    //fitPull( hzp1d_ref, "pullFit_dcb",  Form("ref_%s_%s", res.Data(), name.Data()), binName);


#endif
    //fitPull( hzp1d_ref, "pullFit_cb",   Form("ref_%s_%s", res.Data(), name.Data()), binName);
    fitPull( hzp1d_ref, "pullFit_dcb2", Form("ref_%s_%s", res.Data(), name.Data()), binName);
    //fitMass2D( hzm2d_ref, "massFit_cb",   Form("ref_%s_%s", res.Data(), name.Data()), binName);
    fitMass2D( hzm2d_ref, "massFit_dcb2", Form("ref_%s_%s", res.Data(), name.Data()), binName);
    //fitMass2D( hzm2d_4l,  "massFit_cb",   Form("4l_%s_%s", res.Data(), name.Data()), binName);
    //fitMass2D( hzm2d_4l,  "massFit_dcb2", Form("4l_%s_%s", res.Data(), name.Data()), binName);
    
}

void ebeClosureOnMC(int energy=8, int cha=0, int ijob = 0, int maxev=-1) {
    prefix = ""; 
    initCanvas(); 
    gStyle->SetOptTitle(0);
    gSystem->Load("libHiggsAnalysisCombinedLimit.so");
    TString sRes("ggH126");
    switch (cha) {
        case 0: sRes += "_4mu";   break;
        case 1: sRes += "_4e";    break;
        case 2: sRes += "_2e2mu"; break;
    }
    gYear = (energy == 7 ? 2011 : 2012);
    TString name = (energy == 8 ? "mc53x" : "mc42x");

    initModel();
    TChain *reco = new TChain("zz4lTree/probe_tree");
    if (energy == 7) {
        reco->AddFile(Form("/u3/emanuele/data/hzz4l/HZZ4L_44X_S1_V18_S2_V10/MC/hzzTree_id1126.root"));
    } else {
      reco->AddFile(Form("/u3/emanuele/data/hzz4l/HZZ4L_53X_S1_V18_S2_V10/MC/hzzTree_id1126.root"));
    }
    switch(ijob) {
        case 0:
            w->var("pull_sigma")->setConstant(true); w->saveSnapshot("clean", w->allVars());
            fitPulls(reco, sRes, name, maxev, "noscale");
            break;
        case 1:
            fitPulls(reco, sRes, name, maxev, "bare");
            break;
        case 2:
            if (!sRes.Contains("4mu")) break;
            {
                TFile *fIn = TFile::Open(Form("mcPull_%d_ptEta_ALLFine.root",gYear));
                ebeCorrH = (TH2*) fIn->Get("ALLFine_rms");
                fitPulls(reco, sRes, name, maxev, "mcPull");
            }
            break;
        case 3:
        case 4:
            {   
                TFile *fIn = TFile::Open("finalCorrections.root");
                ebeCorrH   = (TH2*) fIn->Get("mu_"+name);
                ebeCorrH4E = (TH2*) fIn->Get("el_"+name);
            }
            if (ijob == 3) {
                fitPulls(reco, sRes, name, maxev, "GioFinal");
            } else {
                w->var("pull_sigma")->setConstant(true); w->saveSnapshot("clean", w->allVars());
                fitPulls(reco, sRes, name, maxev, "GioFinalNoScale");
            }
            break;
    }

}
