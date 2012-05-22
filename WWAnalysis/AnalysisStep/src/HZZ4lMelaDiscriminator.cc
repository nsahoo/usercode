#include "WWAnalysis/AnalysisStep/interface/HZZ4lMelaDiscriminator.h"
#include <memory>
#include <RooRealVar.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include "AngularPdfFactory.h"
#include "RooXZsZs_5D.h"
#include "RooqqZZ_JHU.h"

HZZ4LMelaDiscriminator::HZZ4LMelaDiscriminator(Model model, const char *qqZZTemplateFileName) :
    model_(model),
    mzz_(new RooRealVar("mzz","mZZ",80,1000)),
    z1mass_(new RooRealVar("z1mass","m_{Z1}",0,180)),
    z2mass_(new RooRealVar("z2mass","m_{Z2}",0,120)), 
    cosThetaStar_(new RooRealVar("costhetastar","cos#theta^{*}",-1,1)), 
    cosTheta1_(new RooRealVar("costheta1","cos#theta_{1}",-1,1)),  
    cosTheta2_(new RooRealVar("costheta2","cos#theta_{2}",-1,1)),
    phi_(new RooRealVar("phi","#Phi",-3.1415,3.1415)),
    phi1_(new RooRealVar("phi1","#Phi^{*}_{1}",-3.1415,3.1415)),
    APF_(0), QQZZH_(0), QQZZP_(0)
{
    if (model_ != QQZZ) {
        APF_ = new AngularPdfFactory(z1mass_, z2mass_, cosTheta1_, cosTheta2_, phi_, mzz_);
        if (model_ == SMHiggs) APF_->makeSMHiggs();
        if (model_ == PSHiggs) APF_->makePSHiggs();
        APF_->makeParamsConst(true);
    } else {
        TFile *templates = TFile::Open(qqZZTemplateFileName);
        QQZZH_ = new QQZZ8DHistos(templates);
        QQZZP_ = new RooqqZZ_JHU("SMZZ","SMZZ",*z1mass_,*z2mass_,*cosTheta1_,*cosTheta2_,*phi_,*cosThetaStar_,*phi1_,*mzz_);
    }

}

HZZ4LMelaDiscriminator::~HZZ4LMelaDiscriminator()
{
    delete APF_;
    delete QQZZP_;
    delete QQZZH_;
    delete mzz_; delete z1mass_; delete z2mass_;
    delete cosThetaStar_; delete cosTheta1_; delete cosTheta2_;
    delete phi_; delete phi1_;
}

HZZ4LMelaDiscriminator::QQZZ8DHistos::QQZZ8DHistos(TDirectory *dir) :
  h_mzz((TH1F*)(dir->Get("h_mzz")->Clone())),
  h_mzzm1m2((TH3F*)(dir->Get("h_mzzm1m2")->Clone())),
  h_mzzcosthetastar((TH2F*)(dir->Get("h_mzzcosthetastar")->Clone())),
  h_mzzcostheta1((TH2F*)(dir->Get("h_mzzcostheta1")->Clone())),
  h_mzzcostheta2((TH2F*)(dir->Get("h_mzzcostheta2")->Clone())),
  h_mzzphi1((TH2F*)(dir->Get("h_mzzphi1")->Clone())),
  h_mzzphi((TH2F*)(dir->Get("h_mzzphi")->Clone()))
{
    h_mzz->SetDirectory(0);
    h_mzzm1m2->SetDirectory(0);
    h_mzzcosthetastar->SetDirectory(0);
    h_mzzcostheta1->SetDirectory(0);
    h_mzzcostheta2->SetDirectory(0);
    h_mzzphi1->SetDirectory(0);
    h_mzzphi->SetDirectory(0);
}

HZZ4LMelaDiscriminator::QQZZ8DHistos::~QQZZ8DHistos() 
{
    delete h_mzz;
    delete h_mzzm1m2;
    delete h_mzzcosthetastar;
    delete h_mzzcostheta1;
    delete h_mzzcostheta2;
    delete h_mzzphi1;
    delete h_mzzphi;
}

double HZZ4LMelaDiscriminator::QQZZ8DHistos::get(bool normalized,double mZZ, double m1, double m2, double costhetastar, double costheta1, double costheta2, double phi, double phi1) const
{
  //multiply the P values
  double n = h_mzz->GetBinContent(h_mzz->FindBin(mZZ));
  double Pmzzm1m2 = h_mzzm1m2->GetBinContent(h_mzzm1m2->FindBin(mZZ,m1,m2));

  // - - - - - - - - - - - - - - - whitbeck
  // if bin has no events: add 1
  // safety feature to prevent LD = 1 as a
  // result of low statistics
  if(Pmzzm1m2==0){
    Pmzzm1m2++;
  }
  // - - - - - - - - - - - - - - - 

  double Pmzzcosthetastar = h_mzzcosthetastar->GetBinContent(h_mzzcosthetastar->FindBin(mZZ,costhetastar));
  double Pmzzcostheta2 = h_mzzcostheta2->GetBinContent(h_mzzcostheta2->FindBin(mZZ,costheta2));
  double Pmzzcostheta1 = h_mzzcostheta1->GetBinContent(h_mzzcostheta1->FindBin(mZZ,costheta1));
  double Pmzzphi1 = h_mzzphi1->GetBinContent(h_mzzphi1->FindBin(mZZ,phi1));
  double Pmzzphi = h_mzzphi->GetBinContent(h_mzzphi->FindBin(mZZ,phi));


  //normalization
  if (normalized) {
      double binwidth_mzzm1m2 = h_mzzm1m2->GetYaxis()->GetBinWidth(1) * h_mzzm1m2->GetZaxis()->GetBinWidth(1);
      double binwidth_mzzcosthetastar = h_mzzcosthetastar->GetYaxis()->GetBinWidth(1);
      double binwidth_mzzcostheta1 = h_mzzcostheta1->GetYaxis()->GetBinWidth(1);
      double binwidth_mzzcostheta2 = h_mzzcostheta1->GetYaxis()->GetBinWidth(1);
      double binwidth_mzzphi1 = h_mzzphi1->GetYaxis()->GetBinWidth(1);
      double binwidth_mzzphi = h_mzzphi->GetYaxis()->GetBinWidth(1);

      double Pmzzm1m2_norm = Pmzzm1m2/(n*binwidth_mzzm1m2); 
      double Pmzzcosthetastar_norm = Pmzzcosthetastar/(n*binwidth_mzzcosthetastar);
      double Pmzzcostheta1_norm = Pmzzcostheta1/(n*binwidth_mzzcostheta1);
      double Pmzzcostheta2_norm = Pmzzcostheta2/(n*binwidth_mzzcostheta2);
      double Pmzzphi1_norm = Pmzzphi1/(n*binwidth_mzzphi1);
      double Pmzzphi_norm = Pmzzphi/(n*binwidth_mzzphi);

      return Pmzzm1m2_norm*Pmzzcosthetastar_norm*Pmzzcostheta1_norm*Pmzzcostheta2_norm*Pmzzphi1_norm*Pmzzphi_norm;
  } else {
      return Pmzzm1m2*Pmzzcosthetastar*Pmzzcostheta1*Pmzzcostheta2*Pmzzphi1*Pmzzphi;
  }

}

double HZZ4LMelaDiscriminator::get(double mZZ, double m1, double m2, double costhetastar, double costheta1, double costheta2, double phi, double phi1) const
{
  checkZorder(m1,m2,costhetastar,costheta1,costheta2,phi,phi1);

  mzz_->setVal(mZZ);
  z1mass_->setVal(m1);  
  z2mass_->setVal(m2);
  cosTheta1_->setVal(costheta1);
  cosTheta2_->setVal(costheta2);
  cosThetaStar_->setVal(costhetastar);
  phi_->setVal(phi);
  phi1_->setVal(phi1);
  if (mZZ <= 80) {
    return 0;
  } else if (mZZ<180) {
    if (model_ == QQZZ) {
        return 5.0 * QQZZH_->get(true, mZZ,  m1,  m2,  costhetastar,  costheta1,  costheta2,  phi,  phi1);
    } else {
        return APF_->getVal(mZZ);
    }
  } else {
      if ( mZZ <= 2*91.188) {
          z1mass_->setVal(mZZ/2. - 1e-9);
          z2mass_->setVal(mZZ/2. - 1e-9);
      } else {
          z1mass_->setVal(91.188);
          z2mass_->setVal(91.188);
      }
      if (model_ == QQZZ) {
        std::auto_ptr<RooAbsReal> integral(QQZZP_->createIntegral(RooArgSet(*cosThetaStar_,*cosTheta1_,*cosTheta2_,*phi_,*phi1_)));
        return QQZZP_->getVal() / integral->getVal() * 10.0;
      } else {
        std::auto_ptr<RooAbsReal> integral(APF_->PDF->createIntegral(RooArgSet(*cosTheta1_,*cosTheta2_,*phi_)));
        return APF_->PDF->getVal()/integral->getVal();
      }
  }
}

void HZZ4LMelaDiscriminator::checkZorder(double& z1mass, double& z2mass,
                 double& costhetastar, double& costheta1,
                 double& costheta2, double& phi,
                 double& phistar1) const {

  double tempZ1mass=z1mass;
  double tempZ2mass=z2mass;
  double tempH1=costheta1;
  double tempH2=costheta2;
  double tempHs=costhetastar;
  double tempPhi1=phistar1;
  double tempPhi=phi;

  if(z2mass>z1mass){

    z1mass=tempZ2mass;
    z2mass=tempZ1mass;
    costhetastar=-tempHs;
    costheta1=tempH2;
    costheta2=tempH1;
    phi=tempPhi;
    phistar1=-tempPhi1-tempPhi;
    if(phistar1>3.1415)
      phistar1=phistar1-2*3.1415;
    if(phistar1<-3.1415)
      phistar1=phistar1+2*3.1415;

  }else
    return;

}

