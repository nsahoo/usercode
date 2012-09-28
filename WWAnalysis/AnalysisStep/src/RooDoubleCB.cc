#include "RooFit.h"

#include "Riostream.h"
#include "Riostream.h"
#include <math.h>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooMath.h"
#include "TMath.h"

#include "RooDoubleCB.h"

ClassImp(RooDoubleCB)

RooDoubleCB::RooDoubleCB(const char *name, const char *title, RooAbsReal& _m, RooAbsReal& _m0,
                                                 RooAbsReal& _sigma, RooAbsReal& _alphaL, RooAbsReal& _alphaR, RooAbsReal& _nL, RooAbsReal& _nR) :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
    m0("m0", "M0", this, _m0),
    sigmaL("sigmaL", "Sigma", this, _sigma),
  sigmaR("sigmaR", "Sigma", this, _sigma),
    alphaL("alphaL", "AlphaL", this, _alphaL),
    alphaR("alphaR", "AlphaR", this, _alphaR),
    nL("nL", "nL", this, _nL),
    nR("nR", "nR", this, _nR)
{
}

RooDoubleCB::RooDoubleCB(const char *name, const char *title, RooAbsReal& _m, RooAbsReal& _m0, RooAbsReal& _sigmaL,
                                                 RooAbsReal& _sigmaR, RooAbsReal& _alphaL, RooAbsReal& _alphaR, RooAbsReal& _nL, RooAbsReal& _nR) :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
    m0("m0", "M0", this, _m0),
    sigmaL("sigmaL", "Sigma", this, _sigmaL),
    sigmaR("sigmaR", "Sigma", this, _sigmaR),
  alphaL("alphaL", "AlphaL", this, _alphaL),
  alphaR("alphaR", "AlphaR", this, _alphaR),
  nL("nL", "nL", this, _nL),
  nR("nR", "nR", this, _nR)
{
}

RooDoubleCB::RooDoubleCB(const RooDoubleCB& other, const char* name) :
    RooAbsPdf(other, name), m("m", this, other.m), m0("m0", this, other.m0),
    sigmaL("sigmaL", this, other.sigmaL), 
  sigmaR("sigmaR", this, other.sigmaR), 
  alphaL("alphaL", this, other.alphaL),
  alphaR("alphaR", this, other.alphaR), 
  nL("nL", this, other.nL),
  nR("nR", this, other.nR)
{
}

Double_t RooDoubleCB::evaluate() const {

  Double_t arg = m - m0;

  if (arg < 0.0)
    {
      Double_t t = (m-m0)/sigmaL; //t < 0

      Double_t absAlpha = fabs((Double_t)alphaL); //slightly redundant since alpha > 0 anyway, but never mind

      if (t >= -absAlpha) { //-absAlpha <= t < 0
        return exp(-0.5*t*t);
      }
      else {
        Double_t a = TMath::Power(nL/absAlpha,nL)*exp(-0.5*absAlpha*absAlpha);
        Double_t b = nL/absAlpha - absAlpha;

        return a/TMath::Power(b - t, nL); //b - t
      }
    }
  else 
    {
      Double_t t = (m-m0)/sigmaR; //t > 0

      Double_t absAlpha = fabs((Double_t)alphaR);

      if (t <= absAlpha) { //0 <= t <= absAlpha
        return exp(-0.5*t*t);
      }
      else {
        Double_t a = TMath::Power(nR/absAlpha,nR)*exp(-0.5*absAlpha*absAlpha);
        Double_t b = nR/absAlpha - absAlpha;

        return a/TMath::Power(b + t, nR); //b + t
      }
    }

}

