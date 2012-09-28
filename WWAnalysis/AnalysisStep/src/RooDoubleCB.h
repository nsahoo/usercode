#ifndef RooDoubleCB_H
#define RooDoubleCB_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"

using namespace RooFit;

class RooDoubleCB : public RooAbsPdf {
public:
  RooDoubleCB() {} ;
  RooDoubleCB(const char *name, const char *title, RooAbsReal& _m,
        RooAbsReal& _m0, 
        RooAbsReal& _sigma,
        RooAbsReal& _alphaL, RooAbsReal& _alphaR,
        RooAbsReal& _nL, RooAbsReal& _nR) ;
   RooDoubleCB(const char *name, const char *title, RooAbsReal& _m,
        RooAbsReal& _m0,
        RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
        RooAbsReal& _alphaL, RooAbsReal& _alphaR,
        RooAbsReal& _nL, RooAbsReal& _nR) ;
 
  RooDoubleCB(const RooDoubleCB& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { 
    return new RooDoubleCB(*this,newname); }

  inline virtual ~RooDoubleCB() { }

protected:

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sigmaL;
  RooRealProxy sigmaR;
  RooRealProxy alphaL;
  RooRealProxy alphaR;
  RooRealProxy nL;
  RooRealProxy nR;

  Double_t evaluate() const;

private:
  
  ClassDef(RooDoubleCB,1)
};



#endif
