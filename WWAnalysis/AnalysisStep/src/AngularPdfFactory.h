class RooXZsZs_5D;
class RooRealVar; 

class AngularPdfFactory {

public:

  RooRealVar* mZ;     
  RooRealVar* mX;     
  RooRealVar* gamZ;   
    
  RooRealVar* a1Val;  
  RooRealVar* phi1Val;
  RooRealVar* a2Val;  
  RooRealVar* phi2Val;
  RooRealVar* a3Val;  
  RooRealVar* phi3Val;
    
  RooRealVar* R1Val;  
  RooRealVar* R2Val;  
  
  RooXZsZs_5D *PDF;

  int modelIndex;  //0 - SM Higgs, 1 - PS Higgs, 2 - Fully Longitudinal Scalar, -1 - Custom

  AngularPdfFactory(RooRealVar* m1,RooRealVar* m2,RooRealVar* h1,RooRealVar* h2,RooRealVar* Phi,RooRealVar* mZZ) ;
  ~AngularPdfFactory() ;
  void makeSMHiggs() ;
  void makePSHiggs() ;
  void makeLGHiggs() ;          
  void makeCustom(double a1, double a2, double a3,  double phi1, double phi2, double phi3) ;
  void makeParamsConst(bool yesNo=true);
  double getVal(double mZZ);
  double getValIntegrOutAngles(RooRealVar* m1,RooRealVar* m2,RooRealVar* h1,RooRealVar* h2,RooRealVar* Phi,RooRealVar* mZZ);
private:
  double Norm[80];
  void initNorm() ;
};
