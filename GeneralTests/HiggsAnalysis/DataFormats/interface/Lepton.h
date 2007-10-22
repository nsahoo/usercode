#ifndef HiggsAnalysis_DataFormats_Lepton_h
#define HiggsAnalysis_DataFormats_Lepton_h


class Lepton 
{
 public:
  
  Lepton() {;}
  Lepton( double px, double py, double pz, double e) :
    px_(px),py_(py),pz_(pz),e_(e) {;}

  virtual ~Lepton() {;}
    
    double px();
    double py();
    double pz();
    double  e();

 private:
    double px_,py_,pz_,e_;
    //		ClassDef(Lepton, 1);
};

#endif
