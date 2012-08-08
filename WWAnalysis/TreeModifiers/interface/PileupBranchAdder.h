#ifndef TreeModifiers_PileupBranchAdder_h
#define TreeModifiers_PileupBranchAdder_h

#include "TH1F.h"

#include "WWAnalysis/TreeModifiers/interface/BranchAdder.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/** This class implements the calculateValue() method for 
 *  adding cross-section*BR*lumi weight.
 *  Lumi is 1/fb. The cross-section and BR values are hardcoded for the moment.
 *  TO DO: extend the class so that the XS and BR values are read from 
 *  a configuration file.
 */


class PileupBranchAdder : public BranchAdder {
 public:
 PileupBranchAdder(const edm::ParameterSet& pset) ;
  
  virtual ~PileupBranchAdder(){
    delete histoWeights_;
  };

  /// add calculate the value to be added to the new branch
  virtual float calculateValue(TTree& tree,int entry, int id, double& value) ;

  
  virtual void initialize(TTree& tree);


 private:
  //add datamembers and private methods here
  void setWeights2011();
  void setWeights2012();

  float weight(float nTrueInt) const ;
  
  TH1F* histoWeights_;
  float nTrueInt_;
};

#endif
