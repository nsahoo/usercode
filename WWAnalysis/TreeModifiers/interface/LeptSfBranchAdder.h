#ifndef TreeModifiers_LeptSfBranchAdder_h
#define TreeModifiers_LeptSfBranchAdder_h

#include "TH1F.h"

#include "WWAnalysis/TreeModifiers/interface/BranchAdder.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/** This class implements the calculateValue() method for 
 *  adding cross-section*BR*lumi weight.
 *  Lumi is 1/fb. The cross-section and BR values are hardcoded for the moment.
 *  TO DO: extend the class so that the XS and BR values are read from 
 *  a configuration file.
 */


class LeptSfBranchAdder : public BranchAdder {
 public:
 LeptSfBranchAdder(const edm::ParameterSet& pset) ;
  
  virtual ~LeptSfBranchAdder(){
    //delete histoWeights_;
  };

  /// add calculate the value to be added to the new branch
  virtual float calculateValue(TTree& tree,int entry, int id, float& value) ;

  
  virtual void initialize(TTree& tree);


 private:
  //add datamembers and private methods here
  
  //TH1F* histoWeights_;
  float pt1_,pt2_,pt3_,pt4_;
  float eta1_,eta2_,eta3_,eta4_;
  float id1_,id2_,id3_,id4_;
};

#endif
