#ifndef TreeModifiers_PileupBranchAdder_h
#define TreeModifiers_PileupBranchAdder_h

#include "TH1F.h"

#include "WWAnalysis/TreeModifiers/interface/BranchAdder.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

/** This class implements the calculateValue() method for 
 *  adding cross-section*BR*lumi weight.
 */


class PileupBranchAdder : public BranchAdder {
 public:
 PileupBranchAdder(const edm::ParameterSet& pset) ;
  
  virtual ~PileupBranchAdder(){
    delete histoWeights_;
  };

  virtual void initialize(TTree& tree);

  /// add calculate the value to be added to the new branch
  virtual float calculateValue(TTree& tree,int k, float& value) ;

 

 private:
  //add datamembers and private methods here
  void setWeights2011();
  void setWeights2012();

  float weight(float nTrueInt) const {
    //if(nTrueInt>50) return 1;
    return histoWeights_->GetBinContent(histoWeights_->FindBin(nTrueInt));
  }
  
  bool isInitialized_;
  std::string pileupConf_;
  TH1F* histoWeights_;
  //float nTrueInt;
};

#endif
