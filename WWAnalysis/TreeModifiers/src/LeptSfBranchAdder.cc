#include "WWAnalysis/TreeModifiers/interface/LeptSfBranchAdder.h"


LeptSfBranchAdder::LeptSfBranchAdder(const edm::ParameterSet& pset) : BranchAdder(pset){
  /*
  std::string pileupConf_     = pset.getParameter<std::string>("pileupConf"); 
  if(pileupConf_ == "2011") setWeights2011();
  else if(pileupConf_ == "2012") setWeights2012();
  else{
    std::cout << "WARNING: you didn't specify the \"pileupConf\" parameter. Using default: 2012" 
	      << std::endl;
    setWeights2012();    
  }
  */
}

void LeptSfBranchAdder::initialize(TTree& tree) {
  tree.SetBranchAddress("l1pt",&pt1_);tree.SetBranchAddress("l1eta",&eta1_);tree.SetBranchAddress("l1pdgId",&id1_);
  tree.SetBranchAddress("l2pt",&pt2_);tree.SetBranchAddress("l2eta",&eta2_);tree.SetBranchAddress("l2pdgId",&id2_);
  tree.SetBranchAddress("l3pt",&pt3_);tree.SetBranchAddress("l3eta",&eta3_);tree.SetBranchAddress("l3pdgId",&id3_);
  tree.SetBranchAddress("l4pt",&pt4_);tree.SetBranchAddress("l4eta",&eta4_);tree.SetBranchAddress("l4pdgId",&id4_);
}


float LeptSfBranchAdder::calculateValue(TTree& tree,int entry, int id, float& value)  {
  tree.GetEntry(entry);
  //value = getSF(pt1_,eta1_,id1_)*getSF(pt2_,eta2_,id2_)*getSF(pt3_,eta3_,id3_)*getSF(pt4_,eta4_,id4_);
  value = 1.0;

  return 0;
}


