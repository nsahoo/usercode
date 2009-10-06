#include "OctoberXTracking/DataFormats/interface/QuarkoniaHypothesis.h"


void QuarkoniaHypothesis::setRecoLeg(const reco::Muon* mu,int i){
  switch(i) { 
  case 1:
    lowerPtLeg_ = *mu;
    break;
    
  case 2:
    higherPtLeg_ = *mu;
    break;
    
  default:
    break;
  }
}
 void QuarkoniaHypothesis::setRecoLegTk(const reco::Track* tk,int i){
  switch(i) { 
  case 1:
    if(tk) lowerPtLegTk_ = *tk; else lowerPtLegTk_ = reco::Track();
    break;
    
  case 2:
    if(tk) higherPtLegTk_ = *tk; else higherPtLegTk_ = reco::Track();
    break;
    
  default:
    break;
  }
}


void QuarkoniaHypothesis::setSimLeg(const LorentzVector& lv,int i){
  switch(i) { 
  case 1:
    lowerPtSimLeg_ = lv;
    break;
    
  case 2:
    higherPtSimLeg_ = lv;
    break;
    
  default:
    break;
  }
}




reco::Muon QuarkoniaHypothesis::recoLeg(int i) const {
  switch(i) { 
  case 1:
    return lowerPtLeg_;
    
  case 2:
    return higherPtLeg_;
    
  default:
    return reco::Muon();
  };
}

reco::Track QuarkoniaHypothesis::recoLegTk(int i) const {
  switch(i) { 
  case 1:
    return lowerPtLegTk_;
    
  case 2:
    return higherPtLegTk_;
    
  default:
    return reco::Track();
  };
}


LorentzVector QuarkoniaHypothesis::simLeg(int i) const {
  switch(i) { 
  case 1:
    return lowerPtSimLeg_;
    
  case 2:
    return higherPtSimLeg_;
    
  default:
    return LorentzVector();
  };
}
