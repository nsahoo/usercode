#ifndef OctoberXTracking_DataFormats_QuarkoniaHypothesis_h
#define OctoberXTracking_DataFormats_QuarkoniaHypothesis_h

#include <vector>
#include <DataFormats/Candidate/interface/LeafCandidate.h>
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

typedef reco::Candidate::LorentzVector LorentzVector;

// a simple class
class QuarkoniaHypothesis
{
 public:
  QuarkoniaHypothesis() { }
  LorentzVector recoJPsi() const {return recoJPsi_;}
  reco::Muon recoLeg(int i=1) const;
  reco::Track recoLegTk(int i=1) const;
  LorentzVector simJPsi() const { return simJPsi_;}
  LorentzVector simLeg(int i=1) const;  
  reco::Vertex recoVertex() const {return vertex_;}

  void setRecoJPsi(const LorentzVector& lv) {recoJPsi_=lv;}
  void setRecoLeg(const reco::Muon* mu,int i=1);
  void setRecoLegTk(const reco::Track* tk,int i=1);
  void setSimJPsi(const LorentzVector& lv) {simJPsi_=lv;}
  void setSimLeg(const LorentzVector& lv,int i=1);
  void setRecoVertex(const TransientVertex& vertex){
    vertex_= vertex;}

 private:
  LorentzVector recoJPsi_;
  reco::Muon lowerPtLeg_,higherPtLeg_;
  reco::Track lowerPtLegTk_,higherPtLegTk_;
  reco::Vertex vertex_;

  LorentzVector simJPsi_;
  LorentzVector lowerPtSimLeg_,higherPtSimLeg_;

  //int value_;
};

// this is our new product, it is simply a 
// collection of SampleProd held in an std::vector
typedef std::vector<QuarkoniaHypothesis> Hypotheses;

#endif
