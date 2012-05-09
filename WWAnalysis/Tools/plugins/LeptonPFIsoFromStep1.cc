#include <memory>
#include "Math/VectorUtil.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "TMath.h"

class LeptonPFIsoFromStep1 : public edm::EDProducer {
public:
  explicit LeptonPFIsoFromStep1(const edm::ParameterSet&);
  ~LeptonPFIsoFromStep1();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  edm::InputTag leptonLabel_;
  edm::InputTag pfLabel_;
  StringCutObjectSelector<reco::Candidate> pfSelection_;
  double deltaR_, deltaRself_;
  double directional_;
  double vetoConeEndcaps_;

};



LeptonPFIsoFromStep1::LeptonPFIsoFromStep1(const edm::ParameterSet& iConfig) :
  leptonLabel_(iConfig.getParameter<edm::InputTag>("leptonLabel")),
  pfLabel_(iConfig.getParameter<edm::InputTag>("pfLabel")),
  pfSelection_(iConfig.getParameter<std::string>("pfSelection"), true),
  deltaR_(iConfig.getParameter<double>("deltaR")),
  deltaRself_(iConfig.getParameter<double>("deltaRself")),
  directional_(iConfig.getParameter<bool>("directional")),
  vetoConeEndcaps_(iConfig.getParameter<double>("vetoConeEndcaps"))
{
  produces<edm::ValueMap<float> >().setBranchAlias("pfMuIso");
}

void LeptonPFIsoFromStep1::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<reco::CandidateView> muH;
  iEvent.getByLabel(leptonLabel_,muH);

  edm::Handle<reco::CandidateView> pfH;
  iEvent.getByLabel(pfLabel_,pfH);

  std::vector<float> isoV;
  std::auto_ptr<edm::ValueMap<float> > isoM(new edm::ValueMap<float> ());
  edm::ValueMap<float>::Filler isoF(*isoM);

  for(size_t i=0, n = muH->size(); i < n; ++i) {
    const reco::Candidate &mu = muH->at(i);

//     Double_t zLepton = 0.0;
//     if(mu.track().isNonnull()) zLepton = mu.track()->dz(vtxH->at(0).position());

    Double_t ptSum =0.;  
    math::XYZVector isoAngleSum;
    std::vector<math::XYZVector> coneParticles;

    for (size_t j=0; j<pfH->size();j++) {   
      const reco::Candidate &pf = pfH->at(j);
        
      if (!pfSelection_(pf)) continue; 

      double dr = deltaR(pf, mu) ;

      if (pf.charge() != 0 && deltaR(pf, mu) < deltaRself_) continue;

      // dR Veto for Gamma: no-one in EB, dR > 0.08 in EE
      if (fabs(mu.eta()) >1.479 && dr < vetoConeEndcaps_) continue;

     // add the pf pt if it is inside the extRadius 
      if ( dr < deltaR_ ) {

        // scalar sum
        ptSum += pf.pt();

        // directional sum
        math::XYZVector transverse( pf.eta() - mu.eta()
                                    , reco::deltaPhi(pf.phi(), mu.phi())
                                    , 0);
        transverse *= pf.pt() / transverse.rho();
        if (transverse.rho() > 0) {
          isoAngleSum += transverse;
          coneParticles.push_back(transverse);
        }
      }

    }
    if (directional_) {
      double directionalPT = 0;
      for (unsigned int iPtcl = 0; iPtcl < coneParticles.size(); ++iPtcl)
        directionalPT += pow(TMath::ACos( coneParticles[iPtcl].Dot(isoAngleSum) / coneParticles[iPtcl].rho() / isoAngleSum.rho() ),2) * coneParticles[iPtcl].rho();
      isoV.push_back(directionalPT);
    } else isoV.push_back(ptSum);
  }

  isoF.insert(muH,isoV.begin(),isoV.end());

  isoF.fill();
  iEvent.put(isoM);

}

LeptonPFIsoFromStep1::~LeptonPFIsoFromStep1() { }
void LeptonPFIsoFromStep1::beginJob() { }
void LeptonPFIsoFromStep1::endJob() { }
DEFINE_FWK_MODULE(LeptonPFIsoFromStep1);
