// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//BDT MuonIso
#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"


#include <vector>
#include <string>

#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaPhi.h"


//
// class declaration
//

class PatMuonBoosterBDTIso : public edm::EDProducer {
    public:
        explicit PatMuonBoosterBDTIso(const edm::ParameterSet&);
        ~PatMuonBoosterBDTIso();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------
        edm::InputTag muonTag_;
        MuonMVAEstimator* muMVANonTrig;
        std::vector<std::string> manualCatNonTrigWeigths;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
PatMuonBoosterBDTIso::PatMuonBoosterBDTIso(const edm::ParameterSet& iConfig) :
  muonTag_(iConfig.getParameter<edm::InputTag>("src"))
{
  produces<pat::MuonCollection>();  

  // ---- Here I initialize the BDT for the mva MuonIso

  const char *base=getenv("CMSSW_BASE");
  std::string baseFolder(base);
  baseFolder += "/src/Muon/MuonAnalysisTools/data/";

  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml");

  muMVANonTrig  = new MuonMVAEstimator();
  muMVANonTrig->initialize("MuonIso_BDTG_IsoRings",MuonMVAEstimator::kIsoRings,true,manualCatNonTrigWeigths);
  muMVANonTrig->SetPrintMVADebug(kFALSE);

  // ---------
}


PatMuonBoosterBDTIso::~PatMuonBoosterBDTIso() {
  delete muMVANonTrig;
}


void PatMuonBoosterBDTIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<reco::Candidate> > muons;
    iEvent.getByLabel(muonTag_,muons);

    std::auto_ptr<pat::MuonCollection> pOut(new pat::MuonCollection);

    // ----- here is the real loop over the muons ----
    for(edm::View<reco::Candidate>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){    
      const pat::MuonRef musRef = edm::RefToBase<reco::Candidate>(muons,mu-muons->begin()).castTo<pat::MuonRef>();
      pat::Muon clone = *edm::RefToBase<reco::Candidate>(muons,mu-muons->begin()).castTo<pat::MuonRef>();
      
      
      // ------ HERE I ADD THE BDT MU ISO VALUE TO THE MUONS
      double mvaValueNonTrig = muMVANonTrig->isoMvaValue(clone.track()->pt(),
                                                         clone.track()->eta(),
                                                         clone.isGlobalMuon(),
                                                         clone.isTrackerMuon(),
                                                         clone.userFloat("rhoMu"),
                                                         MuonEffectiveArea::kMuEAFall11MC,
                                                         clone.userFloat("muonPFIsoChHad01"),
                                                         clone.userFloat("muonPFIsoChHad02") - clone.userFloat("muonPFIsoChHad01"),
                                                         clone.userFloat("muonPFIsoChHad03") - clone.userFloat("muonPFIsoChHad02"),
                                                         clone.userFloat("muonPFIsoChHad04") - clone.userFloat("muonPFIsoChHad03"),
                                                         clone.userFloat("muonPFIsoChHad05") - clone.userFloat("muonPFIsoChHad04"),
                                                         clone.userFloat("muonPFIsoPhoton01"),
                                                         clone.userFloat("muonPFIsoPhoton02") - clone.userFloat("muonPFIsoPhoton01"),
                                                         clone.userFloat("muonPFIsoPhoton03") - clone.userFloat("muonPFIsoPhoton02"),
                                                         clone.userFloat("muonPFIsoPhoton04") - clone.userFloat("muonPFIsoPhoton03"),
                                                         clone.userFloat("muonPFIsoPhoton05") - clone.userFloat("muonPFIsoPhoton04"),
                                                         clone.userFloat("muonPFIsoNHad01"),
                                                         clone.userFloat("muonPFIsoNHad02") - clone.userFloat("muonPFIsoNHad01"),
                                                         clone.userFloat("muonPFIsoNHad03") - clone.userFloat("muonPFIsoNHad02"),
                                                         clone.userFloat("muonPFIsoNHad04") - clone.userFloat("muonPFIsoNHad03"),
                                                         clone.userFloat("muonPFIsoNHad05") - clone.userFloat("muonPFIsoNHad04"),
                                                         false);

      clone.addUserFloat(std::string("bdtisonontrig"),mvaValueNonTrig);

      // -----------------------------
      

      pOut->push_back(clone);
	
    }
    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void PatMuonBoosterBDTIso::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void PatMuonBoosterBDTIso::endJob() { }


//define this as a plug-in
DEFINE_FWK_MODULE(PatMuonBoosterBDTIso);
