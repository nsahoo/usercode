// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include <DataFormats/PatCandidates/interface/Muon.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//BDT MuonID
#include "Muon/MuonAnalysisTools/interface/MuonMVAEstimator.h"


#include <vector>
#include <string>

#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaPhi.h"


//
// class declaration
//

class PatMuonBoosterBDTID : public edm::EDProducer {
    public:
        explicit PatMuonBoosterBDTID(const edm::ParameterSet&);
        ~PatMuonBoosterBDTID();

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
PatMuonBoosterBDTID::PatMuonBoosterBDTID(const edm::ParameterSet& iConfig) :
  muonTag_(iConfig.getParameter<edm::InputTag>("src"))
{
  produces<pat::MuonCollection>();  

  // ---- Here I initialize the BDT for the mva MuonID

  const char *base=getenv("CMSSW_BASE");
  std::string baseFolder(base);
  baseFolder += "/src/Muon/MuonAnalysisTools/data/";

  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIDMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIDMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIDMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIDMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIDMVA_sixie-Tracker_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIDMVA_sixie-Global_V0_BDTG.weights.xml");

  muMVANonTrig  = new MuonMVAEstimator();
  muMVANonTrig->initialize("MuonID_BDTG",MuonMVAEstimator::kID,true,manualCatNonTrigWeigths);
  muMVANonTrig->SetPrintMVADebug(kFALSE);

  // ---------
}


PatMuonBoosterBDTID::~PatMuonBoosterBDTID() {
  delete muMVANonTrig;
}


void PatMuonBoosterBDTID::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
      double mvaValueNonTrig;
      reco::TrackRef muTrk = musRef->track();
      if (muTrk.isNull()) {
        muTrk = musRef->standAloneMuon();
      }
      if (muTrk.isNull()) {
        //if muon is not standalone either, then return -9999
        mvaValueNonTrig = -9999;
      }
     
      double muNchi2 = 0.0; 
      if (musRef->combinedMuon().isNonnull()) { 
        muNchi2 = musRef->combinedMuon()->chi2() / (Double_t)musRef->combinedMuon()->ndof(); 
      } else if (musRef->standAloneMuon().isNonnull()) { 
        muNchi2 = musRef->standAloneMuon()->chi2() / (Double_t)musRef->standAloneMuon()->ndof(); 
      } else if (musRef->track().isNonnull()) { 
        muNchi2 = musRef->track()->chi2() / (Double_t)musRef->track()->ndof(); 
      }

      mvaValueNonTrig = muMVANonTrig->mvaValue_ID(muTrk->pt(),
                                                  muTrk->eta(),
                                                  musRef->isGlobalMuon(),
                                                  musRef->isTrackerMuon(),
                                                  muTrk->chi2() / (Double_t)muTrk->ndof(),
                                                  muNchi2,
                                                  musRef->globalTrack().isNonnull() ? musRef->globalTrack()->hitPattern().numberOfValidMuonHits() : 0,
                                                  muTrk->numberOfValidHits(),
                                                  muTrk->hitPattern().numberOfValidPixelHits(),
                                                  musRef->numberOfMatches(),
                                                  musRef->combinedQuality().trkKink,
                                                  muon::segmentCompatibility(*musRef, reco::Muon::SegmentAndTrackArbitration),
                                                  musRef->caloCompatibility(),
                                                  musRef->calEnergy().had,
                                                  musRef->calEnergy().em,
                                                  musRef->calEnergy().hadS9,
                                                  musRef->calEnergy().emS9,
                                                  false);

      clone.addUserFloat(std::string("bdtidnontrig"),mvaValueNonTrig);

      // -----------------------------
      

      pOut->push_back(clone);
	
    }
    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void PatMuonBoosterBDTID::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void PatMuonBoosterBDTID::endJob() { }


//define this as a plug-in
DEFINE_FWK_MODULE(PatMuonBoosterBDTID);
