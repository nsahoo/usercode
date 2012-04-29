// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Electron.h>
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//BDT ElectronID
//#include "WWAnalysis/AnalysisStep/interface/ElectronIDMVAHZZ.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"


#include <vector>
#include <string>

#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaPhi.h"


//
// class declaration
//

class PatElectronBoosterBDTIso : public edm::EDProducer {
    public:
        explicit PatElectronBoosterBDTIso(const edm::ParameterSet&);
        ~PatElectronBoosterBDTIso();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------
        edm::InputTag electronTag_;
        EGammaMvaEleEstimator* eleMVANonTrig;
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
PatElectronBoosterBDTIso::PatElectronBoosterBDTIso(const edm::ParameterSet& iConfig) :
        electronTag_(iConfig.getParameter<edm::InputTag>("src"))
{
  produces<pat::ElectronCollection>();  

  // ---- Here I initialize the BDT for the mva ElectronIso

  const char *base=getenv("CMSSW_BASE");
  std::string baseFolder(base);
  baseFolder += "/src/WWAnalysis/AnalysisStep/data/ElectronMVAWeights/";

  manualCatNonTrigWeigths.push_back(baseFolder+"/ElectronIso_BDTG_V0_BarrelPt5To10.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/ElectronIso_BDTG_V0_EndcapPt5To10.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/ElectronIso_BDTG_V0_BarrelPt10ToInf.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/ElectronIso_BDTG_V0_EndcapPt10ToInf.weights.xml");  

  eleMVANonTrig  = new EGammaMvaEleEstimator();
  eleMVANonTrig->initialize("EleIso_BDTG_IsoRings",EGammaMvaEleEstimator::kIsoRings,true,manualCatNonTrigWeigths);
  eleMVANonTrig->SetPrintMVADebug(kFALSE);

  // ---------
}


PatElectronBoosterBDTIso::~PatElectronBoosterBDTIso() {
  delete eleMVANonTrig;
}


void PatElectronBoosterBDTIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<reco::Candidate> > electrons;
    iEvent.getByLabel(electronTag_,electrons);

    std::auto_ptr<pat::ElectronCollection> pOut(new pat::ElectronCollection);

    // ----- here is the real loop over the electrons ----
    for(edm::View<reco::Candidate>::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){    
      const pat::ElectronRef elecsRef = edm::RefToBase<reco::Candidate>(electrons,ele-electrons->begin()).castTo<pat::ElectronRef>();
      pat::Electron clone = *edm::RefToBase<reco::Candidate>(electrons,ele-electrons->begin()).castTo<pat::ElectronRef>();
      
      
      // ------ HERE I ADD THE BDT ELE ID VALUE TO THE ELECTRONS
      double mvaValueNonTrig = eleMVANonTrig->isoMvaValue(clone.pt(),
                                                          clone.superCluster()->eta(),
                                                          clone.userFloat("rhoEl"),
                                                          ElectronEffectiveArea::kEleEAData2011,
                                                          clone.userFloat("electronPFIsoChHad01"),
                                                          clone.userFloat("electronPFIsoChHad02") - clone.userFloat("electronPFIsoChHad01"),
                                                          clone.userFloat("electronPFIsoChHad03") - clone.userFloat("electronPFIsoChHad02"),
                                                          clone.userFloat("electronPFIsoChHad04") - clone.userFloat("electronPFIsoChHad03"),
                                                          clone.userFloat("electronPFIsoChHad05") - clone.userFloat("electronPFIsoChHad04"),
                                                          clone.userFloat("electronPFIsoPhoton01"),
                                                          clone.userFloat("electronPFIsoPhoton02") - clone.userFloat("electronPFIsoPhoton01"),
                                                          clone.userFloat("electronPFIsoPhoton03") - clone.userFloat("electronPFIsoPhoton02"),
                                                          clone.userFloat("electronPFIsoPhoton04") - clone.userFloat("electronPFIsoPhoton03"),
                                                          clone.userFloat("electronPFIsoPhoton05") - clone.userFloat("electronPFIsoPhoton04"),
                                                          clone.userFloat("electronPFIsoNHad01"),
                                                          clone.userFloat("electronPFIsoNHad02") - clone.userFloat("electronPFIsoNHad01"),
                                                          clone.userFloat("electronPFIsoNHad03") - clone.userFloat("electronPFIsoNHad02"),
                                                          clone.userFloat("electronPFIsoNHad04") - clone.userFloat("electronPFIsoNHad03"),
                                                          clone.userFloat("electronPFIsoNHad05") - clone.userFloat("electronPFIsoNHad04"),
                                                          false);

      clone.addUserFloat(std::string("bdtisonontrig"),mvaValueNonTrig);

      // -----------------------------
      

      pOut->push_back(clone);
	
    }
    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void PatElectronBoosterBDTIso::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void PatElectronBoosterBDTIso::endJob() { }


//define this as a plug-in
DEFINE_FWK_MODULE(PatElectronBoosterBDTIso);
