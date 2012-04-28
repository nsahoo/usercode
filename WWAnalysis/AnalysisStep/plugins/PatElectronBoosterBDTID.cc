// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <DataFormats/PatCandidates/interface/Electron.h>

//BDT ElectronID
#include "HiggsAnalysis/HiggsToWW2Leptons/interface/ElectronIDMVA.h"


#include<vector>

#include "Math/VectorUtil.h"
#include "DataFormats/Math/interface/deltaPhi.h"


//
// class declaration
//

class PatElectronBoosterBDTID : public edm::EDProducer {
    public:
        explicit PatElectronBoosterBDTID(const edm::ParameterSet&);
        ~PatElectronBoosterBDTID();

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;

        // ----------member data ---------------------------
        edm::InputTag electronTag_;
        ElectronIDMVA* eleMVA;
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
PatElectronBoosterBDTID::PatElectronBoosterBDTID(const edm::ParameterSet& iConfig) :
        electronTag_(iConfig.getParameter<edm::InputTag>("electronTag"))
{
  produces<pat::ElectronCollection>();  

  // ---- Here I initialize the BDT for the mva ElectronID
  eleMVA  = new ElectronIDMVA();
  const char *base=getenv("CMSSW_BASE");
  std::string baseFolder(base);
  baseFolder += "/src/HiggsAnalysis/HiggsToWW2Leptons/data/ElectronMVAWeights/"; 
  eleMVA->Initialize("BDTG method",
		     (baseFolder+"Subdet0LowPt_WithIPInfo_BDTG.weights.xml").c_str(),
		     (baseFolder+"Subdet1LowPt_WithIPInfo_BDTG.weights.xml").c_str(),
		     (baseFolder+"Subdet2LowPt_WithIPInfo_BDTG.weights.xml").c_str(),
		     (baseFolder+"Subdet0HighPt_WithIPInfo_BDTG.weights.xml").c_str(),
		     (baseFolder+"Subdet1HighPt_WithIPInfo_BDTG.weights.xml").c_str(),
		     (baseFolder+"Subdet2HighPt_WithIPInfo_BDTG.weights.xml").c_str(),                
		     ElectronIDMVA::kWithIPInfo);
  // ---------
}


PatElectronBoosterBDTID::~PatElectronBoosterBDTID() {
  delete eleMVA;
}


void PatElectronBoosterBDTID::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<reco::Candidate> > electrons;
    iEvent.getByLabel(electronTag_,electrons);

    std::auto_ptr<pat::ElectronCollection> pOut(new pat::ElectronCollection);



    // ----- here is the real loop over the electrons ----
    for(edm::View<reco::Candidate>::const_iterator ele=electrons->begin(); ele!=electrons->end(); ++ele){    
      const pat::ElectronRef elecsRef = edm::RefToBase<reco::Candidate>(electrons,ele-electrons->begin()).castTo<pat::ElectronRef>();
      pat::Electron clone = *edm::RefToBase<reco::Candidate>(electrons,ele-electrons->begin()).castTo<pat::ElectronRef>();
      
      
      const reco::CandidateBaseRef elecsRef2(electrons,ele-electrons->begin());

      // ------ HERE I ADD THE BDT ELE ID VALUE TO THE ELECTRONS
      double EleOneOverEMinusOneOverP = (1.0/(clone.superCluster()->energy())) - 1.0 / clone.gsfTrack()->p(); 
      double xieSign  = ( (-clone.userFloat("dxyPV")) >=0 )  ? 1: -1;
      double mvaValue = eleMVA->MVAValue(clone.pt() , 
					 clone.superCluster()->eta(), 
					 clone.sigmaIetaIeta(), 
					 clone.deltaEtaSuperClusterTrackAtVtx(),
					 clone.deltaPhiSuperClusterTrackAtVtx(),
					 clone.hcalOverEcal(),
					 -clone.userFloat("dxyPV"),
					 clone.fbrem(), 
					 clone.eSuperClusterOverP(), 
					 clone.eSeedClusterOverPout(),
					 clone.userFloat("sigmaIphiIphi"),
					 clone.basicClustersSize() - 1,
					 EleOneOverEMinusOneOverP,
					 clone.eSeedClusterOverP(),
					 xieSign*clone.userFloat("ip"),
					 xieSign*clone.userFloat("ip")/clone.userFloat("ipErr"));	
      clone.addUserFloat(std::string("bdt"),mvaValue);
      // -----------------------------
      

      pOut->push_back(clone);
	
    }
    iEvent.put(pOut);
}

// ------------ method called once each job just before starting event loop  ------------
void PatElectronBoosterBDTID::beginJob() { }

// ------------ method called once each job just after ending the event loop  ------------
void PatElectronBoosterBDTID::endJob() { }


//define this as a plug-in
DEFINE_FWK_MODULE(PatElectronBoosterBDTID);
