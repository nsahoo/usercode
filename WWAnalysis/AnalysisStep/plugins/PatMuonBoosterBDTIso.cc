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
        edm::InputTag vertexsTag_;
        edm::InputTag pfCandsTag_;
        edm::InputTag rhoTag_;
        double dzCut_;
        std::string outputName_;
        MuonEffectiveArea::MuonEffectiveAreaTarget effAreaTarget_;
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
  muonTag_(iConfig.getParameter<edm::InputTag>("src")),
  vertexsTag_(iConfig.getParameter<edm::InputTag>("vertexs")),
  pfCandsTag_(iConfig.getParameter<edm::InputTag>("pfCands")),
  rhoTag_(iConfig.getParameter<edm::InputTag>("rho")),
  dzCut_(iConfig.getParameter<double>("dzCut")),
  outputName_(iConfig.getParameter<std::string>("outputName"))
{
  produces<pat::MuonCollection>();  

  // ---- Here I initialize the BDT for the mva MuonIso

  const char *base=getenv("CMSSW_BASE");
  std::string baseFolder(base);
  baseFolder += "/src/Muon/MuonAnalysisTools/data/";

  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-BarrelPt5To10_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-EndcapPt5To10_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-BarrelPt10ToInf_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-EndcapPt10ToInf_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-Tracker_V0_BDTG.weights.xml");
  manualCatNonTrigWeigths.push_back(baseFolder+"/MuonIsoMVA_sixie-Global_V0_BDTG.weights.xml");

  muMVANonTrig  = new MuonMVAEstimator();
  muMVANonTrig->initialize("MuonIso_BDTG_IsoRings",MuonMVAEstimator::kIsoRings,true,manualCatNonTrigWeigths);
  muMVANonTrig->SetPrintMVADebug(kFALSE);
  //muMVANonTrig->SetPrintMVADebug(kTRUE);

  // ---------

  std::string eaTarget = iConfig.getParameter<std::string>("effectiveAreaTarget");
  if      (eaTarget == "NoCorr") effAreaTarget_ = MuonEffectiveArea::kMuEANoCorr;
  else if (eaTarget == "Data2011") effAreaTarget_ = MuonEffectiveArea::kMuEAData2011;
  else if (eaTarget == "Data2012") effAreaTarget_ = MuonEffectiveArea::kMuEAData2012;
  else if (eaTarget == "Summer11MC") effAreaTarget_ = MuonEffectiveArea::kMuEASummer11MC;
  else if (eaTarget == "Fall11MC") effAreaTarget_ = MuonEffectiveArea::kMuEAFall11MC;
  else throw cms::Exception("Configuration") << "Unknown effective area " << eaTarget << std::endl;

}


PatMuonBoosterBDTIso::~PatMuonBoosterBDTIso() {
  delete muMVANonTrig;
}


void PatMuonBoosterBDTIso::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace std;

    edm::Handle<edm::View<pat::Muon> > muons;
    iEvent.getByLabel(muonTag_,muons);

    edm::Handle<reco::VertexCollection> vertexs;
    iEvent.getByLabel(vertexsTag_,vertexs);

    Handle<reco::PFCandidateCollection> pfCands;
    iEvent.getByLabel(pfCandsTag_,pfCands);

    edm::Handle<double> hRho;
    iEvent.getByLabel(rhoTag_,hRho);


    std::auto_ptr<pat::MuonCollection> pOut(new pat::MuonCollection);


    const reco::GsfElectronCollection dummyIdentifiedEleCollection;
    const reco::MuonCollection dummyIdentifiedMuCollection;

    // ----- here is the real loop over the muons ----
    for(edm::View<pat::Muon>::const_iterator mu=muons->begin(); mu!=muons->end(); ++mu){    
      pat::Muon clone = *mu;
      
      //cout << "-- output rho: " << *hRho << endl;

      double mvaValueNonTrig = muMVANonTrig->mvaValue(clone,
						      vertexs->front(),
						      *pfCands,
						      *hRho,MuonEffectiveArea::kMuEAFall11MC,
						      dummyIdentifiedEleCollection,
						      dummyIdentifiedMuCollection,
						      dzCut_);     

      clone.addUserFloat(outputName_,mvaValueNonTrig);

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
