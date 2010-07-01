// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include <DataFormats/JetReco/interface/PFJet.h>
#include <DataFormats/JetReco/interface/PFJetCollection.h>

#include "DataFormats/Common/interface/Ptr.h"

// from ROOT
#include "Math/VectorUtil.h"

//
// class declaration
//

class JetVetoFilter : public edm::EDFilter {
   public:
      explicit JetVetoFilter(const edm::ParameterSet&);
      ~JetVetoFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  edm::InputTag inputJets_;
  edm::InputTag inputLeptons_;
  double etaMax_,etMax_,dRcut_;

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
JetVetoFilter::JetVetoFilter(const edm::ParameterSet& iConfig):
  inputJets_(iConfig.getParameter<edm::InputTag>("srcJets")),
  inputLeptons_(iConfig.getParameter<edm::InputTag>("srcLeptons")),
  etaMax_(iConfig.getParameter<double>("etaMax")),
  etMax_(iConfig.getParameter<double>("etMax")),
  dRcut_(iConfig.getParameter<double>("deltaR"))
{
   //now do what ever initialization is needed

}


JetVetoFilter::~JetVetoFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
JetVetoFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace std;
   using namespace edm;
   using namespace reco;
   
   Handle<PFJetCollection> jets;
   iEvent.getByLabel(inputJets_,jets);

   Handle< vector< Ptr<Candidate> > > cands;
   iEvent.getByLabel(inputLeptons_,cands);


   bool foundJet(false);
   for(PFJetCollection::const_iterator it=jets->begin(); it!=jets->end(); ++it){
     if( abs(it->eta()) >= etaMax_ ) continue;
     if( it->et() <= etMax_ ) continue;
     // if here, found a jet passing jet-selection

     Candidate::LorentzVector jetP4 = it->p4();

     bool thisJetIsLepton(false);
     for (vector< Ptr<Candidate> >::const_iterator itLep = cands->begin(), ed = cands->end();
	  itLep != ed; ++itLep){
       double dR = abs(ROOT::Math::VectorUtil::DeltaR(jetP4,(*itLep)->p4()) );
       if(dR < dRcut_){
	 thisJetIsLepton = true;
	 break;
       }
     }
     
     if(!thisJetIsLepton) foundJet = true;
   }

   return !foundJet ;
}

// ------------ method called once each job just before starting event loop  ------------
void 
JetVetoFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JetVetoFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetVetoFilter);
