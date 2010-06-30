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

#include <DataFormats/METReco/interface/MET.h>
#include <DataFormats/METReco/interface/METFwd.h>


//
// class declaration
//

class METFilter : public edm::EDFilter {
   public:
      explicit METFilter(const edm::ParameterSet&);
      ~METFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
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
METFilter::METFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

}


METFilter::~METFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
METFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   
   Handle<METCollection> met;
   iEvent.getByLabel("tcMET",met);

   Handle< View<Candidate> > cands;
   iEvent.getByLabel("",cands);

   for (View<Candidate>::const_iterator it = cands->begin(), ed = cands->end();
	it != ed; ++it){
     

   }


   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
METFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
METFilter::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(METFilter);
