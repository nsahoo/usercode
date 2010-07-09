#include "HiggsAnalysis/TestAnalyzer/interface/TestAnalyzer.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/TestProduct/interface/McInfo.h"

#include <iostream>

using namespace std;
using namespace math;
using namespace edm;


//
// constructors and destructor
//
TestAnalyzer::TestAnalyzer(const edm::ParameterSet& iConfig)
{
   //register your products
  //Examples
  //produces<McInfo>();
  produces<McInfo>().setBranchAlias("SampleCollection");

  /*
  //if do put with a label
  produces<ExampleData2>("label");
  */
  //now do what ever other initialization is needed
  
}


TestAnalyzer::~TestAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TestAnalyzer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
*/

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::auto_ptr<McInfo> pOut(new McInfo);


/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/

   // get HepMCE form Event
   edm::Handle<HepMCProduct> hepMCEvent;
   iEvent.getByLabel("source", hepMCEvent);

   const HepMC::GenEvent * generatedEvent = hepMCEvent->GetEvent();
   
   for(HepMC::GenEvent::particle_const_iterator myParticle = generatedEvent->particles_begin(); 
       myParticle != generatedEvent->particles_end(); ++myParticle)
     {
       /*
       if( ( (abs( (*MCParticle)->pdg_id() ) == 11) ||
	     (abs( (*MCParticle)->pdg_id() ) == 13) ||
	     (abs( (*MCParticle)->pdg_id() ) == 15)    )	  
	   && (*MCParticle)->status() == 1)
	 {
	   
	   HepMC::GenParticle* mc 	= (*MCParticle);
	   HepMC::GenParticle* mcMom 	= GenParticleMother(*MCParticle);
	   HepMC::GenParticle* mcGMom 	= GenParticleGrandMother(*MCParticle);
	   
	 }
       */

       if(   (*myParticle)->pdg_id()  == 25)
	 {
	   HepMC::GenParticle* mc 	= *myParticle;
	   cout << "=== found higgs type,status, pt: " 
		<< mc->pdg_id() << " , " 
		<< mc->status() << " , " 
		<< mc->momentum().perp() << endl;
	   XYZTLorentzVector vect(mc->momentum().px(),
				  mc->momentum().py(),
				  mc->momentum().pz(),
				  mc->momentum().e());
	   
	   pOut->higgsMass_ = vect.mass();
	   pOut->higgsP_ = vect;
	   
	   cout << "=== vector pt: " << vect.pt() << endl; 
	 }
     }

   iEvent.put(pOut);
}



// ------------ method called once each job just before starting event loop  ------------
void 
TestAnalyzer::beginJob(const edm::EventSetup&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TestAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(TestAnalyzer);
