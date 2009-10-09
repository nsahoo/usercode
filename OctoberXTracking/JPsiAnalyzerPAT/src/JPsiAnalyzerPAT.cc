// -*- C++ -*-
//
// Package:    JPsiAnalyzerPAT
// Class:      JPsiAnalyzerPAT
// 
/**\class JPsiAnalyzerPAT JPsiAnalyzerPAT.cc OctoberXTracking/JPsiAnalyzerPAT/src/JPsiAnalyzerPAT.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  
//         Created:  Fri Oct  9 04:59:40 PDT 2009
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class decleration
//
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/PatCandidates/interface/Muon.h>

class JPsiAnalyzerPAT : public edm::EDAnalyzer {
   public:
      explicit JPsiAnalyzerPAT(const edm::ParameterSet&);
      ~JPsiAnalyzerPAT();


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
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
JPsiAnalyzerPAT::JPsiAnalyzerPAT(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

}


JPsiAnalyzerPAT::~JPsiAnalyzerPAT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
JPsiAnalyzerPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;

   Handle<pat::CompositeCandidateCollection > collH;
   iEvent.getByLabel("myJPsiAnalysisPAT",collH);

   bool passedEvent(false);
   cout << "analyze event" << endl;
   for(vector<pat::CompositeCandidate>::const_iterator it=collH->begin();
       it!=collH->end();++it){
     
     bool step1(true),step2(true),step3(true),step4(true),step5(true),step6(true);


      const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(it->daughter("muon1"));
      const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(it->daughter("muon2"));
      
      if(muon1==0 || muon2==0) {cout << "ERROR: dynamic cast failed" << endl; continue;}
      
      if(muon1->charge() == muon2->charge()) step1=false;
      if(!(muon1->isGlobalMuon() && muon2->isGlobalMuon())) step2=false;
      if(!(muon1->innerTrack()->found() >11  && muon2->innerTrack()->found() >=11 )) step3=false;
      if( muon1->innerTrack()->pt()< 2.5 || muon2->innerTrack()->pt()<2.5 ) step4=false;
      if( it->userFloat("vNChi2") > 4.) step5=false;
      if( it->mass()<2.6 || it->mass()>3.5) step6=false;

      cout << "analyze candidate" << endl;
     if(step1) {  
       cout << "passed step1" << endl;
       //recoMass1->Fill(it->mass());   recoMass1b->Fill(it->mass());
      }

      if(step1 && step2) {
	//recoMass2->Fill(it->mass());	recoMass2b->Fill(it->mass());
      }

      if(step1 && step2 && step3) {
	//recoMass3->Fill(it->mass());	recoMass3b->Fill(it->mass());
      }

      if(step1 && step2 && step3 && step4) {
	//recoMass4->Fill(it->mass());	recoMass4b->Fill(it->mass());
      }

      if(step1 && step2 && step3 && step4 && step5){
	//recoMass5->Fill(it->mass());	recoMass5b->Fill(it->mass());
	cout << "passed step5" << endl;
      }

      if(step1 && step2 && step3 && step4 && step5 && step6) {
	//recoMass6->Fill(it->mass());	recoMass6b->Fill(it->mass());
      }

   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
JPsiAnalyzerPAT::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
JPsiAnalyzerPAT::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JPsiAnalyzerPAT);
