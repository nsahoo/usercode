#ifndef JPsiAnalysis_JPsiAnalysisPAT_h
#define JPsiAnalysis_JPsiAnalysisPAT_h


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class decleration
//

class JPsiAnalysisPAT : public edm::EDProducer {
   public:
      explicit JPsiAnalysisPAT(const edm::ParameterSet&);
      ~JPsiAnalysisPAT();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
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

#endif
