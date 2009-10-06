#ifndef JPsiAnalysisMC_JPsiAnalysisMC_h
#define JPsiAnalysisMC_JPsiAnalysisMC_h


// -*- C++ -*-
//
// Package:    JPsiAnalysisMC
// Class:      JPsiAnalysisMC
// 
/**\class OctoberXTracking/JPsiAnalysis/interface/JPsiAnalysisMC.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Boris MANGANO
//         Created:  Tue Aug  4 11:01:49 CEST 2009
// $Id$
//
//


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

class JPsiAnalysisMC : public edm::EDProducer {
   public:
      explicit JPsiAnalysisMC(const edm::ParameterSet&);
      ~JPsiAnalysisMC();

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
