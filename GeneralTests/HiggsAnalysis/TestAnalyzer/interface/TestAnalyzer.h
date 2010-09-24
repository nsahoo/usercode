#ifndef TestAnalyzer_h
#define TestAnalyzer_h

// -*- C++ -*-
//
// Package:    TestAnalyzer
// Class:      TestAnalyzer
// 
/**\class TestAnalyzer TestAnalyzer.cc HiggsAnalysis/TestAnalyzer/src/TestAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Boris Mangano
//         Created:  Wed Oct 24 19:53:34 CDT 2007
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

class TestAnalyzer : public edm::EDProducer {
   public:
      explicit TestAnalyzer(const edm::ParameterSet&);
      ~TestAnalyzer();

   private:
      virtual void beginJob(const edm::EventSetup&) ;
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


#endif
