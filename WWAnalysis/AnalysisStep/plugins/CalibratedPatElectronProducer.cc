// -*- C++ -*-
//
// Package:    EgammaElectronProducers
// Class:      CalibratedGsfElectronProducer
//
/**\class CalibratedGsfElectronProducer RecoEgamma/ElectronProducers/src/CalibratedGsfElectronProducer.cc

 Description: EDProducer of GsfElectron objects

 Implementation:
     <Notes on implementation>
*/

#include "WWAnalysis/AnalysisStep/interface/CalibratedPatElectronProducer.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "WWAnalysis/AnalysisStep/interface/PatElectronEnergyCalibrator.h"

#include <iostream>

using namespace edm ;
using namespace std ;
using namespace reco ;
using namespace pat ;

CalibratedPatElectronProducer::CalibratedPatElectronProducer( const edm::ParameterSet & cfg )
// : PatElectronBaseProducer(cfg)
 {

  produces<ElectronCollection>();

  inputPatElectrons = cfg.getParameter<edm::InputTag>("inputPatElectronsTag");
  dataset = cfg.getParameter<std::string>("inputDataset");
  isMC = cfg.getParameter<bool>("isMC");
  isAOD = cfg.getParameter<bool>("isAOD");
  updateEnergyError = cfg.getParameter<bool>("updateEnergyError");
  debug = cfg.getParameter<bool>("debug");
  
  //basic checks
  if (isMC&&(dataset!="Summer11"&&dataset!="Fall11"&&dataset!="Summer12"))
   { throw cms::Exception("CalibratedgsfElectronProducer|ConfigError")<<"Unknown MC dataset" ; }
  if (!isMC&&(dataset!="Prompt"&&dataset!="ReReco"&&dataset!="Jan16ReReco"&&dataset!="Prompt2012"))
   { throw cms::Exception("CalibratedgsfElectronProducer|ConfigError")<<"Unknown Data dataset" ; }
  cout << "[CalibratedPatElectronProducer] Correcting scale for dataset " << dataset << endl;

 }
 
CalibratedPatElectronProducer::~CalibratedPatElectronProducer()
 {}

void CalibratedPatElectronProducer::produce( edm::Event & event, const edm::EventSetup & setup )
 {

  edm::Handle<ElectronCollection> gsfElectrons ;
  event.getByLabel(inputPatElectrons,gsfElectrons) ;

  std::auto_ptr<ElectronCollection> electrons( new ElectronCollection ) ;
  const ElectronCollection * oldElectrons = gsfElectrons.product() ;

  ElectronCollection::const_iterator electron ;
  ElectronCollection::iterator ele ;
  // first clone the initial collection
  for
   ( electron = oldElectrons->begin() ;
     electron != oldElectrons->end() ;
     ++electron )
   {     
      electrons->push_back(*electron->clone()) ;
   }
  
  ElectronEnergyCalibrator theEnCorrector(dataset, isAOD, isMC, updateEnergyError, debug);

  for
   ( ele = electrons->begin() ;
     ele != electrons->end() ;
     ++ele )
   {     
      // energy calibration for ecalDriven electrons
      if (ele->core()->ecalDrivenSeed()) {        
        theEnCorrector.correct(*ele, event, setup);
      }
      else {
        //std::cout << "[CalibratedPatElectronProducer] is tracker driven only!!" << std::endl;
      }
   }
   event.put(electrons) ;

 }


#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
DEFINE_FWK_MODULE(CalibratedPatElectronProducer);
