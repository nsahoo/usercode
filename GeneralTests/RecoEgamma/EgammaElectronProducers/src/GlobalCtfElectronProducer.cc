// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "RecoEgamma/EgammaElectronProducers/interface/GlobalCtfElectronProducer.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronAlgoA.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h"
#include "DataFormats/EgammaReco/interface/SeedSuperClusterAssociation.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchElectronFwd.h"

#include <iostream>

using namespace reco;
 
GlobalCtfElectronProducer::GlobalCtfElectronProducer(const edm::ParameterSet& iConfig) : conf_(iConfig)
{
  //register your products
  produces<GlobalCtfElectronCollection>();

  //create algo
  algo_ = new
    ElectronAlgoA(iConfig.getParameter<double>("maxEOverPBarrel"),
		  iConfig.getParameter<double>("maxEOverPEndcaps"),
		  iConfig.getParameter<double>("hOverEConeSize"),
		  iConfig.getParameter<double>("maxHOverE"),
		  iConfig.getParameter<double>("maxDeltaEta"),
		  iConfig.getParameter<double>("maxDeltaPhi"),
		  iConfig.getParameter<double>("ptCut"));

}


GlobalCtfElectronProducer::~GlobalCtfElectronProducer()
{
  delete algo_;
}

void GlobalCtfElectronProducer::beginJob(edm::EventSetup const&iSetup) 
{     
  algo_->setupES(iSetup,conf_);  
}

// ------------ method called to produce the data  ------------
void GlobalCtfElectronProducer::produce(edm::Event& e, const edm::EventSetup& iSetup) 
{

  // Create the output collections   
  std::auto_ptr<GlobalCtfElectronCollection> pOutEle(new GlobalCtfElectronCollection);
  
  // invoke algorithm
    algo_->run(e,*pOutEle);

  // put result into the Event
    e.put(pOutEle);
  
}


