#ifndef GlobalCtfElectronProducer_h
#define GlobalCtfElectronProducer_h
  
//
// Package:         RecoEgamma/EgammaElectronProducers
// Class:           GlobalCtfElectronProducer
// 
// Description:   
  
  
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/EDProduct.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <string>

class ElectronAlgoA;

class GlobalCtfElectronProducer : public edm::EDProducer
{
 public:

  explicit GlobalCtfElectronProducer(const edm::ParameterSet& conf);

  virtual ~GlobalCtfElectronProducer();

  virtual void beginJob(edm::EventSetup const&iSetup);
  virtual void produce(edm::Event& e, const edm::EventSetup& c);

 private:

  const edm::ParameterSet conf_;

  ElectronAlgoA* algo_;
  std::string  seedProducer_;
};
#endif
